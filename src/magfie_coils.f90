module magfie_coils_sub

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use field_coils, only : CoilsField, create_coils_field
  use simple_coordinates, only : transform_vmec_to_cart

  implicit none

  class(CoilsField), allocatable :: coils_field_gc

contains

  subroutine init_magfie_coils_from_file(filename)
    character(*), intent(in) :: filename

    if (allocated(coils_field_gc)) then
      deallocate(coils_field_gc)
    end if

    ! For guiding-centre coils mode we only need direct Biot–Savart
    ! evaluation via evaluate_direct; spline precomputation is not required.
    call create_coils_field(filename, coils_field_gc, should_spline = .false.)
  end subroutine init_magfie_coils_from_file


  subroutine magfie_coils(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    ! Guiding-centre magnetic field interface for direct coils mode.
    !
    ! Coordinates x are VMEC reference coordinates:
    !   x(1) = s (normalized toroidal flux)
    !   x(2) = theta (VMEC poloidal angle)
    !   x(3) = phi   (geometrical toroidal angle)
    !
    ! The magnetic field is evaluated from the coils via Biot–Savart in
    ! Cartesian space and then projected onto the VMEC coordinate basis.
    !
    ! All geometric information (metric, Jacobian) comes from the VMEC
    ! transformation, so the guiding-centre equations remain in the
    ! familiar (s, theta, phi) reference frame.
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

    real(dp), parameter :: hs = 1.0d-3
    real(dp), parameter :: ht = hs*2.0d0*3.14159265358979d0
    real(dp), parameter :: hp = ht/5.0d0

    real(dp) :: x_base(3)
    real(dp) :: bmod_base
    real(dp) :: hcov_base(3), hctr_base(3)
    real(dp) :: sqrtg_base

    real(dp) :: x_plus(3), x_minus(3)
    real(dp) :: b_plus, b_minus
    real(dp) :: hcov_plus(3), hcov_minus(3)
    real(dp) :: sqrtg_tmp
    real(dp) :: hctr_tmp(3)
    real(dp) :: dhdq(3,3)
    real(dp) :: delta(3)

    integer :: i, j, k

    if (.not. allocated(coils_field_gc)) then
      print *, 'magfie_coils: coils field not initialized'
      error stop
    end if

    ! Base state
    x_base = x
    call eval_coils_state(x_base, bmod_base, sqrtg_base, hcov_base, hctr_base)

    bmod  = bmod_base
    sqrtg = sqrtg_base
    hcovar = hcov_base
    hctrvr = hctr_base

    ! Step sizes in (s, theta, phi)
    delta(1) = hs
    delta(2) = ht
    delta(3) = hp

    ! Finite-difference derivatives for log(B) and hcovar
    do i = 1, 3
      x_plus  = x_base
      x_minus = x_base

      x_plus(i)  = x_plus(i)  + delta(i)
      x_minus(i) = x_minus(i) - delta(i)

      ! Enforce periodicity in angles
      if (i == 2 .or. i == 3) then
        x_plus(i)  = modulo(x_plus(i),  2.0d0*3.14159265358979d0)
        x_minus(i) = modulo(x_minus(i), 2.0d0*3.14159265358979d0)
      else
        x_plus(i)  = max(0.0d0, min(1.0d0, x_plus(i)))
        x_minus(i) = max(0.0d0, min(1.0d0, x_minus(i)))
      end if

      call eval_coils_state(x_plus,  b_plus,  sqrtg_tmp, hcov_plus,  hctr_tmp)
      call eval_coils_state(x_minus, b_minus, sqrtg_tmp, hcov_minus, hctr_tmp)

      if (b_plus <= 0.0d0 .or. b_minus <= 0.0d0) then
        print *, 'magfie_coils: non-positive B encountered in finite differences'
        error stop
      end if

      bder(i) = (log(b_plus) - log(b_minus))/(2.0d0*delta(i))

      do j = 1, 3
        dhdq(j,i) = (hcov_plus(j) - hcov_minus(j))/(2.0d0*delta(i))
      end do
    end do

    ! Curl of h in contravariant components:
    !   (curl h)^i = (1/sqrtg) * epsilon^{ijk} * d h_k / d q_j
    do i = 1, 3
      hcurl(i) = 0.0d0
      do j = 1, 3
        do k = 1, 3
          hcurl(i) = hcurl(i) + levi_civita(i,j,k)*dhdq(k,j)
        end do
      end do
      hcurl(i) = hcurl(i)/sqrtg_base
    end do
  end subroutine magfie_coils


  subroutine eval_coils_state(x_ref, bmod, sqrtg, hcovar, hctrvr)
    ! Evaluate B, metric and covariant/contravariant components of h at a
    ! single point in VMEC coordinates x_ref = (s, theta, phi).
    real(dp), intent(in) :: x_ref(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), intent(out) :: hcovar(3), hctrvr(3)

    real(dp) :: x_vmec(3)
    real(dp) :: r
    real(dp) :: Acov(3), hcov(3)
    real(dp) :: xcart(3), dxcart_dxvmec(3,3)
    real(dp) :: e1(3), e2(3), e3(3)
    real(dp) :: g(3,3), g_inv(3,3)
    real(dp) :: J
    real(dp) :: Bcov(3), Bctr(3)

    ! Map (s, theta, phi) to the radius r = sqrt(s) used by
    ! CoilsField%evaluate_direct.
    r = sqrt(max(x_ref(1), 0.0d0))
    x_vmec = [r, x_ref(2), x_ref(3)]

    call coils_field_gc%evaluate_direct(x_vmec, Acov, hcov, bmod)

    ! Geometry from VMEC transformation at (s, theta, phi)
    call transform_vmec_to_cart(x_ref, xcart, dxcart_dxvmec)

    e1 = dxcart_dxvmec(:,1)
    e2 = dxcart_dxvmec(:,2)
    e3 = dxcart_dxvmec(:,3)

    g(1,1) = dot_product(e1,e1)
    g(1,2) = dot_product(e1,e2)
    g(1,3) = dot_product(e1,e3)
    g(2,1) = g(1,2)
    g(2,2) = dot_product(e2,e2)
    g(2,3) = dot_product(e2,e3)
    g(3,1) = g(1,3)
    g(3,2) = g(2,3)
    g(3,3) = dot_product(e3,e3)

    ! Jacobian determinant J = e1 · (e2 × e3)
    J = e1(1)*(e2(2)*e3(3) - e2(3)*e3(2)) - &
        e1(2)*(e2(1)*e3(3) - e2(3)*e3(1)) + &
        e1(3)*(e2(1)*e3(2) - e2(2)*e3(1))

    sqrtg = abs(J)

    call invert_metric(g, g_inv)

    ! CoilsField%evaluate_direct returns covariant h components
    ! with hcov(1) already scaled consistently; just use them.
    hcovar = hcov

    ! B_cov = B * h_cov
    Bcov = bmod*hcovar

    ! Contravariant B and h
    Bctr(1) = g_inv(1,1)*Bcov(1) + g_inv(1,2)*Bcov(2) + g_inv(1,3)*Bcov(3)
    Bctr(2) = g_inv(2,1)*Bcov(1) + g_inv(2,2)*Bcov(2) + g_inv(2,3)*Bcov(3)
    Bctr(3) = g_inv(3,1)*Bcov(1) + g_inv(3,2)*Bcov(2) + g_inv(3,3)*Bcov(3)

    hctrvr = Bctr/bmod
  end subroutine eval_coils_state


  subroutine invert_metric(g, g_inv)
    real(dp), intent(in) :: g(3,3)
    real(dp), intent(out) :: g_inv(3,3)
    real(dp) :: det

    det = g(1,1)*(g(2,2)*g(3,3)-g(2,3)*g(3,2)) - &
          g(1,2)*(g(2,1)*g(3,3)-g(2,3)*g(3,1)) + &
          g(1,3)*(g(2,1)*g(3,2)-g(2,2)*g(3,1))

    if (abs(det) < 1.0d-20) then
      print *, 'invert_metric: singular metric tensor'
      error stop
    end if

    g_inv(1,1) = (g(2,2)*g(3,3)-g(2,3)*g(3,2))/det
    g_inv(1,2) = (g(1,3)*g(3,2)-g(1,2)*g(3,3))/det
    g_inv(1,3) = (g(1,2)*g(2,3)-g(1,3)*g(2,2))/det

    g_inv(2,1) = (g(2,3)*g(3,1)-g(2,1)*g(3,3))/det
    g_inv(2,2) = (g(1,1)*g(3,3)-g(1,3)*g(3,1))/det
    g_inv(2,3) = (g(1,3)*g(2,1)-g(1,1)*g(2,3))/det

    g_inv(3,1) = (g(2,1)*g(3,2)-g(2,2)*g(3,1))/det
    g_inv(3,2) = (g(1,2)*g(3,1)-g(1,1)*g(3,2))/det
    g_inv(3,3) = (g(1,1)*g(2,2)-g(1,2)*g(2,1))/det
  end subroutine invert_metric


  pure function levi_civita(i, j, k) result(eps)
    integer, intent(in) :: i, j, k
    real(dp) :: eps

    if (i == j .or. j == k .or. i == k) then
      eps = 0.0d0
    else if ((i == 1 .and. j == 2 .and. k == 3) .or. &
             (i == 2 .and. j == 3 .and. k == 1) .or. &
             (i == 3 .and. j == 1 .and. k == 2)) then
      eps = 1.0d0
    else
      eps = -1.0d0
    end if
  end function levi_civita

end module magfie_coils_sub
