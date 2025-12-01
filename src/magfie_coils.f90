module magfie_coils_sub

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use field_coils, only : CoilsField, create_coils_field
  use neo_biotsavart, only : compute_magnetic_field

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
    ! Coordinates x are Cartesian:
    !   x(1) = X [cm]
    !   x(2) = Y [cm]
    !   x(3) = Z [cm]
    !
    ! The magnetic field is evaluated from the coils via Biot–Savart in
    ! Cartesian space. The guiding-centre equations use a trivial metric
    ! (sqrtg = 1, covariant = contravariant) in these coordinates.
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

    real(dp), parameter :: hfd = 1.0d-3

    real(dp) :: x_base(3)
    real(dp) :: bmod_base
    real(dp) :: hcov_base(3), hctr_base(3)

    real(dp) :: x_plus(3), x_minus(3)
    real(dp) :: b_plus, b_minus
    real(dp) :: hcov_plus(3), hcov_minus(3), hctr_tmp(3)
    real(dp) :: dhdq(3,3)
    real(dp) :: delta(3)

    integer :: i, j, k

    if (.not. allocated(coils_field_gc)) then
      print *, 'magfie_coils: coils field not initialized'
      error stop
    end if

    ! Base state
    x_base = x
    call eval_coils_state(x_base, bmod_base, hcov_base, hctr_base)

    bmod  = bmod_base
    sqrtg = 1.0d0
    hcovar = hcov_base
    hctrvr = hctr_base

    delta = hfd

    ! Finite-difference derivatives for log(B) and hcovar in Cartesian coordinates
    do i = 1, 3
      x_plus  = x_base
      x_minus = x_base

      x_plus(i)  = x_plus(i)  + delta(i)
      x_minus(i) = x_minus(i) - delta(i)

      call eval_coils_state(x_plus,  b_plus,  hcov_plus,  hctr_tmp)
      call eval_coils_state(x_minus, b_minus, hcov_minus, hctr_tmp)

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
      hcurl(i) = hcurl(i)/sqrtg
    end do
  end subroutine magfie_coils


  subroutine eval_coils_state(x_cart, bmod, hcovar, hctrvr)
    ! Evaluate B and unit vector h at a single point in Cartesian coordinates.
    real(dp), intent(in) :: x_cart(3)
    real(dp), intent(out) :: bmod
    real(dp), intent(out) :: hcovar(3), hctrvr(3)

    real(dp) :: Bcart(3), normB

    Bcart = compute_magnetic_field(coils_field_gc%coils, x_cart)

    normB = sqrt(Bcart(1)**2 + Bcart(2)**2 + Bcart(3)**2)
    if (normB <= 0.0d0) then
      print *, 'eval_coils_state: non-positive |B|'
      error stop
    end if

    bmod = normB
    hcovar = Bcart / normB
    hctrvr = hcovar
  end subroutine eval_coils_state

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
