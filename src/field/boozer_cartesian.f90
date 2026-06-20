module boozer_cartesian
  ! Boozer flux <-> Cartesian map and the Cartesian Larmor displacement shared by
  ! the CP seed (guiding center -> particle) and the CP guiding-center
  ! reconstruction (particle -> guiding center).
  !
  ! The Larmor offset is a true Euclidean vector. Adding it directly to the flux
  ! label (s, vartheta_B, varphi_B) is only first order in rho* and diverges near
  ! the axis, where a fixed physical gyroradius maps to an unbounded angle step
  ! and the radial offset can cross s = 0 (a poloidal reflection the linear add
  ! cannot represent). Here the offset is built and applied in Cartesian, where it
  ! is exact, then the displaced point is mapped back to Boozer by a Newton
  ! inversion of the same forward map.
  !
  ! Forward map: Boozer (s, vartheta_B, varphi_B) -> Cartesian (X, Y, Z) through
  ! the Boozer->VMEC angle map and the splined R, Z, with the analytic Jacobian
  ! d(X,Y,Z)/d(s, vartheta_B, varphi_B). Length unit is the VMEC R, Z unit (cm).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: boozer_to_cart, cart_to_boozer, perp_unit_dir_flux, raise_flux, &
            gc_to_particle, particle_to_gc

contains

  ! Forward map and its Jacobian Jc(a,k) = d x^a / d u^k, u = Boozer coordinate.
  !$acc routine seq
  subroutine boozer_to_cart(u, xyz, Jc)
    use spline_vmec_sub, only: splint_vmec_data_d2
    use boozer_sub, only: delthe_delphi_BV_d2
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: xyz(3), Jc(3,3)

    real(dp) :: s, vartheta_B, varphi_B, theta_V, varphi_V
    real(dp) :: del_t, del_p, ddel_t(3), ddel_p(3), d2del_t(6), d2del_p(6)
    real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Zc, alam
    real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
    real(dp) :: d2R(6), d2Z(6), d2l(6)
    real(dp) :: Jm(3,3), cphi, sphi, dRk, dZk, dphik
    integer :: k

    s = u(1); vartheta_B = u(2); varphi_B = u(3)

    call delthe_delphi_BV_d2(s, vartheta_B, varphi_B, del_t, del_p, &
                             ddel_t, ddel_p, d2del_t, d2del_p)
    theta_V = vartheta_B - del_t
    varphi_V = varphi_B - del_p

    call splint_vmec_data_d2(s, theta_V, varphi_V, A_phi, A_theta, &
         dA_phi_ds, dA_theta_ds, aiota, R, Zc, alam, &
         dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
         dl_ds, dl_dt, dl_dp, d2R, d2Z, d2l)

    ! Angle Jacobian d(s, theta_V, varphi_V)/d(s, vartheta_B, varphi_B).
    Jm = 0.0_dp
    Jm(1,1) = 1.0_dp
    Jm(2,1) = -ddel_t(1); Jm(2,2) = 1.0_dp - ddel_t(2); Jm(2,3) = -ddel_t(3)
    Jm(3,1) = -ddel_p(1); Jm(3,2) = -ddel_p(2);         Jm(3,3) = 1.0_dp - ddel_p(3)

    cphi = cos(varphi_V); sphi = sin(varphi_V)
    xyz(1) = R*cphi
    xyz(2) = R*sphi
    xyz(3) = Zc

    do k = 1, 3
      dRk = dR_ds*Jm(1,k) + dR_dt*Jm(2,k) + dR_dp*Jm(3,k)
      dZk = dZ_ds*Jm(1,k) + dZ_dt*Jm(2,k) + dZ_dp*Jm(3,k)
      dphik = Jm(3,k)
      Jc(1,k) = cphi*dRk - R*sphi*dphik
      Jc(2,k) = sphi*dRk + R*cphi*dphik
      Jc(3,k) = dZk
    end do
  end subroutine boozer_to_cart

  ! Newton inversion of boozer_to_cart: find u (Boozer) with boozer_to_cart(u) =
  ! xyz, started from u_guess. The displacement can be a sizeable fraction of the
  ! major radius (large gyroradius), so a plain Newton step overshoots and can
  ! cross the axis; a backtracking line search keeps the Cartesian residual
  ! monotonically decreasing and the iterate inside 0 < s < 1. ierr = 1 on
  ! non-convergence.
  !$acc routine seq
  subroutine cart_to_boozer(xyz, u_guess, u, ierr)
    real(dp), intent(in) :: xyz(3), u_guess(3)
    real(dp), intent(out) :: u(3)
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 100, maxls = 40
    ! tol is the converged Cartesian residual (length unit cm). 1e-7 cm is far
    ! below any physically relevant displacement yet above the spline noise floor
    ! (~1e-7 relative on a metre-scale R). ok_tol accepts a line-search stall that
    ! has already reached sub-micron accuracy; only a genuine non-convergence
    ! (residual still O(cm)) returns ierr/=0.
    real(dp), parameter :: tol = 1.0e-7_dp, ok_tol = 1.0e-3_dp
    real(dp) :: xc(3), Jc(3,3), res(3), du(3), ut(3), rnew, rn, alpha
    integer :: it, ls

    ierr = 1
    u = u_guess
    call boozer_to_cart(u, xc, Jc)
    res = xc - xyz
    rn = maxval(abs(res))
    do it = 1, maxit
      if (rn < tol) then
        ierr = 0
        return
      end if
      call solve3(Jc, -res, du)
      alpha = 1.0_dp
      do ls = 1, maxls
        ut = u + alpha*du
        if (ut(1) <= 0.0_dp) ut(1) = 1.0e-8_dp
        if (ut(1) >= 1.0_dp) ut(1) = 1.0_dp - 1.0e-8_dp
        call boozer_to_cart(ut, xc, Jc)
        rnew = maxval(abs(xc - xyz))
        if (rnew < rn) exit
        alpha = 0.5_dp*alpha
      end do
      if (rnew >= rn) then         ! line search stalled at the residual floor
        if (rn < ok_tol) ierr = 0
        return
      end if
      u = ut
      res = xc - xyz
      rn = rnew
    end do
    if (rn < ok_tol) ierr = 0
  end subroutine cart_to_boozer

  ! Particle position from guiding center. The CP velocity is seeded at the
  ! PARTICLE point (cpp_canon_init builds the perpendicular direction from the
  ! field there), so the offset must be consistent with that: solve the fixed
  ! point x_p with cart(x_p) - rho(x_p) = cart(x_gc), where rho(x_p) uses the
  ! perpendicular direction and field AT x_p (the same the velocity seed uses).
  ! particle_to_gc(x_p) then returns x_gc exactly, and the seed offset and the
  ! reconstruction are true inverses. vperp0 is the perpendicular speed magnitude;
  ! the gyrophase direction is the shared perp_unit_dir_flux at x_p.
  !$acc routine seq
  subroutine gc_to_particle(x_gc, vperp0, mass, qc, x_particle, ierr)
    use boozer_field_metric, only: boozer_field_metric_eval
    real(dp), intent(in) :: x_gc(3), vperp0, mass, qc
    real(dp), intent(out) :: x_particle(3)
    integer, intent(out) :: ierr
    integer, parameter :: maxfp = 50
    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: xyz_gc(3), xyz(3), Jc(3,3), rho(3), xnew(3)
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: eperp(3), vperp_con(3)
    integer :: it

    call boozer_to_cart(x_gc, xyz_gc, Jc)
    x_particle = x_gc
    do it = 1, maxfp
      call boozer_field_metric_eval(x_particle, g, ginv, sqrtg, dg, Acov, dA, &
                                    Bctr, Bcov, Bmod, dBmod, hcov)
      call perp_unit_dir_flux(g, ginv, hcov, eperp)
      vperp_con = vperp0*eperp
      call boozer_to_cart(x_particle, xyz, Jc)
      call larmor_vector_cart(x_particle, vperp_con, mass, qc, Jc, rho)
      call cart_to_boozer(xyz_gc + rho, x_particle, xnew, ierr)
      if (ierr /= 0) return
      if (maxval(abs(xnew - x_particle)) < tol) then
        x_particle = xnew
        return
      end if
      x_particle = xnew
    end do
    ierr = 0
  end subroutine gc_to_particle

  ! Guiding center from particle: x_gc = x_p - rho, rho built from the particle's
  ! perpendicular velocity vperp_con (contravariant flux components) at x_p.
  !$acc routine seq
  subroutine particle_to_gc(x_particle, vperp_con, mass, qc, x_gc, ierr)
    real(dp), intent(in) :: x_particle(3), vperp_con(3), mass, qc
    real(dp), intent(out) :: x_gc(3)
    integer, intent(out) :: ierr
    real(dp) :: xyz(3), Jc(3,3), rho(3)

    call boozer_to_cart(x_particle, xyz, Jc)
    call larmor_vector_cart(x_particle, vperp_con, mass, qc, Jc, rho)
    call cart_to_boozer(xyz - rho, x_particle, x_gc, ierr)
  end subroutine particle_to_gc

  ! Cartesian Larmor vector rho = (mass/(qc |B|)) (b x v_perp), with b and v_perp
  ! pushed to Cartesian by the map Jacobian Jc at the Boozer point u. Same sign
  ! convention as the legacy flux offset (h x v_perp), so the seed and the
  ! reconstruction are inverse: x_p = x_gc + rho, x_gc = x_p - rho.
  !$acc routine seq
  subroutine larmor_vector_cart(u, vperp_con, mass, qc, Jc, rho)
    use boozer_field_metric, only: boozer_field_metric_eval
    real(dp), intent(in) :: u(3), vperp_con(3), mass, qc, Jc(3,3)
    real(dp), intent(out) :: rho(3)
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: Bcart(3), bhat(3), vperp_cart(3), Bnrm, factor
    integer :: a, i

    call boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
                                  Bctr, Bcov, Bmod, dBmod, hcov)

    ! Push the contravariant field and perpendicular velocity to Cartesian.
    do a = 1, 3
      Bcart(a) = Jc(a,1)*Bctr(1) + Jc(a,2)*Bctr(2) + Jc(a,3)*Bctr(3)
      vperp_cart(a) = Jc(a,1)*vperp_con(1) + Jc(a,2)*vperp_con(2) &
                    + Jc(a,3)*vperp_con(3)
    end do
    Bnrm = sqrt(Bcart(1)**2 + Bcart(2)**2 + Bcart(3)**2)
    if (Bnrm <= 0.0_dp) error stop 'larmor_vector_cart: zero |B|'
    do i = 1, 3
      bhat(i) = Bcart(i)/Bnrm
    end do

    factor = mass/(qc*Bmod)
    rho(1) = factor*(bhat(2)*vperp_cart(3) - bhat(3)*vperp_cart(2))
    rho(2) = factor*(bhat(3)*vperp_cart(1) - bhat(1)*vperp_cart(3))
    rho(3) = factor*(bhat(1)*vperp_cart(2) - bhat(2)*vperp_cart(1))
  end subroutine larmor_vector_cart

  ! Unit perpendicular direction in contravariant flux components: take the raised
  ! radial covector e_r = g^{i1}, project out the field-parallel part, normalize in
  ! the metric. Shared by the CP velocity seed and the position offset so both use
  ! the same gyrophase reference.
  !$acc routine seq
  subroutine perp_unit_dir_flux(g, ginv, hcov, eperp)
    real(dp), intent(in) :: g(3,3), ginv(3,3), hcov(3)
    real(dp), intent(out) :: eperp(3)
    real(dp) :: er(3), hcon(3), hpar, nrm
    integer :: i, j

    er = [ginv(1,1), ginv(2,1), ginv(3,1)]
    call raise_flux(ginv, hcov, hcon)
    hpar = hcov(1)*er(1) + hcov(2)*er(2) + hcov(3)*er(3)
    eperp = er - hpar*hcon

    nrm = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        nrm = nrm + g(i,j)*eperp(i)*eperp(j)
      end do
    end do
    if (nrm <= 0.0_dp) error stop 'perp_unit_dir_flux: degenerate direction'
    eperp = eperp/sqrt(nrm)
  end subroutine perp_unit_dir_flux

  !$acc routine seq
  subroutine raise_flux(ginv, vcov, vcon)
    real(dp), intent(in) :: ginv(3,3), vcov(3)
    real(dp), intent(out) :: vcon(3)
    integer :: i
    do i = 1, 3
      vcon(i) = ginv(i,1)*vcov(1) + ginv(i,2)*vcov(2) + ginv(i,3)*vcov(3)
    end do
  end subroutine raise_flux

  ! Solve A x = b for 3x3 A by cofactor inverse.
  !$acc routine seq
  subroutine solve3(A, b, x)
    real(dp), intent(in) :: A(3,3), b(3)
    real(dp), intent(out) :: x(3)
    real(dp) :: inv(3,3), det

    det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
        - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
        + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
    if (det == 0.0_dp) error stop 'solve3: singular Jacobian'
    inv(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2))/det
    inv(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3))/det
    inv(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2))/det
    inv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))/det
    inv(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1))/det
    inv(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3))/det
    inv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))/det
    inv(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2))/det
    inv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/det
    x(1) = inv(1,1)*b(1) + inv(1,2)*b(2) + inv(1,3)*b(3)
    x(2) = inv(2,1)*b(1) + inv(2,2)*b(2) + inv(2,3)*b(3)
    x(3) = inv(3,1)*b(1) + inv(3,2)*b(2) + inv(3,3)*b(3)
  end subroutine solve3

end module boozer_cartesian
