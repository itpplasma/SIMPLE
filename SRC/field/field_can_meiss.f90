module field_can_meiss

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: SplineData3D, construct_splines_3d, &
    evaluate_splines_3d, evaluate_splines_3d_der, evaluate_splines_3d_der2
use util, only: twopi
use field_can_base, only: FieldCan, n_field_evaluations
use simple_magfie, only: MagneticField

implicit none

class(MagneticField), allocatable :: magfie
integer :: n_r=62, n_th=63, n_phi=64
real(dp) :: xmin(3) = [1d-1, 0d0, 0d0]  ! TODO check limits
real(dp) :: xmax(3) = [1d0, twopi, twopi]

real(dp) :: h_r, h_th, h_phi

! For splining covariant vector potential, h=B/Bmod and Bmod over canonical coordinates
type(SplineData3D) :: spl_Bmod, spl_Ath, spl_Aph, spl_hth, spl_hph

! For splining lambda (difference between canonical and toroidal cylinder angle)
! and chi (gauge transformation)
real(dp), dimension(:,:,:), allocatable :: lam_phi, chi_gauge
type(SplineData3D) :: spl_lam_phi, spl_chi_gauge

integer, parameter :: order(3) = [5, 5, 5]
logical, parameter :: periodic(3) = [.False., .True., .True.]

contains

subroutine init_meiss(magfie_, n_r_, n_th_, n_phi_, rmin, rmax, thmin, thmax)
    use new_vmec_stuff_mod, only : nper

    class(MagneticField), intent(in) :: magfie_
    integer, intent(in), optional :: n_r_, n_th_, n_phi_
    real(dp), intent(in), optional :: rmin, rmax, thmin, thmax

    if (allocated(magfie)) deallocate(magfie)
    allocate(magfie, source=magfie_)

    if (present(n_r_)) n_r = n_r_
    if (present(n_th_)) n_th = n_th_
    if (present(n_phi_)) n_phi = n_phi_

    if (present(rmin)) xmin(1) = rmin
    if (present(rmax)) xmax(1) = rmax
    if (present(thmin)) xmin(2) = thmin
    if (present(thmax)) xmax(2) = thmax
    xmax(3) = twopi/nper

    h_r = (xmax(1)-xmin(1))/(n_r-1)
    h_th = (xmax(2)-xmin(2))/(n_th-1)
    h_phi = (xmax(3)-xmin(3))/(n_phi-1)
end subroutine init_meiss


subroutine evaluate_meiss(f, r, th_c, ph_c, mode_secders)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    real(dp) :: x(3)

    n_field_evaluations = n_field_evaluations + 1

    x = [r, th_c, ph_c]

    if (mode_secders > 0) then

        call evaluate_splines_3d_der2(spl_Ath, x, f%Ath, f%dAth, f%d2Ath)
        call evaluate_splines_3d_der2(spl_Aph, x, f%Aph, f%dAph, f%d2Aph)

        call evaluate_splines_3d_der2(spl_hth, x, f%hth, f%dhth, f%d2hth)
        call evaluate_splines_3d_der2(spl_hph, x, f%hph, f%dhph, f%d2hph)

        call evaluate_splines_3d_der2(spl_Bmod, x, f%Bmod, f%dBmod, f%d2Bmod)

        return
    end if

    call evaluate_splines_3d_der(spl_Ath, x, f%Ath, f%dAth)
    call evaluate_splines_3d_der(spl_Aph, x, f%Aph, f%dAph)

    call evaluate_splines_3d_der(spl_hth, x, f%hth, f%dhth)
    call evaluate_splines_3d_der(spl_hph, x, f%hph, f%dhph)

    call evaluate_splines_3d_der(spl_Bmod, x, f%Bmod, f%dBmod)
end subroutine evaluate_meiss


subroutine can_to_ref_meiss(xcan, xref)
    real(dp), intent(in) :: xcan(3)
    real(dp), intent(out) :: xref(3)
    real(dp) :: lam

    call evaluate_splines_3d(spl_lam_phi, xcan, lam)
    xref(1) = xcan(1)
    xref(2) = modulo(xcan(2), twopi)
    xref(3) = modulo(xcan(3) + lam, twopi)
end subroutine can_to_ref_meiss


subroutine ref_to_can_meiss(xref, xcan)
    real(dp), intent(in) :: xref(3)
    real(dp), intent(out) :: xcan(3)

    real(dp), parameter :: TOL = 1d-12
    integer, parameter :: MAX_ITER = 16

    real(dp) :: lam, dlam(3), phi_can_prev
    integer :: i

    xcan(1) = xref(1)
    xcan(2) = modulo(xref(2), twopi)
    xcan(3) = modulo(xref(3), twopi)

    do i=1, MAX_ITER
        call evaluate_splines_3d_der(spl_lam_phi, xcan, lam, dlam)
        phi_can_prev = xcan(3)
        xcan(3) = phi_can_prev - (phi_can_prev + lam - xref(3))/(1d0 + dlam(3))
        if (abs(xcan(3) - phi_can_prev) < TOL) return
    enddo
    print *, 'WARNING: ref_to_can_meiss did not converge after', MAX_ITER, 'iterations'
end subroutine ref_to_can_meiss


subroutine get_meiss_coordinates
    print *, 'field_can_meiss.init_transformation'
    call init_transformation

    print *, 'field_can_meiss.spline_transformation'
    call spline_transformation

    print *, 'field_can_meiss.init_canonical_field_components'
    call init_canonical_field_components
end subroutine get_meiss_coordinates


subroutine init_transformation
    real(dp) :: y(2)
    integer :: i_r, i_th, i_phi, i_ctr

    if (allocated(lam_phi)) deallocate(lam_phi, chi_gauge)
    allocate(lam_phi(n_r, n_th, n_phi), chi_gauge(n_r, n_th, n_phi))

    i_ctr = 0

    !$omp parallel private(i_r, i_th, i_phi, y)
    !$omp do
    do i_phi=1,n_phi
        !$omp critical
        i_ctr = i_ctr + 1
        call print_progress
        !$omp end critical
        do i_th=1,n_th
            lam_phi(1, i_th, i_phi) = 0d0
            chi_gauge(1, i_th, i_phi) = 0d0
            y = 0d0
            do i_r=2,n_r
                call integrate(i_r, i_th, i_phi, y)
            enddo
        enddo
    enddo
    !$omp end do
    !$omp end parallel

    contains

    subroutine print_progress
        write(*,'(A, I4, A, I4)',advance='no') 'integrate ODE: ', i_ctr, ' of ', n_phi
        if (i_ctr < n_phi) then
            write(*, '(A)', advance="no") char(13)
        else
            write(*, *)
        end if
    end subroutine print_progress
end subroutine init_transformation

subroutine integrate(i_r, i_th, i_phi, y)
    use odeint_sub, only: odeint_allroutines

    integer, intent(in) :: i_r, i_th, i_phi
    real(dp), dimension(2), intent(inout) :: y

    real(dp), parameter :: relerr=1d-11
    real(dp) :: r1, r2
    integer :: ndim=2

    r1 = xmin(1) + h_r*(i_r-2)
    r2 = xmin(1) + h_r*(i_r-1)

    call odeint_allroutines(y, ndim, r1, r2, relerr, rh_can_closure)

    lam_phi(i_r, i_th, i_phi) = y(1)
    chi_gauge(i_r, i_th, i_phi) = y(2)

    contains

    subroutine rh_can_closure(r_c, z, dz)
        real(dp), intent(in) :: r_c
        real(dp), dimension(2), intent(in) :: z
        real(dp), dimension(2), intent(inout) :: dz

        call rh_can(r_c, z, dz, i_th, i_phi)
    end subroutine rh_can_closure

end subroutine integrate

subroutine rh_can(r_c, z, dz, i_th, i_phi)
    real(dp), intent(in) :: r_c
    real(dp), dimension(2), intent(in) :: z  ! lam_phi, chi_gauge
    real(dp), dimension(2), intent(inout) :: dz
    integer, intent(in) :: i_th, i_phi
    real(dp) :: phi_c
    real(dp) :: Ar, Ap, hr, hp

    phi_c = xmin(3) + h_phi*(i_phi-1)
    call ah_cov_on_slice(r_c, modulo(phi_c + z(1), twopi), i_th, Ar, Ap, hr, hp)

    dz(1) = -hr/hp
    dz(2) = Ar + Ap*dz(1)
end subroutine rh_can


subroutine spline_transformation
    call construct_splines_3d(xmin, xmax, lam_phi, order, periodic, spl_lam_phi)
    call construct_splines_3d(xmin, xmax, chi_gauge, order, periodic, spl_chi_gauge)
end subroutine spline_transformation


subroutine ah_cov_on_slice(r, phi, i_th, Ar, Ap, hr, hp)
    real(dp), intent(in) :: r, phi
    integer, intent(in) :: i_th
    real(dp), intent(out) :: Ar, Ap, hr, hp

    real(dp), dimension(3) :: Acov, hcov
    real(dp) :: Bmod
    real(dp) :: th

    ! TODO: Make this more efficient with slices
    ! TODO: Support not only VMEC field
    th = xmin(2) + h_th*(i_th-1)
    call magfie%evaluate([r, th, phi], Acov, hcov, Bmod)

    Ar = Acov(1)
    Ap = Acov(3)
    hr = hcov(1)
    hp = hcov(3)
end subroutine ah_cov_on_slice


subroutine init_canonical_field_components
    real(dp) :: xcan(3)
    real(dp), dimension(:,:,:), allocatable :: Ath, Aphi, hth, hphi, Bmod
    real(dp) :: xref(3), Acov(3), hcov(3), lam, dlam(3), chi, dchi(3)
    integer :: i_r, i_th, i_phi

    allocate(Ath(n_r,n_th,n_phi), Aphi(n_r,n_th,n_phi))
    allocate(hth(n_r,n_th,n_phi), hphi(n_r,n_th,n_phi))
    allocate(Bmod(n_r,n_th,n_phi))

    do i_phi=1,n_phi
        do i_th=1,n_th
            do i_r=1,n_r
                xcan = get_grid_point(i_r, i_th, i_phi)

                call evaluate_splines_3d_der(spl_lam_phi, xcan, lam, dlam)
                call evaluate_splines_3d_der(spl_chi_gauge, xcan, chi, dchi)

                xref(1) = xcan(1)
                xref(2) = modulo(xcan(2), twopi)
                xref(3) = modulo(xcan(3) + lam, twopi)

                call magfie%evaluate(xref, Acov, hcov, Bmod(i_r, i_th, i_phi))

                Ath(i_r, i_th, i_phi) = Acov(2) + Acov(3)*dlam(2) - dchi(2)
                Aphi(i_r, i_th, i_phi) = Acov(3)*(1.0d0 + dlam(3)) - dchi(3)
                hth(i_r, i_th, i_phi) = hcov(2) + hcov(3)*dlam(2)
                hphi(i_r, i_th, i_phi) = hcov(3)*(1.0d0 + dlam(3))
            end do
        end do
    end do

    Ath = Ath - Ath(1,1,1)
    Aphi = Aphi - Aphi(1,1,1)

    call construct_splines_3d(xmin, xmax, Ath, order, periodic, spl_Ath)
    call construct_splines_3d(xmin, xmax, Aphi, order, periodic, spl_Aph)
    call construct_splines_3d(xmin, xmax, hth, order, periodic, spl_hth)
    call construct_splines_3d(xmin, xmax, hphi, order, periodic, spl_hph)
    call construct_splines_3d(xmin, xmax, Bmod, order, periodic, spl_Bmod)
end subroutine init_canonical_field_components


pure function get_grid_point(i_r, i_th, i_phi)
    integer, intent(in) :: i_r, i_th, i_phi
    real(dp) :: get_grid_point(3)

    get_grid_point = [ &
        xmin(1) + h_r*dble(i_r-1), &
        xmin(2) + h_th*dble(i_th-1), &
        xmin(3) + h_phi*dble(i_phi-1) &
    ]
end function get_grid_point


subroutine magfie_meiss(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!  Computes magnetic field and derivatives with bmod in units of the magnetic code
!
!  Input parameters:
!            formal:  x               - array of canonicalized coordinates r, th, ph
!  Output parameters:
!            formal:  bmod            - magnetic field module
!                     sqrtg           - metric determinant
!                     bder            - covariant components of (grad B)/B
!                     hcovar          - covariant components of \bB/B
!                     hctrvr          - contravariant components of \bB/B
!                     hcurl           - contravariant components of curl (\bB/B)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

    type(FieldCan) :: f
    real(dp) :: sqrtg_bmod

    call evaluate_meiss(f, x(1), x(2), x(3), 0)

    bmod = f%Bmod

    sqrtg_bmod = f%hph*f%dAth(1) - f%hth*f%dAph(1)
    sqrtg = sqrtg_bmod/bmod
    bder = f%dBmod/bmod

    hcovar(1) = 0d0
    hcovar(2) = f%hth
    hcovar(3) = f%hph

    hctrvr(1) = (f%dAph(2) - f%dAth(3))/sqrtg_bmod
    hctrvr(2) = -f%dAph(1)/sqrtg_bmod
    hctrvr(3) = f%dAth(1)/sqrtg_bmod

    hcurl(1) = (f%dhph(2) - f%dhth(3))/sqrtg
    hcurl(2) = -f%dhph(1)/sqrtg
    hcurl(3) = f%dhth(1)/sqrtg
end subroutine magfie_meiss

end module field_can_meiss
