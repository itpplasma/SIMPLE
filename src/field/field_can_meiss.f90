module field_can_meiss

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: &
    BatchSplineData3D, construct_batch_splines_3d, &
    evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
    evaluate_batch_splines_3d_der2
use util, only: twopi
use field_can_base, only: FieldCan, n_field_evaluations
use field, only: MagneticField
use field_geoflux, only: geoflux_is_analytical

implicit none

type :: grid_indices_t
    integer :: i_th, i_phi
end type grid_indices_t

class(MagneticField), allocatable :: field_noncan
integer :: n_r=62, n_th=63, n_phi=64
real(dp) :: xmin(3) = [1d-3, 0d0, 0d0]  ! Avoid singular behavior at s -> 0
real(dp) :: xmax(3) = [1d0, twopi, twopi]

real(dp) :: h_r, h_th, h_phi
real(dp), parameter :: hr_small = 1.0d-3

! Batch spline for optimized field evaluation (5 components: Ath, Aph, hth, hph, Bmod)
type(BatchSplineData3D) :: spl_field_batch
logical :: batch_splines_initialized = .false.

! For splining lambda (difference between canonical and toroidal cylinder angle)
! and chi (gauge transformation)
real(dp), dimension(:,:,:), allocatable :: lam_phi, chi_gauge
! Batch spline for transformation components (2 components: lam_phi, chi_gauge)
type(BatchSplineData3D) :: spl_transform_batch
logical :: transform_splines_initialized = .false.

integer, parameter :: order(3) = [5, 5, 5]  ! Restored to original 5th order for accuracy
logical, parameter :: periodic(3) = [.False., .True., .True.]

contains

subroutine rh_can_wrapper(r_c, z, dz, context)
    real(dp), intent(in) :: r_c
    real(dp), dimension(:), intent(in) :: z
    real(dp), dimension(:), intent(out) :: dz
    class(*), intent(in) :: context
    
    integer :: i_th, i_phi
    
    select type(context)
    type is (grid_indices_t)
        i_th = context%i_th
        i_phi = context%i_phi
    end select
    
    call rh_can(r_c, z, dz, i_th, i_phi)
end subroutine rh_can_wrapper

subroutine init_meiss(field_noncan_, n_r_, n_th_, n_phi_, rmin, rmax, thmin, thmax)
    use new_vmec_stuff_mod, only : nper

    class(MagneticField), intent(in) :: field_noncan_
    integer, intent(in), optional :: n_r_, n_th_, n_phi_
    real(dp), intent(in), optional :: rmin, rmax, thmin, thmax

    if (allocated(field_noncan)) deallocate(field_noncan)
    allocate(field_noncan, source=field_noncan_)

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


subroutine cleanup_meiss()
    ! Clean up batch splines to prevent memory leaks
    use interpolate, only: destroy_batch_splines_3d
    
    if (batch_splines_initialized) then
        call destroy_batch_splines_3d(spl_field_batch)
        batch_splines_initialized = .false.
    end if
    
    if (transform_splines_initialized) then
        call destroy_batch_splines_3d(spl_transform_batch)
        transform_splines_initialized = .false.
    end if
    
    ! Clean up allocated arrays
    if (allocated(lam_phi)) deallocate(lam_phi)
    if (allocated(chi_gauge)) deallocate(chi_gauge)
    if (allocated(field_noncan)) deallocate(field_noncan)
end subroutine cleanup_meiss


subroutine evaluate_meiss(f, r, th_c, ph_c, mode_secders)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    real(dp) :: x(3)

    n_field_evaluations = n_field_evaluations + 1

    x = [r, th_c, ph_c]

    if (mode_secders > 0) then
        call evaluate_meiss_batch_der2(f, x)
        return
    end if

    call evaluate_meiss_batch_der(f, x)
end subroutine evaluate_meiss


subroutine can_to_ref_meiss(xcan, xref)
    real(dp), intent(in) :: xcan(3)
    real(dp), intent(out) :: xref(3)
    real(dp) :: y_batch(2)  ! lam_phi, chi_gauge

    call evaluate_batch_splines_3d(spl_transform_batch, xcan, y_batch)
    xref(1) = xcan(1)**2
    xref(2) = modulo(xcan(2), twopi)
    xref(3) = modulo(xcan(3) + y_batch(1), twopi)  ! y_batch(1) is lam_phi
end subroutine can_to_ref_meiss


subroutine ref_to_can_meiss(xref, xcan)
    real(dp), intent(in) :: xref(3)
    real(dp), intent(out) :: xcan(3)

    real(dp), parameter :: TOL = 1d-12
    integer, parameter :: MAX_ITER = 16

    real(dp) :: y_batch(2), dy_batch(3, 2), phi_can_prev
    integer :: i

    xcan(1) = sqrt(xref(1))
    xcan(2) = modulo(xref(2), twopi)
    xcan(3) = modulo(xref(3), twopi)

    do i=1, MAX_ITER
        call evaluate_batch_splines_3d_der(spl_transform_batch, xcan, y_batch, dy_batch)
        phi_can_prev = xcan(3)
        ! y_batch(1) is lam_phi, dy_batch(3,1) is d(lam_phi)/d(phi)
        xcan(3) = phi_can_prev - (phi_can_prev + y_batch(1) - xref(3))/(1d0 + dy_batch(3,1))
!print *, abs(xcan(3) - phi_can_prev)
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
    call init_transformation_arrays()
    call compute_transformation()
end subroutine init_transformation

subroutine init_transformation_arrays()
!> Initialize transformation arrays - can be called separately for diagnostics
    if (allocated(lam_phi)) deallocate(lam_phi, chi_gauge)
    allocate(lam_phi(n_r, n_th, n_phi), chi_gauge(n_r, n_th, n_phi))
    
    ! Initialize to zero (boundary conditions)
    lam_phi = 0.0_dp
    chi_gauge = 0.0_dp
end subroutine init_transformation_arrays

subroutine compute_transformation()
!> Compute transformation data via integration (expensive operation)
    integer :: i_ctr

    i_ctr = 0

    !$omp parallel private(i_ctr)
    !$omp do
    do i_ctr = 1, n_phi
        call compute_phi_slice(i_ctr)
    enddo
    !$omp end do
    !$omp end parallel
end subroutine compute_transformation

subroutine compute_phi_slice(i_phi)
!> Compute transformation for single phi slice
    integer, intent(in) :: i_phi
    real(dp) :: y(2)
    integer :: i_th

    call print_progress(i_phi)
    
    do i_th = 1, n_th
        lam_phi(1, i_th, i_phi) = 0d0
        chi_gauge(1, i_th, i_phi) = 0d0
        y = 0d0
        
        if (is_hr_zero_on_slice(i_th, i_phi)) then
            call set_identity_slice(i_th, i_phi)
        else
            call integrate_radial_slice(i_th, i_phi, y)
        endif
    enddo
end subroutine compute_phi_slice

subroutine set_identity_slice(i_th, i_phi)
!> Set identity transformation for radial slice where hr ≈ 0
    integer, intent(in) :: i_th, i_phi
    integer :: i_r
    
    do i_r = 2, n_r
        lam_phi(i_r, i_th, i_phi) = 0d0
        chi_gauge(i_r, i_th, i_phi) = 0d0
    enddo
end subroutine set_identity_slice

subroutine integrate_radial_slice(i_th, i_phi, y)
!> Integrate along radial direction for given (i_th, i_phi)
    integer, intent(in) :: i_th, i_phi
    real(dp), intent(inout) :: y(2)
    integer :: i_r
    
    do i_r = 2, n_r
        call integrate(i_r, i_th, i_phi, y)
    enddo
end subroutine integrate_radial_slice

subroutine print_progress(i_phi)
!> Print integration progress
    integer, intent(in) :: i_phi
    
    !$omp critical
    write(*,'(A, I4, A, I4)',advance='no') 'integrate ODE: ', i_phi, ' of ', n_phi
    if (i_phi < n_phi) then
        write(*, '(A)', advance="no") char(13)
    else
        write(*, *)
    end if
    !$omp end critical
end subroutine print_progress

subroutine integrate(i_r, i_th, i_phi, y)
    use odeint_allroutines_sub, only: odeint_allroutines

    integer, intent(in) :: i_r, i_th, i_phi
    real(dp), dimension(2), intent(inout) :: y

    real(dp), parameter :: relerr=1d-11
    real(dp) :: r1, r2, hr_test, hp_test, phi_c, Ar_test, Ap_test
    real(dp) :: relaxed_relerr
    integer :: ndim=2
    type(grid_indices_t) :: context

    r1 = xmin(1) + h_r*(i_r-2)
    r2 = xmin(1) + h_r*(i_r-1)
    
    ! Check if hr is near zero at the starting point to detect problematic cases
    phi_c = xmin(3) + h_phi*(i_phi-1)
    call ah_cov_on_slice(r1, modulo(phi_c + y(1), twopi), i_th, Ar_test, Ap_test, hr_test, hp_test)
    
    if (abs(hr_test) < hr_small .or. abs(hp_test) < hr_small) then
        lam_phi(i_r, i_th, i_phi) = y(1)
        chi_gauge(i_r, i_th, i_phi) = y(2)
        return
    end if

    relaxed_relerr = relerr
    context = grid_indices_t(i_th, i_phi)
    call odeint_allroutines(y, ndim, context, r1, r2, relaxed_relerr, rh_can_wrapper)

    lam_phi(i_r, i_th, i_phi) = y(1)
    chi_gauge(i_r, i_th, i_phi) = y(2)

end subroutine integrate

subroutine rh_can(r_c, z, dz, i_th, i_phi)
    real(dp), intent(in) :: r_c
    real(dp), dimension(2), intent(in) :: z  ! lam_phi, chi_gauge
    real(dp), dimension(2), intent(out) :: dz
    integer, intent(in) :: i_th, i_phi
    real(dp) :: phi_c
    real(dp) :: Ar, Ap, hr, hp

    phi_c = xmin(3) + h_phi*(i_phi-1)
    call ah_cov_on_slice(r_c, modulo(phi_c + z(1), twopi), i_th, Ar, Ap, hr, hp)

    if (abs(hp) < hr_small) then
        dz(1) = 0.0_dp
        dz(2) = Ar
    else
        dz(1) = -hr/hp
        dz(2) = Ar + Ap*dz(1)
    end if
end subroutine rh_can


subroutine spline_transformation
    real(dp), dimension(:,:,:,:), allocatable :: y_batch
    integer :: dims(3)
    
    ! Construct batch spline for transformation components (2: lam_phi, chi_gauge)
    dims = shape(lam_phi)
    allocate(y_batch(dims(1), dims(2), dims(3), 2))
    
    y_batch(:,:,:,1) = lam_phi
    y_batch(:,:,:,2) = chi_gauge
    
    call construct_batch_splines_3d(xmin, xmax, y_batch, order, periodic, spl_transform_batch)
    transform_splines_initialized = .true.
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
    call field_noncan%evaluate([r, th, phi], Acov, hcov, Bmod)

    Ar = Acov(1)
    Ap = Acov(3)
    hr = hcov(1)
    hp = hcov(3)
end subroutine ah_cov_on_slice

logical function is_hr_zero_on_slice(i_th, i_phi) result(is_zero)
!> Check if hr ≈ 0 along the entire radial direction for given (i_th, i_phi)
!> Sample hr at several radial points to make this determination
    integer, intent(in) :: i_th, i_phi
    
    real(dp), parameter :: hr_slice_threshold = 1d-9
    integer, parameter :: n_sample_points = 32
    real(dp) :: r_sample, phi_c, Ar, Ap, hr, hp
    integer :: i_sample
    logical :: all_zero
    
    phi_c = xmin(3) + h_phi*(i_phi-1)
    all_zero = .true.
    
    ! Sample hr at several radial points
    do i_sample = 1, n_sample_points
        ! Distribute sample points across radial range
        r_sample = xmin(1) + (xmax(1) - xmin(1)) * real(i_sample-1, dp) / real(n_sample_points-1, dp)
        
        ! Evaluate field at this point (using lam_phi=0 for sampling)
        call ah_cov_on_slice(r_sample, phi_c, i_th, Ar, Ap, hr, hp)
        
        if (abs(hr) >= hr_slice_threshold) then
            all_zero = .false.
            exit  ! Found non-zero hr, no need to check further
        endif
    enddo
    
    is_zero = all_zero
end function is_hr_zero_on_slice


subroutine init_canonical_field_components
    real(dp) :: xcan(3)
    real(dp), dimension(:,:,:), allocatable :: Ath, Aphi, hth, hphi, Bmod
    real(dp), dimension(:,:,:,:), allocatable :: y_batch
    real(dp) :: xref(3), Acov(3), hcov(3), lam, dlam(3), chi, dchi(3)
    integer :: i_r, i_th, i_phi, dims(3)

    allocate(Ath(n_r,n_th,n_phi), Aphi(n_r,n_th,n_phi))
    allocate(hth(n_r,n_th,n_phi), hphi(n_r,n_th,n_phi))
    allocate(Bmod(n_r,n_th,n_phi))

    do i_phi=1,n_phi
        do i_th=1,n_th
            do i_r=1,n_r
                xcan = get_grid_point(i_r, i_th, i_phi)

                ! Use batch evaluation for both transformation components
                block
                    real(dp) :: y_trans(2), dy_trans(3,2)
                    call evaluate_batch_splines_3d_der(spl_transform_batch, xcan, y_trans, dy_trans)
                    lam = y_trans(1)    ! lam_phi
                    chi = y_trans(2)    ! chi_gauge  
                    dlam = dy_trans(:,1)  ! derivatives of lam_phi
                    dchi = dy_trans(:,2)  ! derivatives of chi_gauge
                end block

                xref(1) = xcan(1)
                xref(2) = modulo(xcan(2), twopi)
                xref(3) = modulo(xcan(3) + lam, twopi)

                call field_noncan%evaluate(xref, Acov, hcov, Bmod(i_r, i_th, i_phi))

                Ath(i_r, i_th, i_phi) = Acov(2) + Acov(3)*dlam(2) - dchi(2)
                Aphi(i_r, i_th, i_phi) = Acov(3)*(1.0d0 + dlam(3)) - dchi(3)
                hth(i_r, i_th, i_phi) = hcov(2) + hcov(3)*dlam(2)
                hphi(i_r, i_th, i_phi) = hcov(3)*(1.0d0 + dlam(3))
            end do
        end do
    end do

    ! Construct batch spline for all 5 field components: [Ath, Aph, hth, hph, Bmod]
    dims = shape(Ath)
    allocate(y_batch(dims(1), dims(2), dims(3), 5))
    
    y_batch(:,:,:,1) = Ath
    y_batch(:,:,:,2) = Aphi
    y_batch(:,:,:,3) = hth
    y_batch(:,:,:,4) = hphi
    y_batch(:,:,:,5) = Bmod
    
    call construct_batch_splines_3d(xmin, xmax, y_batch, order, periodic, spl_field_batch)
    batch_splines_initialized = .true.
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


! Batch evaluation helper for first derivatives
subroutine evaluate_meiss_batch_der(f, x)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: x(3)
    
    real(dp) :: y_batch(5), dy_batch(3, 5)
    
    call evaluate_batch_splines_3d_der(spl_field_batch, x, y_batch, dy_batch)
    
    ! Unpack results: order is [Ath, Aph, hth, hph, Bmod]
    f%Ath = y_batch(1)
    f%Aph = y_batch(2) 
    f%hth = y_batch(3)
    f%hph = y_batch(4)
    f%Bmod = y_batch(5)
    
    f%dAth = dy_batch(:, 1)
    f%dAph = dy_batch(:, 2)
    f%dhth = dy_batch(:, 3)
    f%dhph = dy_batch(:, 4)
    f%dBmod = dy_batch(:, 5)
end subroutine evaluate_meiss_batch_der


! Batch evaluation helper for second derivatives  
subroutine evaluate_meiss_batch_der2(f, x)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: x(3)
    
    real(dp) :: y_batch(5), dy_batch(3, 5), d2y_batch(6, 5)
    
    call evaluate_batch_splines_3d_der2(spl_field_batch, x, y_batch, dy_batch, d2y_batch)
    
    ! Unpack results: order is [Ath, Aph, hth, hph, Bmod]
    f%Ath = y_batch(1)
    f%Aph = y_batch(2)
    f%hth = y_batch(3) 
    f%hph = y_batch(4)
    f%Bmod = y_batch(5)
    
    f%dAth = dy_batch(:, 1)
    f%dAph = dy_batch(:, 2)
    f%dhth = dy_batch(:, 3)
    f%dhph = dy_batch(:, 4)
    f%dBmod = dy_batch(:, 5)
    
    f%d2Ath = d2y_batch(:, 1)
    f%d2Aph = d2y_batch(:, 2)
    f%d2hth = d2y_batch(:, 3)
    f%d2hph = d2y_batch(:, 4)
    f%d2Bmod = d2y_batch(:, 5)
end subroutine evaluate_meiss_batch_der2

end module field_can_meiss
