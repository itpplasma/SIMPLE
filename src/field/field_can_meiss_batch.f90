module field_can_meiss_batch

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
    evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
    evaluate_batch_splines_3d_der2, destroy_batch_splines_3d
use util, only: twopi
use field_can_base, only: FieldCan, n_field_evaluations
use field, only: MagneticField
use batch_spline_helpers, only: evaluate_batch_splines_3d_single

implicit none

type :: grid_indices_t
    integer :: i_th, i_phi
end type grid_indices_t

class(MagneticField), allocatable :: field_noncan
integer :: n_r=62, n_th=63, n_phi=64
real(dp) :: xmin(3) = [1d-6, 0d0, 0d0]
real(dp) :: xmax(3) = [1d0, twopi, twopi]

real(dp) :: h_r, h_th, h_phi

! Batched splines for field components that are evaluated together
! Components: [1:Ath, 2:Aph, 3:hth, 4:hph, 5:Bmod]
type(BatchSplineData3D) :: spl_field_batch

! Transformation splines still separate as they're used independently
! Components: [1:lam_phi, 2:chi_gauge]
type(BatchSplineData3D) :: spl_transform_batch

integer, parameter :: order(3) = [5, 5, 5]
logical, parameter :: periodic(3) = [.False., .True., .True.]

! Constants for batch indices
integer, parameter :: IDX_ATH = 1
integer, parameter :: IDX_APH = 2
integer, parameter :: IDX_HTH = 3
integer, parameter :: IDX_HPH = 4
integer, parameter :: IDX_BMOD = 5

integer, parameter :: IDX_LAM_PHI = 1
integer, parameter :: IDX_CHI_GAUGE = 2

! Forward declarations for external procedures
interface
    subroutine rh_can(r_c, z, dz, i_th, i_phi)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        real(dp), intent(in) :: r_c
        real(dp), dimension(:), intent(in) :: z
        real(dp), dimension(:), intent(out) :: dz
        integer, intent(in) :: i_th, i_phi
    end subroutine rh_can
    
    subroutine get_grid_point(f, xcan)
        import :: FieldCan, dp
        type(FieldCan), intent(inout) :: f
        real(dp), intent(in) :: xcan(3)
    end subroutine get_grid_point
end interface

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

subroutine init_meiss_batch(field_noncan_, n_r_, n_th_, n_phi_, rmin, rmax, thmin, thmax)
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
end subroutine init_meiss_batch


subroutine evaluate_meiss_batch(f, r, th_c, ph_c, mode_secders)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    real(dp) :: x(3)
    real(dp) :: y_batch(5), dy_batch(3, 5), d2y_batch(6, 5)

    n_field_evaluations = n_field_evaluations + 1

    x = [r, th_c, ph_c]

    if (mode_secders > 0) then
        ! Evaluate all field components with second derivatives
        call evaluate_batch_splines_3d_der2(spl_field_batch, x, &
            y_batch, dy_batch, d2y_batch)
        
        ! Unpack results
        f%Ath = y_batch(IDX_ATH)
        f%dAth = dy_batch(:, IDX_ATH)
        f%d2Ath = d2y_batch(:, IDX_ATH)
        
        f%Aph = y_batch(IDX_APH)
        f%dAph = dy_batch(:, IDX_APH)
        f%d2Aph = d2y_batch(:, IDX_APH)
        
        f%hth = y_batch(IDX_HTH)
        f%dhth = dy_batch(:, IDX_HTH)
        f%d2hth = d2y_batch(:, IDX_HTH)
        
        f%hph = y_batch(IDX_HPH)
        f%dhph = dy_batch(:, IDX_HPH)
        f%d2hph = d2y_batch(:, IDX_HPH)
        
        f%Bmod = y_batch(IDX_BMOD)
        f%dBmod = dy_batch(:, IDX_BMOD)
        f%d2Bmod = d2y_batch(:, IDX_BMOD)
        
        return
    end if

    ! Evaluate all field components with first derivatives only
    call evaluate_batch_splines_3d_der(spl_field_batch, x, y_batch, dy_batch)
    
    ! Unpack results
    f%Ath = y_batch(IDX_ATH)
    f%dAth = dy_batch(:, IDX_ATH)
    
    f%Aph = y_batch(IDX_APH)
    f%dAph = dy_batch(:, IDX_APH)
    
    f%hth = y_batch(IDX_HTH)
    f%dhth = dy_batch(:, IDX_HTH)
    
    f%hph = y_batch(IDX_HPH)
    f%dhph = dy_batch(:, IDX_HPH)
    
    f%Bmod = y_batch(IDX_BMOD)
    f%dBmod = dy_batch(:, IDX_BMOD)
end subroutine evaluate_meiss_batch


subroutine can_to_ref_meiss_batch(xcan, xref)
    real(dp), intent(in) :: xcan(3)
    real(dp), intent(out) :: xref(3)
    real(dp) :: y_batch(2)

    ! Evaluate only lam_phi (first component)
    call evaluate_batch_splines_3d_single(spl_transform_batch, xcan, &
        IDX_LAM_PHI, y_batch(IDX_LAM_PHI))
    
    xref(1) = xcan(1)**2
    xref(2) = modulo(xcan(2), twopi)
    xref(3) = modulo(xcan(3) + y_batch(IDX_LAM_PHI), twopi)
end subroutine can_to_ref_meiss_batch


subroutine ref_to_can_meiss_batch(xref, xcan)
    real(dp), intent(in) :: xref(3)
    real(dp), intent(out) :: xcan(3)

    real(dp), parameter :: TOL = 1d-12
    integer, parameter :: MAX_ITER = 16

    real(dp) :: y_batch(2), dy_batch(3, 2)
    real(dp) :: lam, dlam(3), phi_can_prev
    integer :: i

    xcan(1) = sqrt(xref(1))
    xcan(2) = modulo(xref(2), twopi)
    xcan(3) = modulo(xref(3), twopi)

    do i=1, MAX_ITER
        ! Evaluate lam_phi and its derivatives
        call evaluate_batch_splines_3d_der_single(spl_transform_batch, xcan, &
            IDX_LAM_PHI, y_batch(IDX_LAM_PHI), dy_batch(:, IDX_LAM_PHI))
        
        lam = y_batch(IDX_LAM_PHI)
        dlam = dy_batch(:, IDX_LAM_PHI)
        
        phi_can_prev = xcan(3)
        xcan(3) = phi_can_prev - (phi_can_prev + lam - xref(3))/(1d0 + dlam(3))
        
        if (abs(xcan(3) - phi_can_prev) < TOL) return
    enddo
    print *, 'WARNING: ref_to_can_meiss_batch did not converge after', MAX_ITER, 'iterations'
end subroutine ref_to_can_meiss_batch


subroutine spline_transformation_batch(lam_phi, chi_gauge)
    real(dp), dimension(:,:,:), intent(in) :: lam_phi, chi_gauge
    real(dp), dimension(:,:,:,:), allocatable :: transform_batch
    
    ! Create batched data array
    allocate(transform_batch(n_r, n_th, n_phi, 2))
    transform_batch(:,:,:, IDX_LAM_PHI) = lam_phi
    transform_batch(:,:,:, IDX_CHI_GAUGE) = chi_gauge
    
    ! Construct batch spline for transformation quantities
    call construct_batch_splines_3d(xmin, xmax, transform_batch, order, &
        periodic, spl_transform_batch)
    
    deallocate(transform_batch)
end subroutine spline_transformation_batch


subroutine init_canonical_field_components_batch(Ath, Aphi, hth, hphi, Bmod)
    real(dp), dimension(:,:,:), intent(in) :: Ath, Aphi, hth, hphi, Bmod
    real(dp), dimension(:,:,:,:), allocatable :: field_batch
    
    ! Create batched data array
    allocate(field_batch(n_r, n_th, n_phi, 5))
    field_batch(:,:,:, IDX_ATH) = Ath
    field_batch(:,:,:, IDX_APH) = Aphi
    field_batch(:,:,:, IDX_HTH) = hth
    field_batch(:,:,:, IDX_HPH) = hphi
    field_batch(:,:,:, IDX_BMOD) = Bmod
    
    ! Construct batch spline for field components
    call construct_batch_splines_3d(xmin, xmax, field_batch, order, &
        periodic, spl_field_batch)
    
    deallocate(field_batch)
end subroutine init_canonical_field_components_batch


subroutine cleanup_meiss_batch
    ! Clean up batch splines
    call destroy_batch_splines_3d(spl_field_batch)
    call destroy_batch_splines_3d(spl_transform_batch)
    
    if (allocated(field_noncan)) deallocate(field_noncan)
end subroutine cleanup_meiss_batch


! Wrapper routines for single quantity evaluation from batch splines


subroutine evaluate_batch_splines_3d_der_single(spl, x, iq, y, dy)
    type(BatchSplineData3D), intent(in) :: spl
    real(dp), intent(in) :: x(3)
    integer, intent(in) :: iq
    real(dp), intent(out) :: y
    real(dp), intent(out) :: dy(3)
    
    real(dp) :: y_batch(spl%num_quantities)
    real(dp) :: dy_batch(3, spl%num_quantities)
    
    call evaluate_batch_splines_3d_der(spl, x, y_batch, dy_batch)
    y = y_batch(iq)
    dy = dy_batch(:, iq)
end subroutine evaluate_batch_splines_3d_der_single

end module field_can_meiss_batch