module field_can_albert_batch

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
    evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
    evaluate_batch_splines_3d_der2, destroy_batch_splines_3d, &
    SplineData3D, evaluate_splines_3d
use field_can_base, only: FieldCan, n_field_evaluations
use psi_transform, only: grid_r_to_psi
use batch_spline_helpers, only: evaluate_batch_splines_3d_single

implicit none

! Grid parameters from field_can_meiss
integer :: n_r, n_th, n_phi
real(dp) :: xmin(3), xmax(3)
integer, parameter :: order(3) = [5, 5, 5]
logical, parameter :: periodic(3) = [.False., .True., .True.]
real(dp), parameter :: twopi = 6.283185307179586d0

! For splining psi
real(dp) :: psi_inner, psi_outer
real(dp), dimension(:,:,:), allocatable :: psi_of_x
real(dp), dimension(:), allocatable :: psi_grid
logical :: dpsi_dr_positive

! Batched spline for Albert canonical field components
! Components: [1:r_of_xc, 2:Aphi_of_xc, 3:hth_of_xc, 4:hph_of_xc, 5:Bmod_of_xc]
type(BatchSplineData3D) :: spl_albert_batch

! Ath normalization constant
real(dp) :: Ath_norm

! Constants for batch indices
integer, parameter :: IDX_R = 1
integer, parameter :: IDX_APHI = 2
integer, parameter :: IDX_HTH = 3
integer, parameter :: IDX_HPH = 4
integer, parameter :: IDX_BMOD = 5

contains

subroutine init_albert_batch(n_r_, n_th_, n_phi_, xmin_, xmax_)
    integer, intent(in) :: n_r_, n_th_, n_phi_
    real(dp), intent(in) :: xmin_(3), xmax_(3)
    
    n_r = n_r_
    n_th = n_th_
    n_phi = n_phi_
    xmin = xmin_
    xmax = xmax_
end subroutine init_albert_batch


subroutine evaluate_albert_batch(f, r, th_c, ph_c, mode_secders)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    real(dp) :: x(3)
    real(dp) :: y_batch(5), dy_batch(5, 3), d2y_batch(5, 6)

    n_field_evaluations = n_field_evaluations + 1

    x = [r, th_c, ph_c]

    ! Ath is linear in r for Albert coordinates
    f%Ath = Ath_norm*x(1)
    f%dAth = [Ath_norm, 0d0, 0d0]

    if (mode_secders > 0) then
        f%d2Ath = 0d0
        
        ! Evaluate batch splines for other components with second derivatives
        ! Note: We don't need r_of_xc here, only components 2-5
        if (spl_albert_batch%num_quantities /= 5) then
            error stop 'Albert batch spline must have exactly 5 quantities'
        end if
        call evaluate_batch_splines_3d_der2(spl_albert_batch, x, &
            y_batch, dy_batch, d2y_batch)
        
        f%Aph = y_batch(IDX_APHI)
        f%dAph = dy_batch(IDX_APHI, :)
        f%d2Aph = d2y_batch(IDX_APHI, :)
        
        f%hth = y_batch(IDX_HTH)
        f%dhth = dy_batch(IDX_HTH, :)
        f%d2hth = d2y_batch(IDX_HTH, :)
        
        f%hph = y_batch(IDX_HPH)
        f%dhph = dy_batch(IDX_HPH, :)
        f%d2hph = d2y_batch(IDX_HPH, :)
        
        f%Bmod = y_batch(IDX_BMOD)
        f%dBmod = dy_batch(IDX_BMOD, :)
        f%d2Bmod = d2y_batch(IDX_BMOD, :)
        
        return
    end if

    ! Evaluate batch splines with first derivatives only
    if (spl_albert_batch%num_quantities /= 5) then
        error stop 'Albert batch spline must have exactly 5 quantities'
    end if
    call evaluate_batch_splines_3d_der(spl_albert_batch, x, y_batch, dy_batch)
    
    f%Aph = y_batch(IDX_APHI)
    f%dAph = dy_batch(IDX_APHI, :)
    
    f%hth = y_batch(IDX_HTH)
    f%dhth = dy_batch(IDX_HTH, :)
    
    f%hph = y_batch(IDX_HPH)
    f%dhph = dy_batch(IDX_HPH, :)
    
    f%Bmod = y_batch(IDX_BMOD)
    f%dBmod = dy_batch(IDX_BMOD, :)
end subroutine evaluate_albert_batch


subroutine can_to_ref_albert_batch(xcan, xref)
    use field_can_meiss, only: can_to_ref_meiss
    
    real(dp), intent(in) :: xcan(3)
    real(dp), intent(out) :: xref(3) 
    real(dp) :: xmeiss(3)
    
    ! Get r_of_xc from batch spline (only first component)
    call evaluate_batch_splines_3d_single(spl_albert_batch, xcan, IDX_R, xmeiss(1))
    xmeiss(2:3) = xcan(2:3)
    call can_to_ref_meiss(xmeiss, xref)
end subroutine can_to_ref_albert_batch
    

subroutine ref_to_can_albert_batch(xref, xcan)
    use field_can_meiss, only: ref_to_can_meiss, spl_field_batch
    use interpolate, only: SplineData3D, evaluate_splines_3d
    
    real(dp), intent(in) :: xref(3)
    real(dp), intent(out) :: xcan(3)
 
    real(dp) :: Ath, xmeiss(3), y_batch(5)

    call ref_to_can_meiss(xref, xmeiss)   
    ! Extract Ath from batch spline (component 1)
    call evaluate_batch_splines_3d(spl_field_batch, xmeiss, y_batch)
    Ath = y_batch(1)
    xcan(1) = Ath/Ath_norm
    xcan(2:3) = xmeiss(2:3)
end subroutine ref_to_can_albert_batch


subroutine init_splines_with_psi_batch(Aphi_meiss, hth_meiss, hph_meiss, Bmod_meiss)
    real(dp), dimension(:,:,:), intent(in) :: Aphi_meiss, hth_meiss, hph_meiss, Bmod_meiss
    real(dp), dimension(:,:,:), allocatable :: r_of_xc, &
        Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc
    real(dp), dimension(:,:,:,:), allocatable :: albert_batch
    real(dp) :: x(3)
    integer :: i_r, i_th, i_phi

    allocate( &
        r_of_xc(n_r, n_th, n_phi), &
        Aph_of_xc(n_r, n_th, n_phi), &
        hth_of_xc(n_r, n_th, n_phi), &
        hph_of_xc(n_r, n_th, n_phi), &
        Bmod_of_xc(n_r, n_th, n_phi) &
    )

    call init_psi_grid_batch

    call grid_r_to_psi(xmin(1), xmax(1), psi_inner, psi_outer, psi_of_x, &
        Aphi_meiss, hth_meiss, hph_meiss, Bmod_meiss, &
        r_of_xc, Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc)

    ! Create batched data array
    allocate(albert_batch(n_r, n_th, n_phi, 5))
    albert_batch(:,:,:, IDX_R) = r_of_xc
    albert_batch(:,:,:, IDX_APHI) = Aph_of_xc
    albert_batch(:,:,:, IDX_HTH) = hth_of_xc
    albert_batch(:,:,:, IDX_HPH) = hph_of_xc
    albert_batch(:,:,:, IDX_BMOD) = Bmod_of_xc

    ! Construct batch spline for Albert field components
    call construct_batch_splines_3d([psi_inner, xmin(2), xmin(3)], &
        [psi_outer, xmax(2), xmax(3)], albert_batch, order, periodic, &
        spl_albert_batch)
    
    deallocate(albert_batch)
    deallocate(r_of_xc, Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc)
end subroutine init_splines_with_psi_batch


subroutine init_psi_grid_batch
    use field_can_meiss, only: spl_field_batch
    
    integer :: i_r, i_th, i_phi
    real(dp) :: x(3), y_batch_local(5)
    
    allocate(psi_of_x(n_r, n_th, n_phi), psi_grid(n_r))
    
    ! Get psi values from the existing Meiss Ath spline
    do i_r = 1, n_r
        do i_th = 1, n_th
            do i_phi = 1, n_phi
                x = [xmin(1) + (i_r-1)*(xmax(1)-xmin(1))/(n_r-1), &
                     xmin(2) + (i_th-1)*(xmax(2)-xmin(2))/(n_th-1), &
                     xmin(3) + (i_phi-1)*(xmax(3)-xmin(3))/(n_phi-1)]
                ! Extract Ath from batch spline (component 1) 
                call evaluate_batch_splines_3d(spl_field_batch, x, y_batch_local)
                psi_of_x(i_r, i_th, i_phi) = y_batch_local(1)
            enddo
        enddo
    enddo
    
    ! Normalize psi_of_x
    Ath_norm = sign(maxval(abs(psi_of_x)), psi_of_x(n_r, n_th/2, n_phi/2))
    psi_of_x = psi_of_x / Ath_norm
    
    ! Determine inner and outer psi bounds
    if(psi_of_x(n_r, n_th/2, n_phi/2) > psi_of_x(1, n_th/2, n_phi/2)) then
        dpsi_dr_positive = .true.
        psi_inner = maxval(psi_of_x(1,:,:))
        psi_outer = minval(psi_of_x(n_r,:,:))
    else
        dpsi_dr_positive = .false.
        psi_inner = maxval(psi_of_x(n_r,:,:))
        psi_outer = minval(psi_of_x(1,:,:))
    endif
    
    ! Create uniform psi grid
    do i_r = 1, n_r
        psi_grid(i_r) = psi_inner + (psi_outer - psi_inner) * (i_r - 1) / (n_r - 1)
    end do
end subroutine init_psi_grid_batch


subroutine cleanup_albert_batch
    ! Clean up batch spline
    call destroy_batch_splines_3d(spl_albert_batch)
    
    if (allocated(psi_of_x)) deallocate(psi_of_x)
    if (allocated(psi_grid)) deallocate(psi_grid)
end subroutine cleanup_albert_batch



end module field_can_albert_batch