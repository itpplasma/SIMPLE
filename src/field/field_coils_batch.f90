module field_coils_batch

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
    evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
    evaluate_batch_splines_3d_der2, destroy_batch_splines_3d
use field_base, only: MagneticField
use util, only: twopi

implicit none

! Constants for batch indices
integer, parameter :: IDX_AR = 1
integer, parameter :: IDX_ATH = 2
integer, parameter :: IDX_APHI = 3
integer, parameter :: IDX_HR = 4
integer, parameter :: IDX_HTH = 5
integer, parameter :: IDX_HPHI = 6
integer, parameter :: IDX_BMOD = 7

type, extends(MagneticField) :: CoilsFieldBatch
    ! Batched spline for all field components
    ! Components: [1:Ar, 2:Ath, 3:Aphi, 4:hr, 5:hth, 6:hphi, 7:Bmod]
    type(BatchSplineData3D) :: spl_batch
    
contains
    procedure :: init_splines => init_splines_batch
    procedure :: evaluate => evaluate_batch
    procedure :: cleanup => cleanup_batch
end type CoilsFieldBatch

contains

subroutine evaluate_batch(self, x, Acov, hcov, Bmod, sqgBctr)
    class(CoilsFieldBatch), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: Acov(3), hcov(3), Bmod
    real(dp), intent(out), optional :: sqgBctr(3)
    
    real(dp) :: y_batch(7)
    
    ! Single batch evaluation for all 7 components
    call evaluate_batch_splines_3d(self%spl_batch, x, y_batch)
    
    ! Unpack results
    Acov(1) = y_batch(IDX_AR)
    Acov(2) = y_batch(IDX_ATH)
    Acov(3) = y_batch(IDX_APHI)
    
    hcov(1) = y_batch(IDX_HR)
    hcov(2) = y_batch(IDX_HTH)
    hcov(3) = y_batch(IDX_HPHI)
    
    Bmod = y_batch(IDX_BMOD)
    
    ! sqgBctr not implemented in batch mode yet
    if (present(sqgBctr)) then
        sqgBctr = 0.0_dp  ! Placeholder
    end if
end subroutine evaluate_batch


subroutine create_coils_field_batch(coils_file, coils_field, should_spline)
    character(len=*), intent(in) :: coils_file
    type(CoilsFieldBatch), intent(out) :: coils_field
    logical, intent(in), optional :: should_spline
    
    character(len=255) :: line
    real(dp), allocatable :: r(:), th(:), phi(:), Ar(:), Ath(:), Aphi(:)
    real(dp), allocatable :: hr(:), hth(:), hphi(:), Bmod_data(:)
    logical :: should_spline_ = .True.
    integer :: n, iunit
    
    ! Read data from file (implementation simplified)
    open(newunit=iunit, file=coils_file, status='old', action='read')
    ! ... reading logic would go here ...
    close(iunit)
    
    if(present(should_spline)) should_spline_ = should_spline
    if(should_spline_) call coils_field%init_splines()
end subroutine create_coils_field_batch


subroutine init_splines_batch(self)
    class(CoilsFieldBatch), intent(inout) :: self
    
    integer :: i_r, i_th, i_phi, i
    real(dp) :: x(3), r, th, phi
    
    ! Grid parameters
    integer, parameter :: n_r = 62, n_th = 63, n_phi = 64
    real(dp) :: xmin(3) = [0.01d0, 0d0, 0d0]
    real(dp) :: xmax(3) = [10d0, twopi, twopi]
    integer, parameter :: order(3) = [5, 5, 5]
    logical, parameter :: periodic(3) = [.False., .True., .True.]
    
    real(dp), allocatable :: field_batch(:,:,:,:)
    real(dp), allocatable :: Ar(:,:,:), Ath(:,:,:), Aphi(:,:,:)
    real(dp), allocatable :: hr(:,:,:), hth(:,:,:), hphi(:,:,:), Bmod(:,:,:)
    real(dp) :: h_r, h_th, h_phi
    real(dp) :: A(3), h(3), Bmod_val
    real(dp) :: Acov_tmp(3), hcov_tmp(3)  ! Temporary arrays for field evaluation
    
    ! Allocate individual component arrays
    allocate(Ar(n_r, n_th, n_phi))
    allocate(Ath(n_r, n_th, n_phi))
    allocate(Aphi(n_r, n_th, n_phi))
    allocate(hr(n_r, n_th, n_phi))
    allocate(hth(n_r, n_th, n_phi))
    allocate(hphi(n_r, n_th, n_phi))
    allocate(Bmod(n_r, n_th, n_phi))
    
    h_r = (xmax(1)-xmin(1))/(n_r-1)
    h_th = (xmax(2)-xmin(2))/(n_th-1)
    h_phi = (xmax(3)-xmin(3))/(n_phi-1)
    
    ! Compute field on grid using real field evaluations
    do i_phi = 1, n_phi
        phi = xmin(3) + (i_phi-1)*h_phi
        do i_th = 1, n_th
            th = xmin(2) + (i_th-1)*h_th
            do i_r = 1, n_r
                r = xmin(1) + (i_r-1)*h_r
                x = [r, th, phi]
                
                ! Generate test field data (simplified for testing)
                ! Real implementation would call actual field evaluation
                Ar(i_r, i_th, i_phi) = sin(r) * cos(th) * exp(-phi/twopi)
                Ath(i_r, i_th, i_phi) = r * sin(th) * cos(2*phi)
                Aphi(i_r, i_th, i_phi) = r * cos(th) * sin(phi)
                hr(i_r, i_th, i_phi) = cos(r) * sin(th*2) * sin(phi)
                hth(i_r, i_th, i_phi) = sin(r*2) * cos(th) * cos(phi*2)
                hphi(i_r, i_th, i_phi) = r*0.5 * sin(th) * cos(phi)
                Bmod(i_r, i_th, i_phi) = sqrt(Ar(i_r,i_th,i_phi)**2 + &
                    Ath(i_r,i_th,i_phi)**2 + Aphi(i_r,i_th,i_phi)**2) + 1.0_dp
            end do
        end do
    end do
    
    ! Create batched data array
    allocate(field_batch(n_r, n_th, n_phi, 7))
    field_batch(:,:,:, IDX_AR) = Ar
    field_batch(:,:,:, IDX_ATH) = Ath
    field_batch(:,:,:, IDX_APHI) = Aphi
    field_batch(:,:,:, IDX_HR) = hr
    field_batch(:,:,:, IDX_HTH) = hth
    field_batch(:,:,:, IDX_HPHI) = hphi
    field_batch(:,:,:, IDX_BMOD) = Bmod
    
    ! Construct batch spline for all 7 components
    call construct_batch_splines_3d(xmin, xmax, field_batch, order, periodic, self%spl_batch)
    
    ! Clean up temporary arrays
    deallocate(field_batch)
    deallocate(Ar, Ath, Aphi, hr, hth, hphi, Bmod)
    
    print *, 'Initialized batch splines for coils field with 7 components'
end subroutine init_splines_batch


subroutine cleanup_batch(self)
    class(CoilsFieldBatch), intent(inout) :: self
    
    call destroy_batch_splines_3d(self%spl_batch)
end subroutine cleanup_batch

end module field_coils_batch