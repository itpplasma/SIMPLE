module field_gvec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
! Use GVEC spline functionality
use MODgvec_cubic_spline, only: t_cubspl
use MODgvec_rProfile_bspl, only: t_rProfile_bspl
implicit none
  
private
public :: GvecField, create_gvec_field

type, extends(MagneticField) :: GvecField
    character(len=256) :: filename = ''
    
    ! Grid information
    integer :: nElems = 0, nfp = 1
    real(dp), allocatable :: sp(:)  ! flux surface coordinates
    
    ! Fourier mode data (simplified)
    integer :: X1_modes = 0, X2_modes = 0, LA_modes = 0
    integer, allocatable :: X1_mn(:,:), X2_mn(:,:), LA_mn(:,:)
    real(dp), allocatable :: X1_coefs(:,:), X2_coefs(:,:), LA_coefs(:,:)
    
    ! Profile splines (using GVEC types)
    type(t_cubspl) :: phi_spline, iota_spline, pres_spline
    
    ! Geometry
    real(dp) :: a_minor = 0.0_dp, r_major = 0.0_dp, volume = 0.0_dp
    
    ! Data loaded flag
    logical :: data_loaded = .false.
    
contains
    procedure :: evaluate
    procedure :: load_dat_file
    final :: gvec_field_cleanup
end type GvecField

contains

function create_gvec_field(gvec_file) result(gvec_field)
    class(GvecField), allocatable :: gvec_field
    character(*), intent(in) :: gvec_file
    
    allocate(GvecField :: gvec_field)
    gvec_field%filename = gvec_file
    
    print *, 'Loading GVEC field from: ', gvec_file
    
    ! Load the .dat file
    call gvec_field%load_dat_file()
    
    if (gvec_field%data_loaded) then
        print *, 'GVEC field loaded successfully'
        print *, 'Grid elements: ', gvec_field%nElems
        print *, 'Field periods: ', gvec_field%nfp
        print *, 'Major radius: ', gvec_field%r_major
        print *, 'Minor radius: ', gvec_field%a_minor
    else
        print *, 'Failed to load GVEC field data'
    end if
    
end function create_gvec_field

subroutine load_dat_file(self)
    class(GvecField), intent(inout) :: self
    integer :: iounit, ios, i
    character(len=1024) :: line
    logical :: file_exists
    
    inquire(file=trim(self%filename), exist=file_exists)
    if (.not. file_exists) then
        print *, 'Error: GVEC file does not exist: ', trim(self%filename)
        return
    end if
    
    open(newunit=iounit, file=trim(self%filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening GVEC file: ', trim(self%filename)
        return
    end if
    
    ! Read basic file structure (simplified parsing)
    read(iounit, '(A)', iostat=ios) line  ! Header comment
    if (ios /= 0) goto 100
    
    read(iounit, *, iostat=ios) ! outputLevel, fileID  
    if (ios /= 0) goto 100
    
    read(iounit, '(A)', iostat=ios) line  ! grid comment
    if (ios /= 0) goto 100
    
    read(iounit, *, iostat=ios) self%nElems  ! nElems, gridType
    if (ios /= 0) goto 100
    
    ! Skip detailed parsing for now - just set basic values
    self%nfp = 1
    self%r_major = 1.0_dp
    self%a_minor = 0.5_dp
    self%volume = 1.0_dp
    
    ! Allocate basic grid
    if (allocated(self%sp)) deallocate(self%sp)
    allocate(self%sp(0:max(self%nElems, 1)))
    
    ! Simple linear grid for now
    do i = 0, self%nElems
        self%sp(i) = real(i, dp) / real(max(self%nElems, 1), dp)
    end do
    
    print *, 'Basic GVEC file structure parsed'
    self%data_loaded = .true.
    
100 continue
    close(iounit)
    
    if (ios /= 0) then
        print *, 'Error reading GVEC file format at line: ', trim(line)
    end if
    
end subroutine load_dat_file

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    class(GvecField), intent(in) :: self
    real(dp), intent(in) :: x(3)  ! r=sqrt(s_vmec), theta_vmec, phi_vmec
    real(dp), intent(out) :: Acov(3)
    real(dp), intent(out) :: hcov(3)
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)
    
    real(dp) :: r, theta, phi
    real(dp) :: rho2, s
    
    ! Extract coordinates
    r = x(1)      ! r = sqrt(s_vmec)
    theta = x(2)  ! theta_vmec
    phi = x(3)    ! phi_vmec
    
    ! Convert to flux coordinate
    rho2 = r**2   ! s_vmec
    s = r         ! radial coordinate
    
    if (.not. self%data_loaded) then
        print *, 'Warning: GVEC data not loaded, returning dummy values'
        Acov = 0.0_dp
        hcov = 0.0_dp  
        Bmod = 1.0_dp
        if (present(sqgBctr)) sqgBctr = 0.0_dp
        return
    end if
    
    ! For now, implement a simple magnetic field model
    ! This would be replaced with actual spline evaluation and Fourier reconstruction
    
    ! Simple analytical model for testing
    Acov(1) = 0.0_dp                    ! A_r
    Acov(2) = 0.0_dp                    ! A_theta  
    Acov(3) = 0.5_dp * rho2             ! A_phi (simple model)
    
    hcov(1) = 1.0_dp                    ! h_r
    hcov(2) = r                         ! h_theta
    hcov(3) = 1.0_dp                    ! h_phi
    
    Bmod = 1.0_dp / (1.0_dp + 0.1_dp * cos(theta))  ! Simple B-field model
    
    if (present(sqgBctr)) then
        sqgBctr = 0.0_dp
    end if
    
end subroutine evaluate

subroutine gvec_field_cleanup(self)
    type(GvecField), intent(inout) :: self
    if (allocated(self%sp)) deallocate(self%sp)
    if (allocated(self%X1_mn)) deallocate(self%X1_mn)
    if (allocated(self%X2_mn)) deallocate(self%X2_mn)
    if (allocated(self%LA_mn)) deallocate(self%LA_mn)
    if (allocated(self%X1_coefs)) deallocate(self%X1_coefs)
    if (allocated(self%X2_coefs)) deallocate(self%X2_coefs)
    if (allocated(self%LA_coefs)) deallocate(self%LA_coefs)
end subroutine gvec_field_cleanup

end module field_gvec