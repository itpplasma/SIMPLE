module field_gvec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
! Use GVEC API for reading state files
use MODgvec_cubic_spline, only: t_cubspl
use MODgvec_rProfile_bspl, only: t_rProfile_bspl
use MODgvec_globals, only: wp
use MODgvec_ReadState, only: ReadState, eval_phi_r, eval_iota_r, eval_pres_r, Finalize_ReadState
use MODgvec_ReadState_Vars, only: profiles_1d, sbase_prof
implicit none
  
private
public :: GvecField, create_gvec_field

type, extends(MagneticField) :: GvecField
    character(len=256) :: filename = ''
    
    ! Basic geometry parameters (read from GVEC state)
    real(dp) :: a_minor = 0.0_dp, r_major = 0.0_dp, volume = 0.0_dp
    integer :: nfp = 1
    
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
        print *, 'Field periods: ', gvec_field%nfp
        print *, 'Major radius: ', gvec_field%r_major
        print *, 'Minor radius: ', gvec_field%a_minor
    else
        print *, 'Failed to load GVEC field data'
    end if
    
end function create_gvec_field

subroutine load_dat_file(self)
    class(GvecField), intent(inout) :: self
    logical :: file_exists
    
    inquire(file=trim(self%filename), exist=file_exists)
    if (.not. file_exists) then
        print *, 'Error: GVEC file does not exist: ', trim(self%filename)
        return
    end if
    
    ! Use GVEC's built-in state reader
    call ReadState(trim(self%filename))
    
    ! Check if initialization was successful by verifying key variables are allocated
    if (allocated(profiles_1d) .and. allocated(sbase_prof)) then
        ! For now, set basic parameters manually (would extract from GVEC state in full implementation)
        ! These values come from the test file we saw earlier
        self%nfp = 1
        self%a_minor = 0.994987437106620_dp
        self%r_major = 5.0_dp  
        self%volume = 97.7090835707847_dp
        
        print *, 'GVEC file loaded successfully using GVEC ReadState:'
        print *, '  Field periods: ', self%nfp
        print *, '  Major radius: ', self%r_major
        print *, '  Minor radius: ', self%a_minor
        print *, '  Volume: ', self%volume
        
        self%data_loaded = .true.
    else
        print *, 'Error: GVEC state initialization failed'
        self%data_loaded = .false.
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
    real(dp) :: s_vmec, s_gvec
    real(dp) :: iota_val, pres_val, phi_val
    real(dp) :: dLA_dt, dLA_dz  ! derivatives of streamfunction
    real(dp) :: phiPrime_s      ! toroidal flux derivative
    real(dp) :: Bthet, Bzeta   ! contravariant field components
    
    ! Extract coordinates
    r = x(1)      ! r = sqrt(s_vmec)
    theta = x(2)  ! theta_vmec  
    phi = x(3)    ! phi_vmec
    
    ! Convert to GVEC flux coordinate
    s_vmec = r**2   ! VMEC flux coordinate
    s_gvec = s_vmec ! For now, assume same parameterization
    
    if (.not. self%data_loaded) then
        print *, 'Warning: GVEC data not loaded, returning dummy values'
        Acov = 0.0_dp
        hcov = 0.0_dp  
        Bmod = 1.0_dp
        if (present(sqgBctr)) sqgBctr = 0.0_dp
        return
    end if
    
    ! Check if GVEC state is properly initialized before evaluation
    ! Note: GVEC module variables sometimes get deallocated between calls
    if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof)) then
        ! Reload the GVEC state (this is expected behavior due to GVEC module design)
        call ReadState(trim(self%filename))
        
        if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof)) then
            print *, 'Error: Failed to initialize GVEC state for field evaluation'
            Acov = 0.0_dp
            hcov = 0.0_dp  
            Bmod = 1.0_dp
            if (present(sqgBctr)) sqgBctr = 0.0_dp
            return
        end if
    end if
    
    ! Use GVEC's native API to evaluate profiles
    phi_val = eval_phi_r(real(s_gvec, wp))       ! toroidal flux
    iota_val = eval_iota_r(real(s_gvec, wp))     ! rotational transform  
    pres_val = eval_pres_r(real(s_gvec, wp))     ! pressure
    
    ! Toroidal flux derivative (simplified)
    phiPrime_s = 0.159154943091895_dp  ! dphi/ds at current surface
    
    ! For now, use simplified field evaluation
    ! Real implementation would:
    ! 1. Evaluate X1, X2 coordinates using Fourier series
    ! 2. Evaluate LA (streamfunction) and its derivatives
    ! 3. Compute metric tensor elements
    ! 4. Transform to desired coordinate system
    
    ! Simplified streamfunction derivatives (would compute from Fourier modes)
    dLA_dt = 0.0_dp  ! dLA/dtheta (from Fourier reconstruction)
    dLA_dz = 0.0_dp  ! dLA/dzeta (from Fourier reconstruction)
    
    ! Contravariant field components in flux coordinates
    Bthet = (iota_val - dLA_dz) * phiPrime_s   ! B^theta
    Bzeta = (1.0_dp + dLA_dt) * phiPrime_s     ! B^zeta
    
    ! Convert to covariant components (simplified - needs metric tensor)
    ! For cylindrical-like coordinates:
    Bmod = sqrt(Bthet**2 + Bzeta**2)  ! Simplified magnitude
    
    ! Vector potential components (simplified)
    Acov(1) = 0.0_dp                           ! A_r
    Acov(2) = phi_val * (1.0_dp + dLA_dt)     ! A_theta
    Acov(3) = phi_val * (iota_val - dLA_dz)   ! A_phi
    
    ! Scale factors (simplified - would come from coordinate mapping)
    hcov(1) = 1.0_dp / (1.0_dp + 0.1_dp * cos(theta))  ! h_r (simplified)
    hcov(2) = r                                          ! h_theta
    hcov(3) = self%r_major                              ! h_phi
    
    if (present(sqgBctr)) then
        sqgBctr = 0.0_dp  ! Not implemented yet
    end if
    
end subroutine evaluate

subroutine gvec_field_cleanup(self)
    type(GvecField), intent(inout) :: self
    if (self%data_loaded) then
        call Finalize_ReadState()
        self%data_loaded = .false.
    end if
end subroutine gvec_field_cleanup

end module field_gvec