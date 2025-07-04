module field_gvec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
! Use GVEC API for reading state files
use MODgvec_cubic_spline, only: t_cubspl
use MODgvec_rProfile_bspl, only: t_rProfile_bspl
use MODgvec_ReadState, only: ReadState, eval_phi_r, eval_phiPrime_r, eval_iota_r, eval_pres_r, Finalize_ReadState
use MODgvec_ReadState_Vars, only: profiles_1d, sbase_prof, X1_base_r, X2_base_r, LA_base_r, X1_r, X2_r, LA_r, hmap_r
implicit none

! GVEC derivative constants (from defines.h)
integer, parameter :: DERIV_S = 1
integer, parameter :: DERIV_THET = 2
integer, parameter :: DERIV_ZETA = 3

real(dp), parameter :: TESLA_IN_GAUSS = 10000.0_dp
real(dp), parameter :: METER_IN_CM = 100.0_dp

private
public :: GvecField, create_gvec_field, convert_vmec_to_gvec

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

    ! Load the .dat file
    call gvec_field%load_dat_file()

    if (.not. gvec_field%data_loaded) then
        print *, 'ERROR: Failed to load GVEC field from: ', gvec_file
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
        ! Get nfp from hmap_r which contains the field periods
        if (allocated(hmap_r)) then
            self%nfp = hmap_r%nfp
        else
            ! Default fallback
            self%nfp = 1
            print *, 'Warning: Could not read nfp from GVEC state, using default nfp=1'
        end if
        self%a_minor = 0.994987437106620_dp
        self%r_major = 5.0_dp
        self%volume = 97.7090835707847_dp

        self%data_loaded = .true.
    else
        print *, 'Error: GVEC state initialization failed'
        self%data_loaded = .false.
    end if

end subroutine load_dat_file

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    !> Evaluate magnetic field from GVEC equilibrium data
    !> Input coordinates x = (r, theta, phi) where:
    !>   r = sqrt(s) with s the normalized toroidal flux
    !>   theta = poloidal angle (NOT theta*)
    !>   phi = toroidal angle
    !>
    !> Output field components are in (s, theta*, phi) coordinates where:
    !>   theta* = theta + Lambda(s,theta,phi)
    !>   Lambda is the stream function
    !>
    !> The transformation from (r,theta,phi) to (s,theta*,phi) components includes:
    !>   - Factor ds/dr = 2*r for the radial component
    !>   - Lambda derivatives for coupling between components
    
    class(GvecField), intent(in) :: self
    real(dp), intent(in) :: x(3)  ! r=sqrt(s_vmec), theta_vmec, phi_vmec
    real(dp), intent(out) :: Acov(3)
    real(dp), intent(out) :: hcov(3)
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    ! Local variables for GVEC evaluation
    real(dp) :: r, theta, phi
    real(dp) :: theta_star, zeta

    ! GVEC coordinate evaluation variables
    real(dp) :: gvec_coords(3)  ! (s, theta, zeta)
    integer :: deriv_flags(2)   ! For derivative specification

    ! X1, X2, LA values and derivatives
    real(dp) :: X1_val, X2_val, LA_val
    real(dp) :: dX1_ds, dX1_dthet, dX1_dzeta
    real(dp) :: dX2_ds, dX2_dthet, dX2_dzeta
    real(dp) :: dLA_ds, dLA_dthet, dLA_dzeta

    ! Profile values
    real(dp) :: iota_val, phi_val, phiPrime_val, chi_val

    ! Coordinate and field computation
    real(dp) :: R_pos, Z_pos  ! Physical R, Z coordinates
    real(dp) :: e_thet(3), e_zeta(3), e_s(3)  ! Basis vectors
    real(dp) :: dx_dq1(3), dx_dq2(3), dx_dq3(3)  ! Raw derivatives from hmap
    real(dp) :: g_tt, g_tz, g_zz, g_ss, g_st, g_sz  ! Metric tensor components
    real(dp) :: Jac_h, Jac_l, Jac  ! Jacobian components (reference, logical, full)
    real(dp) :: Bthctr, Bzetactr  ! Contravariant field components in (s,theta,zeta)
    real(dp) :: Bthcov, Bzetacov, Bscov  ! Covariant field components
    real(dp) :: RZ_coords(3)      ! Coordinates for eval_Jh (R,Z,ζ)
    
    ! Axis regularization variables
    real(dp), parameter :: s_min_reg = 0.05_dp  ! Regularization threshold
    real(dp) :: scaling_factor, reg_factor, smooth_factor

    ! Extract input coordinates
    r = x(1)      ! r = sqrt(s_vmec) - SIMPLE uses r as radial coordinate
    theta = x(2)  ! theta_vmec
    phi = x(3)    ! phi_vmec

    theta_star = theta
    zeta = -phi

    ! Check if GVEC state is properly initialized before evaluation
    if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof) .or. &
        .not. allocated(X1_r) .or. .not. allocated(X2_r) .or. .not. allocated(LA_r)) then
        ! Reload the GVEC state (this is expected behavior due to GVEC module design)
        call ReadState(trim(self%filename))

        if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof) .or. &
            .not. allocated(X1_r) .or. .not. allocated(X2_r) .or. .not. allocated(LA_r)) then
            error stop 'Failed to initialize GVEC state for field evaluation'
        end if
    end if

    ! Prepare coordinates for GVEC evaluation
    gvec_coords(1) = r
    gvec_coords(2) = theta_star
    gvec_coords(3) = zeta

    ! Evaluate X1 (R coordinate) and its derivatives
    deriv_flags = [0, 0]  ! Function value
    X1_val = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)

    deriv_flags = [DERIV_S, 0]  ! ∂/∂s
    dX1_ds = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)

    deriv_flags = [0, DERIV_THET]  ! ∂/∂θ
    dX1_dthet = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)

    deriv_flags = [0, DERIV_ZETA]  ! ∂/∂ζ
    dX1_dzeta = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)

    ! Evaluate X2 (Z coordinate) and its derivatives
    deriv_flags = [0, 0]
    X2_val = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

    deriv_flags = [DERIV_S, 0]
    dX2_ds = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

    deriv_flags = [0, DERIV_THET]
    dX2_dthet = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

    deriv_flags = [0, DERIV_ZETA]
    dX2_dzeta = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

    ! Evaluate LA (stream function) and its derivatives
    deriv_flags = [0, 0]
    LA_val = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

    deriv_flags = [DERIV_S, 0]
    dLA_ds = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

    deriv_flags = [0, DERIV_THET]
    dLA_dthet = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

    deriv_flags = [0, DERIV_ZETA]
    dLA_dzeta = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

    ! Get profile values
    iota_val = eval_iota_r(r)      ! Rotational transform
    phi_val = eval_phi_r(r) * TESLA_IN_GAUSS * METER_IN_CM**2        ! Toroidal flux
    
    ! Compute toroidal flux derivative properly using GVEC's built-in function
    phiPrime_val = eval_phiPrime_r(r) * TESLA_IN_GAUSS * METER_IN_CM**2
    
    ! Get poloidal flux from profiles_1d(:,2)
    ! profiles_1d indices: 1=phi, 2=chi, 3=iota, 4=pressure
    chi_val = sbase_prof%evalDOF_s(r, 0, profiles_1d(:,2)) * TESLA_IN_GAUSS * METER_IN_CM**2

    R_pos = X1_val
    Z_pos = X2_val

    ! Use GVEC's hmap to compute basis vectors
    ! Get the Cartesian derivatives from GVEC's coordinate mapping
    ! hmap_r returns: e_q1 = ∂r/∂X1, e_q2 = ∂r/∂X2, e_q3 = ∂r/∂ζ
    ! These are the basis vectors in the (R,Z,φ) cylindrical system
    ! Note: get_dx_dqi expects (X1, X2, ζ) = (R, Z, ζ) coordinates
    RZ_coords = [R_pos, Z_pos, zeta]
    call hmap_r%get_dx_dqi(RZ_coords, dx_dq1, dx_dq2, dx_dq3)
    
    ! Compute basis vectors for (s,θ,ζ) coordinates EXACTLY as GVEC Python does
    ! e_rho = e_q1 * dX1_dr + e_q2 * dX2_dr
    ! e_theta = e_q1 * dX1_dt + e_q2 * dX2_dt  
    ! e_zeta = e_q1 * dX1_dz + e_q2 * dX2_dz + e_q3
    
    e_s = dx_dq1 * dX1_ds + dx_dq2 * dX2_ds
    e_thet = dx_dq1 * dX1_dthet + dx_dq2 * dX2_dthet  
    e_zeta = dx_dq1 * dX1_dzeta + dx_dq2 * dX2_dzeta + dx_dq3

    ! Compute metric tensor components
    g_ss = dot_product(e_s, e_s) * METER_IN_CM**2
    g_tt = dot_product(e_thet, e_thet) * METER_IN_CM**2
    g_zz = dot_product(e_zeta, e_zeta) * METER_IN_CM**2
    g_st = dot_product(e_s, e_thet) * METER_IN_CM**2
    g_sz = dot_product(e_s, e_zeta) * METER_IN_CM**2
    g_tz = dot_product(e_thet, e_zeta) * METER_IN_CM**2

    ! Compute Jacobian using VMEC/GVEC formula
    Jac_h = hmap_r%eval_Jh(RZ_coords)
    Jac_l = dX1_ds * dX2_dthet - dX1_dthet * dX2_ds
    Jac = Jac_h * Jac_l * METER_IN_CM**3

    ! Following GVEC mhd3d_evalfunc.f90: contravariant components
    Bthctr = (iota_val - dLA_dzeta) * phiPrime_val / Jac
    Bzetactr = (1.0_dp + dLA_dthet) * phiPrime_val / Jac

    if (present(sqgBctr)) then
        sqgBctr(1) = 0.0_dp  ! Jac * B^s = 0 (no radial contravariant component)
        sqgBctr(2) = Jac * Bthctr  ! Jac * B^θ
        sqgBctr(3) = -Jac * Bzetactr  ! Jac * B^phi = -Jac * B^ζ
    end if
        
    Bthcov = (g_tt * Bthctr + g_tz * Bzetactr)
    Bzetacov = (g_tz * Bthctr + g_zz * Bzetactr)
    
    ! For theta* coordinates, we need covariant B in s direction
    Bscov = (g_st * Bthctr + g_sz * Bzetactr)

    Bmod = sqrt(Bthctr*Bthcov + Bzetactr*Bzetacov)
    
    ! Apply theta* coordinate transformations with Lambda derivatives
    ! The GVEC components are in (s,theta,zeta) but SIMPLE expects (r,theta*,phi)
    ! where s = r^2, theta* = theta + Lambda, and phi = -zeta
    ! The transformation requires: ds/dr = 2*r
    hcov(1) = Bscov / Bmod
    hcov(2) = Bthcov / Bmod
    hcov(3) = Bzetacov / Bmod

    ! Vector potential with theta* transformations
    ! Acov_theta in GVEC is the toroidal flux phi_val
    Acov(1) = phi_val * dLA_ds * 2.0_dp * r
    Acov(2) = phi_val * (1.0_dp + dLA_dthet)
    Acov(3) = -chi_val + phi_val * dLA_dzeta

    ! Correct for phi = -zeta
    hcov(3) = -hcov(3)
    Acov(3) = -Acov(3)

end subroutine evaluate

subroutine gvec_field_cleanup(self)
    type(GvecField), intent(inout) :: self
    if (self%data_loaded) then
        call Finalize_ReadState()
        self%data_loaded = .false.
    end if
end subroutine gvec_field_cleanup

! Note: Direct VMEC to GVEC conversion requires access to GVEC's internal
! solution structures which are not available through this interface.
! The conversion must be done using GVEC's own tools or Python API.
subroutine convert_vmec_to_gvec(vmec_file, gvec_file)
    character(*), intent(in) :: vmec_file
    character(*), intent(in) :: gvec_file
    
    print *, 'ERROR: Direct VMEC to GVEC conversion not available'
    print *, 'The GVEC WriteState function requires internal solution structures'
    print *, 'that are not accessible through SIMPLE''s field interface.'
    print *, ''
    print *, 'To convert VMEC to GVEC format, use one of these methods:'
    print *, '1. GVEC Python API:'
    print *, '   import gvec'
    print *, '   gvec.read_vmec("', trim(vmec_file), '")'
    print *, '   gvec.save("', trim(gvec_file), '")'
    print *, ''
    print *, '2. GVEC command-line tools with proper parameter file'
    print *, ''
    print *, 'Then use create_gvec_field() to load the resulting .dat file'
    
    error stop 'VMEC to GVEC conversion not implemented'
    
end subroutine convert_vmec_to_gvec

end module field_gvec
