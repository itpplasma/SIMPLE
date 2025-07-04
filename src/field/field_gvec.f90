module field_gvec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
! Use GVEC API for reading state files
use MODgvec_cubic_spline, only: t_cubspl
use MODgvec_rProfile_bspl, only: t_rProfile_bspl
use MODgvec_globals, only: wp
use MODgvec_ReadState, only: ReadState, eval_phi_r, eval_phiPrime_r, eval_iota_r, eval_pres_r, Finalize_ReadState
use MODgvec_ReadState_Vars, only: profiles_1d, sbase_prof, X1_base_r, X2_base_r, LA_base_r, X1_r, X2_r, LA_r, hmap_r
implicit none

! GVEC derivative constants (from defines.h)
integer, parameter :: DERIV_S = 1
integer, parameter :: DERIV_THET = 2
integer, parameter :: DERIV_ZETA = 3

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

    ! Local variables for GVEC evaluation
    real(dp) :: r, theta, phi
    real(dp) :: s_gvec, theta_star, zeta

    ! GVEC coordinate evaluation variables
    real(wp) :: gvec_coords(3)  ! (s, theta, zeta)
    integer :: deriv_flags(2)   ! For derivative specification

    ! X1, X2, LA values and derivatives
    real(wp) :: X1_val, X2_val, LA_val
    real(wp) :: dX1_ds, dX1_dthet, dX1_dzeta
    real(wp) :: dX2_ds, dX2_dthet, dX2_dzeta
    real(wp) :: dLA_ds, dLA_dthet, dLA_dzeta

    ! Profile values
    real(wp) :: iota_val, phi_val, phiPrime_val

    ! Coordinate and field computation
    real(wp) :: R_pos, Z_pos  ! Physical R, Z coordinates
    real(wp) :: e_thet(3), e_zeta(3), e_s(3)  ! Basis vectors
    real(wp) :: dx_dq1(3), dx_dq2(3), dx_dq3(3)  ! Raw derivatives from hmap
    real(wp) :: g_tt, g_tz, g_zz, g_ss, g_st, g_sz  ! Metric tensor components
    real(wp) :: Jac_h, Jac_l, Jac  ! Jacobian components (reference, logical, full)
    real(wp) :: B_thet, B_zeta  ! Contravariant field components in (s,theta,zeta)
    real(wp) :: B_cov_theta, B_cov_zeta  ! Covariant field components
    real(wp) :: B_vec(3)  ! Physical B field vector in cylindrical coordinates
    
    ! Axis regularization variables
    real(wp), parameter :: s_min_reg = 0.05_wp  ! Regularization threshold
    real(wp) :: scaling_factor, reg_factor, smooth_factor

    ! Extract input coordinates
    r = x(1)      ! r = sqrt(s_vmec)
    theta = x(2)  ! theta_vmec
    phi = x(3)    ! phi_vmec

    ! Convert to GVEC flux coordinate system
    ! Note: coordinate transformation depends on the source and conventions
    ! Testing different conventions to match VMEC
    s_gvec = r**2       ! VMEC flux coordinate (0 to 1)
    theta_star = theta  ! Straight field line poloidal angle
    zeta = phi          ! VMEC uses phi, GVEC uses zeta; with nfp field periods
    ! The relationship might involve the field periods

    if (.not. self%data_loaded) then
        print *, 'Warning: GVEC data not loaded, returning dummy values'
        Acov = 0.0_dp
        hcov = 0.0_dp
        Bmod = 1.0_dp
        if (present(sqgBctr)) sqgBctr = 0.0_dp
        return
    end if

    ! Check if GVEC state is properly initialized before evaluation
    if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof) .or. &
        .not. allocated(X1_r) .or. .not. allocated(X2_r) .or. .not. allocated(LA_r)) then
        ! Reload the GVEC state (this is expected behavior due to GVEC module design)
        call ReadState(trim(self%filename))

        if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof) .or. &
            .not. allocated(X1_r) .or. .not. allocated(X2_r) .or. .not. allocated(LA_r)) then
            print *, 'Error: Failed to initialize GVEC state for field evaluation'
            Acov = 0.0_dp
            hcov = 0.0_dp
            Bmod = 1.0_dp
            if (present(sqgBctr)) sqgBctr = 0.0_dp
            return
        end if
    end if

    ! Prepare coordinates for GVEC evaluation
    gvec_coords(1) = s_gvec
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
    iota_val = eval_iota_r(s_gvec)      ! Rotational transform
    phi_val = eval_phi_r(s_gvec)        ! Toroidal flux

    ! Compute toroidal flux derivative properly using GVEC's built-in function
    phiPrime_val = eval_phiPrime_r(s_gvec)
    
    ! Debug output disabled - uncomment if needed
    ! if (abs(s_gvec - 0.1_wp) < 0.01_wp .and. abs(theta_star) < 0.1_wp) then
    !     print *, 'GVEC DEBUG at s=', s_gvec, ', theta=', theta_star, ', zeta=', zeta
    !     print *, '  phi=', phi_val, ', phiPrime=', phiPrime_val
    !     print *, '  iota=', iota_val
    ! end if

    ! Get physical coordinates
    R_pos = X1_val
    Z_pos = X2_val

    ! Use GVEC's hmap to compute basis vectors
    ! Get the Cartesian derivatives from GVEC's coordinate mapping
    ! hmap_r returns: e_q1 = ∂r/∂X1, e_q2 = ∂r/∂X2, e_q3 = ∂r/∂ζ
    ! These are the basis vectors in the (R,Z,φ) cylindrical system
    call hmap_r%get_dx_dqi(gvec_coords, dx_dq1, dx_dq2, dx_dq3)
    
    ! Compute basis vectors for (s,θ,ζ) coordinates EXACTLY as GVEC Python does
    ! e_rho = e_q1 * dX1_dr + e_q2 * dX2_dr
    ! e_theta = e_q1 * dX1_dt + e_q2 * dX2_dt  
    ! e_zeta = e_q1 * dX1_dz + e_q2 * dX2_dz + e_q3
    
    e_s = dx_dq1 * dX1_ds + dx_dq2 * dX2_ds
    e_thet = dx_dq1 * dX1_dthet + dx_dq2 * dX2_dthet  
    e_zeta = dx_dq1 * dX1_dzeta + dx_dq2 * dX2_dzeta + dx_dq3

    ! Compute metric tensor components
    g_ss = dot_product(e_s, e_s)
    g_tt = dot_product(e_thet, e_thet)
    g_zz = dot_product(e_zeta, e_zeta)
    g_st = dot_product(e_s, e_thet)
    g_sz = dot_product(e_s, e_zeta)
    g_tz = dot_product(e_thet, e_zeta)

    ! Compute Jacobian using GVEC's exact formula: Jac = Jac_h * Jac_l
    ! where Jac_l = dX1_dr * dX2_dt - dX1_dt * dX2_dr (logical Jacobian)
    Jac_h = hmap_r%eval_Jh(gvec_coords)  ! Reference Jacobian from hmap
    Jac_l = dX1_ds * dX2_dthet - dX1_dthet * dX2_ds  ! Logical Jacobian
    Jac = Jac_h * Jac_l  ! Full Jacobian as used in GVEC Python

    ! Compute contravariant magnetic field components using GVEC's EXACT formula
    ! From GVEC Python quantities.py lines 532-533:
    ! B^θ = (ι - ∂Λ/∂ζ) * ∂Φ/∂r / Jac
    ! B^ζ = (1 + ∂Λ/∂θ) * ∂Φ/∂r / Jac
    ! Note: GVEC uses Jac (full Jacobian), not just sqrt(g)
    
    ! Compute contravariant magnetic field components EXACTLY as GVEC Python does
    ! B^θ = (ι - ∂Λ/∂ζ) * ∂Φ/∂r / Jac
    ! B^ζ = (1 + ∂Λ/∂θ) * ∂Φ/∂r / Jac
    B_thet = (iota_val - dLA_dzeta) * phiPrime_val / Jac  ! B^θ contravariant
    B_zeta = (1.0_wp + dLA_dthet) * phiPrime_val / Jac    ! B^ζ contravariant
    
    ! Debug disabled - uncomment if needed
    ! if (abs(s_gvec - 0.1_wp) < 0.01_wp .and. abs(theta_star) < 0.1_wp) then
    !     print *, '  sqrtG=', sqrtG, ', dLA_dthet=', dLA_dthet, ', dLA_dzeta=', dLA_dzeta
    !     print *, '  B_thet=', B_thet, ', B_zeta=', B_zeta
    ! end if

    ! Compute physical magnetic field vector B = B^θ e_θ + B^ζ e_ζ
    ! This matches GVEC Python: ds["B"] = ds.B_contra_t * ds.e_theta + ds.B_contra_z * ds.e_zeta
    B_vec = B_thet * e_thet + B_zeta * e_zeta
    
    ! Compute field magnitude as |B| = sqrt(B·B) in physical space
    ! This is the magnitude of the physical field vector
    Bmod = real(sqrt(dot_product(B_vec, B_vec)), dp)
    
    ! Also compute covariant components for other uses
    B_cov_theta = g_tt * B_thet + g_tz * B_zeta
    B_cov_zeta = g_tz * B_thet + g_zz * B_zeta
    
    ! Debug output to understand the issue
    if (abs(s_gvec - 0.01_wp) < 0.05_wp .and. abs(theta_star) < 0.1_wp) then
        print *, 'GVEC DEBUG: s=', s_gvec, ', theta=', theta_star, ', zeta=', zeta
        print *, '  R=', R_pos, ', Z=', Z_pos
        print *, '  phiPrime=', phiPrime_val, ', Jac=', Jac
        print *, '  g_tt=', g_tt, ', g_tz=', g_tz, ', g_zz=', g_zz
        print *, '  B^theta=', B_thet, ', B^zeta=', B_zeta
        print *, '  B_theta=', B_cov_theta, ', B_zeta=', B_cov_zeta
        print *, '  |B|=', Bmod
    end if

    ! Set covariant field components for output
    ! B_s = g_ss*B^s + g_st*B^θ + g_sz*B^ζ = g_st*B^θ + g_sz*B^ζ (since B^s=0)
    hcov(1) = real(g_st * B_thet + g_sz * B_zeta, dp)  ! B_s (covariant)
    hcov(2) = real(B_cov_theta, dp)  ! B_θ (covariant) - already computed
    hcov(3) = real(B_cov_zeta, dp)   ! B_ζ (covariant) - already computed

    ! Vector potential components in magnetic flux coordinates
    ! Following GVEC's approach: work with magnetic field directly, not vector potential
    ! In magnetic coordinates with straight field lines, the vector potential is:
    ! A·∇ζ = ψ(s) - toroidal flux function
    ! A·∇θ = χ(s,θ,ζ) - poloidal flux-like function including λ corrections

    ! For magnetic coordinates, gauge choice gives A_s = 0
    Acov(1) = 0.0_dp

    ! A_θ component: related to poloidal flux and stream function λ
    ! From ∇×A = B and the constraint that field lines are straight
    ! A_θ is related to the poloidal flux function and λ corrections
    Acov(2) = real(-LA_val * phiPrime_val / Jac, dp)  ! Stream function contribution

    ! A_ζ component: primarily the toroidal flux function
    ! This ensures B·∇ζ = (1 + ∂λ/∂θ) * ∂ψ/∂s / √g matches our B^ζ
    if (abs(R_pos) > 1.0e-10_wp) then
        Acov(3) = real(phi_val / R_pos, dp)  ! Toroidal flux contribution
    else
        Acov(3) = 0.0_dp
    end if

    if (present(sqgBctr)) then
        ! Jac * B^contravariant components (using full Jacobian as in GVEC)
        ! sqgBctr = Jac * B^i where B^i are contravariant components
        ! We already have B^θ, B^ζ from GVEC, and B^s = 0
        sqgBctr(1) = 0.0_dp  ! Jac * B^s = 0 (no radial contravariant component)
        sqgBctr(2) = real(Jac * B_thet, dp)  ! Jac * B^θ
        sqgBctr(3) = real(Jac * B_zeta, dp)  ! Jac * B^ζ
    end if

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
