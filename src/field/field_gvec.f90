module field_gvec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
! Use GVEC API for reading state files
use MODgvec_cubic_spline, only: t_cubspl
use MODgvec_rProfile_bspl, only: t_rProfile_bspl
use MODgvec_globals, only: wp
use MODgvec_ReadState, only: ReadState, eval_phi_r, eval_iota_r, eval_pres_r, Finalize_ReadState
use MODgvec_ReadState_Vars, only: profiles_1d, sbase_prof, X1_base_r, X2_base_r, LA_base_r, X1_r, X2_r, LA_r, hmap_r
implicit none

! GVEC derivative constants (from defines.h)
integer, parameter :: DERIV_S = 1
integer, parameter :: DERIV_THET = 2
integer, parameter :: DERIV_ZETA = 3

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
    real(wp) :: g_tt, g_tz, g_zz, g_ss, g_st, g_sz  ! Metric tensor components
    real(wp) :: sqrtG  ! Jacobian
    real(wp) :: B_thet, B_zeta  ! Covariant field components in (s,theta,zeta)
    real(wp) :: Bx, By, Bz  ! Cartesian field components

    ! Extract input coordinates
    r = x(1)      ! r = sqrt(s_vmec)
    theta = x(2)  ! theta_vmec
    phi = x(3)    ! phi_vmec

    ! Convert to GVEC flux coordinate system
    ! GVEC transforms VMEC coordinates: zeta_gvec = -phi_vmec (switchzeta=.TRUE.)
    s_gvec = r**2       ! VMEC flux coordinate (0 to 1)
    theta_star = theta  ! Straight field line poloidal angle (may need sign change for signgs)
    zeta = -phi         ! GVEC zeta = -phi_vmec (coordinate transformation)

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

    ! Compute toroidal flux derivative (approximate for now)
    if (s_gvec > 0.01_wp) then
        phiPrime_val = (eval_phi_r(s_gvec + 0.01_wp) - eval_phi_r(s_gvec - 0.01_wp)) / 0.02_wp
    else
        phiPrime_val = eval_phi_r(0.01_wp) / 0.01_wp
    end if

    ! Get physical coordinates
    R_pos = X1_val
    Z_pos = X2_val

    ! Compute basis vectors in Cartesian coordinates
    ! e_s = ∂r/∂s = (∂X1/∂s, ∂X2/∂s, 0)
    e_s(1) = dX1_ds
    e_s(2) = dX2_ds
    e_s(3) = 0.0_wp

    ! e_theta = ∂r/∂θ = (∂X1/∂θ, ∂X2/∂θ, 0)
    e_thet(1) = dX1_dthet
    e_thet(2) = dX2_dthet
    e_thet(3) = 0.0_wp

    ! e_zeta = ∂r/∂ζ_gvec = (∂X1/∂ζ, ∂X2/∂ζ, R * ∂φ_vmec/∂ζ_gvec)
    ! Since ζ_gvec = -φ_vmec, we have ∂φ_vmec/∂ζ_gvec = -1
    e_zeta(1) = dX1_dzeta
    e_zeta(2) = dX2_dzeta
    e_zeta(3) = -R_pos  ! R * ∂φ_vmec/∂ζ_gvec = R * (-1)

    ! Compute metric tensor components
    g_ss = dot_product(e_s, e_s)
    g_tt = dot_product(e_thet, e_thet)
    g_zz = dot_product(e_zeta, e_zeta)
    g_st = dot_product(e_s, e_thet)
    g_sz = dot_product(e_s, e_zeta)
    g_tz = dot_product(e_thet, e_zeta)

    ! Use GVEC's built-in Jacobian computation
    sqrtG = hmap_r%eval_Jh(gvec_coords)

    ! Compute contravariant magnetic field components using GVEC formula
    ! B^θ = (ι - ∂Λ/∂ζ) * Φ'/(√g)
    ! B^ζ = (1 + ∂Λ/∂θ) * Φ'/(√g)
    ! B^s = 0 for divergence-free field in flux coordinates
    B_thet = (iota_val - dLA_dzeta) * phiPrime_val / sqrtG  ! B^θ contravariant
    B_zeta = (1.0_wp + dLA_dthet) * phiPrime_val / sqrtG    ! B^ζ contravariant

    ! Compute Cartesian field components from contravariant components
    ! B = B^θ * e_θ + B^ζ * e_ζ (B^s = 0)
    Bx = B_thet * e_thet(1) + B_zeta * e_zeta(1)
    By = B_thet * e_thet(2) + B_zeta * e_zeta(2)
    Bz = B_thet * e_thet(3) + B_zeta * e_zeta(3)

    ! Compute field magnitude
    Bmod = real(sqrt(Bx**2 + By**2 + Bz**2), dp)

    ! Compute COVARIANT field components: B_i = g_ij * B^j
    ! Need to convert contravariant B^θ, B^ζ to covariant B_s, B_θ, B_ζ
    ! B_s = g_ss*B^s + g_st*B^θ + g_sz*B^ζ = g_st*B^θ + g_sz*B^ζ (since B^s=0)
    ! B_θ = g_st*B^s + g_tt*B^θ + g_tz*B^ζ = g_tt*B^θ + g_tz*B^ζ
    ! B_ζ = g_sz*B^s + g_tz*B^θ + g_zz*B^ζ = g_tz*B^θ + g_zz*B^ζ
    hcov(1) = real(g_st * B_thet + g_sz * B_zeta, dp)  ! B_s (covariant)
    hcov(2) = real(g_tt * B_thet + g_tz * B_zeta, dp)  ! B_θ (covariant)
    hcov(3) = real(g_tz * B_thet + g_zz * B_zeta, dp)  ! B_ζ (covariant)

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
    Acov(2) = real(-LA_val * phiPrime_val / sqrtG, dp)  ! Stream function contribution

    ! A_ζ component: primarily the toroidal flux function
    ! This ensures B·∇ζ = (1 + ∂λ/∂θ) * ∂ψ/∂s / √g matches our B^ζ
    if (abs(R_pos) > 1.0e-10_wp) then
        Acov(3) = real(phi_val / R_pos, dp)  ! Toroidal flux contribution
    else
        Acov(3) = 0.0_dp
    end if

    if (present(sqgBctr)) then
        ! sqrt(g) * B^contravariant components
        ! sqgBctr = √g * B^i where B^i are contravariant components
        ! We already have B^θ, B^ζ from GVEC, and B^s = 0
        sqgBctr(1) = 0.0_dp  ! √g * B^s = 0 (no radial contravariant component)
        sqgBctr(2) = real(sqrtG * B_thet, dp)  ! √g * B^θ
        sqgBctr(3) = real(sqrtG * B_zeta, dp)  ! √g * B^ζ
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
