module field_gvec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: MagneticField
    use MODgvec_cubic_spline, only: t_cubspl
    use MODgvec_rProfile_bspl, only: t_rProfile_bspl
    use MODgvec_ReadState, only: ReadState, eval_phi_r, eval_phiPrime_r, eval_iota_r, eval_pres_r, Finalize_ReadState
    use MODgvec_ReadState_Vars, only: profiles_1d, sbase_prof, X1_base_r, X2_base_r, LA_base_r, X1_r, X2_r, LA_r, hmap_r
    implicit none

    integer, parameter :: DERIV_R = 1
    integer, parameter :: DERIV_THET = 2
    integer, parameter :: DERIV_ZETA = 3

    real(dp), parameter :: TESLA_IN_GAUSS = 10000.0_dp
    real(dp), parameter :: METER_IN_CM = 100.0_dp

    private
    public :: GvecField, create_gvec_field

    type, extends(MagneticField) :: GvecField
        character(len=256) :: filename = ''

        ! Field periods
        integer :: nfp = 1

        ! Data loaded flag
        logical :: data_loaded = .false.

    contains
        procedure :: evaluate
        procedure :: load_dat_file
        final :: gvec_field_cleanup
    end type GvecField

contains

    subroutine create_gvec_field(gvec_file, gvec_field)
        character(*), intent(in) :: gvec_file
        class(GvecField), allocatable, intent(out) :: gvec_field

        allocate (GvecField :: gvec_field)
        gvec_field%filename = gvec_file

        ! Load the .dat file
        call gvec_field%load_dat_file()

        if (.not. gvec_field%data_loaded) then
            print *, 'ERROR: Failed to load GVEC field from: ', gvec_file
        end if

    end subroutine create_gvec_field

    subroutine load_dat_file(self)
        class(GvecField), intent(inout) :: self
        logical :: file_exists

        inquire (file=trim(self%filename), exist=file_exists)
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
                self%nfp = 1
                print *, 'Warning: Could not read nfp from GVEC state, using default nfp=1'
            end if

            self%data_loaded = .true.
        else
            print *, 'Error: GVEC state initialization failed'
            self%data_loaded = .false.
        end if

    end subroutine load_dat_file

    subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)

        class(GvecField), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3)
        real(dp), intent(out) :: hcov(3)
        real(dp), intent(out) :: Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: r, theta, varphi, zeta
        real(dp) :: gvec_coords(3), RZ_coords(3)
        integer :: deriv_flags(2)
        real(dp) :: X1_val, X2_val, R_pos, Z_pos
        real(dp) :: dX1_ds, dX1_dthet, dX1_dzeta
        real(dp) :: dX2_ds, dX2_dthet, dX2_dzeta
        real(dp) :: dLA_dr, dLA_dthet, dLA_dzeta
        real(dp) :: iota_val, phi_val, phiPrime_val, chi_val
        real(dp) :: e_thet(3), e_zeta(3), e_s(3)
        real(dp) :: dx_dq1(3), dx_dq2(3), dx_dq3(3)
        real(dp) :: g_tt, g_tz, g_zz, g_ss, g_st, g_sz
        real(dp) :: Jac_h, Jac_l, Jac
        real(dp) :: Bthctr, Bzetactr, Bthcov, Bzetacov, Bscov
        r = x(1)
        theta = x(2)
        varphi = x(3)

        zeta = -varphi
        if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof) .or. &
            .not. allocated(X1_r) .or. .not. allocated(X2_r) .or. .not. allocated(LA_r)) then
            call ReadState(trim(self%filename))

            if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof) .or. &
                .not. allocated(X1_r) .or. .not. allocated(X2_r) .or. .not. allocated(LA_r)) then
                error stop 'Failed to initialize GVEC state for field evaluation'
            end if
        end if
        gvec_coords(1) = r
        gvec_coords(2) = theta
        gvec_coords(3) = zeta

        deriv_flags = [0, 0]
        X1_val = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)

        deriv_flags = [DERIV_R, 0]
        dX1_ds = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)

        deriv_flags = [0, DERIV_THET]
        dX1_dthet = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)

        deriv_flags = [0, DERIV_ZETA]
        dX1_dzeta = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)
        deriv_flags = [0, 0]
        X2_val = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

        deriv_flags = [DERIV_R, 0]
        dX2_ds = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

        deriv_flags = [0, DERIV_THET]
        dX2_dthet = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

        deriv_flags = [0, DERIV_ZETA]
        dX2_dzeta = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)
        deriv_flags = [DERIV_R, 0]
        dLA_dr = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

        deriv_flags = [0, DERIV_THET]
        dLA_dthet = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

        deriv_flags = [0, DERIV_ZETA]
        dLA_dzeta = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

        iota_val = eval_iota_r(r)
        phi_val = eval_phi_r(r)*TESLA_IN_GAUSS*METER_IN_CM**2
        phiPrime_val = eval_phiPrime_r(r)*TESLA_IN_GAUSS*METER_IN_CM**2
        chi_val = sbase_prof%evalDOF_s(r, 0, profiles_1d(:, 2))*TESLA_IN_GAUSS*METER_IN_CM**2

        R_pos = X1_val
        Z_pos = X2_val
        RZ_coords = [R_pos, Z_pos, zeta]
        call hmap_r%get_dx_dqi(RZ_coords, dx_dq1, dx_dq2, dx_dq3)

        e_s = dx_dq1*dX1_ds + dx_dq2*dX2_ds
        e_thet = dx_dq1*dX1_dthet + dx_dq2*dX2_dthet
        e_zeta = dx_dq1*dX1_dzeta + dx_dq2*dX2_dzeta + dx_dq3
        g_ss = dot_product(e_s, e_s)*METER_IN_CM**2
        g_tt = dot_product(e_thet, e_thet)*METER_IN_CM**2
        g_zz = dot_product(e_zeta, e_zeta)*METER_IN_CM**2
        g_st = dot_product(e_s, e_thet)*METER_IN_CM**2
        g_sz = dot_product(e_s, e_zeta)*METER_IN_CM**2
        g_tz = dot_product(e_thet, e_zeta)*METER_IN_CM**2
        Jac_h = hmap_r%eval_Jh(RZ_coords)
        Jac_l = dX1_ds*dX2_dthet - dX1_dthet*dX2_ds
        Jac = Jac_h*Jac_l*METER_IN_CM**3
        Bthctr = (iota_val - dLA_dzeta)*phiPrime_val/Jac
        Bzetactr = (1.0_dp + dLA_dthet)*phiPrime_val/Jac

        if (present(sqgBctr)) then
            sqgBctr(1) = 0.0_dp
            sqgBctr(2) = (-Jac)*Bthctr
            sqgBctr(3) = (-Jac)*(-Bzetactr)
        end if

        Bthcov = (g_tt*Bthctr + g_tz*Bzetactr)
        Bzetacov = (g_tz*Bthctr + g_zz*Bzetactr)
        Bscov = (g_st*Bthctr + g_sz*Bzetactr)

        Bmod = sqrt(Bthctr*Bthcov + Bzetactr*Bzetacov)
        hcov(1) = Bscov/Bmod
        hcov(2) = Bthcov/Bmod
        hcov(3) = -Bzetacov/Bmod
        Acov(1) = phi_val*dLA_dr
        Acov(2) = phi_val*(1.0_dp + dLA_dthet)
        Acov(3) = -(-chi_val + phi_val*dLA_dzeta)

    end subroutine evaluate

    subroutine gvec_field_cleanup(self)
        type(GvecField), intent(inout) :: self
        if (self%data_loaded) then
            call Finalize_ReadState()
            self%data_loaded = .false.
        end if
    end subroutine gvec_field_cleanup

end module field_gvec
