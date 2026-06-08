module field_boozer_chartmap
    !> Boozer field from extended chartmap NetCDF file.
    !>
    !> Reads an extended chartmap file containing Boozer coordinate geometry
    !> (X, Y, Z on a Boozer angle grid) plus magnetic field data:
    !>   - A_phi(rho)        : toroidal vector potential (1D)
    !>   - B_theta(rho)      : covariant poloidal B component (1D, surface function)
    !>   - B_phi(rho)        : covariant toroidal B component (1D, surface function)
    !>   - Bmod(rho,theta,zeta) : field magnitude (3D)
    !>   - torflux           : total toroidal flux (scalar attribute)
    !>
    !> This type serves as a magnetic_field_t for the non-canonical path and
    !> provides data to populate the Boozer batch splines for the symplectic path.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    use interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                           construct_batch_splines_1d, construct_batch_splines_3d, &
                           evaluate_batch_splines_1d_der2, &
                           evaluate_batch_splines_3d, &
                           destroy_batch_splines_1d, destroy_batch_splines_3d
    use new_vmec_stuff_mod, only: vmec_B_scale, vmec_RZ_scale, nper, rmajor
    use boozer_chartmap_io, only: boozer_chartmap_data_t, read_boozer_chartmap
    use netcdf
    use scaled_chartmap_coordinates, only: wrap_scaled_chartmap_coordinate_system

    implicit none

    private
    public :: boozer_chartmap_field_t, create_boozer_chartmap_field, &
              is_boozer_chartmap

    type, extends(magnetic_field_t) :: boozer_chartmap_field_t
        type(BatchSplineData1D) :: aphi_spline
        type(BatchSplineData1D) :: bcovar_spline
        type(BatchSplineData3D) :: bmod_spline
        real(dp) :: torflux = 0.0_dp
        integer :: nfp = 1
        integer :: ns = 0
        integer :: ntheta = 0
        integer :: nphi = 0
        real(dp) :: hs = 0.0_dp
        real(dp) :: h_theta = 0.0_dp
        real(dp) :: h_phi = 0.0_dp
        character(len=512) :: filename = ''
        logical :: initialized = .false.
    contains
        procedure :: evaluate => boozer_chartmap_evaluate
        final :: boozer_chartmap_cleanup
    end type boozer_chartmap_field_t

contains

    function is_boozer_chartmap(filename) result(is_bc)
        character(len=*), intent(in) :: filename
        logical :: is_bc
        integer :: ncid, status, ival

        is_bc = .false.
        status = nf90_open(trim(filename), nf90_nowrite, ncid)
        if (status /= nf90_noerr) return

        status = nf90_get_att(ncid, nf90_global, "boozer_field", ival)
        if (status == nf90_noerr .and. ival == 1) is_bc = .true.

        status = nf90_close(ncid)
    end function is_boozer_chartmap

    subroutine create_boozer_chartmap_field(filename, field)
        use libneo_coordinates, only: coordinate_system_t, &
                                      make_chartmap_coordinate_system

        character(len=*), intent(in) :: filename
        type(boozer_chartmap_field_t), allocatable, intent(out) :: field

        type(boozer_chartmap_data_t) :: d
        real(dp), allocatable :: y_aphi(:, :), y_bcovar(:, :), y_bmod(:, :, :, :)
        real(dp) :: torflux_val
        real(dp) :: b_scale, rz_scale, covar_scale, flux_scale
        integer, parameter :: spline_order_1d = 5
        integer, parameter :: spline_order_3d(3) = [5, 5, 5]
        logical, parameter :: periodic_3d(3) = [.false., .true., .true.]
        real(dp) :: x_min_3d(3), x_max_3d(3)
        class(coordinate_system_t), allocatable :: cs

        allocate (field)

        ! Single shared parse: base-unit arrays on the endpoint-included field grid.
        call read_boozer_chartmap(filename, d)

        b_scale = vmec_B_scale
        rz_scale = vmec_RZ_scale
        covar_scale = b_scale*rz_scale
        flux_scale = covar_scale*rz_scale

        d%A_phi = flux_scale*d%A_phi
        d%B_theta = covar_scale*d%B_theta
        d%B_phi = covar_scale*d%B_phi
        d%Bmod = b_scale*d%Bmod
        torflux_val = flux_scale*d%torflux

        ! Build A_phi spline over the file's uniform rho grid, like B_theta/B_phi
        ! and Bmod. A_phi is a function of rho here; evaluate at rho (see below).
        allocate (y_aphi(d%n_rho, 1))
        y_aphi(:, 1) = d%A_phi
        call construct_batch_splines_1d(d%rho_min, d%rho_max, y_aphi, spline_order_1d, &
                                        .false., field%aphi_spline)

        ! Build B_theta, B_phi spline over rho_tor
        allocate (y_bcovar(d%n_rho, 2))
        y_bcovar(:, 1) = d%B_theta
        y_bcovar(:, 2) = d%B_phi
        call construct_batch_splines_1d(d%rho_min, d%rho_max, y_bcovar, spline_order_1d, &
                                        .false., field%bcovar_spline)

        ! Build Bmod 3D spline over (rho_tor, theta_B, phi_B). The field grid is
        ! endpoint-included, so x_max spans the full 2*pi and 2*pi/nfp period
        ! that the periodic spline expects (period = (n-1)*h_step).
        x_min_3d = [d%rho_min, 0.0_dp, 0.0_dp]
        x_max_3d = [d%rho_max, real(d%n_theta - 1, dp)*d%h_theta, &
                    real(d%n_phi - 1, dp)*d%h_phi]

        allocate (y_bmod(d%n_rho, d%n_theta, d%n_phi, 1))
        y_bmod(:, :, :, 1) = d%Bmod
        call construct_batch_splines_3d(x_min_3d, x_max_3d, y_bmod, &
                                        spline_order_3d, periodic_3d, &
                                        field%bmod_spline)

        ! Store metadata
        field%torflux = torflux_val
        field%nfp = d%nfp
        field%ns = d%n_rho
        field%ntheta = d%n_theta
        field%nphi = d%n_phi
        field%hs = d%h_s
        field%h_theta = d%h_theta
        field%h_phi = d%h_phi
        field%filename = filename

        ! Restore equilibrium periods/major radius so stevvo and params_init
        ! produce the correct dphi/dtaumin/fper for a VMEC-free chartmap run.
        nper = d%nfp
        rmajor = d%rmajor*rz_scale

        ! Set up coordinate system from the same chartmap file
        call make_chartmap_coordinate_system(cs, filename)
        call wrap_scaled_chartmap_coordinate_system(cs, rz_scale)
        allocate (field%coords, source=cs)

        field%initialized = .true.

        print *, 'Loaded Boozer chartmap field from ', trim(filename)
        print *, '  nfp=', d%nfp, ' ns=', d%n_rho, ' ntheta=', d%n_theta, &
                 ' nphi=', d%n_phi
        print *, '  torflux=', torflux_val
    end subroutine create_boozer_chartmap_field

    subroutine boozer_chartmap_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        !> Evaluate field at x = (rho, theta_B, phi_B) in chartmap coordinates.
        !> Input x(1) = rho = sqrt(s), matching the chartmap coordinate system.
        !> Returns covariant components in (rho, theta_B, phi_B) coordinates.
        class(boozer_chartmap_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: rho_tor, s, theta_B, phi_B
        real(dp) :: A_phi_val, A_theta_val
        real(dp) :: y1d_aphi(1), dy1d_aphi(1), d2y1d_aphi(1)
        real(dp) :: y1d_bc(2), dy1d_bc(2), d2y1d_bc(2)
        real(dp) :: x_3d(3), y_bmod(1)
        real(dp) :: B_theta_val, B_phi_val, Bmod_val
        real(dp), parameter :: twopi = 8.0_dp * atan(1.0_dp)

        rho_tor = x(1)
        theta_B = x(2)
        phi_B = x(3)
        s = rho_tor**2

        ! A_theta = torflux * s (linear, by definition)
        ! In rho coordinates: A_theta(rho) = torflux * rho^2
        ! Covariant component transforms: A_rho = A_s * ds/drho = A_s * 2*rho
        ! But A_s = 0 (gauge), so A_rho = 0 regardless.
        A_theta_val = self%torflux * s

        ! A_phi from 1D spline at rho (tabulated on the file's rho grid)
        call evaluate_batch_splines_1d_der2(self%aphi_spline, rho_tor, y1d_aphi, &
                                             dy1d_aphi, d2y1d_aphi)
        A_phi_val = y1d_aphi(1)

        ! B_theta, B_phi from 1D spline at rho_tor
        call evaluate_batch_splines_1d_der2(self%bcovar_spline, rho_tor, y1d_bc, &
                                             dy1d_bc, d2y1d_bc)
        B_theta_val = y1d_bc(1)
        B_phi_val = y1d_bc(2)

        ! Bmod from 3D spline at (rho_tor, theta_B, phi_B)
        x_3d(1) = rho_tor
        x_3d(2) = modulo(theta_B, twopi)
        x_3d(3) = modulo(phi_B, twopi / real(self%nfp, dp))

        call evaluate_batch_splines_3d(self%bmod_spline, x_3d, y_bmod)
        Bmod_val = y_bmod(1)

        Bmod = Bmod_val

        ! Covariant vector potential in (rho, theta, phi) coordinates.
        ! A_theta and A_phi are the same in (s, theta, phi) and (rho, theta, phi)
        ! since the angular coordinates are unchanged. A_rho = 0.
        Acov(1) = 0.0_dp
        Acov(2) = A_theta_val
        Acov(3) = A_phi_val

        ! Covariant unit field direction: h = B/|B|
        ! B_theta and B_phi are covariant and independent of the radial
        ! coordinate choice (they refer to the angular basis vectors).
        hcov(1) = 0.0_dp
        hcov(2) = B_theta_val / Bmod_val
        hcov(3) = B_phi_val / Bmod_val

        if (present(sqgBctr)) then
            error stop "sqgBctr not implemented for boozer_chartmap_field_t"
        end if
    end subroutine boozer_chartmap_evaluate

    subroutine boozer_chartmap_cleanup(self)
        type(boozer_chartmap_field_t), intent(inout) :: self

        if (self%initialized) then
            call destroy_batch_splines_1d(self%aphi_spline)
            call destroy_batch_splines_1d(self%bcovar_spline)
            call destroy_batch_splines_3d(self%bmod_spline)
            self%initialized = .false.
        end if
    end subroutine boozer_chartmap_cleanup

end module field_boozer_chartmap
