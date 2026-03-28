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
    use netcdf

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

        integer :: ncid, status, dimid, varid
        integer :: n_rho, n_theta, n_zeta, nfp_file
        real(dp) :: torflux_val
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: A_phi_arr(:), B_theta_arr(:), B_phi_arr(:)
        real(dp), allocatable :: Bmod_arr(:, :, :)
        real(dp), allocatable :: s_grid(:)
        real(dp), allocatable :: y_aphi(:, :), y_bcovar(:, :), y_bmod(:, :, :, :)
        real(dp) :: s_min, s_max, rho_min, rho_max
        real(dp) :: h_s, h_theta_val, h_phi_val
        integer :: i
        integer, parameter :: spline_order_1d = 5
        integer, parameter :: spline_order_3d(3) = [5, 5, 5]
        logical, parameter :: periodic_3d(3) = [.false., .true., .true.]
        real(dp) :: x_min_3d(3), x_max_3d(3)
        class(coordinate_system_t), allocatable :: cs

        allocate (field)

        ! Open file
        status = nf90_open(trim(filename), nf90_nowrite, ncid)
        call check_nc(status, "open " // trim(filename))

        ! Read dimensions
        call check_nc(nf90_inq_dimid(ncid, "rho", dimid), "inq_dim rho")
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=n_rho), "get_dim rho")
        call check_nc(nf90_inq_dimid(ncid, "theta", dimid), "inq_dim theta")
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=n_theta), "get_dim theta")
        call check_nc(nf90_inq_dimid(ncid, "zeta", dimid), "inq_dim zeta")
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=n_zeta), "get_dim zeta")

        ! Read coordinate arrays
        allocate (rho(n_rho), theta(n_theta), zeta(n_zeta))
        call check_nc(nf90_inq_varid(ncid, "rho", varid), "inq_var rho")
        call check_nc(nf90_get_var(ncid, varid, rho), "get_var rho")
        call check_nc(nf90_inq_varid(ncid, "theta", varid), "inq_var theta")
        call check_nc(nf90_get_var(ncid, varid, theta), "get_var theta")
        call check_nc(nf90_inq_varid(ncid, "zeta", varid), "inq_var zeta")
        call check_nc(nf90_get_var(ncid, varid, zeta), "get_var zeta")

        ! Read scalar attributes and variables
        call check_nc(nf90_get_att(ncid, nf90_global, "torflux", torflux_val), &
                       "get_att torflux")
        call check_nc(nf90_inq_varid(ncid, "num_field_periods", varid), &
                       "inq_var num_field_periods")
        call check_nc(nf90_get_var(ncid, varid, nfp_file), &
                       "get_var num_field_periods")

        ! Read 1D profiles
        allocate (A_phi_arr(n_rho), B_theta_arr(n_rho), B_phi_arr(n_rho))

        call check_nc(nf90_inq_varid(ncid, "A_phi", varid), "inq_var A_phi")
        call check_nc(nf90_get_var(ncid, varid, A_phi_arr), "get_var A_phi")

        call check_nc(nf90_inq_varid(ncid, "B_theta", varid), "inq_var B_theta")
        call check_nc(nf90_get_var(ncid, varid, B_theta_arr), "get_var B_theta")

        call check_nc(nf90_inq_varid(ncid, "B_phi", varid), "inq_var B_phi")
        call check_nc(nf90_get_var(ncid, varid, B_phi_arr), "get_var B_phi")

        ! Read 3D Bmod (stored as zeta, theta, rho in NetCDF → read into rho, theta, zeta)
        allocate (Bmod_arr(n_rho, n_theta, n_zeta))
        call check_nc(nf90_inq_varid(ncid, "Bmod", varid), "inq_var Bmod")
        call read_3d_reordered(ncid, varid, n_rho, n_theta, n_zeta, Bmod_arr)

        call check_nc(nf90_close(ncid), "close")

        ! Grid parameters
        rho_min = rho(1)
        rho_max = rho(n_rho)
        h_theta_val = theta(2) - theta(1)
        h_phi_val = zeta(2) - zeta(1)

        ! Build A_phi spline over s = rho^2 (matching VMEC/boozer_converter convention)
        allocate (s_grid(n_rho))
        do i = 1, n_rho
            s_grid(i) = rho(i)**2
        end do
        s_min = s_grid(1)
        s_max = s_grid(n_rho)
        h_s = (s_max - s_min) / real(n_rho - 1, dp)

        allocate (y_aphi(n_rho, 1))
        y_aphi(:, 1) = A_phi_arr
        call construct_batch_splines_1d(s_min, s_max, y_aphi, spline_order_1d, &
                                        .false., field%aphi_spline)

        ! Build B_theta, B_phi spline over rho_tor
        allocate (y_bcovar(n_rho, 2))
        y_bcovar(:, 1) = B_theta_arr
        y_bcovar(:, 2) = B_phi_arr
        call construct_batch_splines_1d(rho_min, rho_max, y_bcovar, spline_order_1d, &
                                        .false., field%bcovar_spline)

        ! Build Bmod 3D spline over (rho_tor, theta_B, phi_B)
        x_min_3d = [rho_min, theta(1), zeta(1)]
        x_max_3d = [rho_max, theta(n_theta), zeta(n_zeta)]

        allocate (y_bmod(n_rho, n_theta, n_zeta, 1))
        y_bmod(:, :, :, 1) = Bmod_arr
        call construct_batch_splines_3d(x_min_3d, x_max_3d, y_bmod, &
                                        spline_order_3d, periodic_3d, &
                                        field%bmod_spline)

        ! Store metadata
        field%torflux = torflux_val
        field%nfp = nfp_file
        field%ns = n_rho
        field%ntheta = n_theta
        field%nphi = n_zeta
        field%hs = h_s
        field%h_theta = h_theta_val
        field%h_phi = h_phi_val
        field%filename = filename

        ! Set up coordinate system from the same chartmap file
        call make_chartmap_coordinate_system(cs, filename)
        allocate (field%coords, source=cs)

        field%initialized = .true.

        print *, 'Loaded Boozer chartmap field from ', trim(filename)
        print *, '  nfp=', nfp_file, ' ns=', n_rho, ' ntheta=', n_theta, &
                 ' nphi=', n_zeta
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

        ! A_phi from 1D spline at s
        call evaluate_batch_splines_1d_der2(self%aphi_spline, s, y1d_aphi, &
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

    subroutine read_3d_reordered(ncid, varid, n1, n2, n3, arr)
        !> Read a NetCDF variable with dims (zeta, theta, rho) into Fortran
        !> array (n1=rho, n2=theta, n3=zeta). NF90 reverses dimension order
        !> so the data lands directly in (rho, theta, zeta) layout.
        integer, intent(in) :: ncid, varid, n1, n2, n3
        real(dp), intent(out) :: arr(n1, n2, n3)

        call check_nc(nf90_get_var(ncid, varid, arr), "get_var 3D")
    end subroutine read_3d_reordered

    subroutine check_nc(status, location)
        integer, intent(in) :: status
        character(len=*), intent(in) :: location

        if (status /= nf90_noerr) then
            print *, "NetCDF error at ", trim(location), ": ", &
                trim(nf90_strerror(status))
            error stop "NetCDF operation failed"
        end if
    end subroutine check_nc

end module field_boozer_chartmap
