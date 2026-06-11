module boozer_chartmap_io
    !> Single reader for extended Boozer chartmap NetCDF files.
    !>
    !> Returns the raw, base-unit field arrays and grid metadata. Both the
    !> object-based field (boozer_chartmap_field_t) and the module-level Boozer
    !> batch splines (load_boozer_from_chartmap) build on this one parse so the
    !> two paths cannot drift apart on grid layout, periodicity, or metadata.
    !>
    !> A_phi is read on a uniform s grid. B_theta, B_phi, Bmod and geometry
    !> use the rho grid. Bmod is stored on endpoint-excluded theta/zeta grids;
    !> this reader appends exact periodic endpoint planes for the spline backend.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use netcdf

    implicit none

    private
    public :: boozer_chartmap_data_t, read_boozer_chartmap

    type :: boozer_chartmap_data_t
        integer :: n_rho = 0
        integer :: n_s = 0
        integer :: n_theta = 0       !< internal field grid, endpoint-included
        integer :: n_phi = 0         !< internal field grid, endpoint-included
        integer :: nfp = 1
        real(dp) :: torflux = 0.0_dp
        real(dp) :: rmajor = 0.0_dp  !< metres, derived from innermost-surface geometry
        real(dp) :: rho_min = 0.0_dp
        real(dp) :: rho_max = 0.0_dp
        real(dp) :: h_rho = 0.0_dp   !< uniform rho step
        real(dp) :: h_s = 0.0_dp
        real(dp) :: h_theta = 0.0_dp !< 2*pi/(n_theta-1)
        real(dp) :: h_phi = 0.0_dp   !< 2*pi/nfp/(n_phi-1)
        real(dp), allocatable :: rho(:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: A_phi(:)
        real(dp), allocatable :: B_theta(:)
        real(dp), allocatable :: B_phi(:)
        real(dp), allocatable :: Bmod(:, :, :)  !< (n_rho, n_theta, n_phi)
    end type boozer_chartmap_data_t

contains

    subroutine read_boozer_chartmap(filename, d)
        character(len=*), intent(in) :: filename
        type(boozer_chartmap_data_t), intent(out) :: d

        integer :: ncid, status, dimid, varid
        integer :: n_rho, n_theta_geom, n_phi_geom
        real(dp), allocatable :: theta_geom(:), zeta_geom(:)
        real(dp), allocatable :: bmod_file(:, :, :)
        real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

        status = nf90_open(trim(filename), nf90_nowrite, ncid)
        call check(status, "open "//trim(filename))

        ! Geometry grid (endpoint-excluded): used only for step sizes.
        call check(nf90_inq_dimid(ncid, "rho", dimid), "inq_dim rho")
        call check(nf90_inquire_dimension(ncid, dimid, len=n_rho), "len rho")
        call check(nf90_inq_dimid(ncid, "theta", dimid), "inq_dim theta")
        call check(nf90_inquire_dimension(ncid, dimid, len=n_theta_geom), "len theta")
        call check(nf90_inq_dimid(ncid, "zeta", dimid), "inq_dim zeta")
        call check(nf90_inquire_dimension(ncid, dimid, len=n_phi_geom), "len zeta")

        d%n_rho = n_rho
        allocate (d%rho(n_rho), theta_geom(n_theta_geom), zeta_geom(n_phi_geom))
        call require_variable_dimensions(ncid, "rho", [character(len=3) :: "rho"])
        call check(nf90_inq_varid(ncid, "rho", varid), "inq_var rho")
        call check(nf90_get_var(ncid, varid, d%rho), "get rho")
        call require_variable_dimensions(ncid, "theta", [character(len=5) :: "theta"])
        call check(nf90_inq_varid(ncid, "theta", varid), "inq_var theta")
        call check(nf90_get_var(ncid, varid, theta_geom), "get theta")
        call require_variable_dimensions(ncid, "zeta", [character(len=4) :: "zeta"])
        call check(nf90_inq_varid(ncid, "zeta", varid), "inq_var zeta")
        call check(nf90_get_var(ncid, varid, zeta_geom), "get zeta")

        call require_min_points("rho", d%rho)
        call require_min_points("theta", theta_geom)
        call require_min_points("zeta", zeta_geom)
        d%rho_min = d%rho(1)
        d%rho_max = d%rho(n_rho)
        d%h_rho = (d%rho_max - d%rho_min)/real(n_rho - 1, dp)
        call require_uniform_grid("rho", d%rho, d%h_rho)
        d%h_theta = theta_geom(2) - theta_geom(1)
        d%h_phi = zeta_geom(2) - zeta_geom(1)
        call require_uniform_grid("theta", theta_geom, d%h_theta)
        call require_uniform_grid("zeta", zeta_geom, d%h_phi)

        ! Scalars.
        call check(nf90_get_att(ncid, nf90_global, "torflux", d%torflux), &
                   "att torflux")
        call require_scalar_variable(ncid, "num_field_periods")
        call check(nf90_inq_varid(ncid, "num_field_periods", varid), &
                   "inq_var num_field_periods")
        call check(nf90_get_var(ncid, varid, d%nfp), "get num_field_periods")
        call require_positive_nfp(d%nfp)
        call require_endpoint_excluded_grid("theta", theta_geom, d%h_theta, twopi)
        call require_endpoint_excluded_grid("zeta", zeta_geom, d%h_phi, &
                                            twopi/real(d%nfp, dp))

        ! Major radius from geometry: (theta,zeta)-average of sqrt(x^2+y^2) on
        ! the innermost rho surface (the chartmap analogue of libneo's
        ! rmajor = rmnc(1,0) axis fallback in vmecinm_m.f90). x,y are stored in
        ! cm at base scale; rmajor is kept in metres like the former global
        ! attribute, so the vmec_RZ_scale applied downstream by the field
        ! loaders stays correct. A leftover "rmajor" attribute in older files
        ! is ignored.
        call derive_rmajor(ncid, n_theta_geom, n_phi_geom, d%rmajor)

        ! 1D profiles. A_phi has its own abscissa; B_theta/B_phi remain on rho.
        call read_aphi_profile(ncid, d)
        allocate (d%B_theta(n_rho), d%B_phi(n_rho))
        call require_variable_dimensions(ncid, "B_theta", [character(len=3) :: "rho"])
        call check(nf90_inq_varid(ncid, "B_theta", varid), "inq_var B_theta")
        call check(nf90_get_var(ncid, varid, d%B_theta), "get B_theta")
        call require_variable_dimensions(ncid, "B_phi", [character(len=3) :: "rho"])
        call check(nf90_inq_varid(ncid, "B_phi", varid), "inq_var B_phi")
        call check(nf90_get_var(ncid, varid, d%B_phi), "get B_phi")

        ! Bmod is stored on the endpoint-excluded file grid. Append exact
        ! periodic endpoint planes internally because the spline backend uses
        ! a full-period grid.
        d%n_theta = n_theta_geom + 1
        d%n_phi = n_phi_geom + 1
        allocate (bmod_file(n_rho, n_theta_geom, n_phi_geom))
        allocate (d%Bmod(n_rho, d%n_theta, d%n_phi))
        call require_variable_dimensions(ncid, "Bmod", &
                                         [character(len=5) :: "rho", "theta", "zeta"])
        call check(nf90_inq_varid(ncid, "Bmod", varid), "inq_var Bmod")
        call check(nf90_get_var(ncid, varid, bmod_file), "get Bmod")
        d%Bmod(:, 1:n_theta_geom, 1:n_phi_geom) = bmod_file
        d%Bmod(:, d%n_theta, 1:n_phi_geom) = bmod_file(:, 1, :)
        d%Bmod(:, 1:n_theta_geom, d%n_phi) = bmod_file(:, :, 1)
        d%Bmod(:, d%n_theta, d%n_phi) = bmod_file(:, 1, 1)

        call check(nf90_close(ncid), "close")
    end subroutine read_boozer_chartmap

    subroutine derive_rmajor(ncid, n_theta_geom, n_phi_geom, rmajor)
        integer, intent(in) :: ncid, n_theta_geom, n_phi_geom
        real(dp), intent(out) :: rmajor

        integer :: varid
        real(dp), allocatable :: x_in(:, :, :), y_in(:, :, :)
        real(dp), parameter :: cm_to_m = 1.0e-2_dp

        call require_variable_dimensions(ncid, "x", &
                                         [character(len=5) :: "rho", "theta", "zeta"])
        call require_variable_dimensions(ncid, "y", &
                                         [character(len=5) :: "rho", "theta", "zeta"])
        call require_variable_dimensions(ncid, "z", &
                                         [character(len=5) :: "rho", "theta", "zeta"])

        allocate (x_in(1, n_theta_geom, n_phi_geom), &
                  y_in(1, n_theta_geom, n_phi_geom))

        call check(nf90_inq_varid(ncid, "x", varid), "inq_var x")
        call check(nf90_get_var(ncid, varid, x_in, start=[1, 1, 1], &
                                count=[1, n_theta_geom, n_phi_geom]), "get x")
        call check(nf90_inq_varid(ncid, "y", varid), "inq_var y")
        call check(nf90_get_var(ncid, varid, y_in, start=[1, 1, 1], &
                                count=[1, n_theta_geom, n_phi_geom]), "get y")

        rmajor = sum(sqrt(x_in**2 + y_in**2))*cm_to_m &
                 /real(n_theta_geom*n_phi_geom, dp)
    end subroutine derive_rmajor

    subroutine read_aphi_profile(ncid, d)
        integer, intent(in) :: ncid
        type(boozer_chartmap_data_t), intent(inout) :: d

        integer :: varid, ndims, dimids(nf90_max_var_dims), var_s, n_aphi
        integer :: status
        character(len=nf90_max_name) :: dim_name
        character(len=32) :: abscissa

        call check(nf90_inq_varid(ncid, "A_phi", varid), "inq_var A_phi")
        call check(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), &
                   "inquire A_phi")
        if (ndims /= 1) then
            print *, "read_boozer_chartmap: A_phi must be one-dimensional"
            error stop "read_boozer_chartmap failed"
        end if
        call check(nf90_inquire_dimension(ncid, dimids(1), name=dim_name, len=n_aphi), &
                   "A_phi dim")

        abscissa = ""
        status = nf90_get_att(ncid, varid, "radial_abscissa", abscissa)
        if (status /= nf90_noerr) then
            print *, "read_boozer_chartmap: A_phi needs radial_abscissa='s'"
            error stop "read_boozer_chartmap failed"
        end if

        if (trim(abscissa) /= "s") then
            print *, "read_boozer_chartmap: unsupported A_phi radial_abscissa=", &
                trim(abscissa)
            error stop "read_boozer_chartmap failed"
        end if
        if (trim(dim_name) /= "s") then
            print *, "read_boozer_chartmap: A_phi radial_abscissa='s' ", &
                "requires dimension s"
            error stop "read_boozer_chartmap failed"
        end if

        d%n_s = n_aphi
        allocate (d%s(n_aphi), d%A_phi(n_aphi))
        call require_variable_dimensions(ncid, "s", [character(len=1) :: "s"])
        call check(nf90_inq_varid(ncid, "s", var_s), "inq_var s")
        call check(nf90_get_var(ncid, var_s, d%s), "get s")
        call require_min_points("s", d%s)
        d%h_s = (d%s(n_aphi) - d%s(1))/real(n_aphi - 1, dp)
        call require_uniform_grid("s", d%s, d%h_s)
        call require_s_range(d)
        call check(nf90_get_var(ncid, varid, d%A_phi), "get A_phi")
    end subroutine read_aphi_profile

    subroutine require_scalar_variable(ncid, var_name)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: var_name

        integer :: varid, ndims

        call check(nf90_inq_varid(ncid, trim(var_name), varid), &
                   "inq_var "//trim(var_name))
        call check(nf90_inquire_variable(ncid, varid, ndims=ndims), &
                   "inquire "//trim(var_name))
        if (ndims /= 0) then
            print *, "read_boozer_chartmap: ", trim(var_name), " must be scalar"
            error stop "read_boozer_chartmap failed"
        end if
    end subroutine require_scalar_variable

    subroutine require_variable_dimensions(ncid, var_name, expected)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: var_name
        character(len=*), intent(in) :: expected(:)

        integer :: varid, ndims, dimids(nf90_max_var_dims), i
        character(len=nf90_max_name) :: dim_name

        call check(nf90_inq_varid(ncid, trim(var_name), varid), &
                   "inq_var "//trim(var_name))
        call check(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), &
                   "inquire "//trim(var_name))
        if (ndims /= size(expected)) then
            print *, "read_boozer_chartmap: ", trim(var_name), " must have ", &
                size(expected), " dimensions"
            error stop "read_boozer_chartmap failed"
        end if
        do i = 1, size(expected)
            call check(nf90_inquire_dimension(ncid, dimids(i), name=dim_name), &
                       "dimension for "//trim(var_name))
            if (trim(dim_name) /= trim(expected(i))) then
                print *, "read_boozer_chartmap: ", trim(var_name), &
                    " dimension ", i, " is ", trim(dim_name), &
                    " but expected ", trim(expected(i))
                error stop "read_boozer_chartmap failed"
            end if
        end do
    end subroutine require_variable_dimensions

    subroutine require_min_points(name, grid)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: grid(:)

        if (size(grid) < 2) then
            print *, "read_boozer_chartmap: ", trim(name), " needs at least two points"
            error stop "read_boozer_chartmap failed"
        end if
        if (any(.not. ieee_is_finite(grid))) then
            print *, "read_boozer_chartmap: ", trim(name), &
                " grid contains nonfinite values"
            error stop "read_boozer_chartmap failed"
        end if
    end subroutine require_min_points

    subroutine require_uniform_grid(name, grid, h)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: grid(:), h

        integer :: i
        real(dp) :: want, tol

        if (.not. ieee_is_finite(h) .or. h <= 0.0_dp) then
            print *, "read_boozer_chartmap: ", trim(name), " must increase"
            error stop "read_boozer_chartmap failed"
        end if

        tol = 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(grid(size(grid))))
        do i = 1, size(grid)
            want = grid(1) + real(i - 1, dp)*h
            if (abs(grid(i) - want) > tol) then
                print *, "read_boozer_chartmap: nonuniform ", trim(name), &
                    " grid at index ", i
                error stop "read_boozer_chartmap failed"
            end if
        end do
    end subroutine require_uniform_grid

    subroutine require_endpoint_excluded_grid(name, grid, h, period)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: grid(:), h, period

        real(dp) :: tol

        if (.not. ieee_is_finite(period) .or. period <= 0.0_dp) then
            print *, "read_boozer_chartmap: ", trim(name), " period must be positive"
            error stop "read_boozer_chartmap failed"
        end if
        tol = 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(period))
        if (abs(grid(1)) > tol) then
            print *, "read_boozer_chartmap: ", trim(name), " must start at zero"
            error stop "read_boozer_chartmap failed"
        end if
        if (abs(real(size(grid), dp)*h - period) > tol) then
            print *, "read_boozer_chartmap: ", trim(name), &
                " must be endpoint-excluded over one period"
            error stop "read_boozer_chartmap failed"
        end if
    end subroutine require_endpoint_excluded_grid

    subroutine require_positive_nfp(nfp)
        integer, intent(in) :: nfp

        if (nfp <= 0) then
            print *, "read_boozer_chartmap: num_field_periods must be positive"
            error stop "read_boozer_chartmap failed"
        end if
    end subroutine require_positive_nfp

    subroutine require_s_range(d)
        type(boozer_chartmap_data_t), intent(in) :: d

        real(dp) :: tol, s_first, s_last

        s_first = d%rho_min**2
        s_last = d%rho_max**2
        tol = 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(s_last))
        if (abs(d%s(1) - s_first) > tol .or. &
            abs(d%s(d%n_s) - s_last) > tol) then
            print *, "read_boozer_chartmap: s must span rho**2"
            error stop "read_boozer_chartmap failed"
        end if
    end subroutine require_s_range

    subroutine check(status, location)
        integer, intent(in) :: status
        character(len=*), intent(in) :: location

        if (status /= nf90_noerr) then
            print *, "read_boozer_chartmap: NetCDF error at ", trim(location), &
                ": ", trim(nf90_strerror(status))
            error stop "read_boozer_chartmap failed"
        end if
    end subroutine check

end module boozer_chartmap_io
