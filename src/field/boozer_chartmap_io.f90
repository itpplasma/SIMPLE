module boozer_chartmap_io
    !> Single reader for extended Boozer chartmap NetCDF files.
    !>
    !> Returns the raw, base-unit field arrays and grid metadata. Both the
    !> object-based field (boozer_chartmap_field_t) and the module-level Boozer
    !> batch splines (load_boozer_from_chartmap) build on this one parse so the
    !> two paths cannot drift apart on grid layout, periodicity, or metadata.
    !>
    !> Bmod is read on the endpoint-included field grid (theta_field/zeta_field),
    !> which spans the full 2*pi (poloidal) and 2*pi/nfp (toroidal) period. The
    !> geometry grid (theta/zeta) is endpoint-excluded and is used here to
    !> recover the angular step sizes and the major radius.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf

    implicit none

    private
    public :: boozer_chartmap_data_t, read_boozer_chartmap

    type :: boozer_chartmap_data_t
        integer :: n_rho = 0
        integer :: n_theta = 0       !< field grid, endpoint-included (period span)
        integer :: n_phi = 0         !< field grid, endpoint-included (period span)
        integer :: nfp = 1
        real(dp) :: torflux = 0.0_dp
        real(dp) :: rmajor = 0.0_dp  !< metres, derived from innermost-surface geometry
        real(dp) :: rho_min = 0.0_dp
        real(dp) :: rho_max = 0.0_dp
        real(dp) :: h_s = 0.0_dp     !< uniform rho step
        real(dp) :: h_theta = 0.0_dp !< 2*pi/(n_theta-1)
        real(dp) :: h_phi = 0.0_dp   !< 2*pi/nfp/(n_phi-1)
        real(dp), allocatable :: rho(:)
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

        status = nf90_open(trim(filename), nf90_nowrite, ncid)
        call check(status, "open " // trim(filename))

        ! Geometry grid (endpoint-excluded): used only for step sizes.
        call check(nf90_inq_dimid(ncid, "rho", dimid), "inq_dim rho")
        call check(nf90_inquire_dimension(ncid, dimid, len=n_rho), "len rho")
        call check(nf90_inq_dimid(ncid, "theta", dimid), "inq_dim theta")
        call check(nf90_inquire_dimension(ncid, dimid, len=n_theta_geom), "len theta")
        call check(nf90_inq_dimid(ncid, "zeta", dimid), "inq_dim zeta")
        call check(nf90_inquire_dimension(ncid, dimid, len=n_phi_geom), "len zeta")

        d%n_rho = n_rho
        allocate (d%rho(n_rho), theta_geom(n_theta_geom), zeta_geom(n_phi_geom))
        call check(nf90_inq_varid(ncid, "rho", varid), "inq_var rho")
        call check(nf90_get_var(ncid, varid, d%rho), "get rho")
        call check(nf90_inq_varid(ncid, "theta", varid), "inq_var theta")
        call check(nf90_get_var(ncid, varid, theta_geom), "get theta")
        call check(nf90_inq_varid(ncid, "zeta", varid), "inq_var zeta")
        call check(nf90_get_var(ncid, varid, zeta_geom), "get zeta")

        d%rho_min = d%rho(1)
        d%rho_max = d%rho(n_rho)
        d%h_s = (d%rho_max - d%rho_min) / real(n_rho - 1, dp)
        d%h_theta = theta_geom(2) - theta_geom(1)
        d%h_phi = zeta_geom(2) - zeta_geom(1)

        ! Scalars.
        call check(nf90_get_att(ncid, nf90_global, "torflux", d%torflux), &
                    "att torflux")
        call check(nf90_inq_varid(ncid, "num_field_periods", varid), &
                    "inq_var num_field_periods")
        call check(nf90_get_var(ncid, varid, d%nfp), "get num_field_periods")

        ! Major radius from geometry: (theta,zeta)-average of sqrt(x^2+y^2) on
        ! the innermost rho surface (the chartmap analogue of libneo's
        ! rmajor = rmnc(1,0) axis fallback in vmecinm_m.f90). x,y are stored in
        ! cm at base scale; rmajor is kept in metres like the former global
        ! attribute, so the vmec_RZ_scale applied downstream by the field
        ! loaders stays correct. A leftover "rmajor" attribute in older files
        ! is ignored.
        call derive_rmajor(ncid, n_theta_geom, n_phi_geom, d%rmajor)

        ! 1D profiles on the rho grid.
        allocate (d%A_phi(n_rho), d%B_theta(n_rho), d%B_phi(n_rho))
        call check(nf90_inq_varid(ncid, "A_phi", varid), "inq_var A_phi")
        call check(nf90_get_var(ncid, varid, d%A_phi), "get A_phi")
        call check(nf90_inq_varid(ncid, "B_theta", varid), "inq_var B_theta")
        call check(nf90_get_var(ncid, varid, d%B_theta), "get B_theta")
        call check(nf90_inq_varid(ncid, "B_phi", varid), "inq_var B_phi")
        call check(nf90_get_var(ncid, varid, d%B_phi), "get B_phi")

        ! Bmod lives on the endpoint-included field grid.
        call check(nf90_inq_dimid(ncid, "theta_field", dimid), "inq_dim theta_field")
        call check(nf90_inquire_dimension(ncid, dimid, len=d%n_theta), &
                    "len theta_field")
        call check(nf90_inq_dimid(ncid, "zeta_field", dimid), "inq_dim zeta_field")
        call check(nf90_inquire_dimension(ncid, dimid, len=d%n_phi), &
                    "len zeta_field")

        allocate (d%Bmod(n_rho, d%n_theta, d%n_phi))
        call check(nf90_inq_varid(ncid, "Bmod", varid), "inq_var Bmod")
        call check(nf90_get_var(ncid, varid, d%Bmod), "get Bmod")

        call check(nf90_close(ncid), "close")
    end subroutine read_boozer_chartmap

    subroutine derive_rmajor(ncid, n_theta_geom, n_phi_geom, rmajor)
        integer, intent(in) :: ncid, n_theta_geom, n_phi_geom
        real(dp), intent(out) :: rmajor

        integer :: varid
        real(dp), allocatable :: x_in(:, :, :), y_in(:, :, :)
        real(dp), parameter :: cm_to_m = 1.0e-2_dp

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
