program test_chartmap_metadata
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use chartmap_metadata, only: chartmap_metadata_t, read_chartmap_metadata
    use test_utils, only: check_close, check_string
    implicit none

    type(chartmap_metadata_t) :: meta
    integer :: errors

    errors = 0
    call write_min_chartmap("mini.chartmap.nc", "cm", 0.8_dp)
    call read_chartmap_metadata("mini.chartmap.nc", meta)
    call check_close(meta%rho_lcfs, 0.8_dp, 1.0e-15_dp, "rho_lcfs mismatch", &
        errors)
    call check_close(meta%cart_scale_to_m, 0.01_dp, 1.0e-15_dp, &
        "scale mismatch", errors)
    call check_string(meta%cart_units, "cm", "units mismatch", errors)
    if (errors /= 0) error stop 1
contains

    subroutine check_nc(status, where)
        integer, intent(in) :: status
        character(len=*), intent(in) :: where

        if (status /= nf90_noerr) then
            print *, "NetCDF error at ", trim(where), ": ", trim(nf90_strerror(status))
            error stop 1
        end if
    end subroutine check_nc

    subroutine write_min_chartmap(path, x_units, rho_lcfs)
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: x_units
        real(dp), intent(in) :: rho_lcfs

        integer :: ncid, status
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x
        integer :: dims(3)
        real(dp) :: rho(2), theta(2), zeta(1)
        real(dp) :: x(2, 2, 1)

        rho = [0.0_dp, 1.0_dp]
        theta = [0.0_dp, 3.141592653589793_dp]
        zeta = [0.0_dp]
        x = 0.0_dp

        status = nf90_create(path, nf90_netcdf4, ncid)
        call check_nc(status, "create")

        status = nf90_def_dim(ncid, "rho", 2, dim_rho)
        call check_nc(status, "def_dim rho")
        status = nf90_def_dim(ncid, "theta", 2, dim_theta)
        call check_nc(status, "def_dim theta")
        status = nf90_def_dim(ncid, "zeta", 1, dim_zeta)
        call check_nc(status, "def_dim zeta")

        status = nf90_def_var(ncid, "rho", nf90_double, [dim_rho], var_rho)
        call check_nc(status, "def_var rho")
        status = nf90_def_var(ncid, "theta", nf90_double, [dim_theta], var_theta)
        call check_nc(status, "def_var theta")
        status = nf90_def_var(ncid, "zeta", nf90_double, [dim_zeta], var_zeta)
        call check_nc(status, "def_var zeta")

        dims = [dim_rho, dim_theta, dim_zeta]
        status = nf90_def_var(ncid, "x", nf90_double, dims, var_x)
        call check_nc(status, "def_var x")
        status = nf90_put_att(ncid, var_x, "units", trim(x_units))
        call check_nc(status, "put_att x.units")

        status = nf90_put_att(ncid, nf90_global, "rho_lcfs", rho_lcfs)
        call check_nc(status, "put_att rho_lcfs")

        status = nf90_enddef(ncid)
        call check_nc(status, "enddef")

        status = nf90_put_var(ncid, var_rho, rho)
        call check_nc(status, "put_var rho")
        status = nf90_put_var(ncid, var_theta, theta)
        call check_nc(status, "put_var theta")
        status = nf90_put_var(ncid, var_zeta, zeta)
        call check_nc(status, "put_var zeta")
        status = nf90_put_var(ncid, var_x, x)
        call check_nc(status, "put_var x")

        status = nf90_close(ncid)
        call check_nc(status, "close")
    end subroutine write_min_chartmap

end program test_chartmap_metadata
