module chartmap_metadata

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    implicit none

    private
    public :: chartmap_metadata_t
    public :: read_chartmap_metadata, read_chartmap_cart_units

    type :: chartmap_metadata_t
        real(dp) :: rho_lcfs = -1.0_dp
        real(dp) :: cart_scale_to_m = -1.0_dp
        character(len=16) :: cart_units = ""
    end type chartmap_metadata_t

contains

    subroutine check_nc(status, location)
        integer, intent(in) :: status
        character(len=*), intent(in) :: location

        if (status /= nf90_noerr) then
            print *, "NetCDF error at ", trim(location), ": ", &
                trim(nf90_strerror(status))
            error stop "NetCDF operation failed"
        end if
    end subroutine check_nc

    subroutine read_chartmap_metadata(filename, meta)
        character(len=*), intent(in) :: filename
        type(chartmap_metadata_t), intent(out) :: meta

        integer :: ncid, status

        meta%rho_lcfs = -1.0_dp
        meta%cart_scale_to_m = -1.0_dp
        meta%cart_units = ""

        status = nf90_open(trim(filename), nf90_nowrite, ncid)
        call check_nc(status, "nf90_open(chartmap)")

        status = nf90_get_att(ncid, nf90_global, "rho_lcfs", meta%rho_lcfs)
        if (status /= nf90_noerr) then
            status = nf90_close(ncid)
            call check_nc(status, "nf90_close(chartmap)")
            error stop "chartmap missing global attribute rho_lcfs"
        end if

        call read_cart_units(ncid, meta%cart_units, meta%cart_scale_to_m)

        status = nf90_close(ncid)
        call check_nc(status, "nf90_close(chartmap)")
    end subroutine read_chartmap_metadata

    subroutine read_chartmap_cart_units(filename, cart_units)
        character(len=*), intent(in) :: filename
        character(len=*), intent(out) :: cart_units

        integer :: ncid, status
        real(dp) :: cart_scale_to_m

        status = nf90_open(trim(filename), nf90_nowrite, ncid)
        call check_nc(status, "nf90_open(chartmap)")
        call read_cart_units(ncid, cart_units, cart_scale_to_m)
        status = nf90_close(ncid)
        call check_nc(status, "nf90_close(chartmap)")
    end subroutine read_chartmap_cart_units

    subroutine read_cart_units(ncid, cart_units, cart_scale_to_m)
        integer, intent(in) :: ncid
        character(len=*), intent(out) :: cart_units
        real(dp), intent(out) :: cart_scale_to_m

        integer :: status, varid
        character(len=128) :: units

        status = nf90_inq_varid(ncid, "x", varid)
        call check_nc(status, "inq_varid x")
        units = ""
        status = nf90_get_att(ncid, varid, "units", units)
        if (status /= nf90_noerr) units = ""

        select case (trim(adjustl(units)))
        case ("cm", "")
            cart_scale_to_m = 0.01_dp
            cart_units = "cm"
        case ("m")
            cart_scale_to_m = 1.0_dp
            cart_units = "m"
        case default
            print *, "chartmap_metadata: unsupported x.units=", trim(units)
            error stop "unsupported chartmap Cartesian units"
        end select
    end subroutine read_cart_units

end module chartmap_metadata
