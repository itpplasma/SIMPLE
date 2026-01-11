module chartmap_metadata

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    implicit none

    private
    public :: chartmap_metadata_t
    public :: read_chartmap_metadata

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

        integer :: ncid, status, varid
        character(len=128) :: units

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

        status = nf90_inq_varid(ncid, "x", varid)
        call check_nc(status, "inq_varid x")

        units = ""
        status = nf90_get_att(ncid, varid, "units", units)
        if (status /= nf90_noerr) then
            units = ""
        end if

        units = adjustl(units)
        select case (trim(units))
        case ("cm")
            meta%cart_scale_to_m = 0.01_dp
            meta%cart_units = "cm"
        case ("m")
            meta%cart_scale_to_m = 1.0_dp
            meta%cart_units = "m"
        case ("")
            meta%cart_scale_to_m = 0.01_dp
            meta%cart_units = "cm"
        case default
            print *, "chartmap_metadata: unsupported x.units=", trim(units)
            status = nf90_close(ncid)
            call check_nc(status, "nf90_close(chartmap)")
            error stop "unsupported chartmap Cartesian units"
        end select

        status = nf90_close(ncid)
        call check_nc(status, "nf90_close(chartmap)")
    end subroutine read_chartmap_metadata

end module chartmap_metadata
