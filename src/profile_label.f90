module profile_label
    !> Tabulated profile label s(R, phi, Z) for the collision profiles.
    !>
    !> Reference coordinates without toroidal-flux meaning (geometric
    !> chartmaps, rho_convention = unknown) cannot give the radius of the
    !> &profiles. A profile label file supplies that radius as a function of
    !> the physical position instead. It is a netCDF file with
    !>
    !>   dimensions  R, Z, phi
    !>   double R(R), Z(Z)       uniform increasing grids, attribute units
    !>                           "m" or "cm"
    !>   double phi(phi)         phi(i) = (i - 1)*2*pi/(nfp*size(phi)), rad
    !>   double label(phi, Z, R) finite profile label (netCDF/C order)
    !>   global attribute nfp    number of field periods (integer >= 1)
    !>
    !> The label is interpolated trilinearly, periodically in phi modulo
    !> 2*pi/nfp. A position outside the R-Z box is an error: the table must
    !> cover the whole tracing domain.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: load_profile_label, profile_label_at, profile_label_loaded

    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

    !> Grid origin and spacing, in cm and rad.
    real(dp), save :: r0 = 0.0_dp, dr = 0.0_dp, z0 = 0.0_dp, dz = 0.0_dp
    real(dp), save :: period = 0.0_dp, dphi = 0.0_dp
    integer, save :: nr = 0, nz = 0, nphi = 0
    !> label(i_R, i_Z, i_phi)
    real(dp), allocatable, save :: table(:, :, :)

contains

    logical function profile_label_loaded()
        profile_label_loaded = allocated(table)
    end function profile_label_loaded

    subroutine load_profile_label(filename)
        use netcdf
        character(*), intent(in) :: filename
        real(dp), allocatable :: rgrid(:), zgrid(:), phigrid(:)
        real(dp) :: scale_r, scale_z
        integer :: ncid, varid, nfp

        call check(filename, nf90_open(trim(filename), nf90_nowrite, ncid), 'open')
        nr = dim_len(filename, ncid, 'R')
        nz = dim_len(filename, ncid, 'Z')
        nphi = dim_len(filename, ncid, 'phi')
        if (nr < 2 .or. nz < 2 .or. nphi < 1) &
            call fail(filename, 'needs at least 2 R, 2 Z and 1 phi points')
        allocate (rgrid(nr), zgrid(nz), phigrid(nphi))
        if (allocated(table)) deallocate (table)
        allocate (table(nr, nz, nphi))
        call check(filename, nf90_inq_varid(ncid, 'R', varid), 'R')
        call check(filename, nf90_get_var(ncid, varid, rgrid), 'R')
        scale_r = length_scale(filename, ncid, varid, 'R')
        call check(filename, nf90_inq_varid(ncid, 'Z', varid), 'Z')
        call check(filename, nf90_get_var(ncid, varid, zgrid), 'Z')
        scale_z = length_scale(filename, ncid, varid, 'Z')
        call check(filename, nf90_inq_varid(ncid, 'phi', varid), 'phi')
        call check(filename, nf90_get_var(ncid, varid, phigrid), 'phi')
        call check(filename, nf90_inq_varid(ncid, 'label', varid), 'label')
        call check(filename, nf90_get_var(ncid, varid, table), 'label')
        call check(filename, nf90_get_att(ncid, nf90_global, 'nfp', nfp), 'nfp')
        call check(filename, nf90_close(ncid), 'close')

        if (nfp < 1) call fail(filename, 'nfp must be >= 1')
        rgrid = scale_r*rgrid
        zgrid = scale_z*zgrid
        r0 = rgrid(1)
        dr = (rgrid(nr) - rgrid(1))/real(nr - 1, dp)
        z0 = zgrid(1)
        dz = (zgrid(nz) - zgrid(1))/real(nz - 1, dp)
        period = twopi/real(nfp, dp)
        dphi = period/real(nphi, dp)
        if (.not. (dr > 0.0_dp .and. dz > 0.0_dp)) &
            call fail(filename, 'R and Z must increase')
        call require_uniform(filename, rgrid, r0, dr, 'R')
        call require_uniform(filename, zgrid, z0, dz, 'Z')
        call require_uniform(filename, phigrid, 0.0_dp, dphi, 'phi')
        if (.not. all(abs(table) <= huge(1.0_dp))) &
            call fail(filename, 'label must be finite')
        print *, 'profile_label: ', trim(filename), ' R,Z,phi points ', nr, nz, &
            nphi, ' nfp ', nfp
    end subroutine load_profile_label

    subroutine fail(filename, message)
        character(*), intent(in) :: filename, message

        print *, 'load_profile_label: ', trim(filename), ': ', message
        error stop 'load_profile_label: invalid profile label file'
    end subroutine fail

    subroutine check(filename, status, what)
        use netcdf, only: nf90_noerr, nf90_strerror
        character(*), intent(in) :: filename
        integer, intent(in) :: status
        character(*), intent(in) :: what

        if (status /= nf90_noerr) &
            call fail(filename, what//': '//trim(nf90_strerror(status)))
    end subroutine check

    integer function dim_len(filename, ncid, name)
        use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
        character(*), intent(in) :: filename
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        integer :: dimid

        call check(filename, nf90_inq_dimid(ncid, name, dimid), 'dimension '//name)
        call check(filename, nf90_inquire_dimension(ncid, dimid, len=dim_len), &
                   'dimension '//name)
    end function dim_len

    real(dp) function length_scale(filename, ncid, varid, name)
        !> Factor from the file's length unit to cm.
        use netcdf, only: nf90_get_att
        character(*), intent(in) :: filename
        integer, intent(in) :: ncid, varid
        character(*), intent(in) :: name
        character(len=16) :: units

        units = ''
        call check(filename, nf90_get_att(ncid, varid, 'units', units), &
                   name//' units')
        select case (trim(units))
        case ('m')
            length_scale = 100.0_dp
        case ('cm')
            length_scale = 1.0_dp
        case default
            call fail(filename, name//' units must be "m" or "cm"')
            length_scale = 0.0_dp
        end select
    end function length_scale

    subroutine require_uniform(filename, x, x0, dx, name)
        character(*), intent(in) :: filename
        real(dp), intent(in) :: x(:), x0, dx
        character(*), intent(in) :: name
        integer :: i

        do i = 1, size(x)
            if (abs(x(i) - (x0 + real(i - 1, dp)*dx)) > &
                1.0e-9_dp*max(abs(dx), abs(x0) + abs(dx))) &
                call fail(filename, name//' grid must be uniform as documented')
        end do
    end subroutine require_uniform

    !> Profile label at the cylindrical position x = (R, phi, Z), R and Z in
    !> cm. Must stay free of IEEE-module procedures (see profile_s).
    function profile_label_at(x) result(s)
        real(dp), intent(in) :: x(3)
        real(dp) :: s
        real(dp) :: qr, qz, qp, wr, wz, wp, c0, c1
        integer :: ir, iz, ip, jp

        qr = (x(1) - r0)/dr
        qz = (x(3) - z0)/dz
        if (.not. (qr >= 0.0_dp .and. qr <= real(nr - 1, dp) .and. &
                   qz >= 0.0_dp .and. qz <= real(nz - 1, dp))) then
            print *, 'profile_label_at: R, Z [cm] = ', x(1), x(3), &
                ' outside the profile label table'
            error stop 'profile_label_at: position outside the profile label table'
        end if
        ir = min(int(qr), nr - 2)
        iz = min(int(qz), nz - 2)
        wr = qr - real(ir, dp)
        wz = qz - real(iz, dp)
        qp = modulo(x(2), period)/dphi
        ip = min(int(qp), nphi - 1)
        wp = qp - real(ip, dp)
        jp = ip + 1
        if (jp == nphi) jp = 0
        c0 = bilinear(ir, iz, ip + 1, wr, wz)
        c1 = bilinear(ir, iz, jp + 1, wr, wz)
        s = (1.0_dp - wp)*c0 + wp*c1
    end function profile_label_at

    real(dp) function bilinear(ir, iz, k, wr, wz)
        integer, intent(in) :: ir, iz, k
        real(dp), intent(in) :: wr, wz

        bilinear = (1.0_dp - wz)*((1.0_dp - wr)*table(ir + 1, iz + 1, k) &
                                  + wr*table(ir + 2, iz + 1, k)) &
                   + wz*((1.0_dp - wr)*table(ir + 1, iz + 2, k) &
                         + wr*table(ir + 2, iz + 2, k))
    end function bilinear

end module profile_label
