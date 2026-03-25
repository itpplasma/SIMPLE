module gvec_export_data
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                           construct_batch_splines_1d, construct_batch_splines_3d, &
                           evaluate_batch_splines_1d_der2, &
                           evaluate_batch_splines_3d_der2
    use nctools_module, only: nc_close, nc_get
    use netcdf, only: NF90_CHAR, NF90_GLOBAL, NF90_NOERR, NF90_NOWRITE, &
                      nf90_get_att, nf90_inquire_attribute, nf90_open

    implicit none

    private

    integer, parameter, public :: gvec_profile_a_theta = 1
    integer, parameter, public :: gvec_profile_a_phi = 2
    integer, parameter, public :: gvec_profile_da_theta_ds = 3
    integer, parameter, public :: gvec_profile_da_phi_ds = 4
    integer, parameter, public :: gvec_profile_iota = 5

    integer, parameter, public :: gvec_geom_r = 1
    integer, parameter, public :: gvec_geom_z = 2
    integer, parameter, public :: gvec_geom_dr_ds = 3
    integer, parameter, public :: gvec_geom_dr_dt = 4
    integer, parameter, public :: gvec_geom_dr_dp = 5
    integer, parameter, public :: gvec_geom_dz_ds = 6
    integer, parameter, public :: gvec_geom_dz_dt = 7
    integer, parameter, public :: gvec_geom_dz_dp = 8
    integer, parameter, public :: gvec_geom_lambda = 9
    integer, parameter, public :: gvec_geom_dlambda_ds = 10
    integer, parameter, public :: gvec_geom_dlambda_dt = 11
    integer, parameter, public :: gvec_geom_dlambda_dp = 12
    integer, parameter, public :: gvec_geom_sqg = 13
    integer, parameter, public :: gvec_geom_bctr_vartheta = 14
    integer, parameter, public :: gvec_geom_bctr_varphi = 15
    integer, parameter, public :: gvec_geom_bcov_s = 16
    integer, parameter, public :: gvec_geom_bcov_vartheta = 17
    integer, parameter, public :: gvec_geom_bcov_varphi = 18

    character(len=*), parameter :: export_file_type = 'gvec_export'
    character(len=*), parameter :: export_version = '1'
    integer, parameter :: num_geom_splines = 3

    type, public :: gvec_export_data_t
        type(BatchSplineData1D) :: profile_splines
        type(BatchSplineData3D) :: geom_splines(num_geom_splines)
        integer :: nfp = 1
        logical :: initialized = .false.
    contains
        procedure :: evaluate_profiles => gvec_export_evaluate_profiles
        procedure :: evaluate_geometry => gvec_export_evaluate_geometry
        procedure :: phi_period => gvec_export_phi_period
    end type gvec_export_data_t

    public :: gvec_export_is_file
    public :: load_gvec_export_data

contains

    logical function gvec_export_is_file(filename)
        character(*), intent(in) :: filename

        integer :: ncid
        integer :: att_type
        integer :: att_len
        character(len=:), allocatable :: value

        gvec_export_is_file = .false.

        if (len_trim(filename) < 4) return
        if (filename(len_trim(filename) - 2:len_trim(filename)) /= '.nc') return

        if (nf90_open(trim(filename), NF90_NOWRITE, ncid) /= NF90_NOERR) return

        if (nf90_inquire_attribute(ncid, NF90_GLOBAL, 'simple_file_type', &
                                   xtype=att_type, len=att_len) == NF90_NOERR .and. &
            att_type == NF90_CHAR .and. att_len > 0) then
            allocate (character(len=att_len) :: value)
            if (nf90_get_att(ncid, NF90_GLOBAL, 'simple_file_type', value) == &
                NF90_NOERR) then
                gvec_export_is_file = trim(value) == export_file_type
            end if
        end if

        call nc_close(ncid)
    end function gvec_export_is_file

    subroutine load_gvec_export_data(filename, data)
        character(*), intent(in) :: filename
        type(gvec_export_data_t), intent(out) :: data

        integer, parameter :: order_1d = 5
        integer, parameter :: order_3d(3) = [5, 5, 5]
        logical, parameter :: periodic_3d(3) = [.false., .true., .true.]

        integer :: ncid
        integer :: version_att_type
        integer :: version_att_len
        integer :: status
        character(len=:), allocatable :: version_value
        real(dp), allocatable :: s(:), theta(:), varphi(:)
        real(dp), allocatable :: profile_batch(:, :)
        real(dp), allocatable :: geom_batch(:, :, :, :)

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            print *, 'load_gvec_export_data: failed to open ', trim(filename)
            error stop
        end if

        if (.not. gvec_export_is_file(filename)) then
            print *, 'load_gvec_export_data: not a SIMPLE GVEC export file: ', &
                trim(filename)
            error stop
        end if

        status = nf90_inquire_attribute(ncid, NF90_GLOBAL, &
                                        'simple_gvec_export_version', &
                                        xtype=version_att_type, len=version_att_len)
        if (status /= NF90_NOERR .or. version_att_type /= NF90_CHAR .or. &
            version_att_len < 1) then
            print *, 'load_gvec_export_data: missing export version attribute in ', &
                trim(filename)
            error stop
        end if

        allocate (character(len=version_att_len) :: version_value)
        status = nf90_get_att(ncid, NF90_GLOBAL, 'simple_gvec_export_version', &
                              version_value)
        if (status /= NF90_NOERR) then
            print *, 'load_gvec_export_data: failed to read export version in ', &
                trim(filename)
            error stop
        end if
        if (trim(version_value) /= export_version) then
            print *, 'load_gvec_export_data: unsupported export version ', &
                trim(version_value)
            error stop
        end if

        call read_coords(ncid, s, theta, varphi)
        call read_profiles(ncid, size(s), profile_batch)
        call read_geometry(ncid, size(s), size(theta), size(varphi), geom_batch)
        status = nf90_get_att(ncid, NF90_GLOBAL, 'nfp', data%nfp)
        if (status /= NF90_NOERR) then
            print *, 'load_gvec_export_data: missing nfp attribute in ', trim(filename)
            error stop
        end if
        call nc_close(ncid)

        call construct_batch_splines_1d(s(1), s(size(s)), profile_batch, order_1d, &
                                        .false., data%profile_splines)
        call construct_batch_splines_3d([s(1), theta(1), varphi(1)], &
                                        [s(size(s)), theta(size(theta)), &
                                         varphi(size(varphi))], &
                                        geom_batch(:, :, :, 1:8), order_3d, &
                                        periodic_3d, data%geom_splines(1))
        call construct_batch_splines_3d([s(1), theta(1), varphi(1)], &
                                        [s(size(s)), theta(size(theta)), &
                                         varphi(size(varphi))], &
                                        geom_batch(:, :, :, 9:16), order_3d, &
                                        periodic_3d, data%geom_splines(2))
        call construct_batch_splines_3d([s(1), theta(1), varphi(1)], &
                                        [s(size(s)), theta(size(theta)), &
                                         varphi(size(varphi))], &
                                        geom_batch(:, :, :, 17:18), order_3d, &
                                        periodic_3d, data%geom_splines(3))

        data%initialized = .true.
    end subroutine load_gvec_export_data

    subroutine gvec_export_evaluate_profiles(self, s, values, dvalues, d2values)
        class(gvec_export_data_t), intent(in) :: self
        real(dp), intent(in) :: s
        real(dp), intent(out) :: values(:)
        real(dp), intent(out), optional :: dvalues(:)
        real(dp), intent(out), optional :: d2values(:)

        real(dp) :: local_d1(size(values))
        real(dp) :: local_d2(size(values))

        call evaluate_batch_splines_1d_der2(self%profile_splines, s, values, &
                                            local_d1, local_d2)

        if (present(dvalues)) dvalues = local_d1
        if (present(d2values)) d2values = local_d2
    end subroutine gvec_export_evaluate_profiles

    subroutine gvec_export_evaluate_geometry(self, x, values, dvalues, d2values)
        class(gvec_export_data_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: values(:)
        real(dp), intent(out), optional :: dvalues(:, :)
        real(dp), intent(out), optional :: d2values(:, :)

        real(dp) :: values_1(8)
        real(dp) :: values_2(8)
        real(dp) :: values_3(2)
        real(dp) :: dvalues_1(3, 8)
        real(dp) :: dvalues_2(3, 8)
        real(dp) :: dvalues_3(3, 2)
        real(dp) :: d2values_1(6, 8)
        real(dp) :: d2values_2(6, 8)
        real(dp) :: d2values_3(6, 2)

        call evaluate_batch_splines_3d_der2(self%geom_splines(1), x, values_1, &
                                            dvalues_1, d2values_1)
        call evaluate_batch_splines_3d_der2(self%geom_splines(2), x, values_2, &
                                            dvalues_2, d2values_2)
        call evaluate_batch_splines_3d_der2(self%geom_splines(3), x, values_3, &
                                            dvalues_3, d2values_3)

        values(1:8) = values_1
        values(9:16) = values_2
        values(17:18) = values_3

        values(gvec_geom_dr_ds) = dvalues_1(1, 1)
        values(gvec_geom_dr_dt) = dvalues_1(2, 1)
        values(gvec_geom_dr_dp) = dvalues_1(3, 1)
        values(gvec_geom_dz_ds) = dvalues_1(1, 2)
        values(gvec_geom_dz_dt) = dvalues_1(2, 2)
        values(gvec_geom_dz_dp) = dvalues_1(3, 2)
        values(gvec_geom_dlambda_ds) = dvalues_2(1, 1)
        values(gvec_geom_dlambda_dt) = dvalues_2(2, 1)
        values(gvec_geom_dlambda_dp) = dvalues_2(3, 1)

        if (present(dvalues)) then
            dvalues(:, 1:8) = dvalues_1
            dvalues(:, 9:16) = dvalues_2
            dvalues(:, 17:18) = dvalues_3
        end if
        if (present(d2values)) then
            d2values(:, 1:8) = d2values_1
            d2values(:, 9:16) = d2values_2
            d2values(:, 17:18) = d2values_3
        end if
    end subroutine gvec_export_evaluate_geometry

    real(dp) function gvec_export_phi_period(self)
        class(gvec_export_data_t), intent(in) :: self

        gvec_export_phi_period = 2.0_dp * acos(-1.0_dp) / real(self%nfp, dp)
    end function gvec_export_phi_period

    subroutine read_coords(ncid, s, theta, varphi)
        integer, intent(in) :: ncid
        real(dp), allocatable, intent(out) :: s(:)
        real(dp), allocatable, intent(out) :: theta(:)
        real(dp), allocatable, intent(out) :: varphi(:)

        integer :: ns
        integer :: ntheta
        integer :: nvarphi

        call read_vector_1d(ncid, 's', ns, s)
        call read_vector_1d(ncid, 'theta', ntheta, theta)
        call read_vector_1d(ncid, 'varphi', nvarphi, varphi)
    end subroutine read_coords

    subroutine read_profiles(ncid, ns, profile_batch)
        integer, intent(in) :: ncid
        integer, intent(in) :: ns
        real(dp), allocatable, intent(out) :: profile_batch(:, :)

        real(dp), allocatable :: temp(:)

        allocate (profile_batch(ns, 5))

        allocate (temp(ns))
        call nc_get(ncid, 'A_theta', temp)
        profile_batch(:, gvec_profile_a_theta) = temp
        call nc_get(ncid, 'A_phi', temp)
        profile_batch(:, gvec_profile_a_phi) = temp
        call nc_get(ncid, 'dA_theta_ds', temp)
        profile_batch(:, gvec_profile_da_theta_ds) = temp
        call nc_get(ncid, 'dA_phi_ds', temp)
        profile_batch(:, gvec_profile_da_phi_ds) = temp
        call nc_get(ncid, 'iota', temp)
        profile_batch(:, gvec_profile_iota) = temp
    end subroutine read_profiles

    subroutine read_geometry(ncid, ns, ntheta, nvarphi, geom_batch)
        integer, intent(in) :: ncid
        integer, intent(in) :: ns
        integer, intent(in) :: ntheta
        integer, intent(in) :: nvarphi
        real(dp), allocatable, intent(out) :: geom_batch(:, :, :, :)

        real(dp), allocatable :: temp(:, :, :)

        allocate (geom_batch(ns, ntheta, nvarphi, 18))
        allocate (temp(ns, ntheta, nvarphi))

        call nc_get(ncid, 'R', temp)
        geom_batch(:, :, :, gvec_geom_r) = temp
        call nc_get(ncid, 'Z', temp)
        geom_batch(:, :, :, gvec_geom_z) = temp
        call nc_get(ncid, 'dR_ds', temp)
        geom_batch(:, :, :, gvec_geom_dr_ds) = temp
        call nc_get(ncid, 'dR_dt', temp)
        geom_batch(:, :, :, gvec_geom_dr_dt) = temp
        call nc_get(ncid, 'dR_dp', temp)
        geom_batch(:, :, :, gvec_geom_dr_dp) = temp
        call nc_get(ncid, 'dZ_ds', temp)
        geom_batch(:, :, :, gvec_geom_dz_ds) = temp
        call nc_get(ncid, 'dZ_dt', temp)
        geom_batch(:, :, :, gvec_geom_dz_dt) = temp
        call nc_get(ncid, 'dZ_dp', temp)
        geom_batch(:, :, :, gvec_geom_dz_dp) = temp
        call nc_get(ncid, 'Lambda', temp)
        geom_batch(:, :, :, gvec_geom_lambda) = temp
        call nc_get(ncid, 'dLambda_ds', temp)
        geom_batch(:, :, :, gvec_geom_dlambda_ds) = temp
        call nc_get(ncid, 'dLambda_dt', temp)
        geom_batch(:, :, :, gvec_geom_dlambda_dt) = temp
        call nc_get(ncid, 'dLambda_dp', temp)
        geom_batch(:, :, :, gvec_geom_dlambda_dp) = temp
        call nc_get(ncid, 'sqg_symflux', temp)
        geom_batch(:, :, :, gvec_geom_sqg) = temp
        call nc_get(ncid, 'Bctr_vartheta', temp)
        geom_batch(:, :, :, gvec_geom_bctr_vartheta) = temp
        call nc_get(ncid, 'Bctr_varphi', temp)
        geom_batch(:, :, :, gvec_geom_bctr_varphi) = temp
        call nc_get(ncid, 'Bcov_s', temp)
        geom_batch(:, :, :, gvec_geom_bcov_s) = temp
        call nc_get(ncid, 'Bcov_vartheta', temp)
        geom_batch(:, :, :, gvec_geom_bcov_vartheta) = temp
        call nc_get(ncid, 'Bcov_varphi', temp)
        geom_batch(:, :, :, gvec_geom_bcov_varphi) = temp
    end subroutine read_geometry

    subroutine read_vector_1d(ncid, name, n, values)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        integer, intent(out) :: n
        real(dp), allocatable, intent(out) :: values(:)

        integer :: dimid

        call inquire_dim_1d(ncid, name, dimid, n)
        allocate (values(n))
        call nc_get(ncid, name, values)
    end subroutine read_vector_1d

    subroutine inquire_dim_1d(ncid, name, dimid, n)
        use netcdf, only: nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension

        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        integer, intent(out) :: dimid
        integer, intent(out) :: n

        integer :: varid
        integer :: dimids(1)
        integer :: status

        status = nf90_inq_varid(ncid, trim(name), varid)
        if (status /= NF90_NOERR) then
            print *, 'gvec_export_data: missing variable ', trim(name)
            error stop
        end if

        status = nf90_inquire_variable(ncid, varid, dimids=dimids)
        if (status /= NF90_NOERR) then
            print *, 'gvec_export_data: failed to inspect variable ', trim(name)
            error stop
        end if

        dimid = dimids(1)
        status = nf90_inquire_dimension(ncid, dimid, len=n)
        if (status /= NF90_NOERR) then
            print *, 'gvec_export_data: failed to inspect dimension for ', trim(name)
            error stop
        end if
    end subroutine inquire_dim_1d

end module gvec_export_data
