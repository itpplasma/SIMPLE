module gvec_export_data
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                          construct_batch_splines_1d, construct_batch_splines_3d, &
                          evaluate_batch_splines_1d_der2, evaluate_batch_splines_3d, &
                          evaluate_batch_splines_3d_der2
   use nctools_module, only: nc_close, nc_get
   use netcdf, only: NF90_CHAR, NF90_GLOBAL, NF90_NOERR, NF90_NOWRITE, &
                     nf90_get_att, nf90_inq_varid, nf90_inquire_attribute, &
                     nf90_open

   implicit none

   private

   integer, parameter, public :: gvec_family_logical = 1
   integer, parameter, public :: gvec_family_boozer = 2
   integer, parameter, public :: gvec_family_pest = 3

   integer, parameter, public :: gvec_profile_a_theta = 1
   integer, parameter, public :: gvec_profile_a_phi = 2
   integer, parameter, public :: gvec_profile_da_theta_ds = 3
   integer, parameter, public :: gvec_profile_da_phi_ds = 4
   integer, parameter, public :: gvec_profile_iota = 5

   integer, parameter, public :: gvec_geom_x = 1
   integer, parameter, public :: gvec_geom_y = 2
   integer, parameter, public :: gvec_geom_z = 3
   integer, parameter, public :: gvec_geom_dx_ds = 4
   integer, parameter, public :: gvec_geom_dx_dt = 5
   integer, parameter, public :: gvec_geom_dx_dp = 6
   integer, parameter, public :: gvec_geom_dy_ds = 7
   integer, parameter, public :: gvec_geom_dy_dt = 8
   integer, parameter, public :: gvec_geom_dy_dp = 9
   integer, parameter, public :: gvec_geom_dz_ds = 10
   integer, parameter, public :: gvec_geom_dz_dt = 11
   integer, parameter, public :: gvec_geom_dz_dp = 12
   integer, parameter, public :: gvec_geom_lambda = 13
   integer, parameter, public :: gvec_geom_dlambda_ds = 14
   integer, parameter, public :: gvec_geom_dlambda_dt = 15
   integer, parameter, public :: gvec_geom_dlambda_dp = 16

   integer, parameter, public :: gvec_field_acov_s = 1
   integer, parameter, public :: gvec_field_acov_t = 2
   integer, parameter, public :: gvec_field_acov_p = 3
   integer, parameter, public :: gvec_field_hcov_s = 4
   integer, parameter, public :: gvec_field_hcov_t = 5
   integer, parameter, public :: gvec_field_hcov_p = 6
   integer, parameter, public :: gvec_field_bmod = 7
   integer, parameter, public :: gvec_field_sqgbctr_t = 8
   integer, parameter, public :: gvec_field_sqgbctr_p = 9

   character(len=*), parameter :: export_file_type = 'gvec_export'
   character(len=*), parameter :: export_version = '1'
   integer, parameter :: num_geom_splines = 2
   integer, parameter :: num_field_splines = 2

   type, public :: gvec_export_data_t
      type(BatchSplineData1D) :: profile_splines
      type(BatchSplineData3D) :: geom_splines(num_geom_splines)
      type(BatchSplineData3D) :: field_splines(num_field_splines)
      integer :: nfp = 1
      integer :: family = gvec_family_logical
      logical :: initialized = .false.
      logical :: has_lambda = .false.
   contains
      procedure :: evaluate_profiles => gvec_export_evaluate_profiles
      procedure :: evaluate_geometry => gvec_export_evaluate_geometry
      procedure :: evaluate_field => gvec_export_evaluate_field
      procedure :: phi_period => gvec_export_phi_period
      procedure :: uses_vmec_adapter => gvec_export_uses_vmec_adapter
      procedure :: same_family => gvec_export_same_family
      procedure :: min_s => gvec_export_min_s
   end type gvec_export_data_t

   public :: gvec_export_is_file
   public :: load_gvec_export_data

contains

   subroutine periodic_bounds(s, theta, varphi, nfp, lower, upper)
      real(dp), intent(in) :: s(:)
      real(dp), intent(in) :: theta(:)
      real(dp), intent(in) :: varphi(:)
      integer, intent(in) :: nfp
      real(dp), intent(out) :: lower(3)
      real(dp), intent(out) :: upper(3)

      lower = [s(1), theta(1), varphi(1)]
      upper(1) = s(size(s))
      upper(2) = theta(1) + 2.0_dp*acos(-1.0_dp)
      upper(3) = varphi(1) + 2.0_dp*acos(-1.0_dp)/real(nfp, dp)
   end subroutine periodic_bounds

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
      real(dp), allocatable :: field_batch(:, :, :, :)
      real(dp) :: lower(3)
      real(dp) :: upper(3)

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

      call read_family(ncid, data%family)
      call read_coords(ncid, s, theta, varphi)
      call read_profiles(ncid, size(s), profile_batch)
      call read_geometry(ncid, size(s), size(theta), size(varphi), geom_batch, &
                         data%has_lambda)
      call read_fields(ncid, size(s), size(theta), size(varphi), field_batch)
      status = nf90_get_att(ncid, NF90_GLOBAL, 'nfp', data%nfp)
      if (status /= NF90_NOERR) then
         print *, 'load_gvec_export_data: missing nfp attribute in ', trim(filename)
         error stop
      end if
      call nc_close(ncid)
      call periodic_bounds(s, theta, varphi, data%nfp, lower, upper)

      call construct_batch_splines_1d(s(1), s(size(s)), profile_batch, order_1d, &
                                      .false., data%profile_splines)
      call construct_batch_splines_3d([s(1), lower(2), lower(3)], &
                                      [s(size(s)), upper(2), upper(3)], &
                                      geom_batch(:, :, :, 1:8), order_3d, &
                                      periodic_3d, data%geom_splines(1))
      call construct_batch_splines_3d([s(1), lower(2), lower(3)], &
                                      [s(size(s)), upper(2), upper(3)], &
                                      geom_batch(:, :, :, 9:16), order_3d, &
                                      periodic_3d, data%geom_splines(2))
      call construct_batch_splines_3d([s(1), lower(2), lower(3)], &
                                      [s(size(s)), upper(2), upper(3)], &
                                      field_batch(:, :, :, 1:8), order_3d, &
                                      periodic_3d, data%field_splines(1))
      call construct_batch_splines_3d([s(1), lower(2), lower(3)], &
                                      [s(size(s)), upper(2), upper(3)], &
                                      field_batch(:, :, :, 9:9), order_3d, &
                                      periodic_3d, data%field_splines(2))

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
      real(dp) :: dvalues_1(3, 8)
      real(dp) :: dvalues_2(3, 8)
      real(dp) :: d2values_1(6, 8)
      real(dp) :: d2values_2(6, 8)

      call evaluate_batch_splines_3d_der2(self%geom_splines(1), x, values_1, &
                                          dvalues_1, d2values_1)
      call evaluate_batch_splines_3d_der2(self%geom_splines(2), x, values_2, &
                                          dvalues_2, d2values_2)

      values(1:8) = values_1
      values(9:16) = values_2

      values(gvec_geom_dx_ds) = dvalues_1(1, gvec_geom_x)
      values(gvec_geom_dx_dt) = dvalues_1(2, gvec_geom_x)
      values(gvec_geom_dx_dp) = dvalues_1(3, gvec_geom_x)
      values(gvec_geom_dy_ds) = dvalues_1(1, gvec_geom_y)
      values(gvec_geom_dy_dt) = dvalues_1(2, gvec_geom_y)
      values(gvec_geom_dy_dp) = dvalues_1(3, gvec_geom_y)
      values(gvec_geom_dz_ds) = dvalues_1(1, gvec_geom_z)
      values(gvec_geom_dz_dt) = dvalues_1(2, gvec_geom_z)
      values(gvec_geom_dz_dp) = dvalues_1(3, gvec_geom_z)
      if (self%has_lambda) then
         values(gvec_geom_dlambda_ds) = dvalues_2(1, 5)
         values(gvec_geom_dlambda_dt) = dvalues_2(2, 5)
         values(gvec_geom_dlambda_dp) = dvalues_2(3, 5)
      else
         values(gvec_geom_dlambda_ds:gvec_geom_dlambda_dp) = 0.0_dp
      end if

      if (present(dvalues)) then
         dvalues(:, 1:8) = dvalues_1
         dvalues(:, 9:16) = dvalues_2
      end if
      if (present(d2values)) then
         d2values(:, 1:8) = d2values_1
         d2values(:, 9:16) = d2values_2
      end if
   end subroutine gvec_export_evaluate_geometry

   subroutine gvec_export_evaluate_field(self, x, values)
      class(gvec_export_data_t), intent(in) :: self
      real(dp), intent(in) :: x(3)
      real(dp), intent(out) :: values(:)

      real(dp) :: values_1(8)
      real(dp) :: values_2(1)

      call evaluate_batch_splines_3d(self%field_splines(1), x, values_1)
      call evaluate_batch_splines_3d(self%field_splines(2), x, values_2)
      values(1:8) = values_1
      values(9) = values_2(1)
   end subroutine gvec_export_evaluate_field

   real(dp) function gvec_export_phi_period(self)
      class(gvec_export_data_t), intent(in) :: self

      gvec_export_phi_period = 2.0_dp*acos(-1.0_dp)/real(self%nfp, dp)
   end function gvec_export_phi_period

   logical function gvec_export_uses_vmec_adapter(self)
      class(gvec_export_data_t), intent(in) :: self

      gvec_export_uses_vmec_adapter = self%family == gvec_family_logical .and. &
                                      self%has_lambda
   end function gvec_export_uses_vmec_adapter

   logical function gvec_export_same_family(self, other)
      class(gvec_export_data_t), intent(in) :: self
      class(gvec_export_data_t), intent(in) :: other

      gvec_export_same_family = self%family == other%family
   end function gvec_export_same_family

   real(dp) function gvec_export_min_s(self)
      class(gvec_export_data_t), intent(in) :: self

      gvec_export_min_s = self%geom_splines(1)%x_min(1)
   end function gvec_export_min_s

   subroutine read_family(ncid, family)
      integer, intent(in) :: ncid
      integer, intent(out) :: family

      integer :: att_type
      integer :: att_len
      integer :: status
      character(len=:), allocatable :: value

      status = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'coordinate_family', &
                                      xtype=att_type, len=att_len)
      if (status /= NF90_NOERR .or. att_type /= NF90_CHAR .or. att_len < 1) then
         print *, 'gvec_export_data: missing coordinate_family attribute'
         error stop
      end if

      allocate (character(len=att_len) :: value)
      status = nf90_get_att(ncid, NF90_GLOBAL, 'coordinate_family', value)
      if (status /= NF90_NOERR) then
         print *, 'gvec_export_data: failed to read coordinate_family attribute'
         error stop
      end if

      select case (trim(value))
      case ('logical')
         family = gvec_family_logical
      case ('boozer')
         family = gvec_family_boozer
      case ('pest')
         family = gvec_family_pest
      case default
         print *, 'gvec_export_data: unsupported coordinate_family ', trim(value)
         error stop
      end select
   end subroutine read_family

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

   subroutine read_geometry(ncid, ns, ntheta, nvarphi, geom_batch, has_lambda)
      integer, intent(in) :: ncid
      integer, intent(in) :: ns
      integer, intent(in) :: ntheta
      integer, intent(in) :: nvarphi
      real(dp), allocatable, intent(out) :: geom_batch(:, :, :, :)
      logical, intent(out) :: has_lambda

      real(dp), allocatable :: temp(:, :, :)

      allocate (geom_batch(ns, ntheta, nvarphi, 16))
      allocate (temp(ns, ntheta, nvarphi))

      call nc_get(ncid, 'X', temp)
      geom_batch(:, :, :, gvec_geom_x) = temp
      call nc_get(ncid, 'Y', temp)
      geom_batch(:, :, :, gvec_geom_y) = temp
      call nc_get(ncid, 'Z', temp)
      geom_batch(:, :, :, gvec_geom_z) = temp
      call nc_get(ncid, 'dX_ds', temp)
      geom_batch(:, :, :, gvec_geom_dx_ds) = temp
      call nc_get(ncid, 'dX_dt', temp)
      geom_batch(:, :, :, gvec_geom_dx_dt) = temp
      call nc_get(ncid, 'dX_dp', temp)
      geom_batch(:, :, :, gvec_geom_dx_dp) = temp
      call nc_get(ncid, 'dY_ds', temp)
      geom_batch(:, :, :, gvec_geom_dy_ds) = temp
      call nc_get(ncid, 'dY_dt', temp)
      geom_batch(:, :, :, gvec_geom_dy_dt) = temp
      call nc_get(ncid, 'dY_dp', temp)
      geom_batch(:, :, :, gvec_geom_dy_dp) = temp
      call nc_get(ncid, 'dZ_ds', temp)
      geom_batch(:, :, :, gvec_geom_dz_ds) = temp
      call nc_get(ncid, 'dZ_dt', temp)
      geom_batch(:, :, :, gvec_geom_dz_dt) = temp
      call nc_get(ncid, 'dZ_dp', temp)
      geom_batch(:, :, :, gvec_geom_dz_dp) = temp

      has_lambda = has_variable(ncid, 'Lambda')
      if (has_lambda) then
         call nc_get(ncid, 'Lambda', temp)
         geom_batch(:, :, :, gvec_geom_lambda) = temp
         call nc_get(ncid, 'dLambda_ds', temp)
         geom_batch(:, :, :, gvec_geom_dlambda_ds) = temp
         call nc_get(ncid, 'dLambda_dt', temp)
         geom_batch(:, :, :, gvec_geom_dlambda_dt) = temp
         call nc_get(ncid, 'dLambda_dp', temp)
         geom_batch(:, :, :, gvec_geom_dlambda_dp) = temp
      else
         geom_batch(:, :, :, gvec_geom_lambda:gvec_geom_dlambda_dp) = 0.0_dp
      end if
   end subroutine read_geometry

   subroutine read_fields(ncid, ns, ntheta, nvarphi, field_batch)
      integer, intent(in) :: ncid
      integer, intent(in) :: ns
      integer, intent(in) :: ntheta
      integer, intent(in) :: nvarphi
      real(dp), allocatable, intent(out) :: field_batch(:, :, :, :)

      real(dp), allocatable :: temp(:, :, :)

      allocate (field_batch(ns, ntheta, nvarphi, 9))
      allocate (temp(ns, ntheta, nvarphi))

      call nc_get(ncid, 'Acov_s', temp)
      field_batch(:, :, :, gvec_field_acov_s) = temp
      call nc_get(ncid, 'Acov_t', temp)
      field_batch(:, :, :, gvec_field_acov_t) = temp
      call nc_get(ncid, 'Acov_p', temp)
      field_batch(:, :, :, gvec_field_acov_p) = temp
      call nc_get(ncid, 'hcov_s', temp)
      field_batch(:, :, :, gvec_field_hcov_s) = temp
      call nc_get(ncid, 'hcov_t', temp)
      field_batch(:, :, :, gvec_field_hcov_t) = temp
      call nc_get(ncid, 'hcov_p', temp)
      field_batch(:, :, :, gvec_field_hcov_p) = temp
      call nc_get(ncid, 'Bmod', temp)
      field_batch(:, :, :, gvec_field_bmod) = temp
      call nc_get(ncid, 'sqgBctr_t', temp)
      field_batch(:, :, :, gvec_field_sqgbctr_t) = temp
      call nc_get(ncid, 'sqgBctr_p', temp)
      field_batch(:, :, :, gvec_field_sqgbctr_p) = temp
   end subroutine read_fields

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
      use netcdf, only: nf90_inquire_dimension, nf90_inquire_variable

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

   logical function has_variable(ncid, name)
      integer, intent(in) :: ncid
      character(*), intent(in) :: name

      integer :: varid
      integer :: status

      status = nf90_inq_varid(ncid, trim(name), varid)
      has_variable = status == NF90_NOERR
   end function has_variable

end module gvec_export_data
