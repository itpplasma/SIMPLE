module field_gvec
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use field_base, only: magnetic_field_t
   use gvec_export_data, only: gvec_export_data_t, load_gvec_export_data, &
                               gvec_family_logical, gvec_profile_a_theta, &
                               gvec_profile_a_phi, gvec_profile_da_theta_ds, &
                               gvec_profile_da_phi_ds, gvec_field_acov_s, &
                               gvec_field_acov_t, gvec_field_acov_p, &
                               gvec_field_hcov_s, gvec_field_hcov_t, &
                               gvec_field_hcov_p, gvec_field_bmod, &
                               gvec_field_sqgbctr_t, gvec_field_sqgbctr_p, &
                               gvec_geom_dx_ds, gvec_geom_dx_dt, gvec_geom_dx_dp, &
                               gvec_geom_dy_ds, gvec_geom_dy_dt, gvec_geom_dy_dp, &
                               gvec_geom_dz_ds, gvec_geom_dz_dt, gvec_geom_dz_dp
   use gvec_reference_coordinates, only: gvec_coordinate_system_t
   use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                          evaluate_batch_splines_3d, evaluate_batch_splines_3d_der2
   use nctools_module, only: nc_get
   use netcdf, only: NF90_NOERR, NF90_NOWRITE, nf90_close, nf90_get_att, &
                     nf90_inq_dimid, nf90_inq_varid, nf90_inquire_dimension, nf90_open
   use spline_vmec_sub, only: compute_field_components

   implicit none

   private
   public :: gvec_field_t, create_gvec_field

   type, extends(magnetic_field_t) :: gvec_field_t
      type(gvec_export_data_t) :: data
      type(BatchSplineData3D) :: logical_geom_spline
      type(BatchSplineData3D) :: bcart_spline
      character(len=256) :: filename = ''
      logical :: has_logical_geom = .false.
      logical :: has_bcart = .false.
   contains
      procedure :: evaluate => gvec_evaluate
   end type gvec_field_t

contains

   subroutine create_gvec_field(filename, field)
      character(*), intent(in) :: filename
      class(gvec_field_t), allocatable, intent(out) :: field

      allocate (gvec_field_t :: field)
      field%filename = filename
      call load_gvec_export_data(filename, field%data)
      if (field%data%family == gvec_family_logical) then
         call load_logical_geometry_spline(filename, field%logical_geom_spline)
         field%has_logical_geom = .true.
      else
         call load_bcart_spline(filename, field%bcart_spline, field%has_bcart)
      end if

      allocate (gvec_coordinate_system_t :: field%coords)
      select type (coords => field%coords)
      type is (gvec_coordinate_system_t)
         coords%data = field%data
         coords%has_logical_geom = field%has_logical_geom
         if (field%has_logical_geom) coords%logical_geom_spline = field%logical_geom_spline
      class default
         error stop 'create_gvec_field: allocation failure'
      end select
   end subroutine create_gvec_field

   subroutine gvec_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
      class(gvec_field_t), intent(in) :: self
      real(dp), intent(in) :: x(3)
      real(dp), intent(out) :: Acov(3)
      real(dp), intent(out) :: hcov(3)
      real(dp), intent(out) :: Bmod
      real(dp), intent(out), optional :: sqgBctr(3)

      real(dp) :: values(9)
      real(dp) :: profiles(5)
      real(dp) :: geom_values(3)
      real(dp) :: geom_derivs(3, 3)
      real(dp) :: geom_der2(6, 3)
      real(dp) :: sqg
      real(dp) :: Bctr_t
      real(dp) :: Bctr_p
      real(dp) :: Bcov_s
      real(dp) :: Bcov_t
      real(dp) :: Bcov_p
      real(dp) :: geometry(16)
      real(dp) :: bcart(3)
      real(dp) :: e_s(3)
      real(dp) :: e_t(3)
      real(dp) :: e_p(3)

      if (self%has_logical_geom) then
         call self%data%evaluate_profiles(x(1), profiles)
         call evaluate_batch_splines_3d_der2(self%logical_geom_spline, x, geom_values, &
                                             geom_derivs, geom_der2)

         call compute_field_components(geom_values(1), geom_derivs(1, 1), &
                                       geom_derivs(2, 1), geom_derivs(3, 1), &
                                       geom_derivs(1, 2), geom_derivs(2, 2), &
                                       geom_derivs(3, 2), &
                                       profiles(gvec_profile_da_theta_ds), &
                                       profiles(gvec_profile_da_phi_ds), &
                                       geom_derivs(1, 3), geom_derivs(2, 3), &
                                       geom_derivs(3, 3), sqg, Bctr_t, Bctr_p, &
                                       Bcov_s, Bcov_t, Bcov_p)

         Acov(1) = profiles(gvec_profile_a_theta)*geom_derivs(1, 3)
         Acov(2) = profiles(gvec_profile_a_theta)*(1.0_dp + geom_derivs(2, 3))
         Acov(3) = profiles(gvec_profile_a_phi) + &
                   profiles(gvec_profile_a_theta)*geom_derivs(3, 3)

         Bmod = sqrt(Bctr_t*Bcov_t + Bctr_p*Bcov_p)
         hcov(1) = (Bcov_s + Bcov_t*geom_derivs(1, 3))/Bmod
         hcov(2) = Bcov_t*(1.0_dp + geom_derivs(2, 3))/Bmod
         hcov(3) = (Bcov_p + Bcov_t*geom_derivs(3, 3))/Bmod

         if (present(sqgBctr)) then
            sqgBctr = [0.0_dp, sqg*Bctr_t, sqg*Bctr_p]
         end if
         return
      end if

      call self%data%evaluate_field(x, values)
      Acov = [values(gvec_field_acov_s), values(gvec_field_acov_t), &
              values(gvec_field_acov_p)]
      if (.not. self%has_bcart) then
         error stop 'gvec_evaluate: missing cartesian field data for straight export'
      end if

      call self%data%evaluate_geometry(x, geometry)
      call evaluate_batch_splines_3d(self%bcart_spline, x, bcart)

      e_s = [geometry(gvec_geom_dx_ds), geometry(gvec_geom_dy_ds), geometry(gvec_geom_dz_ds)]
      e_t = [geometry(gvec_geom_dx_dt), geometry(gvec_geom_dy_dt), geometry(gvec_geom_dz_dt)]
      e_p = [geometry(gvec_geom_dx_dp), geometry(gvec_geom_dy_dp), geometry(gvec_geom_dz_dp)]

      Bcov_s = dot_product(bcart, e_s)
      Bcov_t = dot_product(bcart, e_t)
      Bcov_p = dot_product(bcart, e_p)
      Bmod = sqrt(dot_product(bcart, bcart))
      hcov = [Bcov_s/Bmod, Bcov_t/Bmod, Bcov_p/Bmod]

      if (present(sqgBctr)) then
         sqgBctr = [0.0_dp, values(gvec_field_sqgbctr_t), &
                    values(gvec_field_sqgbctr_p)]
      end if
   end subroutine gvec_evaluate

   subroutine load_bcart_spline(filename, spline, has_bcart)
      character(*), intent(in) :: filename
      type(BatchSplineData3D), intent(out) :: spline
      logical, intent(out) :: has_bcart

      integer, parameter :: order_3d(3) = [5, 5, 5]
      logical, parameter :: periodic_3d(3) = [.false., .true., .true.]

      integer :: ncid
      integer :: status
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: varphi(:)
      real(dp), allocatable :: batch(:, :, :, :)
      real(dp), allocatable :: temp(:, :, :)
      real(dp) :: lower(3)
      real(dp) :: upper(3)
      integer :: nfp

      has_bcart = .false.
      status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
      if (status /= NF90_NOERR) then
         error stop 'load_bcart_spline: failed to open export file'
      end if
      if (nf90_inq_varid(ncid, 'B_x', status) /= NF90_NOERR) then
         call nf90_close_checked(ncid, 'load_bcart_spline')
         return
      end if

      call read_coord_1d(ncid, 's', s)
      call read_coord_1d(ncid, 'theta', theta)
      call read_coord_1d(ncid, 'varphi', varphi)
      if (nf90_get_att(ncid, 0, 'nfp', nfp) /= NF90_NOERR) then
         error stop 'load_bcart_spline: missing nfp attribute'
      end if

      allocate (batch(size(s), size(theta), size(varphi), 3))
      allocate (temp(size(s), size(theta), size(varphi)))

      call nc_get(ncid, 'B_x', temp)
      batch(:, :, :, 1) = temp
      call nc_get(ncid, 'B_y', temp)
      batch(:, :, :, 2) = temp
      call nc_get(ncid, 'B_z', temp)
      batch(:, :, :, 3) = temp

      lower = [s(1), theta(1), varphi(1)]
      upper = [s(size(s)), theta(1) + 2.0_dp*acos(-1.0_dp), &
               varphi(1) + 2.0_dp*acos(-1.0_dp)/real(nfp, dp)]

      call construct_batch_splines_3d(lower, upper, batch, order_3d, periodic_3d, spline)
      has_bcart = .true.
      call nf90_close_checked(ncid, 'load_bcart_spline')
   end subroutine load_bcart_spline

   subroutine load_logical_geometry_spline(filename, spline)
      character(*), intent(in) :: filename
      type(BatchSplineData3D), intent(out) :: spline

      integer, parameter :: order_3d(3) = [5, 5, 5]
      logical, parameter :: periodic_3d(3) = [.false., .true., .true.]

      integer :: ncid
      integer :: status
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: varphi(:)
      real(dp), allocatable :: batch(:, :, :, :)
      real(dp), allocatable :: temp(:, :, :)
      real(dp) :: lower(3)
      real(dp) :: upper(3)
      integer :: nfp

      status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
      if (status /= NF90_NOERR) then
         error stop 'load_logical_geometry_spline: failed to open export file'
      end if

      call read_coord_1d(ncid, 's', s)
      call read_coord_1d(ncid, 'theta', theta)
      call read_coord_1d(ncid, 'varphi', varphi)
      if (nf90_get_att(ncid, 0, 'nfp', nfp) /= NF90_NOERR) then
         error stop 'load_logical_geometry_spline: missing nfp attribute'
      end if

      allocate (batch(size(s), size(theta), size(varphi), 3))
      allocate (temp(size(s), size(theta), size(varphi)))

      call nc_get(ncid, 'R', temp)
      batch(:, :, :, 1) = temp
      call nc_get(ncid, 'Z_cyl', temp)
      batch(:, :, :, 2) = temp
      call nc_get(ncid, 'Lambda', temp)
      batch(:, :, :, 3) = temp

      lower = [s(1), theta(1), varphi(1)]
      upper = [s(size(s)), theta(1) + 2.0_dp*acos(-1.0_dp), &
               varphi(1) + 2.0_dp*acos(-1.0_dp)/real(nfp, dp)]

      call construct_batch_splines_3d(lower, upper, &
                                      batch, order_3d, periodic_3d, spline)
      call nf90_close_checked(ncid, 'load_logical_geometry_spline')
   end subroutine load_logical_geometry_spline

   subroutine nf90_close_checked(ncid, caller)
      integer, intent(in) :: ncid
      character(*), intent(in) :: caller

      if (nf90_close(ncid) /= NF90_NOERR) then
         error stop trim(caller)//': failed to close export file'
      end if
   end subroutine nf90_close_checked

   subroutine read_coord_1d(ncid, name, values)
      integer, intent(in) :: ncid
      character(*), intent(in) :: name
      real(dp), allocatable, intent(out) :: values(:)

      integer :: dimid
      integer :: n

      if (nf90_inq_dimid(ncid, trim(name), dimid) /= NF90_NOERR) then
         error stop 'read_coord_1d: missing coordinate dimension'
      end if
      if (nf90_inquire_dimension(ncid, dimid, len=n) /= NF90_NOERR) then
         error stop 'read_coord_1d: failed to query coordinate dimension'
      end if

      allocate (values(n))
      call nc_get(ncid, name, values)
   end subroutine read_coord_1d

end module field_gvec
