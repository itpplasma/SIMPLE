module gvec_reference_coordinates
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use cylindrical_cartesian, only: cart_to_cyl
   use gvec_export_data, only: gvec_export_data_t, load_gvec_export_data, &
                               gvec_family_logical, &
                               gvec_geom_x, gvec_geom_y, gvec_geom_z, &
                               gvec_geom_dx_ds, gvec_geom_dx_dt, gvec_geom_dx_dp, &
                               gvec_geom_dy_ds, gvec_geom_dy_dt, gvec_geom_dy_dp, &
                               gvec_geom_dz_ds, gvec_geom_dz_dt, gvec_geom_dz_dp
   use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                          evaluate_batch_splines_3d_der2
   use libneo_coordinates, only: coordinate_system_t
   use math_constants, only: TWOPI
   use nctools_module, only: nc_get
   use netcdf, only: NF90_NOERR, NF90_NOWRITE, nf90_close, nf90_get_att, &
                     nf90_inq_dimid, nf90_inquire_dimension, nf90_open

   implicit none

   private

   type, extends(coordinate_system_t), public :: gvec_coordinate_system_t
      type(gvec_export_data_t) :: data
      type(BatchSplineData3D) :: logical_geom_spline
      logical :: has_logical_geom = .false.
   contains
      procedure :: evaluate_cart => gvec_evaluate_cart
      procedure :: evaluate_cyl => gvec_evaluate_cyl
      procedure :: covariant_basis => gvec_covariant_basis
      procedure :: metric_tensor => gvec_metric_tensor
      procedure :: from_cyl => gvec_from_cyl
      procedure :: phi_period => gvec_coordinate_phi_period
   end type gvec_coordinate_system_t

   public :: make_gvec_coordinate_system

contains

   subroutine make_gvec_coordinate_system(cs, filename)
      class(coordinate_system_t), allocatable, intent(out) :: cs
      character(*), intent(in) :: filename

      allocate (gvec_coordinate_system_t :: cs)
      select type (cs)
      type is (gvec_coordinate_system_t)
         call load_gvec_export_data(filename, cs%data)
         if (cs%data%family == gvec_family_logical) then
            call load_logical_geometry_spline(filename, cs%logical_geom_spline)
            cs%has_logical_geom = .true.
         end if
      class default
         error stop 'make_gvec_coordinate_system: allocation failure'
      end select
   end subroutine make_gvec_coordinate_system

   subroutine gvec_evaluate_cart(self, u, x)
      class(gvec_coordinate_system_t), intent(in) :: self
      real(dp), intent(in) :: u(3)
      real(dp), intent(out) :: x(3)

      real(dp) :: values(16)
      real(dp) :: geom_values(2)
      real(dp) :: geom_derivs(3, 2)
      real(dp) :: geom_der2(6, 2)
      real(dp) :: phi

      if (.not. self%has_logical_geom) then
         call self%data%evaluate_geometry(u, values)
         x = [values(gvec_geom_x), values(gvec_geom_y), values(gvec_geom_z)]
         return
      end if

      call evaluate_batch_splines_3d_der2(self%logical_geom_spline, u, geom_values, &
                                          geom_derivs, geom_der2)
      phi = u(3)
      x(1) = geom_values(1)*cos(phi)
      x(2) = geom_values(1)*sin(phi)
      x(3) = geom_values(2)
   end subroutine gvec_evaluate_cart

   subroutine gvec_evaluate_cyl(self, u, x)
      class(gvec_coordinate_system_t), intent(in) :: self
      real(dp), intent(in) :: u(3)
      real(dp), intent(out) :: x(3)

      real(dp) :: xcart(3)

      call self%evaluate_cart(u, xcart)
      call cart_to_cyl(xcart, x)
   end subroutine gvec_evaluate_cyl

   subroutine gvec_covariant_basis(self, u, e_cov)
      class(gvec_coordinate_system_t), intent(in) :: self
      real(dp), intent(in) :: u(3)
      real(dp), intent(out) :: e_cov(3, 3)

      real(dp) :: values(16)
      real(dp) :: geom_values(2)
      real(dp) :: geom_derivs(3, 2)
      real(dp) :: geom_der2(6, 2)
      real(dp) :: r
      real(dp) :: dr_ds
      real(dp) :: dr_dt
      real(dp) :: dr_dp
      real(dp) :: dz_ds
      real(dp) :: dz_dt
      real(dp) :: dz_dp
      real(dp) :: phi

      if (.not. self%has_logical_geom) then
         call self%data%evaluate_geometry(u, values)
         e_cov(:, 1) = [values(gvec_geom_dx_ds), values(gvec_geom_dy_ds), &
                        values(gvec_geom_dz_ds)]
         e_cov(:, 2) = [values(gvec_geom_dx_dt), values(gvec_geom_dy_dt), &
                        values(gvec_geom_dz_dt)]
         e_cov(:, 3) = [values(gvec_geom_dx_dp), values(gvec_geom_dy_dp), &
                        values(gvec_geom_dz_dp)]
         return
      end if

      call evaluate_batch_splines_3d_der2(self%logical_geom_spline, u, geom_values, &
                                          geom_derivs, geom_der2)
      r = geom_values(1)
      dr_ds = geom_derivs(1, 1)
      dr_dt = geom_derivs(2, 1)
      dr_dp = geom_derivs(3, 1)
      dz_ds = geom_derivs(1, 2)
      dz_dt = geom_derivs(2, 2)
      dz_dp = geom_derivs(3, 2)
      phi = u(3)

      e_cov(:, 1) = [dr_ds*cos(phi), dr_ds*sin(phi), dz_ds]
      e_cov(:, 2) = [dr_dt*cos(phi), dr_dt*sin(phi), dz_dt]
      e_cov(:, 3) = [dr_dp*cos(phi) - r*sin(phi), dr_dp*sin(phi) + r*cos(phi), dz_dp]
   end subroutine gvec_covariant_basis

   subroutine gvec_metric_tensor(self, u, g, ginv, sqrtg)
      class(gvec_coordinate_system_t), intent(in) :: self
      real(dp), intent(in) :: u(3)
      real(dp), intent(out) :: g(3, 3)
      real(dp), intent(out) :: ginv(3, 3)
      real(dp), intent(out) :: sqrtg

      real(dp) :: e_cov(3, 3)
      real(dp) :: det
      integer :: i
      integer :: j

      call self%covariant_basis(u, e_cov)

      do j = 1, 3
         do i = 1, 3
            g(i, j) = dot_product(e_cov(:, i), e_cov(:, j))
         end do
      end do

      det = g(1, 1)*(g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2)) - &
            g(1, 2)*(g(2, 1)*g(3, 3) - g(2, 3)*g(3, 1)) + &
            g(1, 3)*(g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))

      sqrtg = sqrt(abs(det))

      ginv(1, 1) = (g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2))/det
      ginv(1, 2) = (g(1, 3)*g(3, 2) - g(1, 2)*g(3, 3))/det
      ginv(1, 3) = (g(1, 2)*g(2, 3) - g(1, 3)*g(2, 2))/det
      ginv(2, 1) = (g(2, 3)*g(3, 1) - g(2, 1)*g(3, 3))/det
      ginv(2, 2) = (g(1, 1)*g(3, 3) - g(1, 3)*g(3, 1))/det
      ginv(2, 3) = (g(1, 3)*g(2, 1) - g(1, 1)*g(2, 3))/det
      ginv(3, 1) = (g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))/det
      ginv(3, 2) = (g(1, 2)*g(3, 1) - g(1, 1)*g(3, 2))/det
      ginv(3, 3) = (g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1))/det
   end subroutine gvec_metric_tensor

   subroutine gvec_from_cyl(self, xcyl, u, ierr)
      class(gvec_coordinate_system_t), intent(in) :: self
      real(dp), intent(in) :: xcyl(3)
      real(dp), intent(out) :: u(3)
      integer, intent(out) :: ierr

      real(dp) :: axis_values(16)
      real(dp) :: boundary_values(16)
      real(dp) :: s
      real(dp) :: theta
      real(dp) :: phi
      real(dp) :: r_axis
      real(dp) :: z_axis
      real(dp) :: r_bnd
      real(dp) :: z_bnd
      real(dp) :: rs
      real(dp) :: zs
      real(dp) :: dist_axis
      real(dp) :: dist_bnd
      real(dp) :: s_seed(4)
      real(dp) :: theta_seed(4)
      real(dp) :: phi_seed(4)
      integer :: is
      integer :: it
      integer :: ip
      integer :: ierr_try

      ierr = 0

      call evaluate_cyl_state_at(self, [self%data%min_s(), 0.0_dp, 0.0_dp], axis_values)
      r_axis = axis_values(1)
      z_axis = axis_values(3)

      rs = xcyl(1) - r_axis
      zs = xcyl(3) - z_axis
      theta = modulo(atan2(zs, rs), TWOPI)
      phi = modulo(xcyl(2), self%phi_period())

      call evaluate_cyl_state_at(self, [1.0_dp, theta, phi], boundary_values)
      r_bnd = boundary_values(1)
      z_bnd = boundary_values(3)

      dist_axis = sqrt(rs**2 + zs**2)
      dist_bnd = sqrt((r_bnd - r_axis)**2 + (z_bnd - z_axis)**2)
      if (dist_bnd > 1.0e-12_dp) then
         s = min(1.0_dp, max(self%data%min_s(), (dist_axis/dist_bnd)**2))
      else
         s = self%data%min_s()
      end if

      s_seed = [s, max(self%data%min_s(), 0.5_dp*s), min(1.0_dp, 1.5_dp*s), 0.5_dp]
      theta_seed = [theta, modulo(theta + 0.5_dp*pi_dp(), TWOPI), &
                    modulo(theta + pi_dp(), TWOPI), &
                    modulo(theta + 1.5_dp*pi_dp(), TWOPI)]
      phi_seed = [phi, modulo(phi + 0.25_dp*self%phi_period(), self%phi_period()), &
                  modulo(phi + 0.5_dp*self%phi_period(), self%phi_period()), &
                  modulo(phi + 0.75_dp*self%phi_period(), self%phi_period())]

      do is = 1, size(s_seed)
         do it = 1, size(theta_seed)
            do ip = 1, size(phi_seed)
               call newton_from_seed(self, xcyl, s_seed(is), theta_seed(it), &
                                     phi_seed(ip), u, ierr_try)
               if (ierr_try == 0) then
                  ierr = 0
                  return
               end if
            end do
         end do
      end do

      ierr = 1
   end subroutine gvec_from_cyl

   subroutine newton_from_seed(self, xcyl, s0, theta0, phi0, u, ierr)
      class(gvec_coordinate_system_t), intent(in) :: self
      real(dp), intent(in) :: xcyl(3)
      real(dp), intent(in) :: s0
      real(dp), intent(in) :: theta0
      real(dp), intent(in) :: phi0
      real(dp), intent(out) :: u(3)
      integer, intent(out) :: ierr

      integer, parameter :: max_iter = 50
      real(dp), parameter :: tol_res = 1.0e-10_dp
      real(dp), parameter :: tol_step = 1.0e-12_dp

      real(dp) :: xcyl_current(3)
      real(dp) :: f(3)
      real(dp) :: jac(3, 3)
      real(dp) :: delta(3)
      real(dp) :: res_norm
      real(dp) :: res_norm_try
      real(dp) :: alpha
      real(dp) :: s
      real(dp) :: theta
      real(dp) :: phi
      real(dp) :: s_try
      real(dp) :: theta_try
      real(dp) :: phi_try
      integer :: iter
      integer :: k

      s = s0
      theta = theta0
      phi = phi0

      do iter = 1, max_iter
         call evaluate_cyl_state_at(self, [s, theta, phi], xcyl_current, jac)
         f(1) = xcyl_current(1) - xcyl(1)
         f(2) = angle_difference(xcyl_current(2), xcyl(2))
         f(3) = xcyl_current(3) - xcyl(3)
         res_norm = sqrt(sum(f**2))
         if (res_norm < tol_res) then
            u = [s, theta, phi]
            ierr = 0
            return
         end if

         call solve_3x3(jac, -f, delta, ierr)
         if (ierr /= 0) then
            ierr = 2
            return
         end if

         if (maxval(abs(delta)) < tol_step) then
            u = [s, theta, phi]
            ierr = 0
            return
         end if

         alpha = 1.0_dp
         do k = 1, 10
            s_try = min(1.0_dp, max(self%data%min_s(), s + alpha*delta(1)))
            theta_try = modulo(theta + alpha*delta(2), TWOPI)
            phi_try = modulo(phi + alpha*delta(3), self%phi_period())
            call evaluate_cyl_state_at(self, [s_try, theta_try, phi_try], xcyl_current)
            f(1) = xcyl_current(1) - xcyl(1)
            f(2) = angle_difference(xcyl_current(2), xcyl(2))
            f(3) = xcyl_current(3) - xcyl(3)
            res_norm_try = sqrt(sum(f**2))
            if (res_norm_try < res_norm) then
               s = s_try
               theta = theta_try
               phi = phi_try
               exit
            end if
            alpha = 0.5_dp*alpha
         end do

         if (k > 10) then
            ierr = 1
            return
         end if
      end do

      ierr = 1
   end subroutine newton_from_seed

   subroutine evaluate_cyl_state_at(self, u, xcyl, jac)
      class(gvec_coordinate_system_t), intent(in) :: self
      real(dp), intent(in) :: u(3)
      real(dp), intent(out) :: xcyl(3)
      real(dp), intent(out), optional :: jac(3, 3)

      real(dp) :: values(16)
      real(dp) :: geom_values(2)
      real(dp) :: geom_derivs(3, 2)
      real(dp) :: geom_der2(6, 2)
      real(dp) :: r
      real(dp) :: phi
      real(dp) :: dr_ds
      real(dp) :: dr_dt
      real(dp) :: dr_dp

      if (.not. self%has_logical_geom) then
         call self%data%evaluate_geometry(u, values)
         call evaluate_cyl_state(values, xcyl, jac)
         return
      end if

      call evaluate_batch_splines_3d_der2(self%logical_geom_spline, u, geom_values, &
                                          geom_derivs, geom_der2)
      r = geom_values(1)
      phi = u(3)

      xcyl(1) = r
      xcyl(2) = modulo(phi, self%phi_period())
      xcyl(3) = geom_values(2)

      if (.not. present(jac)) return

      dr_ds = geom_derivs(1, 1)
      dr_dt = geom_derivs(2, 1)
      dr_dp = geom_derivs(3, 1)

      jac(1, 1) = dr_ds
      jac(1, 2) = dr_dt
      jac(1, 3) = dr_dp
      jac(2, 1) = 0.0_dp
      jac(2, 2) = 0.0_dp
      jac(2, 3) = 1.0_dp
      jac(3, 1) = geom_derivs(1, 2)
      jac(3, 2) = geom_derivs(2, 2)
      jac(3, 3) = geom_derivs(3, 2)
   end subroutine evaluate_cyl_state_at

   subroutine evaluate_cyl_state(values, xcyl, jac)
      real(dp), intent(in) :: values(16)
      real(dp), intent(out) :: xcyl(3)
      real(dp), intent(out), optional :: jac(3, 3)

      real(dp) :: x
      real(dp) :: y
      real(dp) :: r
      real(dp) :: r2

      x = values(gvec_geom_x)
      y = values(gvec_geom_y)
      r2 = x*x + y*y
      r = sqrt(r2)

      xcyl(1) = r
      xcyl(2) = modulo(atan2(y, x), TWOPI)
      xcyl(3) = values(gvec_geom_z)

      if (.not. present(jac)) return

      if (r > 1.0e-14_dp) then
         jac(1, 1) = (x*values(gvec_geom_dx_ds) + y*values(gvec_geom_dy_ds))/r
         jac(1, 2) = (x*values(gvec_geom_dx_dt) + y*values(gvec_geom_dy_dt))/r
         jac(1, 3) = (x*values(gvec_geom_dx_dp) + y*values(gvec_geom_dy_dp))/r
      else
         jac(1, :) = 0.0_dp
      end if

      if (r2 > 1.0e-24_dp) then
         jac(2, 1) = (x*values(gvec_geom_dy_ds) - y*values(gvec_geom_dx_ds))/r2
         jac(2, 2) = (x*values(gvec_geom_dy_dt) - y*values(gvec_geom_dx_dt))/r2
         jac(2, 3) = (x*values(gvec_geom_dy_dp) - y*values(gvec_geom_dx_dp))/r2
      else
         jac(2, :) = 0.0_dp
      end if

      jac(3, 1) = values(gvec_geom_dz_ds)
      jac(3, 2) = values(gvec_geom_dz_dt)
      jac(3, 3) = values(gvec_geom_dz_dp)
   end subroutine evaluate_cyl_state

   real(dp) function angle_difference(angle, target)
      real(dp), intent(in) :: angle
      real(dp), intent(in) :: target

      angle_difference = modulo(angle - target + 0.5_dp*TWOPI, TWOPI) - 0.5_dp*TWOPI
   end function angle_difference

   real(dp) function pi_dp()
      pi_dp = acos(-1.0_dp)
   end function pi_dp

   real(dp) function gvec_coordinate_phi_period(self)
      class(gvec_coordinate_system_t), intent(in) :: self

      gvec_coordinate_phi_period = self%data%phi_period()
   end function gvec_coordinate_phi_period

   subroutine load_logical_geometry_spline(filename, spline)
      character(*), intent(in) :: filename
      type(BatchSplineData3D), intent(out) :: spline

      integer, parameter :: order_3d(3) = [5, 5, 5]
      logical, parameter :: periodic_3d(3) = [.false., .true., .true.]

      integer :: ncid
      integer :: status
      integer :: nfp
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: varphi(:)
      real(dp), allocatable :: batch(:, :, :, :)
      real(dp), allocatable :: temp(:, :, :)
      real(dp) :: lower(3)
      real(dp) :: upper(3)

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

      allocate (batch(size(s), size(theta), size(varphi), 2))
      allocate (temp(size(s), size(theta), size(varphi)))

      call nc_get(ncid, 'R', temp)
      batch(:, :, :, 1) = temp
      call nc_get(ncid, 'Z_cyl', temp)
      batch(:, :, :, 2) = temp

      lower = [s(1), theta(1), varphi(1)]
      upper = [s(size(s)), theta(1) + 2.0_dp*acos(-1.0_dp), &
               varphi(1) + 2.0_dp*acos(-1.0_dp)/real(nfp, dp)]

      call construct_batch_splines_3d(lower, upper, batch, order_3d, periodic_3d, &
                                      spline)
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

   subroutine solve_3x3(a, b, x, ierr)
      real(dp), intent(in) :: a(3, 3)
      real(dp), intent(in) :: b(3)
      real(dp), intent(out) :: x(3)
      integer, intent(out) :: ierr

      real(dp) :: det
      real(dp) :: inv(3, 3)

      det = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) - &
            a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) + &
            a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
      if (abs(det) < 1.0e-14_dp) then
         ierr = 1
         return
      end if

      inv(1, 1) = (a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2))/det
      inv(1, 2) = (a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3))/det
      inv(1, 3) = (a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2))/det
      inv(2, 1) = (a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3))/det
      inv(2, 2) = (a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1))/det
      inv(2, 3) = (a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3))/det
      inv(3, 1) = (a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))/det
      inv(3, 2) = (a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2))/det
      inv(3, 3) = (a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))/det
      x = matmul(inv, b)
      ierr = 0
   end subroutine solve_3x3

end module gvec_reference_coordinates
