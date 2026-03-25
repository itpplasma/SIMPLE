module gvec_reference_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cylindrical_cartesian, only: cyl_to_cart
    use gvec_export_data, only: gvec_export_data_t, load_gvec_export_data, &
                                gvec_geom_r, gvec_geom_z, gvec_geom_dr_ds, &
                                gvec_geom_dr_dt, gvec_geom_dr_dp, gvec_geom_dz_ds, &
                                gvec_geom_dz_dt, gvec_geom_dz_dp
    use libneo_coordinates, only: coordinate_system_t
    use math_constants, only: TWOPI

    implicit none

    private

    type, extends(coordinate_system_t), public :: gvec_coordinate_system_t
        type(gvec_export_data_t) :: data
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
        class default
            error stop 'make_gvec_coordinate_system: allocation failure'
        end select
    end subroutine make_gvec_coordinate_system

    subroutine gvec_evaluate_cyl(self, u, x)
        class(gvec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: values(18)

        call self%data%evaluate_geometry(u, values)
        x = [values(gvec_geom_r), u(3), values(gvec_geom_z)]
    end subroutine gvec_evaluate_cyl

    subroutine gvec_evaluate_cart(self, u, x)
        class(gvec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: xcyl(3)

        call self%evaluate_cyl(u, xcyl)
        call cyl_to_cart(xcyl, x)
    end subroutine gvec_evaluate_cart

    subroutine gvec_covariant_basis(self, u, e_cov)
        class(gvec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        real(dp) :: values(18)
        real(dp) :: r
        real(dp) :: cos_phi
        real(dp) :: sin_phi

        call self%data%evaluate_geometry(u, values)

        r = values(gvec_geom_r)
        cos_phi = cos(u(3))
        sin_phi = sin(u(3))

        e_cov(1, 1) = values(gvec_geom_dr_ds) * cos_phi
        e_cov(2, 1) = values(gvec_geom_dr_ds) * sin_phi
        e_cov(3, 1) = values(gvec_geom_dz_ds)

        e_cov(1, 2) = values(gvec_geom_dr_dt) * cos_phi
        e_cov(2, 2) = values(gvec_geom_dr_dt) * sin_phi
        e_cov(3, 2) = values(gvec_geom_dz_dt)

        e_cov(1, 3) = values(gvec_geom_dr_dp) * cos_phi - r * sin_phi
        e_cov(2, 3) = values(gvec_geom_dr_dp) * sin_phi + r * cos_phi
        e_cov(3, 3) = values(gvec_geom_dz_dp)
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

        det = g(1, 1) * (g(2, 2) * g(3, 3) - g(2, 3) * g(3, 2)) - &
              g(1, 2) * (g(2, 1) * g(3, 3) - g(2, 3) * g(3, 1)) + &
              g(1, 3) * (g(2, 1) * g(3, 2) - g(2, 2) * g(3, 1))

        sqrtg = sqrt(abs(det))

        ginv(1, 1) = (g(2, 2) * g(3, 3) - g(2, 3) * g(3, 2)) / det
        ginv(1, 2) = (g(1, 3) * g(3, 2) - g(1, 2) * g(3, 3)) / det
        ginv(1, 3) = (g(1, 2) * g(2, 3) - g(1, 3) * g(2, 2)) / det
        ginv(2, 1) = (g(2, 3) * g(3, 1) - g(2, 1) * g(3, 3)) / det
        ginv(2, 2) = (g(1, 1) * g(3, 3) - g(1, 3) * g(3, 1)) / det
        ginv(2, 3) = (g(1, 3) * g(2, 1) - g(1, 1) * g(2, 3)) / det
        ginv(3, 1) = (g(2, 1) * g(3, 2) - g(2, 2) * g(3, 1)) / det
        ginv(3, 2) = (g(1, 2) * g(3, 1) - g(1, 1) * g(3, 2)) / det
        ginv(3, 3) = (g(1, 1) * g(2, 2) - g(1, 2) * g(2, 1)) / det
    end subroutine gvec_metric_tensor

    subroutine gvec_from_cyl(self, xcyl, u, ierr)
        class(gvec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        integer, parameter :: max_iter = 50
        real(dp), parameter :: tol_res = 1.0e-10_dp
        real(dp), parameter :: tol_step = 1.0e-12_dp

        real(dp) :: axis_values(18)
        real(dp) :: boundary_values(18)
        real(dp) :: values(18)
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
        real(dp) :: f(2)
        real(dp) :: jac(2, 2)
        real(dp) :: det
        real(dp) :: delta(2)
        real(dp) :: res_norm
        real(dp) :: res_norm_try
        real(dp) :: alpha
        real(dp) :: s_try
        real(dp) :: theta_try
        integer :: iter
        integer :: k

        ierr = 0
        phi = xcyl(2)

        call self%data%evaluate_geometry([self%data%geom_splines(1)%x_min(1), 0.0_dp, phi], &
                                         axis_values)
        r_axis = axis_values(gvec_geom_r)
        z_axis = axis_values(gvec_geom_z)

        rs = xcyl(1) - r_axis
        zs = xcyl(3) - z_axis
        theta = modulo(atan2(zs, rs), TWOPI)

        call self%data%evaluate_geometry([1.0_dp, theta, phi], boundary_values)
        r_bnd = boundary_values(gvec_geom_r)
        z_bnd = boundary_values(gvec_geom_z)

        dist_axis = sqrt(rs**2 + zs**2)
        dist_bnd = sqrt((r_bnd - r_axis)**2 + (z_bnd - z_axis)**2)
        if (dist_bnd > 1.0e-12_dp) then
            s = min(1.0_dp, max(self%data%geom_splines(1)%x_min(1), &
                                (dist_axis / dist_bnd)**2))
        else
            s = self%data%geom_splines(1)%x_min(1)
        end if

        do iter = 1, max_iter
            call self%data%evaluate_geometry([s, theta, phi], values)

            f(1) = values(gvec_geom_r) - xcyl(1)
            f(2) = values(gvec_geom_z) - xcyl(3)
            res_norm = sqrt(f(1)**2 + f(2)**2)
            if (res_norm < tol_res) exit

            jac(1, 1) = values(gvec_geom_dr_ds)
            jac(1, 2) = values(gvec_geom_dr_dt)
            jac(2, 1) = values(gvec_geom_dz_ds)
            jac(2, 2) = values(gvec_geom_dz_dt)
            det = jac(1, 1) * jac(2, 2) - jac(1, 2) * jac(2, 1)
            if (abs(det) < 1.0e-14_dp) then
                ierr = 2
                return
            end if

            delta(1) = (-f(1) * jac(2, 2) + f(2) * jac(1, 2)) / det
            delta(2) = (-jac(1, 1) * f(2) + jac(2, 1) * f(1)) / det

            if (maxval(abs(delta)) < tol_step) exit

            alpha = 1.0_dp
            do k = 1, 10
                s_try = min(1.0_dp, max(self%data%geom_splines(1)%x_min(1), &
                                        s + alpha * delta(1)))
                theta_try = modulo(theta + alpha * delta(2), TWOPI)
                call self%data%evaluate_geometry([s_try, theta_try, phi], values)
                res_norm_try = sqrt((values(gvec_geom_r) - xcyl(1))**2 + &
                                    (values(gvec_geom_z) - xcyl(3))**2)
                if (res_norm_try < res_norm) then
                    s = s_try
                    theta = theta_try
                    exit
                end if
                alpha = 0.5_dp * alpha
            end do

            if (k > 10) then
                ierr = 1
                return
            end if
        end do

        if (iter > max_iter) then
            ierr = 1
            return
        end if

        u = [s, theta, phi]
    end subroutine gvec_from_cyl

    real(dp) function gvec_coordinate_phi_period(self)
        class(gvec_coordinate_system_t), intent(in) :: self

        gvec_coordinate_phi_period = self%data%phi_period()
    end function gvec_coordinate_phi_period

end module gvec_reference_coordinates
