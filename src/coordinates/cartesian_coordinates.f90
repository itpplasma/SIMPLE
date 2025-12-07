module cartesian_coordinates
    !> Cartesian coordinate system implementation.
    !> Single responsibility: identity coordinate system where all transforms
    !> are trivial (covariant basis = identity matrix).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t
    use cylindrical_cartesian, only: cart_to_cyl, cyl_to_cart

    implicit none

    type, extends(coordinate_system_t) :: cartesian_coordinate_system_t
    contains
        procedure :: evaluate_point => cartesian_evaluate_point
        procedure :: covariant_basis => cartesian_covariant_basis
        procedure :: metric_tensor => cartesian_metric_tensor
        procedure :: from_cyl => cartesian_from_cyl
    end type cartesian_coordinate_system_t

contains

    subroutine cartesian_evaluate_point(self, u, x)
        !> Convert Cartesian to cylindrical.
        class(cartesian_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        call cart_to_cyl(u, x)
    end subroutine cartesian_evaluate_point

    subroutine cartesian_covariant_basis(self, u, e_cov)
        !> Covariant basis for Cartesian: identity matrix.
        !> e_cov(i,j) = dx_cart_i / dx_cart_j = delta_ij
        class(cartesian_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        e_cov = 0d0
        e_cov(1, 1) = 1d0
        e_cov(2, 2) = 1d0
        e_cov(3, 3) = 1d0
    end subroutine cartesian_covariant_basis

    subroutine cartesian_metric_tensor(self, u, g, ginv, sqrtg)
        !> Metric tensor for Cartesian: identity.
        class(cartesian_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg

        g = 0d0
        g(1, 1) = 1d0
        g(2, 2) = 1d0
        g(3, 3) = 1d0

        ginv = g
        sqrtg = 1d0
    end subroutine cartesian_metric_tensor

    subroutine cartesian_from_cyl(self, xcyl, u, ierr)
        !> Convert cylindrical to Cartesian.
        class(cartesian_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        call cyl_to_cart(xcyl, u)
        ierr = 0
    end subroutine cartesian_from_cyl

end module cartesian_coordinates
