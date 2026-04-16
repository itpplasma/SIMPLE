module scaled_chartmap_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t

    implicit none

    private
    public :: scaled_chartmap_coordinate_system_t
    public :: wrap_scaled_chartmap_coordinate_system

    type, extends(chartmap_coordinate_system_t) :: scaled_chartmap_coordinate_system_t
        real(dp) :: cart_scale = 1.0_dp
    contains
        procedure :: evaluate_cart => scaled_chartmap_evaluate_cart
        procedure :: evaluate_cyl => scaled_chartmap_evaluate_cyl
        procedure :: covariant_basis => scaled_chartmap_covariant_basis
        procedure :: metric_tensor => scaled_chartmap_metric_tensor
        procedure :: from_cyl => scaled_chartmap_from_cyl
        procedure :: from_cart => scaled_chartmap_from_cart
    end type scaled_chartmap_coordinate_system_t

contains

    subroutine wrap_scaled_chartmap_coordinate_system(coords, cart_scale)
        class(coordinate_system_t), allocatable, intent(inout) :: coords
        real(dp), intent(in) :: cart_scale

        class(coordinate_system_t), allocatable :: wrapped

        if (abs(cart_scale - 1.0_dp) <= epsilon(1.0_dp)) return
        if (cart_scale <= 0.0_dp) error stop "chartmap cart_scale must be positive"

        select type (base => coords)
        type is (scaled_chartmap_coordinate_system_t)
            base%cart_scale = base%cart_scale*cart_scale
        type is (chartmap_coordinate_system_t)
            allocate (scaled_chartmap_coordinate_system_t :: wrapped)
            select type (scaled => wrapped)
            type is (scaled_chartmap_coordinate_system_t)
                scaled%chartmap_coordinate_system_t = base
                scaled%cart_scale = cart_scale
            class default
                error stop "wrap_scaled_chartmap_coordinate_system: allocation failed"
            end select
            call move_alloc(wrapped, coords)
        class default
            error stop "wrap_scaled_chartmap_coordinate_system: expected chartmap coords"
        end select
    end subroutine wrap_scaled_chartmap_coordinate_system

    subroutine scaled_chartmap_evaluate_cart(self, u, x)
        class(scaled_chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        call self%chartmap_coordinate_system_t%evaluate_cart(u, x)
        x = self%cart_scale*x
    end subroutine scaled_chartmap_evaluate_cart

    subroutine scaled_chartmap_evaluate_cyl(self, u, x)
        class(scaled_chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        call self%chartmap_coordinate_system_t%evaluate_cyl(u, x)
        x(1) = self%cart_scale*x(1)
        x(3) = self%cart_scale*x(3)
    end subroutine scaled_chartmap_evaluate_cyl

    subroutine scaled_chartmap_covariant_basis(self, u, e_cov)
        class(scaled_chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        call self%chartmap_coordinate_system_t%covariant_basis(u, e_cov)
        e_cov = self%cart_scale*e_cov
    end subroutine scaled_chartmap_covariant_basis

    subroutine scaled_chartmap_metric_tensor(self, u, g, ginv, sqrtg)
        class(scaled_chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg

        real(dp) :: scale2

        call self%chartmap_coordinate_system_t%metric_tensor(u, g, ginv, sqrtg)
        scale2 = self%cart_scale*self%cart_scale
        g = scale2*g
        ginv = ginv/scale2
        sqrtg = self%cart_scale*scale2*sqrtg
    end subroutine scaled_chartmap_metric_tensor

    subroutine scaled_chartmap_from_cyl(self, xcyl, u, ierr)
        class(scaled_chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp) :: xcyl_base(3)

        xcyl_base = xcyl
        xcyl_base(1) = xcyl_base(1)/self%cart_scale
        xcyl_base(3) = xcyl_base(3)/self%cart_scale
        call self%chartmap_coordinate_system_t%from_cyl(xcyl_base, u, ierr)
    end subroutine scaled_chartmap_from_cyl

    subroutine scaled_chartmap_from_cart(self, x, u, ierr)
        class(scaled_chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp) :: x_base(3)

        x_base = x/self%cart_scale
        call self%chartmap_coordinate_system_t%from_cart(x_base, u, ierr)
    end subroutine scaled_chartmap_from_cart

end module scaled_chartmap_coordinates
