module coordinate_scaling
    !> Bijective coordinate scaling transforms.
    !> Single responsibility: transform coordinates before/after operations
    !> like splining, with explicit forward and inverse transforms.

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type, abstract :: coordinate_scaling_t
    contains
        procedure(transform_iface), deferred :: transform
        procedure(transform_iface), deferred :: inverse
    end type coordinate_scaling_t

    abstract interface
        pure subroutine transform_iface(self, x_in, x_out, jacobian)
            !> Transform coordinates with optional Jacobian diagonal.
            !> For diagonal scaling, jacobian(i) = dx_out(i)/dx_in(i).
            import :: coordinate_scaling_t, dp
            class(coordinate_scaling_t), intent(in) :: self
            real(dp), intent(in) :: x_in(3)
            real(dp), intent(out) :: x_out(3)
            real(dp), intent(out), optional :: jacobian(3)
        end subroutine transform_iface
    end interface

    type, extends(coordinate_scaling_t) :: sqrt_s_scaling_t
        !> Default scaling: r = sqrt(s), angles unchanged.
        !> Improves spline accuracy near magnetic axis.
    contains
        procedure :: transform => sqrt_s_transform
        procedure :: inverse => sqrt_s_inverse
    end type sqrt_s_scaling_t

    type, extends(coordinate_scaling_t) :: identity_scaling_t
        !> Identity scaling: no transformation.
        !> Use when coordinate system already uses rho = sqrt(s) (RHO_TOR).
    contains
        procedure :: transform => identity_transform
        procedure :: inverse => identity_inverse
    end type identity_scaling_t

contains

    pure subroutine sqrt_s_transform(self, x_in, x_out, jacobian)
        !> Transform (s, theta, phi) -> (r=sqrt(s), theta, phi).
        !> Jacobian: dr/ds = 1/(2*sqrt(s)) = 1/(2r).
        class(sqrt_s_scaling_t), intent(in) :: self
        real(dp), intent(in) :: x_in(3)
        real(dp), intent(out) :: x_out(3)
        real(dp), intent(out), optional :: jacobian(3)

        real(dp) :: r

        r = sqrt(x_in(1))
        x_out(1) = r
        x_out(2) = x_in(2)
        x_out(3) = x_in(3)

        if (present(jacobian)) then
            if (r > 0d0) then
                jacobian(1) = 0.5d0 / r
            else
                jacobian(1) = 0d0
            end if
            jacobian(2) = 1d0
            jacobian(3) = 1d0
        end if
    end subroutine sqrt_s_transform

    pure subroutine sqrt_s_inverse(self, x_in, x_out, jacobian)
        !> Transform (r, theta, phi) -> (s=r^2, theta, phi).
        !> Jacobian: ds/dr = 2r.
        class(sqrt_s_scaling_t), intent(in) :: self
        real(dp), intent(in) :: x_in(3)
        real(dp), intent(out) :: x_out(3)
        real(dp), intent(out), optional :: jacobian(3)

        x_out(1) = x_in(1)**2
        x_out(2) = x_in(2)
        x_out(3) = x_in(3)

        if (present(jacobian)) then
            jacobian(1) = 2d0 * x_in(1)
            jacobian(2) = 1d0
            jacobian(3) = 1d0
        end if
    end subroutine sqrt_s_inverse

    pure subroutine identity_transform(self, x_in, x_out, jacobian)
        !> Identity transform: x_out = x_in.
        class(identity_scaling_t), intent(in) :: self
        real(dp), intent(in) :: x_in(3)
        real(dp), intent(out) :: x_out(3)
        real(dp), intent(out), optional :: jacobian(3)

        x_out = x_in
        if (present(jacobian)) jacobian = 1d0
    end subroutine identity_transform

    pure subroutine identity_inverse(self, x_in, x_out, jacobian)
        !> Identity inverse: x_out = x_in.
        class(identity_scaling_t), intent(in) :: self
        real(dp), intent(in) :: x_in(3)
        real(dp), intent(out) :: x_out(3)
        real(dp), intent(out), optional :: jacobian(3)

        x_out = x_in
        if (present(jacobian)) jacobian = 1d0
    end subroutine identity_inverse

end module coordinate_scaling
