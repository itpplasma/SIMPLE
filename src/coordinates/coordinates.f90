module simple_coordinates

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vmec_coordinates, only: vmec_to_cyl_lib => vmec_to_cyl, &
                                 vmec_to_cart_lib => vmec_to_cart, &
                                 cyl_to_cart_lib  => cyl_to_cart

    implicit none

    abstract interface
        subroutine transform_i(xfrom, xto, dxto_dxfrom)
            import :: dp
            real(dp), intent(in) :: xfrom(3)
            real(dp), intent(out) :: xto(3)
            real(dp), intent(out), optional :: dxto_dxfrom(3,3)
        end subroutine transform_i
    end interface

contains

subroutine transform_vmec_to_cyl(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call vmec_to_cyl_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_vmec_to_cyl


subroutine transform_vmec_to_cart(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call vmec_to_cart_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_vmec_to_cart


subroutine transform_cyl_to_cart(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call cyl_to_cart_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_cyl_to_cart


subroutine transform_cart_to_vmec(xfrom, xto, dxto_dxfrom)
    ! Inverse mapping of vmec_to_cart: find (s, theta, phi) such that
    ! vmec_to_cart_lib((s, theta, phi)) = xfrom using Newton iterations.
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    real(dp), parameter :: tol = 1.0d-10
    integer, parameter :: max_iter = 20

    real(dp) :: q(3), x_eval(3), J(3,3), rhs(3), Jinv(3,3), delta(3)
    real(dp) :: res_norm, rxy
    integer :: iter

    ! Initial guess for (s, theta, phi)
    rxy = sqrt(xfrom(1)**2 + xfrom(2)**2)
    q(1) = 0.5d0
    q(2) = atan2(xfrom(3), max(rxy, 1.0d-8))
    q(3) = atan2(xfrom(2), xfrom(1))

    do iter = 1, max_iter
        call vmec_to_cart_lib(q, x_eval, J)

        rhs = x_eval - xfrom
        res_norm = sqrt(rhs(1)**2 + rhs(2)**2 + rhs(3)**2)
        if (res_norm < tol) exit

        call invert_3x3(J, Jinv)
        delta = matmul(Jinv, rhs)
        q = q - delta

        q(1) = max(0.0d0, min(1.0d0, q(1)))
        q(2) = modulo(q(2), 2.0d0*3.14159265358979d0)
        q(3) = modulo(q(3), 2.0d0*3.14159265358979d0)
    end do

    xto = q

    if (present(dxto_dxfrom)) then
        dxto_dxfrom = Jinv
    end if
end subroutine transform_cart_to_vmec


subroutine invert_3x3(a, a_inv)
    real(dp), intent(in) :: a(3,3)
    real(dp), intent(out) :: a_inv(3,3)
    real(dp) :: det

    det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - &
          a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) + &
          a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

    if (abs(det) < 1.0d-20) then
        print *, "invert_3x3: singular matrix"
        error stop
    end if

    a_inv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))/det
    a_inv(1,2) = (a(1,3)*a(3,2)-a(1,2)*a(3,3))/det
    a_inv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))/det

    a_inv(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))/det
    a_inv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))/det
    a_inv(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))/det

    a_inv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))/det
    a_inv(3,2) = (a(1,2)*a(3,1)-a(1,1)*a(3,2))/det
    a_inv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))/det
end subroutine invert_3x3

function get_transform(from, to)
    procedure(transform_i), pointer :: get_transform
    character(*), intent(in) :: from, to

    get_transform => null()

    select case (trim(from))
    case('cyl')
        select case (trim(to))
        case ('cart')
            get_transform => transform_cyl_to_cart
        case default
            call handle_transform_error(from, to)
        end select
    case ('vmec')
        select case (trim(to))
        case ('cart')
            get_transform => transform_vmec_to_cart
        case ('cyl')
            get_transform => transform_vmec_to_cyl
        case default
            call handle_transform_error(from, to)
        end select
    case default
        print *, "get_transform: Unknown transform from ", from
        error stop
    end select
end function get_transform


subroutine handle_transform_error(from, to)
    character(*), intent(in) :: from
    character(*), intent(in), optional :: to

    if (present(to)) then
        print *, "Unknown transform from ", from, " to ", to
    else
        print *, "Unknown transform from ", from
    end if
    error stop
end subroutine handle_transform_error

end module simple_coordinates
