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


subroutine transform_cyl_to_vmec(x_cyl, x_vmec, ierr)
    use spline_vmec_sub, only: splint_vmec_data

    real(dp), intent(in) :: x_cyl(3)
    real(dp), intent(out) :: x_vmec(3)
    integer, intent(out) :: ierr

    real(dp) :: R_target, phi_target, Z_target
    real(dp) :: s, theta, phi
    real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
    real(dp) :: R_eval, Z_eval, alam
    real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
    real(dp) :: dl_ds, dl_dt, dl_dp
    real(dp) :: f(2), jac(2,2), delta(2), det
    real(dp) :: res_norm
    integer :: iter
    integer, parameter :: MAX_ITER = 50
    real(dp), parameter :: TOL = 1d-10

    ierr = 0

    R_target = x_cyl(1)
    phi_target = x_cyl(2)
    Z_target = x_cyl(3)

    phi = phi_target
    s = 0.5d0
    theta = 0d0

    do iter = 1, MAX_ITER
        call splint_vmec_data(s, theta, phi, A_phi, A_theta, dA_phi_ds, &
            dA_theta_ds, aiota, R_eval, Z_eval, alam, dR_ds, dR_dt, dR_dp, &
            dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

        f(1) = R_eval - R_target
        f(2) = Z_eval - Z_target

        res_norm = sqrt(f(1)**2 + f(2)**2)
        if (res_norm < TOL) exit

        jac(1,1) = dR_ds
        jac(1,2) = dR_dt
        jac(2,1) = dZ_ds
        jac(2,2) = dZ_dt

        det = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)
        if (abs(det) < 1d-30) then
            ierr = 2
            return
        endif

        delta(1) = (jac(2,2) * f(1) - jac(1,2) * f(2)) / det
        delta(2) = (-jac(2,1) * f(1) + jac(1,1) * f(2)) / det

        s = s - delta(1)
        theta = theta - delta(2)

        s = max(0d0, min(1d0, s))
        theta = modulo(theta, 2d0 * 3.141592653589793d0)
    enddo

    if (iter > MAX_ITER) then
        ierr = 1
        return
    endif

    x_vmec(1) = sqrt(s)
    x_vmec(2) = theta
    x_vmec(3) = phi

end subroutine transform_cyl_to_vmec


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
