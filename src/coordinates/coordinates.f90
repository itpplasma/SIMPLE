module simple_coordinates

use, intrinsic :: iso_fortran_env, only: dp => real64

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
    use spline_vmec_sub, only: splint_vmec_data
    use spline_vmec_libneo, only: splint_vmec_data_libneo

    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                dl_ds, dl_dt, dl_dp
    
    logical, parameter :: use_libneo = .true.

    if (use_libneo) then
        call splint_vmec_data_libneo(xfrom(1), xfrom(2), xfrom(3), &
            A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, &
            dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
    else
        call splint_vmec_data(xfrom(1), xfrom(2), xfrom(3), &
            A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, &
            dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
    end if

    xto(1) = R
    xto(2) = xfrom(3)
    xto(3) = Z

    if (present(dxto_dxfrom)) then
        dxto_dxfrom(1,1) = dR_ds
        dxto_dxfrom(1,2) = dR_dt
        dxto_dxfrom(1,3) = dR_dp
        dxto_dxfrom(2,1) = 0d0
        dxto_dxfrom(2,2) = 0d0
        dxto_dxfrom(2,3) = 1d0
        dxto_dxfrom(3,1) = dZ_ds
        dxto_dxfrom(3,2) = dZ_dt
        dxto_dxfrom(3,3) = dZ_dp
    end if
end subroutine transform_vmec_to_cyl


subroutine transform_vmec_to_cart(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    real(dp) :: xcyl(3), dxcyl_dxvmec(3,3), dxcart_dxcyl(3,3)

    if (present(dxto_dxfrom)) then
        call transform_vmec_to_cyl(xfrom, xcyl, dxcyl_dxvmec)
        call transform_cyl_to_cart(xcyl, xto, dxcart_dxcyl)
        dxto_dxfrom = matmul(dxcart_dxcyl, dxcyl_dxvmec)
    else
        call transform_vmec_to_cyl(xfrom, xcyl)
        call transform_cyl_to_cart(xcyl, xto)
    end if
end subroutine transform_vmec_to_cart


subroutine transform_cyl_to_cart(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    xto(1) = xfrom(1)*cos(xfrom(2))
    xto(2) = xfrom(1)*sin(xfrom(2))
    xto(3) = xfrom(3)

    if (present(dxto_dxfrom)) then
        dxto_dxfrom(1,1) = cos(xfrom(2))
        dxto_dxfrom(1,2) = -xfrom(1)*sin(xfrom(2))
        dxto_dxfrom(1,3) = 0d0
        dxto_dxfrom(2,1) = sin(xfrom(2))
        dxto_dxfrom(2,2) = xfrom(1)*cos(xfrom(2))
        dxto_dxfrom(2,3) = 0d0
        dxto_dxfrom(3,1) = 0d0
        dxto_dxfrom(3,2) = 0d0
        dxto_dxfrom(3,3) = 1d0
    end if
end subroutine transform_cyl_to_cart


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
