program test_magfie_coils
    !> Test coils field evaluation comparing VMEC field, raw Biot-Savart,
    !> and splined field.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use magfie_sub, only: VMEC
    use velo_mod, only: isw_field_type
    use field_base, only: magnetic_field_t
    use field_vmec, only: vmec_field_t
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use magfie_sub, only: magfie_vmec
    use util, only: twopi
    use cylindrical_cartesian, only: cyl_to_cart

    implicit none

    type(vmec_field_t) :: vmec_field
    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: splined_coils
    real(dp) :: dummy, x(3), Acov(3), hcov(3), Bmod

    logical :: wout_exists, coils_exists

    isw_field_type = VMEC

    inquire(file='wout.nc', exist=wout_exists)
    inquire(file='coils.simple', exist=coils_exists)
    if (.not. wout_exists .or. .not. coils_exists) then
        print *, 'SKIP: Missing wout.nc or coils.simple'
        stop 0
    end if

    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates('wout.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, splined_coils)

    x = [0.3d0, 0.2d0, 0.1d0]

    print *, 'x = ', x

    call vmec_field%evaluate(x, Acov, hcov, Bmod)
    print *, 'vmec_field%evaluate'
    print *, 'A = ', Acov
    print *, 'h = ', hcov
    print *, 'B = ', Bmod

    call evaluate_raw_at_ref(raw_coils, x, Acov, hcov, Bmod)
    print *, 'raw_coils (Biot-Savart at ref point)'
    print *, 'A = ', Acov
    print *, 'h = ', hcov
    print *, 'B = ', Bmod

    call splined_coils%evaluate(x, Acov, hcov, Bmod)
    print *, 'splined_coils%evaluate'
    print *, 'A = ', Acov
    print *, 'h = ', hcov
    print *, 'B = ', Bmod

    call test_curve
    call test_magfie

contains

    subroutine evaluate_raw_at_ref(field, x_spline, Acov, hcov, Bmod)
        !> Helper to evaluate raw coils at a spline coordinate point.
        !> x_spline is in (r, theta, phi) where r = sqrt(s).
        !> ref_coords expects VMEC coords (s, theta, phi).
        type(coils_field_t), intent(in) :: field
        real(dp), intent(in) :: x_spline(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod

        real(dp) :: x_vmec(3), x_cyl(3), x_cart(3)

        ! Convert spline coords (r, theta, phi) to VMEC coords (s, theta, phi)
        x_vmec = [x_spline(1)**2, x_spline(2), x_spline(3)]
        call ref_coords%evaluate_point(x_vmec, x_cyl)
        call cyl_to_cart(x_cyl, x_cart)
        call field%evaluate(x_cart, Acov, hcov, Bmod)
    end subroutine evaluate_raw_at_ref


    subroutine test_curve
        real(dp) :: x(3), Acov(3), hcov(3), Bmod
        integer :: i, N = 1000

        x = [0.3d0, 0.2d0, 0.0d0]

        do i = 0, N
            x(3) = x(3) + twopi / N
            call vmec_field%evaluate(x, Acov, hcov, Bmod)
            write(1, *) x, Acov, hcov, Bmod
            call splined_coils%evaluate(x, Acov, hcov, Bmod)
            write(2, *) x, Acov, hcov, Bmod
        end do
    end subroutine test_curve


    subroutine test_magfie
        real(dp) :: bmod, sqrtg
        real(dp), dimension(3) :: bder, hcovar, hctrvr, hcurl

        call magfie_vmec(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        print *, 'magfie_vmec'
        print *, 'B = ', bmod
        print *, 'sqrtg = ', sqrtg
        print *, 'Bder = ', bder
        print *, 'hcovar = ', hcovar
        print *, 'hctrvr = ', hctrvr
        print *, 'hcurl = ', hcurl
    end subroutine test_magfie

end program test_magfie_coils
