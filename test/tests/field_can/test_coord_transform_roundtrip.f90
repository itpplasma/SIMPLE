program test_coord_transform_roundtrip
    !> Test coordinate transform roundtrips for Meiss and Albert canonical coordinates.
    !>
    !> Verifies that:
    !>   x_ref == integ_to_ref(ref_to_integ(x_ref))
    !> for multiple test points across the domain.
    !>
    !> The transforms handle:
    !>   - s (flux) <-> r = sqrt(s) scaling for radial coordinate
    !>   - theta modulo 2*pi
    !>   - phi + lambda (canonical angle shift for Meiss)

    use, intrinsic :: iso_fortran_env, only: dp => real64

    use magfie_sub, only: MEISS, ALBERT
    use field_can_mod, only: init_field_can, integ_to_ref, ref_to_integ
    use field_can_meiss, only: integ_to_ref_meiss, ref_to_integ_meiss, twopi
    use field_can_albert, only: integ_to_ref_albert, ref_to_integ_albert
    use simple, only: init_vmec
    use field, only: vmec_field_t

    implicit none

    integer :: n_failed
    real(dp) :: fper
    class(vmec_field_t), allocatable :: magfie

    n_failed = 0

    print *, 'Initializing VMEC...'
    call init_vmec('wout.nc', 5, 5, 5, fper)

    magfie = vmec_field_t()

    print *, ''
    print *, '=== Testing Meiss coordinate transforms ==='
    call test_meiss_roundtrip(n_failed)

    print *, ''
    print *, '=== Testing Albert coordinate transforms ==='
    call test_albert_roundtrip(n_failed)

    print *, ''
    if (n_failed == 0) then
        print *, 'PASSED: All coordinate transform roundtrip tests passed'
    else
        print *, 'FAILED:', n_failed, 'test(s) failed'
        error stop 'Coordinate transform roundtrip tests failed'
    end if

contains

    subroutine test_meiss_roundtrip(n_failed)
        use field_can_meiss, only: init_meiss, get_meiss_coordinates, cleanup_meiss

        integer, intent(inout) :: n_failed
        real(dp), parameter :: TOL = 1d-10
        real(dp) :: test_points(3, 5)
        real(dp) :: x_ref(3), x_integ(3), x_ref_back(3)
        real(dp) :: err(3), max_err
        integer :: i
        character(len=32) :: point_name

        print *, 'Initializing Meiss coordinates...'
        call init_meiss(magfie, 32, 32, 32, 0.01d0, 0.99d0, 0.0d0, twopi)
        call get_meiss_coordinates()

        ! Test points: (s, theta, phi) in reference coordinates
        ! s ranges from 0 to 1 (normalized flux)
        test_points(:, 1) = [0.1d0, 0.5d0, 0.3d0]       ! Near axis
        test_points(:, 2) = [0.25d0, 1.0d0, 1.0d0]      ! Inner region
        test_points(:, 3) = [0.5d0, 3.14d0, 0.5d0]      ! Mid-radius
        test_points(:, 4) = [0.7d0, 5.0d0, 1.2d0]       ! Outer region
        test_points(:, 5) = [0.9d0, 0.1d0, 0.1d0]       ! Near edge

        do i = 1, 5
            select case (i)
            case (1); point_name = 'near axis'
            case (2); point_name = 'inner region'
            case (3); point_name = 'mid-radius'
            case (4); point_name = 'outer region'
            case (5); point_name = 'near edge'
            end select

            x_ref = test_points(:, i)

            call ref_to_integ_meiss(x_ref, x_integ)
            call integ_to_ref_meiss(x_integ, x_ref_back)

            ! For s coordinate, direct comparison
            err(1) = abs(x_ref(1) - x_ref_back(1))
            ! For angular coordinates, use modulo 2*pi
            err(2) = abs(mod(x_ref(2) - x_ref_back(2) + twopi, twopi))
            if (err(2) > twopi/2) err(2) = twopi - err(2)
            err(3) = abs(mod(x_ref(3) - x_ref_back(3) + twopi, twopi))
            if (err(3) > twopi/2) err(3) = twopi - err(3)

            max_err = maxval(err)

            if (max_err > TOL) then
                print *, '  FAILED at ', trim(point_name), ':'
                print *, '    x_ref      =', x_ref
                print *, '    x_integ    =', x_integ
                print *, '    x_ref_back =', x_ref_back
                print *, '    error      =', err
                n_failed = n_failed + 1
            else
                print *, '  PASSED:', trim(point_name), ' (max_err =', max_err, ')'
            end if
        end do

        ! Test s <-> r = sqrt(s) conversion explicitly
        call test_meiss_s_r_conversion(n_failed)

        call cleanup_meiss()
    end subroutine test_meiss_roundtrip


    subroutine test_meiss_s_r_conversion(n_failed)
        !> Verify the s <-> r = sqrt(s) scaling is applied correctly.
        integer, intent(inout) :: n_failed
        real(dp), parameter :: TOL = 1d-12
        real(dp) :: s_test, r_expected, r_actual
        real(dp) :: x_ref(3), x_integ(3)

        print *, '  Testing s <-> r = sqrt(s) conversion...'

        ! Test: s = 0.25 should give r = 0.5
        s_test = 0.25d0
        r_expected = 0.5d0
        x_ref = [s_test, 1.0d0, 0.5d0]

        call ref_to_integ_meiss(x_ref, x_integ)
        r_actual = x_integ(1)

        if (abs(r_actual - r_expected) > TOL) then
            print *, '    FAILED: s=0.25 should give r=0.5'
            print *, '      Expected r =', r_expected
            print *, '      Got r      =', r_actual
            n_failed = n_failed + 1
        else
            print *, '    PASSED: s=0.25 -> r=0.5 conversion correct'
        end if

        ! Test: s = 0.81 should give r = 0.9
        s_test = 0.81d0
        r_expected = 0.9d0
        x_ref = [s_test, 1.0d0, 0.5d0]

        call ref_to_integ_meiss(x_ref, x_integ)
        r_actual = x_integ(1)

        if (abs(r_actual - r_expected) > TOL) then
            print *, '    FAILED: s=0.81 should give r=0.9'
            print *, '      Expected r =', r_expected
            print *, '      Got r      =', r_actual
            n_failed = n_failed + 1
        else
            print *, '    PASSED: s=0.81 -> r=0.9 conversion correct'
        end if
    end subroutine test_meiss_s_r_conversion


    subroutine test_albert_roundtrip(n_failed)
        !> Albert coordinates use spline interpolation for the psi <-> r mapping:
        !>   ref_to_integ: s -> r -> Ath(r,th,ph) -> psi = Ath/Ath_norm
        !>   integ_to_ref: psi -> r(psi,th,ph) via spline -> s = r^2
        !> This introduces inherent numerical error from spline interpolation.
        !> A tolerance of 1e-3 is appropriate for the 32^3 grid resolution.
        use field_can_meiss, only: init_meiss, cleanup_meiss
        use field_can_albert, only: get_albert_coordinates

        integer, intent(inout) :: n_failed
        real(dp), parameter :: TOL = 1d-3  ! Relaxed for spline interpolation
        real(dp) :: test_points(3, 5)
        real(dp) :: x_ref(3), x_integ(3), x_ref_back(3)
        real(dp) :: err(3), max_err
        integer :: i
        character(len=32) :: point_name

        print *, 'Initializing Albert coordinates...'
        call init_meiss(magfie, 32, 32, 32, 0.01d0, 0.99d0, 0.0d0, twopi)
        call get_albert_coordinates()

        ! Test points: (s, theta, phi) in reference coordinates
        test_points(:, 1) = [0.1d0, 0.5d0, 0.3d0]       ! Near axis
        test_points(:, 2) = [0.25d0, 1.0d0, 1.0d0]      ! Inner region
        test_points(:, 3) = [0.5d0, 3.14d0, 0.5d0]      ! Mid-radius
        test_points(:, 4) = [0.7d0, 5.0d0, 1.2d0]       ! Outer region
        test_points(:, 5) = [0.9d0, 0.1d0, 0.1d0]       ! Near edge

        do i = 1, 5
            select case (i)
            case (1); point_name = 'near axis'
            case (2); point_name = 'inner region'
            case (3); point_name = 'mid-radius'
            case (4); point_name = 'outer region'
            case (5); point_name = 'near edge'
            end select

            x_ref = test_points(:, i)

            call ref_to_integ_albert(x_ref, x_integ)
            call integ_to_ref_albert(x_integ, x_ref_back)

            ! For s coordinate, direct comparison
            err(1) = abs(x_ref(1) - x_ref_back(1))
            ! For angular coordinates, use modulo 2*pi
            err(2) = abs(mod(x_ref(2) - x_ref_back(2) + twopi, twopi))
            if (err(2) > twopi/2) err(2) = twopi - err(2)
            err(3) = abs(mod(x_ref(3) - x_ref_back(3) + twopi, twopi))
            if (err(3) > twopi/2) err(3) = twopi - err(3)

            max_err = maxval(err)

            if (max_err > TOL) then
                print *, '  FAILED at ', trim(point_name), ':'
                print *, '    x_ref      =', x_ref
                print *, '    x_integ    =', x_integ
                print *, '    x_ref_back =', x_ref_back
                print *, '    error      =', err
                n_failed = n_failed + 1
            else
                print *, '  PASSED:', trim(point_name), ' (max_err =', max_err, ')'
            end if
        end do

        call cleanup_meiss()
    end subroutine test_albert_roundtrip

end program test_coord_transform_roundtrip
