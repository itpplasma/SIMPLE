program test_scaling_abstraction
    !> Test that Meiss and Albert coordinates work with different scaling options.
    !>
    !> Verifies:
    !>   1. Default (sqrt_s_scaling_t) works as before
    !>   2. identity_scaling_t works when reference coords already use r (not s)
    !>   3. Transform roundtrips are correct for both scalings

    use, intrinsic :: iso_fortran_env, only: dp => real64

    use coordinate_scaling, only: sqrt_s_scaling_t, identity_scaling_t
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, cleanup_meiss, &
        integ_to_ref_meiss, ref_to_integ_meiss, twopi
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
    print *, '=== Testing default sqrt_s_scaling_t ==='
    call test_sqrt_s_scaling(n_failed)

    print *, ''
    print *, '=== Testing identity_scaling_t ==='
    call test_identity_scaling(n_failed)

    print *, ''
    if (n_failed == 0) then
        print *, 'PASSED: All scaling abstraction tests passed'
    else
        print *, 'FAILED:', n_failed, 'test(s) failed'
        error stop 'Scaling abstraction tests failed'
    end if

contains

    subroutine test_sqrt_s_scaling(n_failed)
        !> Test default sqrt_s_scaling_t produces r = sqrt(s).
        integer, intent(inout) :: n_failed
        real(dp), parameter :: TOL = 1d-10
        type(sqrt_s_scaling_t) :: scaling
        real(dp) :: x_ref(3), x_integ(3), x_ref_back(3)
        real(dp) :: err(3), max_err

        print *, 'Initializing Meiss with sqrt_s_scaling_t...'
        call init_meiss(magfie, 32, 32, 32, 0.01d0, 0.99d0, 0.0d0, twopi, scaling)
        call get_meiss_coordinates()

        ! Test point: s=0.25 should give r=0.5
        x_ref = [0.25d0, 1.0d0, 0.5d0]
        call ref_to_integ_meiss(x_ref, x_integ)

        if (abs(x_integ(1) - 0.5d0) > TOL) then
            print *, '  FAILED: s=0.25 should give r=0.5'
            print *, '    Expected r = 0.5'
            print *, '    Got r =', x_integ(1)
            n_failed = n_failed + 1
        else
            print *, '  PASSED: s=0.25 -> r=0.5 conversion correct'
        end if

        ! Test roundtrip
        call integ_to_ref_meiss(x_integ, x_ref_back)
        err(1) = abs(x_ref(1) - x_ref_back(1))
        err(2) = abs(mod(x_ref(2) - x_ref_back(2) + twopi, twopi))
        if (err(2) > twopi/2) err(2) = twopi - err(2)
        err(3) = abs(mod(x_ref(3) - x_ref_back(3) + twopi, twopi))
        if (err(3) > twopi/2) err(3) = twopi - err(3)
        max_err = maxval(err)

        if (max_err > TOL) then
            print *, '  FAILED: roundtrip error too large'
            print *, '    x_ref      =', x_ref
            print *, '    x_ref_back =', x_ref_back
            print *, '    error      =', err
            n_failed = n_failed + 1
        else
            print *, '  PASSED: roundtrip with sqrt_s_scaling (max_err =', max_err, ')'
        end if

        call cleanup_meiss()
    end subroutine test_sqrt_s_scaling


    subroutine test_identity_scaling(n_failed)
        !> Test identity_scaling_t produces r = s (no scaling).
        !> This would be used when reference coordinates already use rho_tor.
        integer, intent(inout) :: n_failed
        real(dp), parameter :: TOL = 1d-10
        type(identity_scaling_t) :: scaling
        real(dp) :: x_ref(3), x_integ(3), x_ref_back(3)
        real(dp) :: err(3), max_err

        print *, 'Initializing Meiss with identity_scaling_t...'
        call init_meiss(magfie, 32, 32, 32, 0.01d0, 0.99d0, 0.0d0, twopi, scaling)
        call get_meiss_coordinates()

        ! Test point: s=0.25 should give r=0.25 (no scaling)
        x_ref = [0.25d0, 1.0d0, 0.5d0]
        call ref_to_integ_meiss(x_ref, x_integ)

        if (abs(x_integ(1) - 0.25d0) > TOL) then
            print *, '  FAILED: s=0.25 should give r=0.25 (identity)'
            print *, '    Expected r = 0.25'
            print *, '    Got r =', x_integ(1)
            n_failed = n_failed + 1
        else
            print *, '  PASSED: s=0.25 -> r=0.25 identity conversion correct'
        end if

        ! Test roundtrip
        call integ_to_ref_meiss(x_integ, x_ref_back)
        err(1) = abs(x_ref(1) - x_ref_back(1))
        err(2) = abs(mod(x_ref(2) - x_ref_back(2) + twopi, twopi))
        if (err(2) > twopi/2) err(2) = twopi - err(2)
        err(3) = abs(mod(x_ref(3) - x_ref_back(3) + twopi, twopi))
        if (err(3) > twopi/2) err(3) = twopi - err(3)
        max_err = maxval(err)

        if (max_err > TOL) then
            print *, '  FAILED: roundtrip error too large'
            print *, '    x_ref      =', x_ref
            print *, '    x_ref_back =', x_ref_back
            print *, '    error      =', err
            n_failed = n_failed + 1
        else
            print *, '  PASSED: roundtrip with identity_scaling (max_err =', max_err, ')'
        end if

        call cleanup_meiss()
    end subroutine test_identity_scaling

end program test_scaling_abstraction
