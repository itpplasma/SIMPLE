program test_coordinate_scaling
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use coordinate_scaling, only: sqrt_s_scaling_t
    implicit none

    integer :: errors

    errors = 0

    call test_sqrt_s_roundtrip(errors)
    call test_sqrt_s_jacobian_product(errors)
    call test_sqrt_s_specific_values(errors)

    if (errors == 0) then
        print *, "All coordinate_scaling tests passed!"
    else
        print *, "ERROR: ", errors, " test(s) failed!"
        error stop 1
    end if

contains

    subroutine test_sqrt_s_roundtrip(errors)
        integer, intent(inout) :: errors
        type(sqrt_s_scaling_t) :: scaling
        real(dp) :: x_orig(3), x_scaled(3), x_back(3)
        real(dp), parameter :: tol = 1.0d-14
        integer :: i

        print *, "Testing sqrt_s_scaling roundtrip..."

        do i = 1, 5
            x_orig = [0.2d0 * i, 1.5d0, 0.8d0]

            call scaling%transform(x_orig, x_scaled)
            call scaling%inverse(x_scaled, x_back)

            if (abs(x_back(1) - x_orig(1)) > tol) then
                print *, "ERROR: Roundtrip failed for s =", x_orig(1)
                print *, "  Original:", x_orig(1), "Got:", x_back(1)
                errors = errors + 1
            end if

            if (abs(x_back(2) - x_orig(2)) > tol) then
                print *, "ERROR: Theta not preserved in roundtrip"
                errors = errors + 1
            end if

            if (abs(x_back(3) - x_orig(3)) > tol) then
                print *, "ERROR: Phi not preserved in roundtrip"
                errors = errors + 1
            end if
        end do

        if (errors == 0) then
            print *, "  Roundtrip test PASSED"
        end if
    end subroutine test_sqrt_s_roundtrip

    subroutine test_sqrt_s_jacobian_product(errors)
        integer, intent(inout) :: errors
        type(sqrt_s_scaling_t) :: scaling
        real(dp) :: x_s(3), x_r(3), x_back(3)
        real(dp) :: jac_fwd(3), jac_inv(3)
        real(dp), parameter :: tol = 1.0d-14

        print *, "Testing sqrt_s_scaling Jacobian chain rule..."

        x_s = [0.64d0, 1.0d0, 2.0d0]

        call scaling%transform(x_s, x_r, jac_fwd)
        call scaling%inverse(x_r, x_back, jac_inv)

        if (abs(jac_fwd(1) * jac_inv(1) - 1.0d0) > tol) then
            print *, "ERROR: Jacobian product should be 1"
            print *, "  dr/ds =", jac_fwd(1), "ds/dr =", jac_inv(1)
            print *, "  Product:", jac_fwd(1) * jac_inv(1)
            errors = errors + 1
        end if

        if (abs(jac_fwd(2) - 1.0d0) > tol .or. abs(jac_inv(2) - 1.0d0) > tol) then
            print *, "ERROR: Angular Jacobians should be 1"
            errors = errors + 1
        end if

        if (abs(jac_fwd(3) - 1.0d0) > tol .or. abs(jac_inv(3) - 1.0d0) > tol) then
            print *, "ERROR: Angular Jacobians should be 1"
            errors = errors + 1
        end if

        if (errors == 0) then
            print *, "  Jacobian chain rule test PASSED"
        end if
    end subroutine test_sqrt_s_jacobian_product

    subroutine test_sqrt_s_specific_values(errors)
        integer, intent(inout) :: errors
        type(sqrt_s_scaling_t) :: scaling
        real(dp) :: x_in(3), x_out(3), jac(3)
        real(dp), parameter :: tol = 1.0d-14

        print *, "Testing sqrt_s_scaling specific values..."

        x_in = [0.25d0, 1.0d0, 2.0d0]
        call scaling%transform(x_in, x_out, jac)

        if (abs(x_out(1) - 0.5d0) > tol) then
            print *, "ERROR: sqrt(0.25) should be 0.5"
            print *, "  Got:", x_out(1)
            errors = errors + 1
        end if

        if (abs(jac(1) - 1.0d0) > tol) then
            print *, "ERROR: dr/ds at s=0.25 should be 1/(2*0.5) = 1"
            print *, "  Got:", jac(1)
            errors = errors + 1
        end if

        x_in = [0.5d0, 1.0d0, 2.0d0]
        call scaling%inverse(x_in, x_out, jac)

        if (abs(x_out(1) - 0.25d0) > tol) then
            print *, "ERROR: 0.5^2 should be 0.25"
            print *, "  Got:", x_out(1)
            errors = errors + 1
        end if

        if (abs(jac(1) - 1.0d0) > tol) then
            print *, "ERROR: ds/dr at r=0.5 should be 2*0.5 = 1"
            print *, "  Got:", jac(1)
            errors = errors + 1
        end if

        if (errors == 0) then
            print *, "  Specific values test PASSED"
        end if
    end subroutine test_sqrt_s_specific_values

end program test_coordinate_scaling
