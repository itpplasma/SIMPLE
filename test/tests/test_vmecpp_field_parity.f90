program test_vmecpp_field_parity
    use field_vmecpp, only: evaluate_vmecpp_geometry, magfie_vmecpp, &
        open_vmecpp_geometry, vmecpp_geometry_point_t
    use magfie_sub, only: magfie_vmec
    use simple, only: init_vmec
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: relative_tolerance = 2.5e-1_dp
    character(1024) :: input_file, wout_file
    real(dp) :: fper, x(3), direct(14), reference(14), error
    type(vmecpp_geometry_point_t) :: point
    integer :: i

    call get_command_argument(1, wout_file)
    call get_command_argument(2, input_file)
    if (len_trim(wout_file) == 0 .or. len_trim(input_file) == 0) &
        error stop 'expected wout and VMEC++ input paths'

    call init_vmec(trim(wout_file), 5, 5, 3, fper)
    call open_vmecpp_geometry(trim(input_file))
    do i = 1, 3
        x = [0.2_dp+0.2_dp*i, 0.31_dp*i, 0.17_dp*i]
        call evaluate_vmecpp_geometry(x(1), x(2), x(3), point)
        x(2) = x(2) + point%lambda(1)
        call evaluate(magfie_vmec, x, reference)
        x(2) = x(2) - point%lambda(1)
        call evaluate(magfie_vmecpp, x, direct)
        error = maxval(abs(direct(1:2)-reference(1:2))/ &
            max(abs(reference(1:2)), 1.0_dp))
        error = max(error, maxval(abs(direct(10:11)-reference(10:11))/ &
            max(abs(reference(10:11)), 1.0e-3_dp)))
        if (error > relative_tolerance) then
            print *, 'relative error:', error
            print *, 'wout:  ', reference
            print *, 'direct:', direct
            error stop 'VMEC++ direct field does not match wout field'
        end if
    end do

contains

    subroutine evaluate(field_function, position, values)
        interface
            subroutine field_function(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
                import :: dp
                real(dp), intent(in) :: x(3)
                real(dp), intent(out) :: bmod, sqrtg
                real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
            end subroutine field_function
        end interface
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: values(14)

        call field_function(position, values(1), values(2), values(3:5), &
            values(6:8), values(9:11), values(12:14))
    end subroutine evaluate

end program test_vmecpp_field_parity
