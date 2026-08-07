program test_start_integration_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_mod, only: integ_to_ref
    use samplers, only: save_starting_points, sample_read_integration, START_FILE

    implicit none

    integer, parameter :: particle_count = 3
    real(dp) :: input(5, particle_count), actual(5, particle_count)
    real(dp) :: expected(5, particle_count)
    integer :: particle

    do particle = 1, particle_count
        input(:, particle) = [0.1_dp*particle, 0.2_dp*particle, &
            -0.3_dp*particle, 1.0_dp, -0.5_dp + 0.2_dp*particle]
        call known_transform(input(1:3, particle), expected(1:3, particle))
        expected(4:5, particle) = input(4:5, particle)
    end do

    integ_to_ref => known_transform
    call save_starting_points(input)
    call sample_read_integration(actual, START_FILE)

    if (maxval(abs(actual - expected)) > 1.0e-14_dp) then
        print *, 'integration-coordinate start conversion mismatch'
        error stop 1
    end if

contains

    subroutine known_transform(xinteg, xref)
        real(dp), intent(in) :: xinteg(3)
        real(dp), intent(out) :: xref(3)

        xref = [2.0_dp*xinteg(1), xinteg(2) + 1.0_dp, xinteg(3) - 1.0_dp]
    end subroutine known_transform

end program test_start_integration_coordinates
