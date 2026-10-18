program test_chartmap_meiss_debug
    !> A chartmap whose radial coordinate is already rho = sqrt(s) must select
    !> the identity Meiss scaling; sqrt_s scaling would transform it twice.
    use libneo_coordinates, only: coordinate_system_t, make_chartmap_coordinate_system
    use field_splined, only: splined_field_t
    use field_can_meiss, only: choose_default_scaling
    use coordinate_scaling, only: coordinate_scaling_t, sqrt_s_scaling_t, &
                                  identity_scaling_t
    implicit none

    character(len=512) :: chartmap_file
    class(coordinate_system_t), allocatable :: coords
    integer :: n_tests_failed

    n_tests_failed = 0

    call get_command_argument(1, chartmap_file)
    if (len_trim(chartmap_file) == 0) then
        print *, "Usage: test_chartmap_meiss_debug <chartmap.nc>"
        stop 1
    end if

    print *, "Testing chartmap: ", trim(chartmap_file)
    call make_chartmap_coordinate_system(coords, trim(chartmap_file))

    call test_meiss_scaling_selection()

    if (n_tests_failed > 0) stop 1
    print *, "PASS: chartmap Meiss scaling selection"

contains

    subroutine test_meiss_scaling_selection()
        class(coordinate_scaling_t), allocatable :: scaling
        type(splined_field_t) :: dummy_field

        allocate (dummy_field%coords, source=coords)

        call choose_default_scaling(dummy_field, scaling)

        select type (s => scaling)
        type is (identity_scaling_t)
            print *, "  identity scaling selected (correct for rho=sqrt(s))"
        type is (sqrt_s_scaling_t)
            print *, "  FAIL: sqrt_s scaling would double-transform"
            n_tests_failed = n_tests_failed + 1
        class default
            print *, "  FAIL: unknown scaling type"
            n_tests_failed = n_tests_failed + 1
        end select
    end subroutine test_meiss_scaling_selection

end program test_chartmap_meiss_debug
