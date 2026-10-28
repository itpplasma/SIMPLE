! Share immutable field/wall initialization; run every case through the same
! production routines as simple.x, rebuilding parameters and collisions per case.
program test_collisional_stl_wall_fast
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use simple, only: tracer_t
    use simple_main, only: load_config, init_run, init_sampling_field, &
        init_collisions, sample_start, init_tracing_field, &
        init_counters, trace_parallel, write_output
    use params, only: params_init, swcoll, ntestpart
    use simple_main, only: wall
    use netcdf_results_output, only: write_results_netcdf
    use diag_counters, only: diag_counters_init
    use timing, only: init_timer, get_wtime
    implicit none

    interface
        function change_directory(path) bind(C, name='chdir') result(status)
            import c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function change_directory
    end interface

    type(tracer_t) :: tracer
    character(len=4096) :: directory
    integer :: i
    real(dp) :: start
    real(dp), allocatable :: weights(:)

    call init_timer
    call load_config('simple.in')
    call init_run(tracer)
    call init_sampling_field
    do i = 1, command_argument_count()
        call get_command_argument(i, directory)
        if (change_directory(trim(directory)//c_null_char) /= 0) &
            error stop 'cannot enter wall-test case directory'
        start = get_wtime()
        call load_config('simple.in')
        call params_init
        if (swcoll) call init_collisions
        call sample_start
        call init_tracing_field
        call init_counters
        call diag_counters_init
        call trace_parallel(tracer)
        call write_output
        if (allocated(weights)) deallocate (weights)
        allocate (weights(ntestpart))
        weights = 1.0_dp
        weights(1) = 2.0_dp
        call write_results_netcdf('weighted-results.nc', wall%facet_count, weights)
        print *, 'Case ', trim(directory), ' tracing/output seconds: ', &
            get_wtime() - start
    end do
end program test_collisional_stl_wall_fast
