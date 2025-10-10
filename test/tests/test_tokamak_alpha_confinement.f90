program test_tokamak_alpha_confinement
    !> System test: Verify alpha particle confinement in analytical tokamak
    !>
    !> This test runs the full SIMPLE workflow with:
    !> - Analytical GS field via geoflux coordinates
    !> - Meiss canonical coordinates
    !> - 128 alpha particles starting at s=0.3
    !> - 1 ms integration time
    !>
    !> Expected: Zero particles lost (perfect confinement without ripple)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple
    use simple_main, only: init_field
    use params, only: read_config, notrace_passing, ntestpart, trace_time, sbeg
    implicit none

    type(Tracer) :: tracer
    character(*), parameter :: config_file = &
        '../../examples/tokamak_alpha_confinement/simple.in'
    integer :: n_lost, n_confined, i
    logical :: test_passed
    real(dp) :: s_min, s_max

    print *, '=== Tokamak Alpha Confinement System Test ==='
    print *, ''

    print *, '1. Loading configuration from example...'
    call read_config(config_file)
    print '(A,I0)', '   Particles: ', ntestpart
    print '(A,F6.3)', '   Starting s: ', sbeg
    print '(A,E10.3,A)', '   Duration: ', trace_time, ' s'
    print *, ''

    print *, '2. Initializing analytical field via geoflux...'
    ! Note: This will fail until field.F90 is updated with analytical support
    ! TODO: Uncomment when analytical field integration is complete
    ! call init_field(tracer, 'analytical', 5, 5, 7, 1)
    print *, '   TODO: Analytical field initialization not yet implemented'
    print *, '   Skipping test until field.F90 is updated'
    print *, ''
    error stop 'Test not yet implemented - waiting for field.F90 integration'

    ! The following code will be activated once integration is complete:

    ! print *, '3. Running particle tracing...'
    ! ! Run SIMPLE with the configuration
    ! ! This would call the main simulation loop
    ! print *, '   (trace complete)'
    ! print *, ''

    ! print *, '4. Checking results...'
    ! ! Parse output to check particle losses
    ! ! For now, this is a placeholder
    ! n_lost = 0
    ! n_confined = ntestpart
    ! s_min = 0.25_dp
    ! s_max = 0.35_dp

    ! print '(A,I0,A,I0)', '   Particles lost: ', n_lost, ' / ', ntestpart
    ! print '(A,F6.3,A,F6.3,A)', '   Flux range: s ∈ [', s_min, ', ', s_max, ']'
    ! print *, ''

    ! test_passed = (n_lost == 0) .and. (s_min >= 0.2_dp) .and. (s_max <= 0.4_dp)

    ! if (test_passed) then
    !     print *, '=== TEST PASSED ==='
    !     print *, 'Perfect confinement achieved in axisymmetric tokamak'
    ! else
    !     print *, '=== TEST FAILED ==='
    !     if (n_lost > 0) then
    !         print *, 'ERROR: Particles lost in axisymmetric configuration'
    !     end if
    !     if (s_min < 0.2_dp .or. s_max > 0.4_dp) then
    !         print *, 'ERROR: Particles drifted far from starting flux surface'
    !     end if
    !     error stop 'Alpha confinement test failed'
    ! end if

end program test_tokamak_alpha_confinement
