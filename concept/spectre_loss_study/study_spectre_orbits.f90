program study_spectre_orbits
    ! Production GC/FO driver. Both solver revisions use this same measurement
    ! source and accepted-FO-clock probe. Physical losses and solver stops are
    ! separate outputs; a numerical failure never counts as physical loss.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_wtime, omp_get_max_threads, omp_set_num_threads
    use simple_api, only: simple_read_config, simple_init, simple_sample_volume
    use simple, only: init_sympl
    use params, only: ntestpart, trace_time, v0, dtaumin, relerr, integmode, &
        crossing_level
    use magfie_sub, only: init_magfie, SPECTRE, spectre_field, M_TO_CM
    use field_can_mod, only: field_can_t, ref_to_integ, integ_to_ref
    use field_can_spectre, only: spectre_mvol, set_spectre_volume_lock
    use parmot_mod, only: ro0
    use spectre_fo_hybrid, only: spectre_fo_enter, SPECTRE_FO_OK
    use spectre_fo_field, only: spectre_fo_field_local, spectre_fo_field_locate, &
        spectre_fo_local_t
    use orbit_symplectic_base, only: symplectic_integrator_t
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
        orbit_microstep_sympl_spectre, sympl_landing_stats_reset, sympl_fo_stats, &
        SYMPL_SPECTRE_OK, SYMPL_SPECTRE_LOSS
    use util, only: sqrt2, twopi
    implicit none

    character(len=1024) :: config_file, prefix, starts_file, fixture
    character(len=256) :: message, sampling, line
    character(len=16) :: study_mode
    real(dp), allocatable :: starts(:, :), records(:, :, :), result(:, :)
    integer, allocatable :: samples(:)
    real(dp) :: setup0, setup_s, trace0, trace_s, dt, begin_marker, elapsed
    real(dp) :: phase_offset
    real(dp) :: y(5), z(5), fraction, worst_shell, bmag, a(3), hc(3)
    real(dp) :: g(3, 3), gi(3, 3), jac, radial_larmor, u(3)
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: field
    type(sympl_spectre_state_t) :: state
    integer :: status, m, k, nsteps, stride, capacity, rec, unit, j, nfail
    integer :: entries, exits, losses, failures, counts(5), terminal, iface, side
    integer :: requested_samples

    call get_command_argument(1, fixture)
    call get_environment_variable('SIMPLE_SPECTRE_STUDY_INPUT', config_file)
    call get_environment_variable('SIMPLE_SPECTRE_STUDY_OUTPUT', prefix)
    call get_environment_variable('SIMPLE_SPECTRE_STUDY_STARTS', starts_file)
    call get_environment_variable('SIMPLE_SPECTRE_STUDY_SAMPLING', sampling)
    call get_environment_variable('SIMPLE_SPECTRE_STUDY_MODE', study_mode)
    phase_offset = 0.0_dp
    requested_samples = 65
    call get_environment_variable('SIMPLE_SPECTRE_STUDY_PHASE', message)
    if (len_trim(message) > 0) read (message, *) phase_offset
    call get_environment_variable('SIMPLE_SPECTRE_STUDY_SAMPLES', message)
    if (len_trim(message) > 0) read (message, *) requested_samples
    call get_environment_variable('OMP_NUM_THREADS', message)
    if (len_trim(message) == 0) call omp_set_num_threads(8)
    if (len_trim(prefix) == 0) prefix = 'study'
    if (len_trim(sampling) == 0) sampling = 'volume'
    if (len_trim(config_file) == 0) then
        config_file = 'study_smoke.in'
        open (newunit=unit, file=trim(config_file), status='replace')
        write (unit, '(a)') '&config'
        write (unit, '(a)') "field_input='"//trim(fixture)//"'"
        write (unit, '(a)') 'integ_coords=6, integmode=3, crossing_level=3'
        write (unit, '(a)') 'ntestpart=8, ntimstep=17, npoiper2=1024'
        write (unit, '(a)') 'facE_al=1000, trace_time=2d-7, relerr=1d-12'
        write (unit, '(a)') 'deterministic=.true., ran_seed=12345 /'
        close (unit)
    end if
    setup0 = omp_get_wtime()
    call simple_read_config(trim(config_file), status, message)
    if (status /= 0) then
        print *, trim(message)
        error stop 'study configuration'
    end if
    call simple_init
    call init_magfie(SPECTRE)
    setup_s = omp_get_wtime() - setup0
    allocate (starts(5, ntestpart))
    if (len_trim(starts_file) > 0) then
        open (newunit=unit, file=trim(starts_file), status='old')
        read (unit, '(a)') line
        do m = 1, ntestpart
            read (unit, *) j, starts(:, m)
            if (j /= m) error stop 'study marker identity'
        end do
        close (unit)
    else if (trim(sampling) == 'volume') then
        call simple_sample_volume(ntestpart, 0.1_dp, &
            real(spectre_mvol, dp) - 0.05_dp, starts)
    else
        ! Identical raw-field seeds on both revisions. Half-gyroradius offsets
        ! on both sides of every sheet ensure actual interface encounters.
        do m = 1, ntestpart
            iface = 1 + mod((m - 1)/2, spectre_mvol - 1)
            side = merge(1, -1, mod(m, 2) == 0)
            u = [real(iface, dp) + side*1.0e-10_dp, &
                modulo(2.399963229728653_dp*m, twopi), &
                modulo(1.414213562373095_dp*m, twopi)]
            y = [u, 1.0_dp, -0.9_dp + &
                1.8_dp*real(m - 1, dp)/real(max(1, ntestpart - 1), dp)]
            call spectre_field%evaluate(u, a, hc, bmag)
            call spectre_field%coords%metric_tensor(u, g, gi, jac)
            ! ro0 in G cm; raw B and metric use T and m.
            radial_larmor = ro0*sqrt(1 - y(5)**2)/(bmag*1.0e4_dp) &
                *sqrt(gi(1, 1))/M_TO_CM
            y(1) = real(iface, dp) + side*0.5_dp*radial_larmor
            starts(:, m) = y
        end do
    end if
    open (newunit=unit, file=trim(prefix)//'_starts.csv', status='replace')
    write (unit, '(a)') 'marker,rho,theta,zeta,pitch_speed,lambda'
    do m = 1, ntestpart
        write (unit, '(i0,5(",",es24.16))') m, starts(:, m)
    end do
    close (unit)
    dt = dtaumin/v0
    nsteps = nint(trace_time/dt)
    stride = max(1, ceiling(real(nsteps, dp)/max(1, requested_samples - 1)))
    capacity = nsteps/stride + 3
    allocate (records(12, capacity, ntestpart), result(7, ntestpart))
    allocate (samples(ntestpart))
    records = 0.0_dp
    result = 0.0_dp
    call sympl_landing_stats_reset
    trace0 = omp_get_wtime()
    !$omp parallel do private(m, k, y, z, si, field, state, status, fraction, &
    !$omp begin_marker, elapsed, worst_shell, rec, terminal)
    do m = 1, ntestpart
        begin_marker = omp_get_wtime()
        y = starts(:, m)
        ! Match production's per-marker dispatch reset before initialization.
        call set_spectre_volume_lock(0)
        call ref_to_integ(y(1:3), z(1:3))
        z(4:5) = y(4:5)
        call init_sympl(si, field, z, dtaumin, dtaumin, relerr, integmode)
        call sympl_spectre_reset(state, si, spectre_mvol, integmode, crossing_level)
        if (trim(study_mode) == 'fo') then
            call spectre_fo_enter(state%fo, y, &
                min(spectre_mvol, int(y(1)) + 1), nint(y(1)), &
                ro0/sqrt2, status)
            if (status /= SPECTRE_FO_OK) error stop 'FO reference launch'
            if (trim(sampling) == 'phase') then
                call rotate_launch(state, y, phase_offset + &
                    twopi*real(m - 1, dp)/real(ntestpart, dp))
            else if (phase_offset /= 0.0_dp) then
                call rotate_launch(state, y, phase_offset)
            end if
            state%fo%exit_rho = 1.0e6_dp
        end if
        rec = 1
        elapsed = 0.0_dp
        worst_shell = 0.0_dp
        terminal = SYMPL_SPECTRE_OK
        call sample_state(state, si, field, elapsed, records(:, rec, m))
        do k = 1, nsteps
            call orbit_microstep_sympl_spectre(state, si, field, status, fraction)
            elapsed = elapsed + dt*fraction
            if (state%fo%active) then
                worst_shell = max(worst_shell, &
                    abs(sum(state%fo%v**2)/2 - starts(4, m)**2))
            else
                worst_shell = max(worst_shell, &
                    abs(field%H - starts(4, m)**2))
            end if
            if (mod(k, stride) == 0 .or. status /= SYMPL_SPECTRE_OK .or. &
                k == nsteps) then
                rec = rec + 1
                call sample_state(state, si, field, elapsed, records(:, rec, m))
            end if
            if (status /= SYMPL_SPECTRE_OK) then
                terminal = status
                exit
            end if
        end do
        samples(m) = rec
        result(:, m) = [elapsed, real(terminal, dp), &
            state%fo_elapsed*sqrt2/v0/max(elapsed, tiny(1.0_dp)), &
            worst_shell, omp_get_wtime() - begin_marker, &
            real(state%warning_hold_count, dp), real(min(k, nsteps), dp)]
        call set_spectre_volume_lock(0)
    end do
    !$omp end parallel do
    trace_s = omp_get_wtime() - trace0
    call sympl_fo_stats(entries, exits, losses, failures, counts)
    nfail = count(nint(result(2, :)) /= SYMPL_SPECTRE_OK .and. &
        nint(result(2, :)) /= SYMPL_SPECTRE_LOSS)
    open (newunit=unit, file=trim(prefix)//'_markers.csv', status='replace')
    write (unit, '(a)') 'marker,time_s,status,fo_fraction,shell_error,'// &
        'wall_s,warning_hold_count,microsteps'
    do m = 1, ntestpart
        write (unit, '(i0,7(",",es24.16))') m, result(:, m)
    end do
    close (unit)
    open (newunit=unit, file=trim(prefix)//'_trace.csv', status='replace')
    write (unit, '(a)') 'marker,time_s,rho,theta,zeta,energy,lambda,'// &
        'fo,owner,x_cm,y_cm,z_cm,fo_time_s'
    do m = 1, ntestpart
        do j = 1, samples(m)
            write (unit, '(i0,12(",",es24.16))') m, records(:, j, m)
        end do
    end do
    close (unit)
    open (newunit=unit, file=trim(prefix)//'_summary.json', status='replace')
    write (unit, '(a)') '{'
    write (unit, '(a,i0,a)') '"markers":', ntestpart, ','
    write (unit, '(a,i0,a)') '"threads":', omp_get_max_threads(), ','
    write (unit, '(a,es24.16,a)') '"duration_s":', trace_time, ','
    write (unit, '(a,es24.16,a)') '"setup_s":', setup_s, ','
    write (unit, '(a,es24.16,a)') '"trace_s":', trace_s, ','
    write (unit, '(a,es24.16,a)') '"v0_cm_s":', v0, ','
    write (unit, '(a,i0,a)') '"physical_losses":', &
        count(nint(result(2, :)) == SYMPL_SPECTRE_LOSS), ','
    write (unit, '(a,i0,a)') '"numerical_stops":', nfail, ','
    write (unit, '(a,i0,a)') '"fo_entries":', entries, ','
    write (unit, '(a,i0,a)') '"fo_exits":', exits, ','
    write (unit, '(a,i0,a)') '"fo_failures":', failures, ','
    write (unit, '(a,i0)') '"warning_hold_count":', nint(sum(result(6, :)))
    write (unit, '(a)') '}'
    close (unit)
    print '(a,2f12.4,a,3i6)', 'study setup/trace s ', setup_s, trace_s, &
        ' entries/exits/stops ', entries, exits, nfail
    if (.not. all(ieee_is_finite(result))) error stop 'nonfinite study result'
    if (nfail /= 0) error stop 'study recorded numerical stops'
    if (failures /= 0 .or. any(result(6, :) > 0)) &
        error stop 'study contains unresolved intervals'

contains

    subroutine rotate_launch(st, y_gc, phase)
        ! Research-only gyrophase control at the same prescribed GC start.
        ! Rotate v_perp around B and apply the same leading-order Larmor map
        ! as production entry; the orbit integrator itself is unchanged.
        type(sympl_spectre_state_t), intent(inout) :: st
        real(dp), intent(in) :: y_gc(5), phase
        type(spectre_fo_local_t) :: loc
        real(dp) :: vp(3), rotated(3), bxv(3), parallel
        integer :: ierr

        call spectre_fo_field_local(y_gc(1:3), st%fo%entry_owner, loc)
        parallel = dot_product(st%fo%v, loc%b)
        vp = st%fo%v - parallel*loc%b
        bxv = cross_product(loc%b, vp)
        rotated = cos(phase)*vp + sin(phase)*bxv
        st%fo%v = parallel*loc%b + rotated
        st%fo%x = loc%x + ro0/sqrt2/loc%bmod*cross_product(loc%b, rotated)
        st%fo%u = y_gc(1:3)
        st%fo%owner = st%fo%entry_owner
        call spectre_fo_field_locate(st%fo%x, st%fo%u, st%fo%owner, ierr)
        if (ierr /= 0) error stop 'phase-control launch outside field'
    end subroutine rotate_launch

    pure function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross_product

    subroutine sample_state(st, integrator, f, time, row)
        type(sympl_spectre_state_t), intent(in) :: st
        type(symplectic_integrator_t), intent(in) :: integrator
        type(field_can_t), intent(in) :: f
        real(dp), intent(in) :: time
        real(dp), intent(out) :: row(12)
        real(dp) :: coordinate(3), pos(3), speed, lambda
        integer :: owner

        if (st%fo%active) then
            coordinate = st%fo%u
            pos = st%fo%x
            speed = sum(st%fo%v**2)/2
            lambda = st%fo%last_y(5)
            owner = st%fo%owner
        else
            call integ_to_ref(integrator%z(1:3), coordinate)
            call spectre_field%coords%evaluate_cart(coordinate, pos)
            pos = pos*M_TO_CM
            speed = f%H
            lambda = f%vpar/sqrt(2*f%H)
            owner = min(spectre_mvol, int(coordinate(1)) + 1)
        end if
        row = [time, coordinate, speed, lambda, &
            real(merge(1, 0, st%fo%active), dp), real(owner, dp), &
            pos, st%fo_elapsed*sqrt2/v0]
    end subroutine sample_state
end program study_spectre_orbits
