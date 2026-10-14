program bench_spectre_fo_step
    ! Single-thread FO cost microbenchmark (not a test). FO everywhere from GC
    ! states near an interface; reports CPU us per accepted RK step and per
    ! gyro-period, and the speed envelope from the returned GC y(4).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: c, e_charge, p_mass, ev
    use magfie_sub, only: set_magfie_spectre_field, spectre_field
    use field_spectre, only: spectre_field_t, create_spectre_field
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_configure, &
        spectre_fo_enter, spectre_fo_advance_until_exit, SPECTRE_FO_OK, &
        SPECTRE_FO_LOSS, spectre_fo_step_count
    implicit none
    type(spectre_field_t) :: sf
    type(spectre_fo_state_t) :: st
    character(len=1024) :: h5file, arg
    real(dp) :: ro0, robar, rho0, y(5), yo(5), t0, t1, dt, used, v0, tgyro
    real(dp) :: emin, emax, total_cpu, cpu(16), steps(16), bscale
    integer :: status, owner, m, rep, nrep, nmark, ngyro, k, chunks
    logical :: exited
    integer(8) :: nsteps

    call get_command_argument(1, h5file)
    call get_command_argument(2, arg); read (arg, *) rho0
    call get_command_argument(3, arg); read (arg, *) nrep
    call get_command_argument(4, arg); read (arg, *) ngyro
    nmark = 8
    call create_spectre_field(sf, trim(h5file), status)
    if (status /= 0) error stop 'field'
    call set_magfie_spectre_field(sf)
    call cpu_time(t0)
    call spectre_fo_configure(1.0e-9_dp)
    call cpu_time(t1)
    print '(a,f8.2,a)', 'configure cpu ', t1 - t0, ' s'
    v0 = sqrt(2.0_dp*3.5e6_dp/1000.0_dp*ev/(4.0_dp*p_mass))
    ro0 = v0*4.0_dp*p_mass*c/(2.0_dp*e_charge)
    robar = ro0/sqrt(2.0_dp)
    bscale = 1.0e4_dp   ! ~1 T reference for the gyro-period unit
    tgyro = 8.0_dp*atan(1.0_dp)*robar/bscale
    do rep = 1, nrep
        nsteps = 0
        emin = huge(1.0_dp); emax = -huge(1.0_dp)
        call cpu_time(t0)
        do m = 1, nmark
            y = [rho0, 0.7_dp*real(m, dp), 0.3_dp*real(m, dp), 1.0_dp, &
                -0.9_dp + 1.8_dp*real(m - 1, dp)/real(nmark - 1, dp)]
            call spectre_fo_enter(st, y, int(rho0) + 1, nint(rho0), robar, status)
            if (status /= SPECTRE_FO_OK) then
                print *, 'enter failed', m, status
                cycle
            end if
            st%exit_rho = 1.0e6_dp
            spectre_fo_step_count = 0
            chunks = 10*ngyro
            dt = tgyro/10.0_dp
            do k = 1, chunks
                call spectre_fo_advance_until_exit(st, dt, robar, used, yo, &
                    owner, exited, status)
                if (status /= SPECTRE_FO_OK) exit
                emin = min(emin, yo(4)); emax = max(emax, yo(4))
            end do
            if (status /= SPECTRE_FO_OK .and. status /= SPECTRE_FO_LOSS) &
                print *, 'advance status', m, status, k
            nsteps = nsteps + spectre_fo_step_count
        end do
        call cpu_time(t1)
        cpu(rep) = t1 - t0
        steps(rep) = real(nsteps, dp)
        print '(a,i2,a,f9.3,a,i10,a,f9.3,a,es10.3)', 'rep ', rep, ' cpu ', &
            cpu(rep), ' s steps ', nsteps, ' us/step ', 1.0e6_dp*cpu(rep)/ &
            steps(rep), ' rel speed envelope ', (emax - emin)/emax
    end do
    total_cpu = sum(cpu(1:nrep))
    print '(a,f9.3,a,f9.3,a,f9.3)', 'us/step median-ish mean ', &
        1.0e6_dp*total_cpu/sum(steps(1:nrep)), ' min ', &
        1.0e6_dp*minval(cpu(1:nrep)/steps(1:nrep)), ' max ', &
        1.0e6_dp*maxval(cpu(1:nrep)/steps(1:nrep))
    print '(a,f10.3)', 'steps per gyro-period (1 T ref) ', &
        sum(steps(1:nrep))/real(nrep*nmark*ngyro, dp)
end program bench_spectre_fo_step
