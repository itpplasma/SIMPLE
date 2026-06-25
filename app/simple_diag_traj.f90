program diag_traj_main
!> Stream the trajectory of one particle through the PRODUCTION integrator
!> selected by integmode (RK orbit_timestep_axis for integmode<=0, symplectic
!> orbit_timestep_sympl otherwise). Unlike diag_orbit.x it does not hardcode a
!> solver and does not preallocate ntimstep*ntau points: it writes s, theta, phi,
!> pphi, H every `stride` substeps and reports min(s), loss time, and the number
!> of axis crossings (EVT_R_NEGATIVE). Built to localize the symplectic-vs-RK
!> non-prompt loss discrepancy on the internal Boozer field (#370).
!>
!> Usage: ./diag_traj.x [config_file] particle_number [stride]

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, integmode, &
             params_init, dtaumin, relerr, ntestpart, zstart, startmode, grid_density, &
                      special_ants_file, reuse_batch, num_surf, sbeg, trace_time, v0, &
                      orbit_model, ORBIT_FULL_ORBIT
    use simple, only: tracer_t, init_sympl, init_fo
    use simple_main, only: init_field
    use magfie_sub, only: init_magfie, VMEC
    use samplers, only: init_starting_surf, sample, START_FILE
    use field_can_mod, only: field_can_t, get_val, eval_field => evaluate, ref_to_integ
    use orbit_symplectic, only: orbit_timestep_sympl
    use orbit_symplectic_base, only: symplectic_integrator_t
    use alpha_lifetime_sub, only: orbit_timestep_axis
    use orbit_fo_boris, only: fo_step, fo_to_gc, fo_energy, FO_OK, FO_LOSS
    use diag_counters, only: diag_counters_init, diag_counters_reset, &
                             diag_counters_total, EVT_R_NEGATIVE
    use util, only: twopi

    implicit none

    character(256) :: config_file, arg
    type(tracer_t) :: norb
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    integer :: pnum, stride, nargs
    real(dp), dimension(5) :: z
    real(dp) :: t, s, smin, H
    integer(8) :: kt
    integer :: ierr, unit
    character(64) :: fname

    config_file = 'simple.in'
    stride = 1
    nargs = command_argument_count()
    if (nargs == 1) then
        call get_command_argument(1, arg); read (arg, *) pnum
    elseif (nargs == 2) then
        call get_command_argument(1, arg)
        ! second arg is either config (non-numeric) or stride; assume "pnum stride"
        read (arg, *) pnum
        call get_command_argument(2, arg); read (arg, *) stride
    elseif (nargs == 3) then
        call get_command_argument(1, config_file)
        call get_command_argument(2, arg); read (arg, *) pnum
        call get_command_argument(3, arg); read (arg, *) stride
    else
        print *, 'Usage: ./diag_traj.x [config] particle_number [stride]'
        stop
    end if

    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call params_init
    call init_magfie(VMEC)
    call init_starting_surf
    if (startmode == 2) then
        call sample(zstart, START_FILE)
    elseif (startmode == 5) then
        if (num_surf == 1) then
            call sample(zstart, 0.0d0, sbeg(1))
        else
            call sample(zstart, sbeg(1), sbeg(num_surf))
        end if
    else
        call sample(zstart, START_FILE)
    end if

    if (pnum < 1 .or. pnum > ntestpart) then
        print *, 'particle out of range', pnum, ntestpart; error stop
    end if

    call diag_counters_init()
    call diag_counters_reset()

    call ref_to_integ(zstart(1:3, pnum), z(1:3))
    z(4:5) = zstart(4:5, pnum)

    if (orbit_model == ORBIT_FULL_ORBIT) then
        call trace_full_orbit(norb, z, pnum, stride)
        stop
    end if

    if (integmode > 0) call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)

    write (fname, '(A,I0,A,I0,A)') 'traj_p', pnum, '_im', integmode, '.dat'
    open (newunit=unit, file=trim(fname), status='replace')
    write (unit, '(A)') '# t[s]    s    theta    phi    pphi/z4    H'

    kt = 0; t = 0.0_dp; smin = z(1); ierr = 0
    do
        if (integmode <= 0) then
            call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
            s = z(1)
            H = -1.0_dp
        else
            call orbit_timestep_sympl(si, f, ierr)
            z(1:4) = si%z
            s = si%z(1)
            call get_val(f, si%z(4))
            H = f%H
        end if
        kt = kt + 1
        t = kt*dtaumin/v0
        smin = min(smin, s)
        if (mod(kt, int(stride, 8)) == 0_8 .or. ierr /= 0) then
            write (unit, '(6ES16.8)') t, s, mod(z(2), twopi), mod(z(3), twopi), z(4), H
        end if
        if (ierr /= 0) exit
        if (t >= trace_time) exit
    end do
    close (unit)

    print '(A,I0,A,I0)', 'particle ', pnum, '  integmode ', integmode
    print '(A,L1,A,ES12.5)', 'lost = ', (ierr /= 0), '   t_end[s] = ', t
    print '(A,ES12.5)', 'min(s) reached       = ', smin
    print '(A,I0)', 'axis crossings (R<0) = ', int(diag_counters_total(EVT_R_NEGATIVE))
    print '(A,A)', 'trajectory written   : ', trim(fname)

contains

    ! Stream the full-orbit (orbit_model=7) trajectory of one particle, recording
    ! the guiding-centre s and the particle s each step. A row's status flags 0 OK,
    ! 1 guiding-centre loss (s>=1), 2 particle-field inversion fault. The point is
    ! to see whether a terminal fault happens while the guiding centre is still
    ! inside the plasma (a confined orbit miscounted) or genuinely at the edge.
    subroutine trace_full_orbit(norb, z, pnum, stride)
        type(tracer_t), intent(inout) :: norb
        real(dp), intent(in) :: z(5)
        integer, intent(in) :: pnum, stride
        real(dp) :: zz(5), tt, sgc, spart, the, phi, vp, E0, sgcmin
        integer(8) :: k
        integer :: u, st, gst, code
        character(64) :: fn

        zz = z
        call init_fo(norb%fo, zz, dtaumin)
        E0 = max(fo_energy(norb%fo), 1.0e-300_dp)
        sgc = zz(1); spart = norb%fo%u(1)**2; the = zz(2); phi = zz(3); sgcmin = sgc

        write (fn, '(A,I0,A)') 'traj_fo_p', pnum, '.dat'
        open (newunit=u, file=trim(fn), status='replace')
        write (u, '(A)') '# t[s]    s_gc    s_part    theta    phi    E/E0    status'
        write (u, '(6ES16.8,I3)') 0.0_dp, sgc, spart, the, phi, 1.0_dp, 0

        k = 0; gst = FO_OK
        do
            call fo_step(norb%fo, st)
            k = k + 1
            tt = k*dtaumin/v0
            if (st /= FO_OK) then            ! inversion fault
                write (u, '(6ES16.8,I3)') tt, sgc, spart, the, phi, &
                    fo_energy(norb%fo)/E0, 2
                exit
            end if
            spart = norb%fo%u(1)**2
            call fo_to_gc(norb%fo, sgc, the, phi, vp, gst)
            sgcmin = min(sgcmin, sgc)
            code = 0
            if (gst == FO_LOSS) code = 1
            if (gst /= FO_OK .and. gst /= FO_LOSS) code = 2
            if (mod(k, int(stride, 8)) == 0_8 .or. gst /= FO_OK) &
                write (u, '(6ES16.8,I3)') tt, sgc, spart, mod(the, twopi), &
                mod(phi, twopi), fo_energy(norb%fo)/E0, code
            if (gst /= FO_OK) exit
            if (tt >= trace_time) exit
        end do
        close (u)

        print '(A,I0,A)', 'particle ', pnum, '  orbit_model = 7 (full orbit)'
        print '(A,ES12.5)', 't_end[s]          = ', tt
        if (st /= FO_OK) then
            print '(A)', 'exit: particle inversion FAULT (not a physical loss)'
        else if (gst == FO_LOSS) then
            print '(A)', 'exit: guiding-centre loss (s>=1)'
        else
            print '(A)', 'exit: reached trace_time (confined)'
        end if
        print '(A,ES12.5,A,ES12.5)', 'min(s_gc) = ', sgcmin, '   last s_part = ', spart
        print '(A,A)', 'trajectory written : ', trim(fn)
    end subroutine trace_full_orbit

end program diag_traj_main
