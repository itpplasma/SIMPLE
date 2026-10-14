program compare
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use omp_lib, only: omp_get_max_threads
    use candidate_field, only: model_t, local_t, geometry, side_of, frame, &
        radius, canonical_pphi
    use candidate_orbits, only: options_t, orbit_t, initialize, advance, particle
    use transforms, only: kernel_available
    use fo_sheet, only: poincare_params_t, poincare_gauss2_step
    use sheet_field, only: region_of, psi_value, hamiltonian
    implicit none
    integer, parameter :: nmarker = 256, nreverse = 32, ncol = 19
    type(model_t) :: pm
    type(options_t) :: opt
    type(orbit_t) :: endings(nmarker)
    real(dp) :: initial(6, nmarker), data(ncol, nmarker), speed, tend, lossr, width
    real(dp) :: randoms(4, nmarker), dt, aout, wall, revwall
    character(256) :: arg, filename, case_name, method_name
    integer(int64) :: clk0, clk1, rate, rng
    integer :: i, j, unit
    if (command_argument_count() < 12) error stop &
        'usage: METHOD CASE AOUT SPEED TEND TOL DIN ORDER GL2_STEP GL2_WIDTH SAMPLE FILE'
    if (omp_get_max_threads() > 16) error stop 'OMP_NUM_THREADS must be <=16'
    call get_command_argument(1, method_name)
    call get_command_argument(2, case_name)
    call get_command_argument(3, arg)
    read (arg, *) aout
    call get_command_argument(4, arg)
    read (arg, *) speed
    call get_command_argument(5, arg)
    read (arg, *) tend
    call get_command_argument(6, arg)
    read (arg, *) opt%tol
    call get_command_argument(7, arg)
    read (arg, *) opt%din
    opt%dout = opt%din + 1.0_dp
    call get_command_argument(8, arg)
    read (arg, *) opt%order
    call get_command_argument(9, arg)
    read (arg, *) dt
    call get_command_argument(10, arg)
    read (arg, *) width
    call get_command_argument(11, arg)
    read (arg, *) opt%sample
    call get_command_argument(12, filename)
    select case (trim(method_name))
    case ('FO')
        opt%method = 0
    case ('C3')
        opt%method = 3
    case ('GL2')
        opt%method = 6
        if (trim(case_name) /= 'A') error stop 'GL2 prototype is axisymmetric'
    case default
        error stop 'unknown method'
    end select
    if (opt%method == 3) then
        if (.not. kernel_available(opt%order)) &
            error stop 'requested transform/RHS/energy order unavailable'
    end if
    pm%torus%a_in = 0.25_dp
    pm%torus%a_out = aout
    if (trim(case_name) == 'B') pm%delta = 0.1_dp
    lossr = 0.55_dp
    rng = 20261014_int64
    do i = 1, nmarker
        do j = 1, 4
            rng = modulo(16807_int64*rng, 2147483647_int64)
            randoms(j, i) = real(rng, dp)/2147483647.0_dp
        end do
        call marker(randoms(:, i), initial(:, i))
    end do
    ! Initial states always use the sharp target field, even for layer controls.
    if (opt%method == 6) then
        pm%torus%width = width*speed
        pm%torus%analytic = opt%method == 6
    end if
    call system_clock(clk0, rate)
    !$omp parallel do default(shared) private(i) schedule(dynamic)
    do i = 1, nmarker
        call trace(i, data(:, i), endings(i))
    end do
    !$omp end parallel do
    call system_clock(clk1)
    wall = real(clk1 - clk0, dp)/real(rate, dp)
    call system_clock(clk0)
    !$omp parallel do default(shared) private(i) schedule(dynamic)
    do i = 1, nreverse
        call reverse(i, endings(i), data(10:11, i))
    end do
    !$omp end parallel do
    call system_clock(clk1)
    revwall = real(clk1 - clk0, dp)/real(rate, dp)
    open (newunit=unit, file=trim(filename), status='replace')
    write (unit, '(a,es24.16,a,es24.16,a,i0)') '# wall=', wall, &
        ',reverse_wall=', revwall, ',threads=', omp_get_max_threads()
    write (unit, '(a)') 'marker,lost,failed,emax,pmax,efinal,pfinal,rmin,rmax,' &
        // 'excursion,rev_position,rev_velocity,switches,crossings,root_failures,' &
        // 'steps,evaluations,fo_fraction,switch_energy,speed'
    do i = 1, nmarker
        write (unit, '(i0,19(a,es24.16))') i, (',', data(j, i), j=1, ncol)
    end do
    close (unit)
    print '(5a,f6.3,a,f7.3,a,f9.2,a,i0,a,i0)', trim(method_name), &
        ' T-', trim(case_name), ' v=', ' ', speed, ' aout=', aout, &
        ' wall=', wall, ' lost=', int(sum(data(1, :))), &
        ' failed=', int(sum(data(2, :)))
contains

    subroutine marker(rnd, y)
        real(dp), intent(in) :: rnd(4)
        real(dp), intent(out) :: y(6)
        real(dp) :: theta, phi, phase, pitch, rr, r, pi, e1(3), e2(3)
        type(local_t) :: lf
        pi = acos(-1.0_dp)
        theta = 2*pi*rnd(1)
        phi = 2*pi*rnd(2)
        phase = 2*pi*rnd(3)
        pitch = 2*rnd(4) - 1
        r = pm%torus%r_sheet - 0.35_dp*speed
        rr = pm%torus%r_major + r*cos(theta)
        y(1:3) = [rr*cos(phi), rr*sin(phi), r*sin(theta)]
        call geometry(pm, y(1:3), 1, lf)
        call frame(lf%b, e1, e2)
        y(4:6) = speed*(pitch*lf%b + sqrt(1 - pitch*pitch)* &
            (cos(phase)*e1 + sin(phase)*e2))
    end subroutine

    subroutine trace(im, row, ending)
        integer, intent(in) :: im
        real(dp), intent(out) :: row(ncol)
        type(orbit_t), intent(out) :: ending
        type(orbit_t) :: orbit
        real(dp) :: t, span, y(6), e0, p0, r, rinit, emin, pmax, emax
        real(dp) :: rmin, rmax, efinal, pfinal, glz(4), glh, gle, gltime
        type(model_t) :: local_pm
        type(poincare_params_t) :: pp
        logical :: good
        integer :: ncross, steps
        local_pm = pm
        y = initial(:, im)
        e0 = 0.5_dp*sum(y(4:6)**2)
        p0 = canonical_pphi(pm, y)
        rinit = radius(pm, y(1:3))
        call initialize(pm, opt, y, orbit)
        rmin = rinit
        rmax = rinit
        emax = 0.0_dp
        pmax = 0.0_dp
        t = 0.0_dp
        ncross = 0
        steps = 0
        pp%g_min = 0.1_dp
        pp%length = 0.2_dp*speed
        if (opt%method == 6) then
            local_pm%torus%p_phi = p0
            call to_reduced(local_pm, y, glz)
            gle = hamiltonian(local_pm%torus, 1, glz)
        end if
        do while (t < tend - 1.0e-10_dp)
            span = min(opt%sample, tend - t)
            if (opt%method == 6) then
                gltime = 0.0_dp
                do while (gltime < span - 1.0e-12_dp)
                    glh = dt
                    call gl_advance(local_pm, pp, glz, gle, glh, &
                        span - gltime, gltime, good)
                    steps = steps + 1
                    if (.not. good) exit
                end do
                orbit%failed = .not. good
                call from_reduced(local_pm, glz, y)
                orbit%y = y
                orbit%steps = steps
                orbit%fo_time = t + span
            else
                call advance(pm, opt, orbit, span)
                call particle(pm, opt, orbit, y, good)
                if (.not. good) orbit%failed = .true.
            end if
            if (orbit%failed) exit
            t = t + span
            r = radius(pm, y(1:3))
            rmin = min(rmin, r)
            rmax = max(rmax, r)
            emin = 0.5_dp*sum(y(4:6)**2)/e0 - 1
            emax = max(emax, abs(emin))
            if (abs(pm%delta) < tiny(1.0_dp)) &
                pmax = max(pmax, abs(canonical_pphi(pm, y) - p0))
        end do
        efinal = 0.5_dp*sum(y(4:6)**2)/e0 - 1
        pfinal = canonical_pphi(pm, y) - p0
        if (abs(pm%delta) > tiny(1.0_dp)) pmax = -1.0_dp
        row = [merge(1.0_dp, 0.0_dp, rmax >= lossr), &
            merge(1.0_dp, 0.0_dp, orbit%failed), emax, pmax, efinal, pfinal, &
            rmin, rmax, rmax - rmin, -1.0_dp, -1.0_dp, &
            real(orbit%switches, dp), real(orbit%crossings, dp), &
            real(orbit%roots_failed, dp), real(orbit%steps, dp), &
            real(orbit%evaluations, dp), orbit%fo_time/max(t, 1.0e-30_dp), &
            orbit%switch_energy, speed]
        ending = orbit
    end subroutine

    subroutine reverse(im, orbit, errors)
        integer, intent(in) :: im
        type(orbit_t), intent(inout) :: orbit
        real(dp), intent(out) :: errors(2)
        real(dp) :: t, span, y(6), z(4), gltime, gle, glz_for_reverse(4)
        type(model_t) :: local_pm
        type(poincare_params_t) :: pp
        logical :: good
        errors = -1.0_dp
        if (orbit%failed) return
        t = 0.0_dp
        if (opt%method == 6) then
            local_pm = pm
            local_pm%torus%p_phi = orbit%p0
            call to_reduced(local_pm, orbit%y, z)
            gle = orbit%energy0
            pp%g_min = 0.1_dp
            pp%length = 0.2_dp*speed
        end if
        do while (t < tend - 1.0e-10_dp)
            span = min(opt%sample, tend - t)
            if (opt%method == 6) then
                gltime = 0.0_dp
                do while (gltime < span - 1.0e-12_dp)
                    call gl_advance(local_pm, pp, z, gle, -dt, &
                        span - gltime, gltime, good)
                    if (.not. good) exit
                end do
                call from_reduced(local_pm, z, y)
            else
                call advance(pm, opt, orbit, -span)
                call particle(pm, opt, orbit, y, good)
            end if
            if (.not. good .or. orbit%failed) return
            t = t + span
        end do
        ! GL2 reduction quotients out the ignorable toroidal angle: compare
        ! R,Z,v_R,v_phi,v_Z, not an unintegrated azimuthal angle.
        if (opt%method == 6) then
            call to_reduced(local_pm, initial(:, im), glz_for_reverse)
            errors(1) = norm2(z(1:2) - glz_for_reverse(1:2))/speed
            errors(2) = norm2(z(3:4) - glz_for_reverse(3:4))/speed
        else
            errors(1) = norm2(y(1:3) - initial(1:3, im))/speed
            errors(2) = norm2(y(4:6) - initial(4:6, im))/speed
        end if
    end subroutine

    subroutine to_reduced(local_pm, y, z)
        type(model_t), intent(in) :: local_pm
        real(dp), intent(in) :: y(6)
        real(dp), intent(out) :: z(4)
        real(dp) :: rr, vr
        rr = hypot(y(1), y(2))
        vr = dot_product(y(1:2), y(4:5))/rr
        z = [rr, y(3), vr, y(6) - local_pm%torus%f_tor*log(rr)]
    end subroutine

    subroutine from_reduced(local_pm, z, y)
        type(model_t), intent(in) :: local_pm
        real(dp), intent(in) :: z(4)
        real(dp), intent(out) :: y(6)
        real(dp) :: psi
        psi = psi_value(local_pm%torus, region_of(local_pm%torus, z(1:2)), z(1:2))
        y = [z(1), 0.0_dp, z(2), z(3), (local_pm%torus%p_phi - psi)/z(1), &
            z(4) + local_pm%torus%f_tor*log(z(1))]
    end subroutine

    subroutine gl_advance(local_pm, pp, z, e, h0, remain, elapsed, good)
        type(model_t), intent(in) :: local_pm
        type(poincare_params_t), intent(in) :: pp
        real(dp), intent(inout) :: z(4), elapsed
        real(dp), intent(in) :: e, h0, remain
        logical, intent(out) :: good
        real(dp) :: zn(4), h, used, a, b, ta, tb, target, direction
        integer :: it
        h = h0
        direction = sign(1.0_dp, h0)
        call poincare_gauss2_step(local_pm%torus, pp, z, e, h, zn, used, good)
        if (.not. good) return
        if (abs(used) > remain) then
            ! Endpoint alignment affects exact fixed-fictitious-time
            ! symplecticity; do not attribute it to this sampled map.
            a = 0.0_dp
            b = abs(h)
            ta = 0.0_dp
            tb = abs(used)
            target = remain
            do it = 1, 60
                h = a + (b - a)*(target - ta)/(tb - ta)
                call poincare_gauss2_step(local_pm%torus, pp, z, e, &
                    direction*h, zn, used, good)
                if (.not. good) return
                if (abs(abs(used) - target) < 1.0e-13_dp) exit
                if (abs(used) < target) then
                    a = h
                    ta = abs(used)
                else
                    b = h
                    tb = abs(used)
                end if
            end do
        end if
        z = zn
        elapsed = elapsed + abs(used)
    end subroutine

end program
