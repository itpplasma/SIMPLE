program bench_secders_cost
    ! Measures rho = cost(eval_field mode_secders=2) / cost(eval_field mode_secders=0)
    ! on the production Boozer chart.
    !
    ! Why this number decides a design question. A two-derivative Runge-Kutta
    ! (TDRK) step on the first-order guiding-center system costs one evaluation
    ! of F (first derivatives only, mode_secders=0) plus s evaluations of
    ! G = F'F (needs second derivatives, mode_secders=2). A classical explicit
    ! Runge-Kutta step of the same order costs s_RK evaluations of F. So TDRK is
    ! cheaper at matched order iff
    !
    !     1 + s * rho  <  s_RK      i.e.   rho < (s_RK - 1) / s
    !
    ! The order-4 spike gives s = 2, s_RK = 4, so the break-even is rho = 1.5.
    ! If the measured rho exceeds the break-even, TDRK loses on this field
    ! representation and the SIMPLE-side TDRK work should not be built.
    !
    ! Timing method: alternate the two modes in an interleaved loop so cache
    ! state and CPU frequency drift affect both equally, and evaluate at
    ! pseudo-random points on a fixed grid so the spline branch pattern is
    ! representative rather than a single hot cell.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_mod, only: eval_field => evaluate, field_can_t
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: coord_input, field_input
    use timing, only: get_wtime

    implicit none

    integer, parameter :: NPT = 20000
    integer, parameter :: NREP = 20

    type(tracer_t)   :: norb
    type(field_can_t) :: f
    real(dp), allocatable :: r(:), th(:), ph(:)
    real(dp) :: t0, t_sec0, t_sec2, rho, breakeven
    real(dp) :: sink
    integer  :: i, rep
    character(len=256) :: wout

    call get_wout_path(wout)
    field_input = trim(wout)
    coord_input = trim(wout)
    call init_field(norb, trim(wout), 5, 5, 3, 2)

    call build_points(r, th, ph)

    ! Warm up both paths so neither pays first-touch or page-fault cost.
    sink = 0.0_dp
    do i = 1, 200
        call eval_field(f, r(i), th(i), ph(i), 0)
        sink = sink + f%Bmod
        call eval_field(f, r(i), th(i), ph(i), 2)
        sink = sink + f%Bmod
    end do

    ! Optional isolation mode: "0" or "2" as argv(1) times only that path, so an
    ! in-process interaction between the two loops can be ruled out.
    call handle_isolated_mode(f, r, th, ph, sink)

    ! Alternate which mode goes first on each repetition. If the two totals
    ! disagree with a fixed ordering but agree here, the difference was an
    ! ordering artefact (frequency ramp, cache warming) rather than real work.
    t_sec0 = 0.0_dp
    t_sec2 = 0.0_dp

    do rep = 1, NREP
        if (mod(rep, 2) == 1) then
            t0 = get_wtime()
            do i = 1, NPT
                call eval_field(f, r(i), th(i), ph(i), 0)
                sink = sink + f%Bmod
            end do
            t_sec0 = t_sec0 + (get_wtime() - t0)

            t0 = get_wtime()
            do i = 1, NPT
                call eval_field(f, r(i), th(i), ph(i), 2)
                sink = sink + f%Bmod
            end do
            t_sec2 = t_sec2 + (get_wtime() - t0)
        else
            t0 = get_wtime()
            do i = 1, NPT
                call eval_field(f, r(i), th(i), ph(i), 2)
                sink = sink + f%Bmod
            end do
            t_sec2 = t_sec2 + (get_wtime() - t0)

            t0 = get_wtime()
            do i = 1, NPT
                call eval_field(f, r(i), th(i), ph(i), 0)
                sink = sink + f%Bmod
            end do
            t_sec0 = t_sec0 + (get_wtime() - t0)
        end if
    end do

    rho = t_sec2 / t_sec0
    breakeven = 1.5_dp

    write (*, "(a)") "=========================================================="
    write (*, "(a)") " eval_field cost ratio, production Boozer chart"
    write (*, "(a)") "=========================================================="
    write (*, "(a,i0,a,i0)") " points per pass ", NPT, ", passes ", NREP
    write (*, "(a,es12.4,a)") " mode_secders=0 total ", t_sec0, " s"
    write (*, "(a,es12.4,a)") " mode_secders=2 total ", t_sec2, " s"
    write (*, "(a,f8.3)")     " rho = cost(secders=2)/cost(secders=0) = ", rho
    write (*, "(a)") ""
    write (*, "(a,f6.2,a)") " order-4 TDRK break-even rho = ", breakeven, &
        "  (1F+2G vs RK4 4F)"
    if (rho < breakeven) then
        write (*, "(a)") " VERDICT: TDRK is cheaper than RK4 at matched order."
    else
        write (*, "(a)") " VERDICT: TDRK is NOT cheaper than RK4 at matched order"
        write (*, "(a)") "          on this field representation."
    end if
    write (*, "(a,es12.4)") " (sink, ignore) ", sink

contains

    ! Time a single mode in isolation and stop. Used to check that the
    ! interleaved measurement is not distorted by an in-process interaction
    ! between the two code paths.
    subroutine handle_isolated_mode(f, r, th, ph, sink)
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: r(:), th(:), ph(:)
        real(dp), intent(inout) :: sink

        character(len=8) :: arg
        integer :: mode, i, rep, nargs
        real(dp) :: t0, total

        nargs = command_argument_count()
        if (nargs < 1) return
        call get_command_argument(1, arg)
        if (trim(arg) == '0') then
            mode = 0
        else if (trim(arg) == '2') then
            mode = 2
        else
            return
        end if

        total = 0.0_dp
        do rep = 1, NREP
            t0 = get_wtime()
            do i = 1, NPT
                call eval_field(f, r(i), th(i), ph(i), mode)
                sink = sink + f%Bmod
            end do
            total = total + (get_wtime() - t0)
        end do

        write (*, "(a,i0,a,es12.4,a,f9.1,a)") "ISOLATED mode_secders=", mode, &
            "  total ", total, " s   per call ", &
            1.0e9_dp*total/real(NREP*NPT, dp), " ns"
        write (*, "(a,es12.4)") " (sink, ignore) ", sink
        stop 0
    end subroutine handle_isolated_mode

    ! Orbit-like sampling: small increments along a smooth path, which is how
    ! the integrator actually calls the field. Scattered sampling would measure
    ! worst-case cache-miss behaviour instead, and is not representative of the
    ! access pattern the cost model is about. Deterministic, so reproducible.
    subroutine build_points(r, th, ph)
        real(dp), allocatable, intent(out) :: r(:), th(:), ph(:)

        real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
        integer :: i
        real(dp) :: s

        allocate (r(NPT), th(NPT), ph(NPT))
        do i = 1, NPT
            s = real(i - 1, dp)
            ! Step sizes comparable to an orbit substep: ~1e-2 rad per call,
            ! with a slow radial excursion that stays inside the chart.
            r(i) = 0.5_dp + 0.35_dp*sin(1.0e-3_dp*s)
            th(i) = mod(1.7e-2_dp*s, twopi)
            ph(i) = mod(0.9e-2_dp*s, twopi)
        end do
    end subroutine build_points

    subroutine get_wout_path(wout)
        character(len=*), intent(out) :: wout

        logical :: found

        wout = 'wout.nc'
        inquire (file=trim(wout), exist=found)
        if (found) return

        wout = 'test/test_data/wout.nc'
        inquire (file=trim(wout), exist=found)
        if (found) return

        error stop 'bench_secders_cost: no wout.nc found (run from project root)'
    end subroutine get_wout_path

end program bench_secders_cost
