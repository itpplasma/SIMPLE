program test_spectre_sympl_crossing
    !> Long-run invariant witnesses for the symplectic multi-volume crossing
    !> pipeline (#441) on the axisymmetric two-volume tok2vol fixture.
    !>
    !> Eight passing markers (lambda = +-0.9) seeded near the interface run
    !> under the implicit midpoint scheme with exact-landing substeps, the
    !> crossing map, and per-volume re-canonicalization. Witnesses:
    !>
    !>   * Interface-crossing markers: the linear-fit slope of H(t) through
    !>     thousands of flow-jump-flow compositions is consistent with zero
    !>     (|slope|*T/H0 below SLOPE_TOL) -- the symplecticity witness. Their
    !>     peak |dH/H| stays below the in-volume scheme floor with margin.
    !>   * Control markers that never cross: canonical p_phi = si%z(4) is
    !>     conserved to round-off (exact axisymmetric momentum map), bounding
    !>     the in-volume floor that H comparisons are measured against.
    !>   * Crossing markers: between events p_phi is constant to round-off;
    !>     every change beyond round-off is a logged crossing kick. The TOTAL
    !>     p_phi walk is reported, not asserted: the F5 mirror-fallback rescale
    !>     (DOC/spectre-interface-crossing.md) conserves H and mu but not
    !>     p_phi, so the map itself moves p_phi at each fallback crossing.
    !>   * Every landing satisfies |rho_g - k| < 1e-8 and no CROSS_STOP occurs.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t, init_sympl
    use simple_main, only: init_field
    use params, only: read_config, params_init, netcdffile, ns_s, ns_tp, &
                      multharm, dtaumin, relerr, v0, crossing_level, integmode
    use magfie_sub, only: init_magfie, SPECTRE
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_val, &
                             ref_to_integ
    use field_can_spectre, only: spectre_mvol, set_spectre_volume_lock
    use orbit_symplectic_base, only: symplectic_integrator_t
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
                                   orbit_microstep_sympl_spectre, recanon_pphi, &
                                   sympl_landing_stats_reset, &
                                   sympl_landing_stats, SYMPL_SPECTRE_OK
    use parmot_mod, only: ro0
    use interface_crossing, only: crossing_log_reset, crossing_log_count_type, &
                                  CROSS_CROSSING, CROSS_STOP
    use util, only: twopi, sqrt2

    implicit none

    integer, parameter :: NMARKER = 8, NSTEP = 60000, NSAMPLE = 16
    integer, parameter :: NREC = NSTEP/NSAMPLE, MIN_REC = 100
    real(dp), parameter :: RHO0 = 0.97_dp, LAM0 = 0.9_dp
    !> Zero-consistency bound on the fitted relative H slope of crossing
    !> markers; measured values sit at 1e-8..1e-7 while the in-volume floor
    !> (control markers, no crossings) drifts at ~1e-5.
    real(dp), parameter :: SLOPE_TOL = 1.0d-6
    real(dp), parameter :: H_DRIFT_TOL = 5.0d-5
    real(dp), parameter :: PPHI_CONTROL_TOL = 1.0d-12
    real(dp), parameter :: PPHI_SEGMENT_TOL = 1.0d-6
    real(dp), parameter :: PPHI_JUMP_EPS = 1.0d-9
    real(dp), parameter :: LANDING_TOL = 1.0d-8
    integer, parameter :: MIN_CROSSINGS = 2000

    character(len=1024) :: h5file
    character(len=256) :: config_file
    type(tracer_t) :: norb
    real(dp) :: h_series(NREC), p_series(NREC), t_series(NREC)
    real(dp) :: slope, serr, drift, seg_drift, pphi_walk, max_resid
    integer :: im, nrec_got, n_cross, n_stop, landings, stops, n_crossers
    logical :: failed, crosser

    if (command_argument_count() < 1) then
        print *, 'usage: test_spectre_sympl_crossing <tok2vol.h5>'
        error stop 1
    end if
    call get_command_argument(1, h5file)

    call write_input(trim(h5file))
    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call params_init
    call init_magfie(SPECTRE)
    call crossing_log_reset(1024)
    call sympl_landing_stats_reset

    failed = .false.
    n_crossers = 0

    call check_recanon_energy(failed)

    do im = 1, NMARKER
        call trace_marker(im, nrec_got, h_series, p_series, t_series)
        if (nrec_got < MIN_REC) then
            print '(A,I0,A,I0,A)', 'FAIL: marker ', im, ' ended after ', &
                nrec_got, ' records'
            failed = .true.
            cycle
        end if

        call segment_stats(p_series(1:nrec_got), crosser, seg_drift, pphi_walk)
        call fit_series(t_series(1:nrec_got), h_series(1:nrec_got), slope, serr, drift)

        if (crosser) then
            n_crossers = n_crossers + 1
            print '(A,I0,A,I0,A,ES11.3,A,ES11.3,A,ES11.3)', 'crosser ', im, &
                ' (', nrec_got, ' rec): H slope*T/H0 = ', slope, ' +- ', &
                3.0_dp*serr, '  max|dH/H| = ', drift
            print '(A,ES11.3,A,ES11.3)', '    pphi: in-segment drift = ', &
                seg_drift, '  total walk (map kicks) = ', pphi_walk
            if (abs(slope) > max(3.0_dp*serr, SLOPE_TOL)) then
                print '(A,I0)', 'FAIL: secular H drift for crossing marker ', im
                failed = .true.
            end if
            if (drift >= H_DRIFT_TOL) then
                print '(A,I0)', 'FAIL: H excursion above floor for marker ', im
                failed = .true.
            end if
            if (seg_drift >= PPHI_SEGMENT_TOL) then
                print '(A,I0)', 'FAIL: p_phi drifts between crossings, marker ', &
                    im
                failed = .true.
            end if
        else
            print '(A,I0,A,I0,A,ES11.3,A,ES11.3)', 'control ', im, ' (', nrec_got, &
                ' rec): max|dH/H| = ', drift, '  max|dpphi/pphi| = ', pphi_walk
            if (pphi_walk >= PPHI_CONTROL_TOL) then
                print '(A,I0)', 'FAIL: in-volume p_phi not conserved, marker ', &
                    im
                failed = .true.
            end if
            if (drift >= H_DRIFT_TOL) then
                print '(A,I0)', 'FAIL: in-volume H floor exceeded, marker ', im
                failed = .true.
            end if
        end if
    end do

    n_cross = crossing_log_count_type(CROSS_CROSSING)
    n_stop = crossing_log_count_type(CROSS_STOP)
    call sympl_landing_stats(landings, max_resid, stops)

    print '(A,I0,A,I0,A,I0,A,I0)', 'events: crossings=', n_cross, &
        ' landings=', landings, ' cross_stop=', n_stop, &
        ' crossing_markers=', n_crossers
    print '(A,ES12.4)', 'landing: max |rho_g - k| = ', max_resid

    if (n_crossers < 2) then
        print '(A,I0)', 'FAIL: fewer than 2 crossing markers: ', n_crossers
        failed = .true.
    end if
    if (n_cross < MIN_CROSSINGS) then
        print '(A,I0,A,I0)', 'FAIL: crossings ', n_cross, ' < ', MIN_CROSSINGS
        failed = .true.
    end if
    if (n_stop /= 0 .or. stops /= 0) then
        print '(A,I0)', 'FAIL: CROSS_STOP events: ', max(n_stop, stops)
        failed = .true.
    end if
    if (max_resid >= LANDING_TOL) then
        print '(A,ES12.4)', 'FAIL: landing residual above 1e-8: ', max_resid
        failed = .true.
    end if

    if (failed) error stop 1
    print *, 'symplectic multi-volume invariants PASS'

contains

    subroutine write_input(h5)
        character(*), intent(in) :: h5
        integer :: unit

        open (newunit=unit, file='simple.in', status='replace', action='write')
        write (unit, '(A)') '&config'
        write (unit, '(A)') "  field_input = '"//h5//"'"
        write (unit, '(A)') '  integ_coords = 6'
        write (unit, '(A)') '  integmode = 3'
        write (unit, '(A)') '  ntestpart = 8'
        write (unit, '(A)') '  ntimstep = 100'
        write (unit, '(A)') '  npoiper2 = 4096'
        write (unit, '(A)') '  relerr = 1d-13'
        write (unit, '(A)') '  facE_al = 500.0d0'
        write (unit, '(A)') '  trace_time = 4.0d-5'
        write (unit, '(A)') '  sbeg = 0.97d0'
        write (unit, '(A)') '/'
        close (unit)
    end subroutine write_input

    subroutine check_recanon_energy(failed)
        !> Consistency guard on the canonical-momentum reconstruction: for a
        !> physical state (p, lambda) the p_phi from recanon_pphi, decoded by
        !> get_val in the same volume, must return H = p^2 (mu and vpar both
        !> built from (p, lambda), so vpar^2/2 + mu*B collapses to p^2). This
        !> pins the hph and sqrt2 factors of the reconstruction. It also prints
        !> Aph/ro0: on this axisymmetric Meiss fixture that gauge term is ~0
        !> (canonical A_phi is gauged away, and zeta is ignorable), so the
        !> substantive re-canonicalization content is the target-volume field
        !> re-evaluation, red-proofed by the crossing-marker H-drift witnesses
        !> below (skipping it drives their H slope from ~1e-8 to ~1e-2).
        logical, intent(inout) :: failed

        type(field_can_t) :: f
        real(dp), parameter :: TOL = 1.0d-11
        real(dp) :: p, lambda, rho, pphi, relerr_h
        integer :: lvol

        p = 0.8_dp
        lambda = 0.4_dp
        do lvol = 1, spectre_mvol
            call set_spectre_volume_lock(lvol)
            rho = real(lvol, dp) - 0.5_dp
            call eval_field(f, rho, 0.7_dp, 0.3_dp, 0)
            f%ro0 = ro0/sqrt2
            f%mu = p**2*(1.0_dp - lambda**2)/f%Bmod
            pphi = recanon_pphi(f, p, lambda)
            call get_val(f, pphi)
            relerr_h = abs(f%H - p**2)/p**2
            print '(A,I0,A,ES11.3,A,ES11.3,A,ES11.3)', 'recanon vol ', lvol, &
                ': |H - p^2|/p^2 = ', relerr_h, '  Aph/ro0 = ', f%Aph/f%ro0, &
                '  vpar = ', f%vpar
            if (relerr_h > TOL) then
                print '(A,I0)', 'FAIL: recanon energy identity broken, vol ', lvol
                failed = .true.
            end if
        end do
        call set_spectre_volume_lock(0)
    end subroutine check_recanon_energy

    subroutine trace_marker(im, nrec_got, h_series, p_series, t_series)
        integer, intent(in) :: im
        integer, intent(out) :: nrec_got
        real(dp), intent(out) :: h_series(NREC), p_series(NREC), t_series(NREC)

        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f, ftmp
        type(sympl_spectre_state_t) :: state
        real(dp) :: zphys(5), z(5), lam, th0, t_frac
        integer :: k, ierr

        lam = LAM0
        if (mod(im, 2) == 0) lam = -LAM0
        th0 = twopi*real((im - 1)/2, dp)/real(NMARKER/2, dp)

        zphys = [RHO0, th0, 0.0_dp, 1.0_dp, lam]
        call ref_to_integ(zphys(1:3), z(1:3))
        z(4:5) = zphys(4:5)
        call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)
        call sympl_spectre_reset(state, si, spectre_mvol, integmode, &
                                 crossing_level)

        nrec_got = 0
        do k = 1, NSTEP
            call orbit_microstep_sympl_spectre(state, si, f, im, &
                                               real(k - 1, dp)*dtaumin/v0, &
                                               dtaumin/v0, ierr, t_frac)
            if (ierr /= SYMPL_SPECTRE_OK) exit
            if (mod(k, NSAMPLE) == 0) then
                ftmp = f
                call eval_field(ftmp, si%z(1), si%z(2), si%z(3), 0)
                call get_val(ftmp, si%z(4))
                nrec_got = nrec_got + 1
                h_series(nrec_got) = ftmp%H
                p_series(nrec_got) = si%z(4)
                t_series(nrec_got) = real(k, dp)*dtaumin/v0
            end if
        end do
        call set_spectre_volume_lock(0)
    end subroutine trace_marker

    subroutine segment_stats(p, crosser, seg_drift, walk)
        !> Split the p_phi series into constant segments at every relative jump
        !> above PPHI_JUMP_EPS (crossing kicks and chatter pairs; in-volume
        !> conservation is exact to round-off by axisymmetry). seg_drift is the
        !> peak relative deviation within any segment, walk the peak relative
        !> excursion of the whole series.
        real(dp), intent(in) :: p(:)
        logical, intent(out) :: crosser
        real(dp), intent(out) :: seg_drift, walk

        real(dp) :: p0, seg_ref
        integer :: k

        p0 = abs(p(1))
        seg_ref = p(1)
        seg_drift = 0.0_dp
        crosser = .false.
        do k = 2, size(p)
            if (abs(p(k) - p(k - 1)) > PPHI_JUMP_EPS*p0) then
                crosser = .true.
                seg_ref = p(k)
            else
                seg_drift = max(seg_drift, abs(p(k) - seg_ref)/p0)
            end if
        end do
        walk = maxval(abs(p - p(1)))/p0
    end subroutine segment_stats

    subroutine fit_series(t, y, slope_rel, serr_rel, drift)
        !> Least-squares slope of y(t) scaled to relative drift over the span
        !> (slope*T/|y(1)|), its 1-sigma error from the residual scatter, and
        !> the peak relative excursion from the initial value.
        real(dp), intent(in) :: t(:), y(:)
        real(dp), intent(out) :: slope_rel, serr_rel, drift

        real(dp) :: tm, ym, stt, sty, ssr, res, span, slope
        integer :: k, n

        n = size(t)
        tm = sum(t)/n
        ym = sum(y)/n
        stt = 0.0_dp
        sty = 0.0_dp
        do k = 1, n
            stt = stt + (t(k) - tm)**2
            sty = sty + (t(k) - tm)*(y(k) - ym)
        end do
        slope = sty/stt

        ssr = 0.0_dp
        do k = 1, n
            res = y(k) - ym - slope*(t(k) - tm)
            ssr = ssr + res**2
        end do
        span = t(n) - t(1)
        slope_rel = slope*span/abs(y(1))
        serr_rel = sqrt(ssr/real(n - 2, dp)/stt)*span/abs(y(1))
        drift = maxval(abs(y - y(1)))/abs(y(1))
    end subroutine fit_series

end program test_spectre_sympl_crossing
