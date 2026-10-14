program test_spectre_fo_interfaces
    !> Full orbit near SPECTRE interfaces (crossing_level = 3) on the
    !> axisymmetric current-sheet fixture tok2vol.
    !>
    !> Markers start next to the interface at 3.5 keV, where the gyroradius is a
    !> few percent of the volume width, so the guiding centre is valid away from
    !> the interface and the orbit switches between guiding centre and full
    !> orbit many times. Witnesses:
    !>
    !>   * the full orbit is entered and left, without failures or losses;
    !>   * the momentum shell is kept: |p| in the guiding-centre phases and
    !>     |v|^2/2 = p^2 in the full-orbit phases;
    !>   * the canonical toroidal momentum, exactly conserved physically, stays
    !>     within the integrator's pzeta_tol (relative) over all full-orbit
    !>     phases of a marker in the same volume. It is sampled in SPECTRE's gauge
    !>     of that volume, so its spread measures the switching and guiding-
    !>     centre errors accumulated between the phases.
    !>
    !> usage: test_spectre_fo_interfaces <tok2vol.h5>
    !> RK45 on the sharp field with exact interface crossings.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use simple, only: init_sympl
    use simple_api, only: simple_read_config, simple_init
    use params, only: dtaumin, relerr, v0, crossing_level, integmode
    use magfie_sub, only: init_magfie, SPECTRE
    use field_can_mod, only: field_can_t, ref_to_integ
    use field_can_spectre, only: spectre_mvol
    use orbit_symplectic_base, only: symplectic_integrator_t
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
        orbit_microstep_sympl_spectre, sympl_landing_stats_reset, sympl_fo_stats, &
        SYMPL_SPECTRE_OK
    use spectre_fo_hybrid, only: spectre_fo_canonical_pzeta, SPECTRE_FO_OK
    use parmot_mod, only: ro0
    use util, only: twopi, sqrt2

    implicit none

    integer, parameter :: NMARKER = 8, NSTEP = 20000, NSAMPLE = 4
    real(dp), parameter :: RHO0 = 0.97_dp, LAM0 = 0.5_dp
    !> RK45 bounds include adaptive FO energy drift and finite-order GC/FO
    !> conversion error; sharp-sheet events also have an independent RK4 test.
    real(dp), parameter :: shell_tol = 1.0e-4_dp, pzeta_tol = 1.0e-5_dp

    character(len=1024) :: h5file
    character(len=256) :: config_file
    character(len=256) :: message
    integer :: config_status
    real(dp) :: spread, worst_spread, worst_shell, shell
    integer :: im, entries, exits, losses, failures, status(5), phases
    logical :: failed

    if (command_argument_count() < 1) then
        print *, 'usage: test_spectre_fo_interfaces <tok2vol.h5>'
        error stop 1
    end if
    call get_command_argument(1, h5file)

    call write_input(trim(h5file))
    config_file = 'simple.in'
    call simple_read_config(config_file, config_status, message)
    if (config_status /= 0) error stop 'library configuration failed'
    call simple_init
    call init_magfie(SPECTRE)
    call sympl_landing_stats_reset

    failed = .false.
    worst_spread = 0.0_dp
    worst_shell = 0.0_dp
    do im = 1, NMARKER
        call trace_marker(im, spread, shell, phases)
        print '(A,I0,A,I0,A,ES10.3,A,ES10.3)', 'marker ', im, ': FO phases ', &
            phases, '  rel p_zeta spread ', spread, '  shell defect ', shell
        worst_spread = max(worst_spread, spread)
        worst_shell = max(worst_shell, shell)
    end do
    call sympl_fo_stats(entries, exits, losses, failures, status)
    print '(A,4(I0,1X),A,5(I0,1X))', 'FO entries exits losses failures: ', &
        entries, exits, losses, failures, ' by status: ', status

    if (entries < NMARKER .or. exits < NMARKER) then
        print '(A)', 'FAIL: markers did not switch between GC and FO'
        failed = .true.
    end if
    if (failures /= 0 .or. losses /= 0) then
        print '(A)', 'FAIL: full-orbit failures or losses'
        failed = .true.
    end if
    if (worst_shell >= shell_tol) then
        print '(A)', 'FAIL: momentum shell not kept'
        failed = .true.
    end if
    if (worst_spread >= pzeta_tol) then
        print '(A)', 'FAIL: canonical toroidal momentum not conserved'
        failed = .true.
    end if
    if (failed) error stop 1
    print '(A)', 'PASS'

contains

    subroutine write_input(h5)
        character(*), intent(in) :: h5
        integer :: unit

        open (newunit=unit, file='simple.in', status='replace', action='write')
        write (unit, '(A)') '&config'
        write (unit, '(A)') "  field_input = '"//h5//"'"
        write (unit, '(A)') '  integ_coords = 6'
        write (unit, '(A)') '  integmode = 3'
        write (unit, '(A)') '  crossing_level = 3'
        ! Exercise RK45 through the public library path.
        write (unit, '(A)') '  spectre_fo_tolerance = 1d-9'
        write (unit, '(A)') '  spectre_ncon_phi = 32'
        write (unit, '(A)') '  ntestpart = 8'
        write (unit, '(A)') '  ntimstep = 100'
        write (unit, '(A)') '  npoiper2 = 1024'
        write (unit, '(A)') '  relerr = 1d-13'
        write (unit, '(A)') '  facE_al = 1000.0d0'
        write (unit, '(A)') '  trace_time = 1.0d-3'
        write (unit, '(A)') '/'
        close (unit)
    end subroutine write_input

    subroutine trace_marker(im, spread, shell, phases)
        !> Relative spread of the full-orbit p_zeta per owning volume, the
        !> worst momentum-shell defect, and the number of full-orbit phases.
        integer, intent(in) :: im
        real(dp), intent(out) :: spread, shell
        integer, intent(out) :: phases

        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f
        type(sympl_spectre_state_t) :: state
        real(dp) :: zphys(5), z(5), lam, th0, t_frac, pzeta
        real(dp) :: pz_min(2), pz_max(2)
        integer :: k, ierr, pz_status, owner, samples(2)
        integer :: entries0, exits0, entries1, exits1, losses, failures, counts(5)
        logical :: was_fo, visited(2)

        lam = merge(-LAM0, LAM0, mod(im, 2) == 0)
        th0 = twopi*(real((im - 1)/2, dp) + 0.25_dp)/real(NMARKER/2, dp)
        zphys = [RHO0, th0, 0.0_dp, 1.0_dp, lam]
        call ref_to_integ(zphys(1:3), z(1:3))
        z(4:5) = zphys(4:5)
        call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)
        call sympl_spectre_reset(state, si, spectre_mvol, integmode, &
                                 crossing_level)

        pz_min = huge(1.0_dp)
        pz_max = -huge(1.0_dp)
        shell = 0.0_dp
        phases = 0
        was_fo = .false.
        samples = 0
        visited = .false.
        call sympl_fo_stats(entries0, exits0, losses, failures, counts)
        do k = 1, NSTEP
            call orbit_microstep_sympl_spectre(state, si, f, im, &
                real(k - 1, dp)*dtaumin/v0, dtaumin/v0, ierr, t_frac)
            if (ierr /= SYMPL_SPECTRE_OK) then
                print '(A,I0,A,I0)', 'marker ', im, ' stopped with ', ierr
                error stop 'incomplete marker trace'
            end if
            if (.not. ieee_is_finite(t_frac)) error stop 'invalid trace fraction'
            if (abs(t_frac - 1.0_dp) > 1.0e-12_dp) &
                error stop 'incomplete microstep'
            if (state%fo%active .and. .not. was_fo) phases = phases + 1
            was_fo = state%fo%active
            if (state%fo%active) then
                owner = state%fo%owner
                if (owner < 1 .or. owner > 2) error stop 'invalid FO owner'
                visited(owner) = .true.
            end if
            if (mod(k, NSAMPLE) /= 0) cycle
            if (state%fo%active) then
                if (.not. all(ieee_is_finite(state%fo%v))) &
                    error stop 'nonfinite full-orbit shell sample'
                shell = max(shell, abs(0.5_dp*dot_product(state%fo%v, &
                    state%fo%v) - 1.0_dp))
                call spectre_fo_canonical_pzeta(state%fo, ro0/sqrt2, pzeta, &
                    pz_status)
                if (pz_status /= SPECTRE_FO_OK) error stop 'invariant sample failed'
                if (.not. ieee_is_finite(pzeta)) error stop 'nonfinite invariant'
                owner = state%fo%owner
                if (owner < 1 .or. owner > 2) error stop 'invalid sample owner'
                visited(owner) = .true.
                samples(owner) = samples(owner) + 1
                pz_min(owner) = min(pz_min(owner), pzeta)
                pz_max(owner) = max(pz_max(owner), pzeta)
            else
                if (.not. ieee_is_finite(si%pabs)) &
                    error stop 'nonfinite guiding-centre shell sample'
                shell = max(shell, abs(si%pabs - 1.0_dp))
            end if
        end do
        if (k /= NSTEP + 1) error stop 'marker did not complete trace'
        call sympl_fo_stats(entries1, exits1, losses, failures, counts)
        print '(A,I0,A,2(I0,1X))', 'marker ', im, ' entries/exits: ', &
            entries1 - entries0, exits1 - exits0
        if (entries1 <= entries0) error stop 'marker did not switch to full orbit'
        ! Sheet-confined markers can remain in FO for the whole interval.
        ! The first pair are round-trip witnesses; every marker must enter FO.
        if (im <= 2 .and. exits1 <= exits0) &
            error stop 'round-trip marker did not return to guiding centre'
        if (sum(samples) < 2) error stop 'missing invariant samples'
        if (any(visited .and. samples < 2)) &
            error stop 'missing visited-volume invariant samples'
        if (.not. ieee_is_finite(shell)) error stop 'nonfinite shell defect'
        spread = 0.0_dp
        do owner = 1, 2
            if (.not. visited(owner)) cycle
            spread = max(spread, (pz_max(owner) - pz_min(owner)) &
                /max(abs(pz_max(owner)), abs(pz_min(owner))))
        end do
        if (.not. ieee_is_finite(spread)) error stop 'invalid invariant spread'
    end subroutine trace_marker

end program test_spectre_fo_interfaces
