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
    !>     within PZETA_TOL (relative) over all full-orbit phases of a marker in
    !>     the same volume. It is sampled in the full orbit in SPECTRE's gauge
    !>     of that volume, so its spread measures the switching and guiding-
    !>     centre errors accumulated between the phases.
    !>
    !> usage: test_spectre_fo_interfaces <tok2vol.h5> <boris|rk45>
    !> Boris runs on a layer of 0.2 gyroradii with the step bound 0.0125, RK45
    !> on the sharp field with exact interface crossings.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t, init_sympl
    use simple_main, only: init_field
    use params, only: read_config, params_init, netcdffile, ns_s, ns_tp, &
        multharm, dtaumin, relerr, v0, crossing_level, integmode, &
        spectre_fo_rotation, spectre_fo_layer_larmor, spectre_fo_integrator, &
        spectre_fo_tolerance
    use magfie_sub, only: init_magfie, SPECTRE
    use field_can_mod, only: field_can_t, ref_to_integ
    use field_can_spectre, only: spectre_mvol
    use orbit_symplectic_base, only: symplectic_integrator_t
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
        orbit_microstep_sympl_spectre, sympl_landing_stats_reset, sympl_fo_stats, &
        SYMPL_SPECTRE_OK
    use spectre_fo_hybrid, only: spectre_fo_canonical_pzeta, spectre_fo_configure, &
        SPECTRE_FO_OK, FO_BORIS, FO_RK45
    use interface_crossing, only: crossing_log_reset
    use parmot_mod, only: ro0
    use util, only: twopi, sqrt2

    implicit none

    integer, parameter :: NMARKER = 8, NSTEP = 20000, NSAMPLE = 4
    real(dp), parameter :: RHO0 = 0.97_dp, LAM0 = 0.5_dp
    !> Boris: measured p_zeta spread 0.93 % (the sharp field at the default
    !> bound 0.05 reaches 14 %), shell 2e-13. RK45 at tolerance 1e-9: p_zeta
    !> spread 6e-8, 1.2e-6 for a marker with a tangential touch of the
    !> interface (taken with one volume's field), shell drift 1.6e-5.
    real(dp) :: shell_tol, pzeta_tol
    character(len=16) :: integrator

    character(len=1024) :: h5file
    character(len=256) :: config_file
    type(tracer_t) :: norb
    real(dp) :: spread, worst_spread, worst_shell, shell
    integer :: im, entries, exits, losses, failures, status(5), phases
    integer :: worst_phases
    logical :: failed

    if (command_argument_count() < 2) then
        print *, 'usage: test_spectre_fo_interfaces <tok2vol.h5> <boris|rk45>'
        error stop 1
    end if
    call get_command_argument(1, h5file)
    call get_command_argument(2, integrator)
    if (trim(integrator) == 'rk45') then
        shell_tol = 1.0e-4_dp
        pzeta_tol = 1.0e-5_dp
    else
        shell_tol = 1.0e-9_dp
        pzeta_tol = 2.0e-2_dp
    end if

    call write_input(trim(h5file), trim(integrator))
    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call params_init
    call init_magfie(SPECTRE)
    call crossing_log_reset(1024)
    call sympl_landing_stats_reset
    call spectre_fo_configure(spectre_fo_rotation, spectre_fo_layer_larmor, &
        merge(FO_RK45, FO_BORIS, spectre_fo_integrator == 'rk45'), &
        spectre_fo_tolerance)

    worst_spread = 0.0_dp
    worst_shell = 0.0_dp
    worst_phases = huge(1)
    do im = 1, NMARKER
        call trace_marker(im, spread, shell, phases)
        print '(A,I0,A,I0,A,ES10.3,A,ES10.3)', 'marker ', im, ': FO phases ', &
            phases, '  rel p_zeta spread ', spread, '  shell defect ', shell
        worst_spread = max(worst_spread, spread)
        worst_shell = max(worst_shell, shell)
        worst_phases = min(worst_phases, phases)
    end do
    call sympl_fo_stats(entries, exits, losses, failures, status)
    print '(A,4(I0,1X),A,5(I0,1X))', 'FO entries exits losses failures: ', &
        entries, exits, losses, failures, ' by status: ', status

    failed = .false.
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

    subroutine write_input(h5, integrator)
        character(*), intent(in) :: h5, integrator
        integer :: unit

        open (newunit=unit, file='simple.in', status='replace', action='write')
        write (unit, '(A)') '&config'
        write (unit, '(A)') "  field_input = '"//h5//"'"
        write (unit, '(A)') '  integ_coords = 6'
        write (unit, '(A)') '  integmode = 3'
        write (unit, '(A)') '  crossing_level = 3'
        if (integrator == 'rk45') then
            write (unit, '(A)') "  spectre_fo_integrator = 'rk45'"
            write (unit, '(A)') '  spectre_fo_tolerance = 1d-9'
        else
            write (unit, '(A)') '  spectre_fo_layer_larmor = 0.2d0'
            write (unit, '(A)') '  spectre_fo_rotation = 0.0125d0'
        end if
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
        integer :: k, ierr, pz_status, owner
        logical :: was_fo

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
        do k = 1, NSTEP
            call orbit_microstep_sympl_spectre(state, si, f, im, &
                real(k - 1, dp)*dtaumin/v0, dtaumin/v0, ierr, t_frac)
            if (ierr /= SYMPL_SPECTRE_OK) then
                print '(A,I0,A,I0)', 'marker ', im, ' stopped with ', ierr
                exit
            end if
            if (state%fo%active .and. .not. was_fo) phases = phases + 1
            was_fo = state%fo%active
            if (mod(k, NSAMPLE) /= 0) cycle
            if (state%fo%active) then
                shell = max(shell, abs(0.5_dp*dot_product(state%fo%v, &
                    state%fo%v) - 1.0_dp))
                call spectre_fo_canonical_pzeta(state%fo, ro0/sqrt2, pzeta, &
                    pz_status)
                if (pz_status /= SPECTRE_FO_OK) cycle
                owner = min(state%fo%owner, 2)
                pz_min(owner) = min(pz_min(owner), pzeta)
                pz_max(owner) = max(pz_max(owner), pzeta)
            else
                shell = max(shell, abs(si%pabs - 1.0_dp))
            end if
        end do
        spread = 0.0_dp
        do owner = 1, 2
            if (pz_max(owner) < pz_min(owner)) cycle
            spread = max(spread, (pz_max(owner) - pz_min(owner)) &
                /max(abs(pz_max(owner)), abs(pz_min(owner))))
        end do
    end subroutine trace_marker

end program test_spectre_fo_interfaces
