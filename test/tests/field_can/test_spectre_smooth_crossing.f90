program test_spectre_smooth_crossing
    !> Direct guiding-centre crossing of an interface whose field jump is below
    !> spectre_fo_jump_tol. On tok2vol (largest jump 4.5e-2) with the tolerance
    !> raised to 0.1, markers started next to the interface must cross it in
    !> both directions without any full-orbit phase, keep |p| to the landing
    !> solver tolerance and complete.
    !>
    !> usage: test_spectre_smooth_crossing <tok2vol.h5>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use simple, only: init_sympl
    use simple_api, only: simple_read_config, simple_init
    use params, only: dtaumin, relerr, crossing_level, integmode
    use magfie_sub, only: init_magfie, SPECTRE
    use field_can_mod, only: field_can_t, ref_to_integ
    use field_can_spectre, only: spectre_mvol
    use orbit_symplectic_base, only: symplectic_integrator_t
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
                                   orbit_microstep_sympl_spectre, sympl_fo_stats, &
                                   SYMPL_SPECTRE_OK
    use spectre_fo_hybrid, only: spectre_fo_smooth_interface, interface_jump
    use util, only: twopi

    implicit none

    integer, parameter :: NMARKER = 8, NSTEP = 20000
    real(dp), parameter :: RHO0 = 0.97_dp, LAM0 = 0.5_dp, SHELL_TOL = 1.0e-6_dp

    character(len=1024) :: h5file
    character(len=256) :: message
    integer :: status, im, k, ierr, home, crossings, total_crossings
    integer :: entries, exits, losses, failures, counts(5)
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    type(sympl_spectre_state_t) :: state
    real(dp) :: zphys(5), z(5), t_frac, shell

    if (command_argument_count() < 1) error stop 'usage: <tok2vol.h5>'
    call get_command_argument(1, h5file)
    call write_input(trim(h5file))
    call simple_read_config('simple.in', status, message)
    if (status /= 0) error stop 'library configuration failed'
    call simple_init
    call init_magfie(SPECTRE)

    print '(A,ES10.3)', 'interface jump: ', interface_jump(1)
    if (.not. spectre_fo_smooth_interface(1)) &
        error stop 'interface not treated as smooth at tolerance 0.1'

    total_crossings = 0
    shell = 0.0_dp
    do im = 1, NMARKER
        zphys = [RHO0, twopi*(real((im - 1)/2, dp) + 0.25_dp)/real(NMARKER/2, dp), &
                 0.0_dp, 1.0_dp, merge(-LAM0, LAM0, mod(im, 2) == 0)]
        call ref_to_integ(zphys(1:3), z(1:3))
        z(4:5) = zphys(4:5)
        call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)
        call sympl_spectre_reset(state, si, spectre_mvol, integmode, crossing_level)
        home = nint(state%home_hi)
        crossings = 0
        do k = 1, NSTEP
            call orbit_microstep_sympl_spectre(state, si, f, ierr, t_frac)
            if (ierr /= SYMPL_SPECTRE_OK) error stop 'marker stopped'
            if (abs(t_frac - 1.0_dp) > 1.0e-12_dp) error stop 'incomplete microstep'
            if (.not. ieee_is_finite(si%pabs)) error stop 'nonfinite momentum'
            shell = max(shell, abs(si%pabs - 1.0_dp))
            if (nint(state%home_hi) /= home) then
                crossings = crossings + 1
                home = nint(state%home_hi)
            end if
        end do
        call sympl_fo_stats(entries, exits, losses, failures, counts)
        print '(A,I0,A,I0,A,I0)', 'marker ', im, ' crossings: ', crossings, &
            '  full-orbit entries so far: ', entries
        total_crossings = total_crossings + crossings
    end do
    call sympl_fo_stats(entries, exits, losses, failures, counts)
    print '(A,I0,A,ES10.3)', 'total crossings: ', total_crossings, &
        '  shell defect: ', shell
    if (entries /= 0) error stop 'full-orbit entries at a smooth interface'
    if (total_crossings < NMARKER) error stop 'too few direct crossings'
    if (shell > SHELL_TOL) error stop 'momentum shell not kept'
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
        write (unit, '(A)') '  spectre_fo_jump_tol = 0.1d0'
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

end program test_spectre_smooth_crossing
