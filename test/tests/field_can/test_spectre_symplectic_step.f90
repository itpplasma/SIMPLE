program test_spectre_symplectic_step
    !> Symplecticity of one interface-crossing microstep (crossing_level 2).
    !>
    !> On the pure-magnitude current sheet of tok2vol_sheet the canonical
    !> charts of the two volumes join presymplectically, so the crossing step
    !> (midpoint landing, flow along the interface field X_s, remainder
    !> midpoint, flow length from matched discrete energies) is an exactly
    !> symplectic map from (theta, phi, p_theta, p_phi) of the home chart to
    !> those of the target chart. Its central-difference Jacobian M satisfies
    !> M^T J M = J up to the O(eps^2) difference error. The Level-1 map on the
    !> same states (identity transfer, exact-energy jump) is not symplectic and
    !> serves as the control: its defect must stay far above that bound.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t, init_sympl
    use simple_main, only: init_field
    use params, only: read_config, params_init, netcdffile, ns_s, ns_tp, &
                      multharm, dtaumin, relerr, v0, integmode
    use magfie_sub, only: init_magfie, SPECTRE
    use field_can_mod, only: field_can_t, eval_field => evaluate, &
                             get_derivatives
    use field_can_spectre, only: spectre_mvol, set_spectre_volume_lock
    use orbit_symplectic_base, only: symplectic_integrator_t
    use orbit_symplectic, only: orbit_sympl_init
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
                                   orbit_microstep_sympl_spectre, &
                                   sympl_symplectic_stats, sympl_landing_stats_reset, &
                                   SYMPL_SPECTRE_OK
    use interface_crossing, only: crossing_log_reset, &
                                  CROSSING_SYMPLECTIC
    use util, only: twopi

    implicit none

    integer, parameter :: NSEARCH = 20000, CROSSING_LEGACY = 1
    real(dp), parameter :: RHO0 = 0.99_dp, LAM0 = 0.9_dp
    !> Difference size in units of the perturbation scales. The O(eps^2)
    !> difference error lies below the solver round-off floor, which grows
    !> as 1/eps (measured 4.9e-9 at eps = 1e-2, 8.4e-8 at 1e-3).
    real(dp), parameter :: EPS = 1.0d-2
    !> Bound on the Level-2 defect and minimum Level-1 defect (measured 4.6).
    real(dp), parameter :: SYMPLECTIC_TOL = 1.0d-7
    real(dp), parameter :: CONTROL_MIN = 1.0d-1

    character(len=1024) :: h5file
    character(len=256) :: config_file
    type(tracer_t) :: norb
    type(symplectic_integrator_t) :: si_cross
    type(field_can_t) :: f_cross
    real(dp) :: u0(4), a_scale, p_scale, defect, defect_legacy
    integer :: home
    logical :: failed

    if (command_argument_count() < 1) then
        print *, 'usage: test_spectre_symplectic_step <tok2vol_sheet.h5>'
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

    call find_crossing_state(si_cross, f_cross, home)
    call canonical_state(si_cross%z, f_cross, home, u0)
    call perturbation_scales(si_cross%z, f_cross, home, a_scale, p_scale)
    print '(A,I0,A,4ES14.6)', 'crossing step from volume ', home, &
        ': (theta, phi, p_theta, p_phi) = ', u0

    defect = symplectic_defect(CROSSING_SYMPLECTIC, EPS)
    defect_legacy = symplectic_defect(CROSSING_LEGACY, EPS)
    print '(A,ES11.3)', 'level 2: max|M^T J M - J| = ', defect
    print '(A,ES11.3)', 'level 1: max|M^T J M - J| = ', defect_legacy

    failed = .false.
    if (defect >= SYMPLECTIC_TOL) then
        print '(A)', 'FAIL: crossing step of level 2 is not symplectic'
        failed = .true.
    end if
    if (defect_legacy < CONTROL_MIN) then
        print '(A)', 'FAIL: control map of level 1 shows no symplecticity defect'
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
        write (unit, '(A)') '  spectre_ncon_phi = 32'
        write (unit, '(A)') '  ntestpart = 1'
        write (unit, '(A)') '  npoiper2 = 256'
        write (unit, '(A)') '  relerr = 1d-13'
        write (unit, '(A)') '  facE_al = 500.0d0'
        write (unit, '(A)') '  trace_time = 1.0d-5'
        write (unit, '(A)') '/'
        close (unit)
    end subroutine write_input

    subroutine find_crossing_state(si, f, home)
        !> Trace a passing marker from inside the axis volume until a microstep
        !> transmits it through the interface with a single crossing event;
        !> return the state before that microstep.
        type(symplectic_integrator_t), intent(out) :: si
        type(field_can_t), intent(out) :: f
        integer, intent(out) :: home

        type(sympl_spectre_state_t) :: state
        type(symplectic_integrator_t) :: si_prev
        type(field_can_t) :: f_prev
        real(dp) :: z(5), t_frac
        integer :: k, ierr, home_prev, trans0, refl0, fall0, trans, refl, fall

        z = [RHO0, 0.75_dp*twopi, 0.0_dp, 1.0_dp, LAM0]
        call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)
        call sympl_spectre_reset(state, si, spectre_mvol, integmode, &
                                 CROSSING_SYMPLECTIC)
        do k = 1, NSEARCH
            si_prev = si
            f_prev = f
            home_prev = nint(state%home_hi)
            call sympl_symplectic_stats(trans0, refl0, fall0)
            call orbit_microstep_sympl_spectre(state, si, f, 1, 0.0_dp, &
                dtaumin/v0, ierr, t_frac)
            if (ierr /= SYMPL_SPECTRE_OK) error stop 'search marker stopped'
            call sympl_symplectic_stats(trans, refl, fall)
            if (trans - trans0 == 1 .and. refl == refl0 .and. fall == fall0) then
                si = si_prev
                f = f_prev
                home = home_prev
                return
            end if
        end do
        error stop 'no transmitting microstep found'
    end subroutine find_crossing_state

    subroutine canonical_state(z, f, lvol, u)
        !> (theta, phi, p_theta, p_phi) of the integrator state z in volume lvol.
        real(dp), intent(in) :: z(4)
        type(field_can_t), intent(in) :: f
        integer, intent(in) :: lvol
        real(dp), intent(out) :: u(4)

        type(field_can_t) :: fe

        fe = f
        call set_spectre_volume_lock(lvol)
        call eval_field(fe, z(1), z(2), z(3), 0)
        call get_derivatives(fe, z(4))
        u = [z(2), z(3), fe%pth, z(4)]
    end subroutine canonical_state

    subroutine perturbation_scales(z, f, lvol, a_scale, p_scale)
        !> Perturbation units that move the start radius by at most the gap d
        !> to the interface: p_scale = |dp_theta/dr| d for momenta, a_scale for
        !> angles from the angular slope of p_theta (capped at one radian).
        real(dp), intent(in) :: z(4)
        type(field_can_t), intent(in) :: f
        integer, intent(in) :: lvol
        real(dp), intent(out) :: a_scale, p_scale

        type(field_can_t) :: fe

        fe = f
        call set_spectre_volume_lock(lvol)
        call eval_field(fe, z(1), z(2), z(3), 0)
        call get_derivatives(fe, z(4))
        p_scale = abs(fe%dpth(1))*min(abs(z(1) - real(lvol, dp)), &
                                      abs(z(1) - real(lvol - 1, dp)))
        a_scale = min(1.0_dp, p_scale/max(abs(fe%dpth(2)), abs(fe%dpth(3)), &
                                          tiny(1.0_dp)))
    end subroutine perturbation_scales

    subroutine integrator_from_canonical(u, f, lvol, z)
        !> Invert p_theta = P(r, theta, phi, p_phi) for r in volume lvol.
        real(dp), intent(in) :: u(4)
        type(field_can_t), intent(inout) :: f
        integer, intent(in) :: lvol
        real(dp), intent(inout) :: z(4)

        real(dp) :: residual
        integer :: it

        call set_spectre_volume_lock(lvol)
        z(2:4) = [u(1), u(2), u(4)]
        do it = 1, 60
            call eval_field(f, z(1), z(2), z(3), 0)
            call get_derivatives(f, z(4))
            residual = f%pth - u(3)
            z(1) = z(1) - residual/f%dpth(1)
            if (abs(residual) <= 4.0_dp*epsilon(1.0_dp)*max(abs(u(3)), 1.0_dp)) &
                return
        end do
        error stop 'radius inversion did not converge'
    end subroutine integrator_from_canonical

    subroutine crossing_map(level, u, v, ok)
        !> One microstep from canonical home-chart state u to the canonical
        !> target-chart state v; ok if it is a single transmission.
        integer, intent(in) :: level
        real(dp), intent(in) :: u(4)
        real(dp), intent(out) :: v(4)
        logical, intent(out) :: ok

        type(sympl_spectre_state_t) :: state
        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f
        real(dp) :: z(4), t_frac
        integer :: ierr, trans0, refl0, fall0, trans, refl, fall

        si = si_cross
        f = f_cross
        z = si_cross%z
        call integrator_from_canonical(u, f, home, z)
        call orbit_sympl_init(si, f, z, si_cross%dt, 1, si_cross%rtol, integmode)
        call sympl_spectre_reset(state, si, spectre_mvol, integmode, level)
        call sympl_symplectic_stats(trans0, refl0, fall0)
        call orbit_microstep_sympl_spectre(state, si, f, 1, 0.0_dp, &
            dtaumin/v0, ierr, t_frac)
        call sympl_symplectic_stats(trans, refl, fall)
        ok = ierr == SYMPL_SPECTRE_OK .and. nint(state%home_hi) /= home
        if (level == CROSSING_SYMPLECTIC) ok = ok .and. trans - trans0 == 1 &
            .and. refl == refl0 .and. fall == fall0
        call canonical_state(si%z, f, nint(state%home_hi), v)
    end subroutine crossing_map

    real(dp) function symplectic_defect(level, eps) result(defect)
        !> max |M^T J M - J| of the central-difference Jacobian M in angles
        !> scaled by a_scale and momenta by p_scale; the same scaling on both
        !> sides leaves the symplectic condition invariant.
        integer, intent(in) :: level
        real(dp), intent(in) :: eps

        real(dp) :: m(4, 4), jm(4, 4), u(4), v_plus(4), v_minus(4), v0(4)
        real(dp) :: scale_in(4), scale_out(4)
        integer :: i
        logical :: ok

        jm = 0.0_dp
        jm(1, 3) = 1.0_dp
        jm(2, 4) = 1.0_dp
        jm(3, 1) = -1.0_dp
        jm(4, 2) = -1.0_dp
        call crossing_map(level, u0, v0, ok)
        if (.not. ok) error stop 'reference state does not transmit'
        scale_in = [a_scale, a_scale, p_scale, p_scale]
        scale_out = scale_in
        do i = 1, 4
            u = u0
            u(i) = u0(i) + eps*scale_in(i)
            call crossing_map(level, u, v_plus, ok)
            if (.not. ok) error stop 'perturbed state left the transmission branch'
            u(i) = u0(i) - eps*scale_in(i)
            call crossing_map(level, u, v_minus, ok)
            if (.not. ok) error stop 'perturbed state left the transmission branch'
            m(:, i) = (v_plus - v_minus)/(2.0_dp*eps)/scale_out
        end do
        defect = maxval(abs(matmul(transpose(m), matmul(jm, m)) - jm))
    end function symplectic_defect

end program test_spectre_symplectic_step
