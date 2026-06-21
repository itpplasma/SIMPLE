program diag_banana
  ! Banana-orbit and invariant-conservation comparison: trace ONE trapped alpha
  ! with the guiding centre (orbit_timestep_axis) and the explicit full-orbit Boris
  ! CP (orbit_cpp_boris) from the same start, on a reactor-scale field. Dumps the
  ! poloidal orbit (R, Z) and the invariants (energy, magnetic moment mu) along
  ! each trajectory so the FLR banana width and the energy/mu conservation can be
  ! seen side by side. Config via environment:
  !   BANANA_WOUT   path to the VMEC wout (default wout.nc)
  !   BANANA_BSCALE vmec_B_scale  (default 1)
  !   BANANA_RZSCALE vmec_RZ_scale (default 1)
  !   BANANA_TAG    output prefix  (default run) -> banana_<tag>_{gc,cp}.dat
  !   BANANA_TRACE  trace time [s] (default 5e-5)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use parmot_mod, only: ro0
  use new_vmec_stuff_mod, only: vmec_B_scale, vmec_RZ_scale
  use simple, only: init_params, init_sympl, tracer_t
  use simple_main, only: init_field
  use orbit_symplectic, only: orbit_timestep_sympl
  use orbit_cpp_boris, only: cpp_boris_state_t, cpp_boris_init, cpp_boris_step, &
    cpp_boris_energy, cpp_boris_mu, cpp_boris_to_gc
  use boozer_cartesian, only: boozer_to_cart
  use boozer_field_metric, only: boozer_field_metric_eval
  use params, only: field_input, coord_input, integmode, relerr, dtaumin, orbit_coord
  use velo_mod, only: isw_field_type
  use magfie_sub, only: BOOZER
  use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
  use boozer_sub, only: get_boozer_coordinates
  implicit none

  ! Field resolution MUST match the loss campaign (template.in: ns_s=3, ns_tp=3,
  ! multharm=5) so the banana is the same physics as the loss runs. Quintic
  ! (ns_s=ns_tp=5) on a high-mode wout (W7-X) blows the canonical-field grid to
  ! tens of GB; cubic is the production setting.
  integer, parameter :: ans_s = 3, ans_tp = 3, amultharm = 5
  type(tracer_t) :: norb
  real(dp) :: z0(5), ro0_bar, mu0, vpar_bar, vperp0, bmod, trace_s, v0
  real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3), Acov(3), dA(3,3)
  real(dp) :: Bctr(3), Bcov(3), dBmod(3), hcov(3)
  character(len=256) :: woutfile, tag, sval
  integer :: nstep_gc, nstep_cp, np2, ilam
  real(dp) :: lam, sbeg

  woutfile = getenv_def('BANANA_WOUT', 'wout.nc')
  tag = getenv_def('BANANA_TAG', 'run')
  call get_environment_variable('BANANA_BSCALE', sval)
  if (len_trim(sval) > 0) then; read(sval,*) vmec_B_scale; else; vmec_B_scale = 1.0_dp; end if
  call get_environment_variable('BANANA_RZSCALE', sval)
  if (len_trim(sval) > 0) then; read(sval,*) vmec_RZ_scale; else; vmec_RZ_scale = 1.0_dp; end if
  call get_environment_variable('BANANA_TRACE', sval)
  if (len_trim(sval) > 0) then; read(sval,*) trace_s; else; trace_s = 5.0e-5_dp; end if
  call get_environment_variable('BANANA_NP', sval)
  if (len_trim(sval) > 0) then; read(sval,*) np2; else; np2 = 16384; end if
  call get_environment_variable('BANANA_LAMBDA', sval)
  if (len_trim(sval) > 0) then; read(sval,*) lam; else; lam = 0.2_dp; end if
  call get_environment_variable('BANANA_SBEG', sval)
  if (len_trim(sval) > 0) then; read(sval,*) sbeg; else; sbeg = 0.5_dp; end if

  isw_field_type = BOOZER
  field_input = trim(woutfile); coord_input = trim(woutfile)
  orbit_coord = 1; integmode = 1; relerr = 1.0d-10
  call init_field(norb, trim(woutfile), ans_s, ans_tp, amultharm, integmode)
  use_B_r = .true.; use_del_tp_B = .true.
  call get_boozer_coordinates
  call init_params(norb, 2, 4, 3.5e6_dp, np2, 1, 1.0d-10)
  dtaumin = norb%dtaumin
  v0 = norb%v0
  ro0_bar = ro0/sqrt(2.0_dp)

  ! one trapped particle at the requested surface and pitch
  z0 = [sbeg, 0.0_dp, 0.0_dp, 1.0_dp, lam]
  call boozer_field_metric_eval(z0(1:3), g, ginv, sqrtg, dg, Acov, dA, &
       Bctr, Bcov, bmod, dBmod, hcov)
  mu0 = 0.5_dp*z0(4)**2*(1.0_dp - z0(5)**2)/bmod*2.0_dp
  vpar_bar = z0(4)*z0(5)*sqrt(2.0_dp)
  vperp0 = sqrt(max(2.0_dp*mu0*bmod, 0.0_dp))

  ! same physical trace time for both; GC steps at dtaumin, CP at dtaumin/sqrt2
  nstep_gc = nint(trace_s*v0/dtaumin)
  nstep_cp = nint(trace_s*v0/(dtaumin/sqrt(2.0_dp)))

  print '(A,A)',    '  config tag      = ', trim(tag)
  print '(A,A)',    '  wout            = ', trim(woutfile)
  print '(A,2F10.5)', '  B_scale RZ_scale= ', vmec_B_scale, vmec_RZ_scale
  print '(A,ES12.4)', '  |B| at start (G)= ', bmod
  print '(A,ES12.4)', '  ro0 (cm)        = ', ro0
  print '(A,2I9)',  '  nstep gc / cp   = ', nstep_gc, nstep_cp

  call trace_gc(z0, mu0, nstep_gc, trim(tag))
  call trace_cp(z0, vpar_bar, vperp0, mu0, ro0_bar, nstep_cp, trim(tag))

contains

  function getenv_def(name, def) result(val)
    character(*), intent(in) :: name, def
    character(len=256) :: val
    call get_environment_variable(name, val)
    if (len_trim(val) == 0) val = def
  end function getenv_def

  subroutine trace_gc(z0, mu0, nstep, tag)
    real(dp), intent(in) :: z0(5), mu0
    integer, intent(in) :: nstep
    character(*), intent(in) :: tag
    real(dp) :: zinit(5), xyz(3), Jc(3,3), R, Zc, E, E0, Emax, mu_now
    integer :: it, ierr, u, sub

    ! production symplectic GC (orbit_timestep_sympl); state in norb%si%z, the GC
    ! invariants in norb%f (mu fixed parameter, Bmod and vpar updated each step).
    zinit = z0
    call init_sympl(norb%si, norb%f, zinit, dtaumin, dtaumin, relerr, integmode)
    open(newunit=u, file='banana_'//tag//'_gc.dat', status='replace')
    write(u,'(A)') '# t_tau s theta R Z energy mu'
    E0 = -1.0_dp; Emax = 0.0_dp; sub = max(1, nstep/4000)
    do it = 1, nstep
      call orbit_timestep_sympl(norb%si, norb%f, ierr)
      if (ierr /= 0) exit
      E = norb%f%mu*norb%f%Bmod + 0.5_dp*norb%f%vpar**2   ! GC energy
      mu_now = norb%f%mu                                  ! fixed GC invariant
      if (E0 < 0.0_dp) E0 = E
      Emax = max(Emax, abs((E - E0)/E0))
      if (mod(it, sub) == 0) then
        call boozer_to_cart(norb%si%z(1:3), xyz, Jc)
        R = sqrt(xyz(1)**2 + xyz(2)**2); Zc = xyz(3)
        write(u,'(7ES16.8)') it*dtaumin, norb%si%z(1), norb%si%z(2), R, Zc, E, mu_now
      end if
    end do
    close(u)
    print '(A,ES10.2,A)', '  GC  max|dE/E0| = ', Emax, '   (mu exact: GC parameter)'
  end subroutine trace_gc

  subroutine trace_cp(z0, vpar_bar, vperp0, mu0, ro0_bar, nstep, tag)
    real(dp), intent(in) :: z0(5), vpar_bar, vperp0, mu0, ro0_bar
    integer, intent(in) :: nstep
    character(*), intent(in) :: tag
    type(cpp_boris_state_t) :: st
    real(dp) :: R, Zc, E, E0, Emax, mu_now, mumin, mumax, s, th, ph, vpar
    real(dp) :: gw, gsum, gn, mu_first, mu_last, mu_drift
    integer :: it, ierr, u, sub, nwin

    call cpp_boris_init(st, .false., z0(1:3), vpar_bar, vperp0, mu0, 1.0_dp, &
                        1.0_dp, dtaumin/sqrt(2.0_dp), ro0_bar, z0(4))
    open(newunit=u, file='banana_'//tag//'_cp.dat', status='replace')
    write(u,'(A)') '# t_tau s theta R Z energy mu'
    E0 = cpp_boris_energy(st); Emax = 0.0_dp
    mu_now = cpp_boris_mu(st); mumin = mu_now; mumax = mu_now
    sub = max(1, nstep/4000)
    ! gyro-averaged secular drift: average mu over a ~1-gyration window at the
    ! start and at the end. The gyroperiod ~ 2 pi ro0_bar/|B| ~ a few dt; window
    ! = max(50, nstep/200) steps spans several gyrations, averaging out the
    ! gyro-phase oscillation to isolate the true (adiabatic) mu change.
    nwin = max(50, nstep/200)
    mu_first = 0.0_dp; mu_last = 0.0_dp; gn = 0.0_dp
    do it = 1, nstep
      call cpp_boris_step(st, ierr)
      if (ierr /= 0) exit
      E = cpp_boris_energy(st); mu_now = cpp_boris_mu(st)
      Emax = max(Emax, abs((E - E0)/E0))
      mumin = min(mumin, mu_now); mumax = max(mumax, mu_now)
      if (it <= nwin) mu_first = mu_first + mu_now
      if (it > nstep - nwin) then; mu_last = mu_last + mu_now; gn = gn + 1.0_dp; end if
      if (mod(it, sub) == 0) then
        call cpp_boris_to_gc(st, s, th, ph, vpar, ierr)
        R = sqrt(st%x(1)**2 + st%x(2)**2); Zc = st%x(3)
        write(u,'(7ES16.8)') it*dtaumin/sqrt(2.0_dp), s, th, R, Zc, E, mu_now
      end if
    end do
    close(u)
    mu_drift = -1.0_dp
    if (gn > 0.0_dp) mu_drift = abs(mu_last/gn - mu_first/nwin)/(mu_first/nwin)
    print '(A,ES10.2,A,F6.3,A,ES10.2)', '  CP  max|dE/E0| = ', Emax, &
      '   mu gyro-osc band = ', (mumax-mumin)/mu0, &
      '   gyro-avg secular drift = ', mu_drift
  end subroutine trace_cp

end program diag_banana
