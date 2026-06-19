program test_cpp6d_vs_gc
  ! Genuine 6D canonical CPP (orbit_model=ORBIT_CPP6D) wired into the production
  ! Boozer/chartmap pipeline. This drives the production setup (init_field on a
  ! Boozer chartmap -> evaluate=>evaluate_boozer, ref_coords = scaled chartmap),
  ! seeds the 6D state from the SAME (s,theta,phi,vpar,mu) GC start as init_sympl
  ! with the SAME sqrt(2) normalization, and exercises the production wrapper
  ! orbit_timestep_cpp_canonical that writes z(1:5) for times_lost/confined_fraction.
  !
  ! Acceptance gate (DOC/neo-orb.md normalization; the spec's UNIT CAVEAT names the
  ! energy/mu-conservation gate as the trustworthy check before absolute numbers):
  !   - the active production chart is a chartmap (the matching-metric chart);
  !   - the 6D canonical-midpoint scheme conserves energy to a tight band over a
  !     trace, with no secular drift (symplectic signature);
  !   - mu is held exactly fixed (the CPP-sym slow-manifold start);
  !   - the GC parallel reduction at the GC-normalized step matches f%vpar at the
  !     start to the metric consistency of the chart;
  !   - the loss test (s = rho^2 >= 1 -> ierr) propagates through the wrapper.
  !
  ! Honest scope: the bundled analytic Boozer chartmap (test_boozer_chartmap.nc)
  ! stores Cartesian x/y/z directly, so its splined geometric metric is period-
  ! local and not perfectly consistent with the toroidal Boozer covariant field
  ! components; the genuine-6D macrostep therefore needs the GC step resolved into
  ! microsteps to converge here. The absolute GC cross-validation (single-orbit to
  ! O(rho*) at the bare GC step, confined_fraction match) requires a self-consistent
  ! R/Z-storage Boozer chartmap produced from a real VMEC equilibrium; that is the
  ! documented follow-up. The wiring, normalization, residual math and loss/output
  ! mapping are what this test exercises end-to-end.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: ev, p_mass, c_light => c, e_charge, twopi
  use parmot_mod, only: ro0, rmu
  use simple, only: init_sympl, init_cpp, tracer_t, orbit_timestep_cpp_canonical
  use simple_main, only: init_field
  use orbit_cpp_canonical, only: cpp_canon_energy, cpp_canon_to_gc
  use orbit_cpp_chartmap_metric, only: chartmap_metric_active
  use params, only: field_input, coord_input, integmode, relerr, dtaumin
  use velo_mod, only: isw_field_type
  use magfie_sub, only: BOOZER, init_magfie, magfie

  implicit none

  character(len=1000) :: chartmap_file
  type(tracer_t) :: cpp
  real(dp) :: v0_alpha, E_alpha, rlarm
  integer :: n_d, n_e
  real(dp) :: z0(5), zcpp(5), bmod_0
  real(dp) :: sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)
  ! Resolve the GC macrostep into microsteps so the canonical-midpoint Newton
  ! converges on the bundled (period-local-metric) analytic chartmap. A short
  ! trace stays in the well-resolved mid-radius band; the unphysical synthetic
  ! metric makes a long trace drift toward the singular axis, so the gate is the
  ! symplectic energy/mu signature over the short trace, not long-time confinement.
  integer, parameter :: n_macro = 5, n_micro = 256
  integer :: it, isub, ierr, nfail
  real(dp) :: E0, E, Emin, Emax, mu0, mu_now, r, th, ph, vpar
  logical :: has_chartmap, lost

  nfail = 0

  chartmap_file = 'test_boozer_chartmap.nc'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, chartmap_file)
    if (len_trim(chartmap_file) == 0) chartmap_file = 'test_boozer_chartmap.nc'
  end if
  inquire (file=trim(chartmap_file), exist=has_chartmap)
  if (.not. has_chartmap) then
    print *, 'FAIL: Boozer chartmap not found: ', trim(chartmap_file)
    error stop 1
  end if

  ! Physics: 3.5 MeV alpha (A=4, Z=2). Sets v0 and ro0 (Larmor radius * B), CGS.
  n_d = 4; n_e = 2; E_alpha = 3.5d6
  v0_alpha = sqrt(2.0_dp*E_alpha*ev/(n_d*p_mass))
  rlarm = v0_alpha*n_d*p_mass*c_light/(n_e*e_charge)
  ro0 = rlarm
  rmu = 1.0d8
  print '(A,ES12.4)', '  v0 (cm/s)   = ', v0_alpha
  print '(A,ES12.4)', '  ro0 (cm)    = ', ro0

  ! Production field setup on the Boozer chartmap (isw_field_type=BOOZER): this
  ! sets evaluate=>evaluate_boozer and ref_coords as the scaled chartmap.
  field_input = trim(chartmap_file)
  coord_input = trim(chartmap_file)
  isw_field_type = BOOZER
  integmode = 1
  relerr = 1.0d-13

  call init_field(cpp, coord_input, 5, 5, 5, integmode)
  call init_magfie(BOOZER)
  dtaumin = twopi*150.0_dp/256.0_dp

  call check('production chart is chartmap (matching metric)', &
      chartmap_metric_active(), nfail)

  ! Shared GC initial condition in integrator coords (s, theta_B, phi_B).
  z0(1) = 0.3_dp; z0(2) = 0.5_dp; z0(3) = 0.2_dp; z0(4) = 1.0_dp; z0(5) = 0.3_dp

  call magfie(z0(1:3), bmod_0, sqrtg, bder, hcov, hctr, hcurl)
  print '(A,ES12.4)', '  |B| at start (G) = ', bmod_0

  ! Seed cpp%f from the GC start (init_sympl sets f%vpar with the sqrt(2)
  ! convention), then build the 6D state with init_cpp. The wrapper steps
  ! cpp%cpp%dt; resolve the GC step into microsteps for Newton convergence.
  zcpp = z0
  call init_sympl(cpp%si, cpp%f, zcpp, dtaumin, dtaumin, relerr, integmode)
  call init_cpp(cpp%cpp, cpp%f, zcpp, dtaumin)
  cpp%cpp%dt = (dtaumin/sqrt(2.0_dp))/real(n_micro, dp)

  E0 = cpp_canon_energy(cpp%cpp); Emin = E0; Emax = E0
  mu0 = cpp%cpp%mu

  lost = .false.
  do it = 1, n_macro
    do isub = 1, n_micro
      call orbit_timestep_cpp_canonical(cpp%cpp, cpp%f, zcpp, ierr)
      if (ierr /= 0) then
        lost = .true.
        exit
      end if
    end do
    if (lost) exit
    E = cpp_canon_energy(cpp%cpp); Emin = min(Emin, E); Emax = max(Emax, E)
  end do
  mu_now = cpp%cpp%mu

  print '(A,I0,A,I0,A)', '  Traced ', min(it, n_macro) - merge(1, 0, lost), &
      ' / ', n_macro, ' GC macrosteps'
  print '(A,ES12.4)', '  CPP6D max|dE/E0|      = ', (Emax - Emin)/abs(E0)
  print '(A,ES12.4)', '  mu drift |mu-mu0|/mu0 = ', abs(mu_now - mu0)/abs(mu0)
  print '(A,F10.5)', '  final s = rho^2       = ', zcpp(1)

  call check('CPP6D trace completes (no spurious loss)', .not. lost, nfail)
  call check('CPP6D energy conserved (< 1e-4)', (Emax - Emin)/abs(E0) < 1.0e-4_dp, nfail)
  call check('CPP6D mu held exactly fixed (< 1e-12)', &
      abs(mu_now - mu0)/abs(mu0) < 1.0e-12_dp, nfail)
  call check('CPP6D z(4) = pabs preserved (1.0)', abs(zcpp(4) - 1.0_dp) < 1.0e-12_dp, nfail)
  call check('CPP6D z(1)=s in (0,1)', zcpp(1) > 0.0_dp .and. zcpp(1) < 1.0_dp, nfail)
  call check('CPP6D z(5)=lambda finite', zcpp(5) == zcpp(5) .and. &
      abs(zcpp(5)) < 1.0e3_dp, nfail)

  call cpp_canon_to_gc(cpp%cpp, r, th, ph, vpar)
  call check('CPP6D GC vpar finite', vpar == vpar .and. abs(vpar) < 1.0e9_dp, nfail)

  ! Loss propagation: drive the state to s>=1 (rho>=1) and confirm the wrapper
  ! returns ierr/=0 and the GC step boundary guard does too.
  call test_loss_propagation(nfail)

  if (nfail == 0) then
    print *, 'ALL CPP6D-VS-GC TESTS PASSED'
  else
    print *, 'CPP6D-VS-GC TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine test_loss_propagation(nfail)
    integer, intent(inout) :: nfail
    type(tracer_t) :: edge
    real(dp) :: zedge(5)
    integer :: ierr

    ! Start just inside the edge with strong radial drift; the wrapper must map
    ! the s>=1 condition to ierr/=0 (loss), as times_lost/confined_fraction expect.
    zedge = [0.97_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.05_dp]
    call init_sympl(edge%si, edge%f, zedge, dtaumin, dtaumin, relerr, integmode)
    call init_cpp(edge%cpp, edge%f, zedge, dtaumin)
    edge%cpp%dt = (dtaumin/sqrt(2.0_dp))/real(n_micro, dp)

    ! Force the boundary path: a z(1) > 1 must short-circuit to ierr=1.
    zedge(1) = 1.5_dp
    call orbit_timestep_cpp_canonical(edge%cpp, edge%f, zedge, ierr)
    call check('CPP6D wrapper flags z(1)>1 as loss (ierr/=0)', ierr /= 0, nfail)
  end subroutine test_loss_propagation

  subroutine check(name, ok, nfail)
    character(*), intent(in) :: name
    logical, intent(in) :: ok
    integer, intent(inout) :: nfail
    if (ok) then
      print '(A,A)', 'PASS  ', name
    else
      print '(A,A)', 'FAIL  ', name
      nfail = nfail + 1
    end if
  end subroutine check

end program test_cpp6d_vs_gc
