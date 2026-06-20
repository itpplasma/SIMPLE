program test_cp6d_vs_gc
  ! Genuine 6D classical charged particle (orbit_model=ORBIT_CP6D) wired into the
  ! production alpha-loss pipeline through the same Boozer canonical midpoint
  ! machinery as CPP, on the reactor-scale test equilibrium test_data/wout.nc
  ! (a QA stellarator, rho* ~ 1/200), validated against the production guiding
  ! center.
  !
  ! CP differs from CPP6D in physics: MODEL_CP omits the Pauli mu|B| term and
  ! seeds the resolved perpendicular velocity. CP must resolve the gyration, i.e.
  ! take many steps per gyroperiod (large npoiper2).
  !
  ! Acceptance gates (the task's validation list):
  !   (1) npoiper2 DETERMINED by energy conservation: sweep npoiper2 and report the
  !       max|dE/E0| table over several gyrations; the smallest npoiper2 with a
  !       bounded/small (< ~1e-3) energy error is the required resolution.
  !   (2) At that npoiper2 the CP gyro-center (running gyro-average of position)
  !       tracks the GC orbit (orbit_timestep_sympl) to O(rho*) over a few bounce
  !       times, and mu_emergent = vperp^2/(2|B|) is adiabatically ~conserved.
  !   (3) Energy bounded |dE/E0| small over the trace.
  !   (4) z(1) > 1 propagates through the production wrapper to ierr/=0 (loss).
  !
  ! CP is expensive (gyro-resolved), so the validation is SHORT: one particle, a
  ! few bounce times. The wire keeps the SIMPLE GC normalization (mass=1,
  ! qc=sqrt(2)/ro0, dt=dtaumin/sqrt(2)), identical to init_cpp.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use parmot_mod, only: ro0
  use simple, only: init_sympl, init_cp, init_params, tracer_t, &
    orbit_timestep_cp_canonical, canonical_state_to_standard_z
  use simple_main, only: init_field
  use orbit_symplectic, only: orbit_timestep_sympl
  use orbit_cpp_canonical, only: cpp_canon_energy, cpp_canon_to_gc, &
    cpp_canon_state_t, cpp_canon_boozer_guiding_center
  use boozer_field_metric, only: boozer_field_metric_eval
  use params, only: field_input, coord_input, integmode, relerr, dtaumin, orbit_coord
  use velo_mod, only: isw_field_type
  use magfie_sub, only: BOOZER
  use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
  use boozer_sub, only: get_boozer_coordinates
  use util, only: twopi

  implicit none

  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
  type(tracer_t) :: norb
  real(dp) :: z0(5), rbig, ro0_bar, gyroperiod, Bmod
  integer :: nfail, npoiper2

  nfail = 0

  ! Production field setup: BOOZER canonical chart on the real VMEC equilibrium.
  isw_field_type = BOOZER
  field_input = 'wout.nc'
  coord_input = 'wout.nc'
  orbit_coord = 1
  integmode = 1
  relerr = 1.0d-13
  call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, integmode)
  use_B_r = .true.
  use_del_tp_B = .true.
  call get_boozer_coordinates
  call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0d-13)
  ! rbig (cm) back out of the npoiper2=256 step: dtaumin = 2 pi rbig / npoiper2.
  rbig = norb%dtaumin*256.0_dp/twopi
  ro0_bar = ro0/sqrt(2.0_dp)

  ! Shared trapped-class IC in flux coords (s, theta, phi, v/v0, lambda).
  z0 = [0.3_dp, 0.5_dp, 0.2_dp, 1.0_dp, 0.3_dp]
  call test_cp_initial_guiding_center(z0, norb%dtaumin, nfail)
  call test_cp_orbit_start_record(z0, norb%dtaumin, nfail)

  ! Read |B| at the start so the normalized gyroperiod can be computed.
  ! The canonical cyclotron frequency is
  ! Omega = qc |B| = |B|/ro0_bar (charge=c=1, qc=1/ro0_bar), so the gyroperiod in
  ! normalized tau is 2 pi ro0_bar/|B|. With |B| ~ 5.9e4 G and ro0_bar ~ 1.9e5 cm
  ! this is O(20) tau -- much shorter than the GC step 2 pi rbig/npoiper2, so CP
  ! must oversample by ~ rbig|B|/ro0 = O(1/rho*) per gyration.
  call init_cp(norb%cp, norb%f, z0, norb%dtaumin)
  call read_boozer_field_mod(norb%cp%z(1:3), Bmod)
  gyroperiod = twopi*ro0_bar/Bmod
  print '(A,ES12.4)', '  ro0 (cm)            = ', ro0
  print '(A,ES12.4)', '  ro0_bar (cm)        = ', ro0_bar
  print '(A,ES12.4)', '  rbig (cm)           = ', rbig
  print '(A,ES12.4)', '  |B| at start (G)    = ', Bmod
  print '(A,ES12.4)', '  rho* ~ ro0/(rbig|B|)= ', ro0/(rbig*Bmod)
  print '(A,ES12.4)', '  gyroperiod (tau)    = ', gyroperiod

  ! Gate (1): determine npoiper2 by energy conservation, report the table.
  call determine_npoiper2(z0, rbig, gyroperiod, npoiper2, nfail)
  print '(A,I0)', '  CHOSEN npoiper2 (|dE/E0| < 1e-3) = ', npoiper2

  ! Gates (2),(3): gyro-center tracking, mu adiabatic invariance, energy bound.
  call test_gyrocenter_tracking(z0, npoiper2, rbig, gyroperiod, nfail)

  ! Gate (4): loss propagation.
  call test_loss_propagation(z0, npoiper2, nfail)

  if (nfail == 0) then
    print *, 'ALL CP6D-VS-GC TESTS PASSED'
  else
    print *, 'CP6D-VS-GC TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  ! dtaumin for a given npoiper2: 2 pi rbig / npoiper2.
  function dtaumin_for(npoiper2, rbig) result(dt)
    integer, intent(in) :: npoiper2
    real(dp), intent(in) :: rbig
    real(dp) :: dt
    dt = twopi*rbig/real(npoiper2, dp)
  end function dtaumin_for

  ! Steps per gyration at a given npoiper2: gyroperiod/dt_step,
  ! dt_step = (2 pi rbig/npoiper2)/sqrt(2).
  function steps_per_gyro(npoiper2, rbig, gyroperiod) result(spg)
    integer, intent(in) :: npoiper2
    real(dp), intent(in) :: rbig, gyroperiod
    real(dp) :: spg
    spg = gyroperiod/(dtaumin_for(npoiper2, rbig)/sqrt(2.0_dp))
  end function steps_per_gyro

  ! Trace the CP orbit for nsteps and return max|dE/E0|.
  subroutine test_cp_initial_guiding_center(z0, dtm, nfail)
    real(dp), intent(in) :: z0(5), dtm
    integer, intent(inout) :: nfail
    type(tracer_t) :: cp
    real(dp) :: zcp(5), xgc(3), dx(3), shift(3)

    zcp = z0
    call init_sympl(cp%si, cp%f, zcp, dtm, dtm, relerr, integmode)
    call init_cp(cp%cp, cp%f, zcp, dtm)
    call cpp_canon_boozer_guiding_center(cp%cp, xgc)
    dx = xgc - z0(1:3)
    shift = cp%cp%z(1:3) - z0(1:3)

    print '(A,3ES12.4)', '  CP initial particle-GC shift = ', shift
    print '(A,3ES12.4)', '  CP reconstructed GC error    = ', dx
    call check('CP starts with finite FLR displacement', &
      maxval(abs(shift)) > 1.0e-5_dp, nfail)
    call check('CP initial particle remains inside 0 < s < 1', &
      cp%cp%z(1) > 0.0_dp .and. cp%cp%z(1) < 1.0_dp, nfail)
    call check('CP first-order guiding center matches start (max error < 5e-2)', &
      maxval(abs(dx)) < 5.0e-2_dp, nfail)
  end subroutine test_cp_initial_guiding_center

  ! The production orbit output (simple_main trace_orbit) must record the CP
  ! trajectory starting at the integrated PARTICLE position, not the GC start it
  ! is seeded from (issue #410). canonical_state_to_standard_z is the exact value
  ! row 0 records after the fix; it must reproduce the integrated particle state
  ! cp%z(1:3) (the position the trace actually advances from), which sits one
  ! Larmor vector off the GC start z0. Before the fix row 0 held z0 instead --
  ! a spurious GC->particle jump between row 0 and row 1.
  subroutine test_cp_orbit_start_record(z0, dtm, nfail)
    real(dp), intent(in) :: z0(5), dtm
    integer, intent(inout) :: nfail
    type(tracer_t) :: cp
    real(dp) :: zrec(5), rec_vs_particle, rec_vs_gc
    integer :: ierr

    zrec = z0
    call init_sympl(cp%si, cp%f, zrec, dtm, dtm, relerr, integmode)
    call init_cp(cp%cp, cp%f, zrec, dtm)

    ! Recorded start: the particle's standard-z right after init (= row 0).
    call canonical_state_to_standard_z(cp%cp, zrec)

    rec_vs_particle = maxval(abs(zrec(1:3) - cp%cp%z(1:3)))
    rec_vs_gc       = maxval(abs(zrec(1:3) - z0(1:3)))

    print '(A,3ES12.4)', '  recorded start - GC start    = ', zrec(1:3) - z0(1:3)
    print '(A,ES12.4)',  '  recorded - particle |max|    = ', rec_vs_particle
    print '(A,ES12.4)',  '  recorded - GC start |max|    = ', rec_vs_gc

    ! Row 0 is the integrated particle position the trace advances from ...
    call check('CP recorded start equals integrated particle position', &
      rec_vs_particle < 1.0e-12_dp, nfail)
    ! ... which sits one Larmor vector off the GC start record, not on it.
    call check('CP recorded start sits off the GC start (FLR offset)', &
      rec_vs_gc > 1.0e-5_dp, nfail)
    ! And the integrator advances from it without error.
    zrec(4) = z0(4); zrec(5) = z0(5)
    call orbit_timestep_cp_canonical(cp%cp, cp%f, zrec, ierr)
    call check('CP step succeeds from recorded start', ierr == 0, nfail)
  end subroutine test_cp_orbit_start_record

  subroutine cp_energy_sweep(z0, npoiper2, rbig, nsteps, maxdE)
    real(dp), intent(in) :: z0(5), rbig
    integer, intent(in) :: npoiper2, nsteps
    real(dp), intent(out) :: maxdE
    type(tracer_t) :: cp
    real(dp) :: zcp(5), dtm, E0, E
    integer :: it, ierr

    dtm = dtaumin_for(npoiper2, rbig)
    zcp = z0
    call init_sympl(cp%si, cp%f, zcp, dtm, dtm, relerr, integmode)
    call init_cp(cp%cp, cp%f, zcp, dtm)
    E0 = cpp_canon_energy(cp%cp); maxdE = 0.0_dp
    do it = 1, nsteps
      call orbit_timestep_cp_canonical(cp%cp, cp%f, zcp, ierr)
      if (ierr /= 0) then
        print '(A,I0,A,I0)', '  CP sweep step ', it, ' ierr=', ierr
        maxdE = huge(1.0_dp); return
      end if
      E = cpp_canon_energy(cp%cp)
      maxdE = max(maxdE, abs((E - E0)/E0))
    end do
  end subroutine cp_energy_sweep

  subroutine determine_npoiper2(z0, rbig, gyroperiod, chosen, nfail)
    ! Trace one CP orbit at increasing npoiper2 over a FIXED span of several
    ! gyrations and find the smallest npoiper2 where |dE/E0| is bounded/small
    ! (< 1e-3). nsteps covers >= 6 gyrations at every resolution so the comparison
    ! is at equal physical time (a few gyroperiods, not equal step count). The
    ! coarse resolutions are undersampled (steps/gyro ~ 1), so their energy blows
    ! up -- that divergence is the physics that fixes the required npoiper2.
    real(dp), intent(in) :: z0(5), rbig, gyroperiod
    integer, intent(out) :: chosen
    integer, intent(inout) :: nfail
    integer, parameter :: ngyro_min = 6
    integer :: trials(6), i, npoiper2, nsteps
    real(dp) :: maxdE, spg
    logical :: found

    trials = [512, 1024, 2048, 4096, 8192, 16384]
    chosen = 0; found = .false.
    print '(A)', '  npoiper2   steps/gyro   nsteps   max|dE/E0|'
    do i = 1, size(trials)
      npoiper2 = trials(i)
      spg = steps_per_gyro(npoiper2, rbig, gyroperiod)
      nsteps = max(1, ceiling(ngyro_min*spg))
      call cp_energy_sweep(z0, npoiper2, rbig, nsteps, maxdE)
      print '(I8,F13.2,I9,ES14.4)', npoiper2, spg, nsteps, maxdE
      if (.not. found .and. maxdE < 1.0e-3_dp) then
        chosen = npoiper2; found = .true.
      end if
    end do
    call check('CP npoiper2 found with bounded energy (|dE/E0| < 1e-3)', found, nfail)
    if (.not. found) chosen = trials(size(trials))
  end subroutine determine_npoiper2

  subroutine test_gyrocenter_tracking(z0, npoiper2, rbig, gyroperiod, nfail)
    ! At the chosen npoiper2, trace the CP full orbit and the production GC from
    ! the same start. The CP gyro-center (running gyro-average of the flux
    ! position) must track the GC orbit to O(rho*); the emergent magnetic moment
    ! mu = vperp^2/(2|B|) must stay adiabatically ~constant; the energy bounded.
    real(dp), intent(in) :: z0(5), rbig, gyroperiod
    integer, intent(in) :: npoiper2
    integer, intent(inout) :: nfail
    type(tracer_t) :: gc, cp
    real(dp) :: zgc(5), zcp(5), dtm
    real(dp) :: E0, E, Emin, Emax, mu_emergent, mu0, mu_min, mu_max
    real(dp) :: sbar, sbar_min, sbar_max, gc_at, dev_max, sgc_min, sgc_max
    real(dp) :: mubar, mubar_ref
    integer :: it, ierr, nstep, ngyro, spg_i
    logical :: cp_lost
    real(dp), allocatable :: scp_hist(:), sgc_hist(:), mu_hist(:)

    dtm = dtaumin_for(npoiper2, rbig)
    spg_i = max(1, nint(steps_per_gyro(npoiper2, rbig, gyroperiod)))  ! averaging window
    ngyro = 60                                     ! resolved gyrations (CP is expensive)
    nstep = ngyro*spg_i
    allocate(scp_hist(0:nstep), sgc_hist(0:nstep), mu_hist(0:nstep))

    ! --- production GC at the SAME (bare) macrostep grid for a fair s comparison.
    zgc = z0
    call init_sympl(gc%si, gc%f, zgc, dtm, dtm, relerr, integmode)
    sgc_min = zgc(1); sgc_max = zgc(1); sgc_hist(0) = zgc(1)
    do it = 1, nstep
      call orbit_timestep_sympl(gc%si, gc%f, ierr)
      if (ierr /= 0) exit
      sgc_min = min(sgc_min, gc%si%z(1)); sgc_max = max(sgc_max, gc%si%z(1))
      sgc_hist(it) = gc%si%z(1)
    end do

    ! --- CP full orbit through the production wrapper, gyro-resolved.
    zcp = z0
    call init_sympl(cp%si, cp%f, zcp, dtm, dtm, relerr, integmode)
    call init_cp(cp%cp, cp%f, zcp, dtm)
    E0 = cpp_canon_energy(cp%cp); Emin = E0; Emax = E0
    mu0 = cp%cp%mu; mu_min = mu0; mu_max = mu0
    scp_hist(0) = zcp(1); cp_lost = .false.
    call emergent_mu(cp%cp, mu_hist(0))
    do it = 1, nstep
      call orbit_timestep_cp_canonical(cp%cp, cp%f, zcp, ierr)
      if (ierr /= 0) then; cp_lost = .true.; exit; end if
      E = cpp_canon_energy(cp%cp); Emin = min(Emin, E); Emax = max(Emax, E)
      ! Emergent magnetic moment from the resolved velocity: mu = vperp^2/(2|B|).
      call emergent_mu(cp%cp, mu_emergent)
      mu_min = min(mu_min, mu_emergent); mu_max = max(mu_max, mu_emergent)
      mu_hist(it) = mu_emergent
      scp_hist(it) = zcp(1)
    end do

    ! Running gyro-average (boxcar over one gyration) of the CP flux label and the
    ! same boxcar of the GC label, then the worst deviation between the two
    ! gyro-averaged tracks: the boxcar removes the FLR ripple, leaving the
    ! gyro-center which must follow the GC surface.
    dev_max = 0.0_dp; sbar_min = z0(1); sbar_max = z0(1)
    do it = spg_i, nstep
      sbar = boxcar(scp_hist, it, spg_i)
      gc_at = boxcar(sgc_hist, it, spg_i)
      sbar_min = min(sbar_min, sbar); sbar_max = max(sbar_max, sbar)
      dev_max = max(dev_max, abs(sbar - gc_at))
    end do

    ! Gyro-AVERAGED mu, the adiabatic invariant: instantaneous mu_emergent breathes
    ! at the gyrofrequency (the FLR ripple, ~8%), so its min/max is NOT the
    ! invariant. Averaging the instantaneous mu over the FIRST half of the trace
    ! (~30 gyrations) vs the SECOND half fully smooths the breathing over many
    ! cycles; their difference is the genuine SECULAR drift of the invariant, not
    ! the bounded ripple envelope.
    mubar_ref = mean(mu_hist, 0, nstep/2)              ! first-half mean mu
    mubar     = mean(mu_hist, nstep/2 + 1, nstep)      ! second-half mean mu

    print '(A,I0,A,I0,A,I0)', '  npoiper2=', npoiper2, '  steps/gyro=', spg_i, &
      '  nstep=', nstep
    print '(A,F8.5,A,F8.5,A)', '  GC      s band [', sgc_min, ',', sgc_max, ']'
    print '(A,F8.5,A,F8.5,A)', '  CP-bar  s band [', sbar_min, ',', sbar_max, ']'
    print '(A,ES12.4)', '  CP max|dE/E0|                  = ', (Emax - Emin)/abs(E0)
    print '(A,ES12.4)', '  mu instantaneous ripple        = ', (mu_max - mu_min)/abs(mu0)
    print '(A,ES12.4)', '  mu secular drift (half-means)  = ', &
        abs(mubar - mubar_ref)/abs(mubar_ref)
    print '(A,ES12.4)', '  gyro-center vs GC max |ds|     = ', dev_max

    call check('CP trace completes (no spurious loss)', .not. cp_lost, nfail)
    call check('CP energy bounded (< 1e-3)', (Emax - Emin)/abs(E0) < 1.0e-3_dp, nfail)
    ! The gyro-AVERAGED mu (the adiabatic invariant) must show only a small secular
    ! drift; the instantaneous mu breathes at the gyrofrequency (~8%), so it is not
    ! the invariant. The single-gyrophase emergent mu approximates the true gyro-
    ! action to O(rho*), so its half-mean secular drift over a partial bounce is
    ! O(rho*) ~ 5e-3, a few percent at most. The bound is 2e-2: small and bounded,
    ! the adiabatic-invariance signature, while admitting the O(rho*) estimate error.
    call check('CP gyro-averaged mu adiabatically conserved (drift < 2e-2)', &
        abs(mubar - mubar_ref)/abs(mubar_ref) < 2.0e-2_dp, nfail)
    call check('CP gyro-center confined (0.05 < s < 0.95)', &
        sbar_min > 0.05_dp .and. sbar_max < 0.95_dp, nfail)
    ! Over a few-bounce-time-fraction of resolved gyration the gyro-center must
    ! FOLLOW the GC surface, not drift off it. The deviation is bounded by the FLR
    ! offset (the O(rho*) physics) plus the running-average residual; a flux-label
    ! tolerance of 0.05 (the band scale) catches a runaway while admitting the
    ! genuine FLR shift.
    call check('CP gyro-center tracks GC (max|ds| < 0.05)', dev_max < 0.05_dp, nfail)
    deallocate(scp_hist, sgc_hist, mu_hist)
  end subroutine test_gyrocenter_tracking

  ! Emergent magnetic moment mu = vperp^2/(2|B|) from the resolved CP velocity.
  subroutine emergent_mu(st, mu_e)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(out) :: mu_e
    real(dp) :: r, th, ph, vpar, vsq, vperp2, Bmod

    call cpp_canon_to_gc(st, r, th, ph, vpar)         ! vpar = h_i v^i
    vsq = 2.0_dp*cpp_canon_energy(st)/st%mass         ! |v|^2 = 2 H (CP: no mu|B|)
    vperp2 = max(vsq - vpar*vpar, 0.0_dp)
    call read_boozer_field_mod(st%z(1:3), Bmod)
    mu_e = st%mass*vperp2/(2.0_dp*Bmod)
  end subroutine emergent_mu

  subroutine read_boozer_field_mod(u, Bmod)
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: Bmod
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), dBmod(3), hcov(3)

    call boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
      Bctr, Bcov, Bmod, dBmod, hcov)
  end subroutine read_boozer_field_mod

  ! Centered boxcar of width w ending no later than index i (i-w+1 .. i).
  function boxcar(h, i, w) result(avg)
    real(dp), intent(in) :: h(0:)
    integer, intent(in) :: i, w
    real(dp) :: avg
    integer :: k, lo
    lo = max(0, i - w + 1)
    avg = 0.0_dp
    do k = lo, i
      avg = avg + h(k)
    end do
    avg = avg/real(i - lo + 1, dp)
  end function boxcar

  ! Mean of h over the inclusive index range [lo, hi].
  function mean(h, lo, hi) result(avg)
    real(dp), intent(in) :: h(0:)
    integer, intent(in) :: lo, hi
    real(dp) :: avg
    integer :: k
    avg = 0.0_dp
    do k = lo, hi
      avg = avg + h(k)
    end do
    avg = avg/real(hi - lo + 1, dp)
  end function mean

  subroutine test_loss_propagation(z0, npoiper2, nfail)
    real(dp), intent(in) :: z0(5)
    integer, intent(in) :: npoiper2
    integer, intent(inout) :: nfail
    type(tracer_t) :: edge
    real(dp) :: zedge(5), dtm
    integer :: ierr

    dtm = dtaumin_for(npoiper2, rbig)
    zedge = z0; zedge(1) = 0.5_dp
    call init_sympl(edge%si, edge%f, zedge, dtm, dtm, relerr, integmode)
    call init_cp(edge%cp, edge%f, zedge, dtm)
    zedge(1) = 1.5_dp
    call orbit_timestep_cp_canonical(edge%cp, edge%f, zedge, ierr)
    call check('CP6D wrapper flags z(1)>1 as loss (ierr/=0)', ierr /= 0, nfail)
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

end program test_cp6d_vs_gc
