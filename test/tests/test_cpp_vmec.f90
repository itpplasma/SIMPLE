program test_cpp_vmec
  ! 6D canonical-midpoint integrator on REAL VMEC flux coordinates
  ! (test/test_data/wout.nc, nfp=2 stellarator). The same generalized full-metric
  ! residual that reproduces the analytic-tokamak oracle (test_cpp_canonical) here
  ! runs on the libneo metric/Christoffel (#322) + SIMPLE native VMEC field.
  !
  ! Assertions:
  !   CP  (gyro-resolved, small dt): symplectic energy bound -- |dE/E0| stays
  !       below a tight band with no secular drift over a long run.
  !   CPP (Pauli, BIG dt): the guiding-center reduction stays on a bounded radial
  !       band (the banana/passing band, not lost to the s=1 edge or the axis),
  !       i.e. the big-step CPP reproduces the GC confinement signature.
  !
  ! Honest limitation: this is a 2-field-period stellarator, so the toroidal
  ! canonical momentum is NOT a conserved quantity (no axisymmetry); we assert the
  ! Hamiltonian energy and the radial band, not p_phi. Near the magnetic axis
  ! (s -> 0) the flux-coordinate metric is singular and the central-difference
  ! field gradients lose accuracy; the test starts at mid-radius s ~ 0.3.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use orbit_cpp_vmec_metric, only: vmec_metric_init, vmec_metric_ready
  use orbit_cpp_canonical, only: cpp_canon_state_t, cpp_canon_init, cpp_canon_step, &
       cpp_canon_energy, cpp_canon_to_gc, MODEL_CP, MODEL_CPP_SYM, COORD_VMEC
  implicit none

  integer :: nfail
  real(dp), parameter :: x0(3) = [0.3_dp, 0.6_dp, 0.2_dp]  ! (s, vartheta, varphi)
  real(dp), parameter :: mass = 1.0_dp, charge = 1.0_dp

  nfail = 0

  call vmec_metric_init('wout.nc')
  if (.not. vmec_metric_ready()) then
    print *, 'FAIL  VMEC metric not initialized'
    error stop 1
  end if
  print *, 'VMEC metric/field initialized from wout.nc'

  call test_metric_sane(nfail)
  call test_cp_energy(nfail)
  call test_cpp_banana(nfail)

  if (nfail == 0) then
    print *, 'ALL VMEC 6D TESTS PASSED'
  else
    print *, 'VMEC 6D TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine test_metric_sane(nfail)
    ! The libneo metric must be symmetric positive (g_ii > 0), g g^-1 = I, and the
    ! field modulus must be the VMEC-scale |B| (CGS Gauss, ~5e4 G for B0~5 T).
    use orbit_cpp_vmec_metric, only: vmec_eval_metric, vmec_eval_field
    integer, intent(inout) :: nfail
    real(dp) :: g(3,3), ginv(3,3), dg(3,3,3), prod(3,3)
    real(dp) :: Acov(3), Bmod, dBmod(3), hcov(3)
    integer :: i, j, k
    real(dp) :: offdiag

    call vmec_eval_metric(x0, g, ginv, dg)
    call vmec_eval_field(x0, Acov, Bmod, dBmod, hcov)

    prod = matmul(g, ginv)
    offdiag = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        if (i == j) then
          offdiag = max(offdiag, abs(prod(i,j) - 1.0_dp))
        else
          offdiag = max(offdiag, abs(prod(i,j)))
        end if
      end do
    end do
    print '(A,ES12.4)', '  |g g^-1 - I| = ', offdiag
    call check('metric g g^-1 = I', offdiag < 1.0e-10_dp, nfail)
    call check('metric g_11 > 0', g(1,1) > 0.0_dp, nfail)
    call check('metric g_22 > 0', g(2,2) > 0.0_dp, nfail)
    call check('metric g_33 > 0', g(3,3) > 0.0_dp, nfail)
    ! Metric is symmetric by construction (g_ij = e_i . e_j); compare relative to
    ! the metric scale, which is O(R^2) in cm^2 (~1e4) for a VMEC equilibrium.
    offdiag = (abs(g(1,2)-g(2,1)) + abs(g(1,3)-g(3,1)) + abs(g(2,3)-g(3,2))) &
              / max(abs(g(1,1)) + abs(g(2,2)) + abs(g(3,3)), 1.0_dp)
    print '(A,ES12.4)', '  metric relative asymmetry = ', offdiag
    call check('metric symmetric', offdiag < 1.0e-12_dp, nfail)
    print '(A,ES12.4)', '  |B| at start (Gauss) = ', Bmod
    call check('|B| at VMEC scale (1e4..1e5 G)', Bmod > 1.0e4_dp .and. Bmod < 1.0e5_dp, nfail)
  end subroutine test_metric_sane

  subroutine test_cp_energy(nfail)
    ! Full charged particle, gyro-resolved. dt small enough to resolve the gyro
    ! orbit (Larmor scale); symplectic energy stays bounded with no secular drift.
    integer, intent(inout) :: nfail
    type(cpp_canon_state_t) :: st
    real(dp) :: vperp0, mu, dt, E0, E, Emin, Emax, Eend, drift
    integer :: it, ierr, nsteps

    ! Pick mu and dt for a resolved gyro-orbit at the VMEC (CGS) scale. vperp0 is a
    ! thermal-ish speed; mu = m vperp^2 / (2|B|) follows in cpp_canon_init.
    vperp0 = 3.0e5_dp        ! cm/s scale velocity
    mu = 0.0_dp              ! CP derives mu from vperp0
    dt = 2.0e-8_dp           ! s; resolves the gyro period at |B|~5e4 G
    nsteps = 2000

    call cpp_canon_init(st, MODEL_CP, COORD_VMEC, x0, 0.0_dp, vperp0, mu, &
                        mass, charge, dt)
    E0 = cpp_canon_energy(st); Emin = E0; Emax = E0; Eend = E0
    do it = 1, nsteps
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) then
        print '(A,I0,A,I0)', '  CP step ', it, ' ierr=', ierr
        call check('CP VMEC run completes', .false., nfail); return
      end if
      E = cpp_canon_energy(st)
      Emin = min(Emin, E); Emax = max(Emax, E); Eend = E
    end do
    drift = (Eend - E0)/abs(E0)
    print '(A,ES12.4,A,ES12.4)', '  CP VMEC max|dE/E0| = ', (Emax-Emin)/abs(E0), &
         '  end-drift = ', drift
    call check('CP VMEC energy bounded (<1e-2)', (Emax-Emin)/abs(E0) < 1.0e-2_dp, nfail)
    call check('CP VMEC no secular drift (<5e-3)', abs(drift) < 5.0e-3_dp, nfail)
  end subroutine test_cp_energy

  subroutine test_cpp_banana(nfail)
    ! Pauli CPP at a BIG (guiding-center-sized) dt. The orbit must stay on a
    ! bounded radial band -- the GC banana/passing confinement signature -- and
    ! keep the GC parallel reduction finite, not run to the s=1 edge or the axis.
    integer, intent(inout) :: nfail
    type(cpp_canon_state_t) :: st
    real(dp) :: mu, dt, smin, smax, vpar0
    real(dp) :: sg, thg, phg, vparg, E0, E, Emax, Emin
    real(dp) :: sprev, sdir, sdir_prev
    integer :: it, ierr, nsteps, nturns

    mu = 1.0e-3_dp           ! magnetic moment (CGS) for the Pauli term
    vpar0 = 1.0e5_dp         ! parallel start speed (cm/s)
    dt = 5.0e-7_dp           ! big GC-scale step
    nsteps = 1000

    call cpp_canon_init(st, MODEL_CPP_SYM, COORD_VMEC, x0, vpar0, 0.0_dp, mu, &
                        mass, charge, dt)
    smin = st%z(1); smax = st%z(1)
    E0 = cpp_canon_energy(st); Emin = E0; Emax = E0
    sprev = st%z(1); sdir_prev = 0.0_dp; nturns = 0
    do it = 1, nsteps
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) then
        print '(A,I0,A,I0)', '  CPP step ', it, ' ierr=', ierr
        call check('CPP VMEC banana run completes', .false., nfail); return
      end if
      smin = min(smin, st%z(1)); smax = max(smax, st%z(1))
      E = cpp_canon_energy(st); Emin = min(Emin, E); Emax = max(Emax, E)
      ! Count radial turning points: sign flips of ds between steps. A banana/
      ! drift orbit reverses radially; a lost orbit drifts monotonically out.
      sdir = sign(1.0_dp, st%z(1) - sprev)
      if (sdir_prev /= 0.0_dp .and. sdir /= sdir_prev) nturns = nturns + 1
      sdir_prev = sdir; sprev = st%z(1)
    end do
    call cpp_canon_to_gc(st, sg, thg, phg, vparg)
    print '(A,ES12.4,A,ES12.4)', '  CPP banana s band = ', smax - smin, &
         '  max|dE/E0| = ', (Emax-Emin)/abs(E0)
    print '(A,F8.4,A,F8.4,A,I0)', '  CPP banana s in [', smin, ',', smax, &
         ']  radial turning points = ', nturns
    call check('CPP banana confined (0.05<s<0.95)', smin > 0.05_dp .and. smax < 0.95_dp, nfail)
    call check('CPP banana s oscillates (band > 1e-4)', smax - smin > 1.0e-4_dp, nfail)
    call check('CPP banana bounces (radial turning points > 2)', nturns > 2, nfail)
    call check('CPP banana energy bounded (<5e-2)', (Emax-Emin)/abs(E0) < 5.0e-2_dp, nfail)
    call check('CPP GC vpar finite', abs(vparg) < 1.0e8_dp .and. vparg == vparg, nfail)
  end subroutine test_cpp_banana

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

end program test_cpp_vmec
