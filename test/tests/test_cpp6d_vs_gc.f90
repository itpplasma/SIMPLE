program test_cpp6d_vs_gc
  ! Genuine 6D canonical CPP (orbit_model=ORBIT_CPP6D) wired into the production
  ! alpha-loss pipeline through REAL VMEC flux coordinates (COORD_VMEC) on the
  ! real reactor-scale equilibrium test_data/wout.nc (a QA stellarator,
  ! rho* ~ 1/200), validated against the production guiding center.
  !
  ! WHY COORD_VMEC, not the Boozer chartmap: the Cartesian-storage Boozer chartmap
  ! was diagnosed inconsistent. libneo splines the chartmap Cartesian x/y/z with a
  ! PERIODIC fit over one field period, but for nfp>1 the Cartesian x,y are not
  ! field-period-periodic (they rotate by 2pi/nfp), so the periodic spline
  ! destroys the secular toroidal rotation: the analytic spline e_phi loses its ~R
  ! magnitude and the geometric metric gives h_i g^ij h_j ~ nfp^2 instead of 1
  ! (the covariant unit-field invariant |h|^2). The defect is upstream in libneo's
  ! Cartesian-storage path and in the storage convention itself; it cannot be
  ! repaired in the SIMPLE metric post-processor. The VMEC flux metric from libneo
  ! is consistent (test_cpp_vmec: |g g^-1 - I| < 1e-10), so the production loss
  ! path runs there. See DOC/coordinates-and-fields.md, "6D canonical CPP".
  !
  ! Acceptance gates:
  !   (a) METRIC CONSISTENCY -- the exact check the chartmap failed: on the
  !       production COORD_VMEC chart h_i g^ij h_j = |h|^2 = 1 to central-
  !       difference (Christoffel-from-FD) accuracy. The chartmap gave 228..472 at
  !       the same kind of point; COORD_VMEC gives ~1.
  !   (b) The 6D canonical-midpoint scheme conserves energy and holds mu fixed over
  !       a short resolved trace (the symplectic / fixed-mu signature).
  !   (c) The 6D->GC reduction stays on a bounded flux band overlapping the GC band
  !       (the chart-independent s label), i.e. the 6D orbit follows the GC surface.
  !   (d) The s>=1 loss propagates through the production wrapper to ierr/=0.
  !
  ! The wire keeps the SIMPLE GC normalization (mass=1, qc=sqrt(2)/ro0,
  ! dt=dtaumin/sqrt(2)). The consistent VMEC metric (|h|^2=1) makes the 6D
  ! Hamiltonian reduce to the GC one exactly, so the bare production GC macrostep
  ! runs without sub-cycling. The FD-Jacobian host path uses an FD-matched Newton
  ! step tolerance (a central-difference Jacobian cannot reach the analytic-path
  ! 1e-12 floor); see orbit_cpp_canonical.cpp_canon_step.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use parmot_mod, only: ro0
  use simple, only: init_sympl, init_cpp, init_params, tracer_t, &
    orbit_timestep_cpp_canonical
  use simple_main, only: init_field
  use orbit_symplectic, only: orbit_timestep_sympl
  use orbit_cpp_canonical, only: cpp_canon_energy
  use orbit_cpp_vmec_metric, only: vmec_eval_metric, vmec_eval_field, &
    vmec_metric_ready
  use params, only: field_input, coord_input, integmode, relerr, dtaumin
  use velo_mod, only: isw_field_type
  use magfie_sub, only: BOOZER

  implicit none

  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
  type(tracer_t) :: norb
  real(dp) :: z0(5)
  integer :: nfail

  nfail = 0

  ! Production field setup: BOOZER canonical chart on the real VMEC equilibrium.
  ! init_field splines wout.nc; init_params sets v0/ro0 (3.5 MeV alpha) and the
  ! production normalized step dtaumin.
  isw_field_type = BOOZER
  field_input = 'wout.nc'
  coord_input = 'wout.nc'
  integmode = 1
  relerr = 1.0d-13
  call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, integmode)
  call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0d-13)
  dtaumin = norb%dtaumin
  print '(A,ES12.4)', '  ro0 (cm)    = ', ro0
  print '(A,ES12.4)', '  dtaumin     = ', dtaumin

  ! Shared trapped-class IC in flux coords (s, theta, phi, v/v0, lambda).
  z0 = [0.3_dp, 0.5_dp, 0.2_dp, 1.0_dp, 0.3_dp]

  call test_metric_consistency(z0, nfail)
  call test_trace_and_tracking(norb, z0, nfail)
  call test_loss_propagation(nfail)

  if (nfail == 0) then
    print *, 'ALL CPP6D-VS-GC TESTS PASSED'
  else
    print *, 'CPP6D-VS-GC TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine test_metric_consistency(z0, nfail)
    ! The defect the chartmap had: h_i g^ij h_j must be 1 (h is the covariant unit
    ! field; g^ij raises it to h^i, so h_i g^ij h_j = |h|^2 = 1). On the production
    ! COORD_VMEC chart it holds to central-difference (Christoffel) accuracy; the
    ! chartmap gave O(nfp^2) = hundreds.
    real(dp), intent(in) :: z0(5)
    integer, intent(inout) :: nfail
    real(dp) :: u(3), g(3,3), ginv(3,3), dg(3,3,3)
    real(dp) :: Acov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: hgh, hcon(3)
    integer :: i, j

    if (.not. vmec_metric_ready()) call init_cpp(norb%cpp, norb%f, z0, dtaumin)
    u = [z0(1), z0(2), z0(3)]
    call vmec_eval_metric(u, g, ginv, dg)
    call vmec_eval_field(u, Acov, Bmod, dBmod, hcov)

    do i = 1, 3
      hcon(i) = 0.0_dp
      do j = 1, 3
        hcon(i) = hcon(i) + ginv(i,j)*hcov(j)
      end do
    end do
    hgh = 0.0_dp
    do i = 1, 3
      hgh = hgh + hcov(i)*hcon(i)
    end do
    print '(A,F12.8)', '  h_i g^ij h_j (must be ~1) = ', hgh
    print '(A,ES12.4)', '  |B| (Gauss)               = ', Bmod
    ! Central-difference Christoffel -> FD-level accuracy (~1e-2), per the
    ! diagnosis (0.998, 1.008, 0.946 at s=0.3,0.5,0.7). NOT the chartmap's 228+.
    call check('COORD_VMEC metric consistent (|h_i g^ij h_j - 1| < 3e-2)', &
        abs(hgh - 1.0_dp) < 3.0e-2_dp, nfail)
    call check('COORD_VMEC NOT the broken chartmap (h_i g^ij h_j < 2)', &
        hgh < 2.0_dp, nfail)
  end subroutine test_metric_consistency

  subroutine test_trace_and_tracking(norb, z0, nfail)
    ! Drive the production GC (orbit_timestep_sympl on the BOOZER chart) and the
    ! genuine 6D CPP through the PRODUCTION wrapper (COORD_VMEC) at the BARE GC
    ! macrostep -- no sub-cycling -- from the SAME (s,theta,phi,v,lambda) start.
    ! Both must stay confined and conserve their invariants; the 6D s band must
    ! overlap the GC band. s is the chart-independent flux label, so the comparison
    ! is fair across the Boozer (GC) and VMEC (6D) angle conventions.
    type(tracer_t), intent(inout) :: norb
    real(dp), intent(in) :: z0(5)
    integer, intent(inout) :: nfail
    type(tracer_t) :: gc, cpp
    real(dp) :: zgc(5), zcpp(5)
    real(dp) :: sgc_min, sgc_max, scpp_min, scpp_max
    real(dp) :: E0, E, Emin, Emax, mu0, mu_now
    integer :: it, ierr, nstep
    logical :: gc_lost, cpp_lost

    nstep = 2000

    ! --- production GC ---
    zgc = z0
    call init_sympl(gc%si, gc%f, zgc, dtaumin, dtaumin, relerr, integmode)
    sgc_min = zgc(1); sgc_max = zgc(1); gc_lost = .false.
    do it = 1, nstep
      call orbit_timestep_sympl(gc%si, gc%f, ierr)
      if (ierr /= 0) then; gc_lost = .true.; exit; end if
      sgc_min = min(sgc_min, gc%si%z(1)); sgc_max = max(sgc_max, gc%si%z(1))
    end do

    ! --- genuine 6D CPP through the production wrapper at the bare GC step ---
    zcpp = z0
    call init_sympl(cpp%si, cpp%f, zcpp, dtaumin, dtaumin, relerr, integmode)
    call init_cpp(cpp%cpp, cpp%f, zcpp, dtaumin)
    E0 = cpp_canon_energy(cpp%cpp); Emin = E0; Emax = E0; mu0 = cpp%cpp%mu
    scpp_min = zcpp(1); scpp_max = zcpp(1); cpp_lost = .false.
    do it = 1, nstep
      call orbit_timestep_cpp_canonical(cpp%cpp, cpp%f, zcpp, ierr)
      if (ierr /= 0) then; cpp_lost = .true.; exit; end if
      scpp_min = min(scpp_min, zcpp(1)); scpp_max = max(scpp_max, zcpp(1))
      E = cpp_canon_energy(cpp%cpp); Emin = min(Emin, E); Emax = max(Emax, E)
    end do
    mu_now = cpp%cpp%mu

    print '(A,F8.5,A,F8.5,A)', '  GC    s band [', sgc_min, ',', sgc_max, ']'
    print '(A,F8.5,A,F8.5,A)', '  CPP6D s band [', scpp_min, ',', scpp_max, ']'
    print '(A,ES12.4)', '  CPP6D max|dE/E0|      = ', (Emax - Emin)/abs(E0)
    print '(A,ES12.4)', '  mu drift |mu-mu0|/mu0 = ', abs(mu_now - mu0)/abs(mu0)

    call check('GC confined over trace', .not. gc_lost, nfail)
    call check('CPP6D trace completes at bare step (no spurious loss)', &
        .not. cpp_lost, nfail)
    call check('CPP6D energy conserved (< 1e-3)', &
        (Emax - Emin)/abs(E0) < 1.0e-3_dp, nfail)
    call check('CPP6D mu held exactly fixed (< 1e-12)', &
        abs(mu_now - mu0)/abs(mu0) < 1.0e-12_dp, nfail)
    call check('CPP6D z(4) = pabs preserved (1.0)', &
        abs(zcpp(4) - 1.0_dp) < 1.0e-12_dp, nfail)
    ! Both orbits stay on the same flux band: the 6D reduction follows the GC
    ! surface. The bands need not coincide bit-for-bit (different angles), but
    ! they must overlap and neither may eject.
    call check('GC stays confined (0.05 < s < 0.95)', &
        sgc_min > 0.05_dp .and. sgc_max < 0.95_dp, nfail)
    call check('CPP6D stays confined (0.05 < s < 0.95)', &
        scpp_min > 0.05_dp .and. scpp_max < 0.95_dp, nfail)
    call check('CPP6D radial band tracks GC band (overlap, edges within 0.1)', &
        abs(scpp_min - sgc_min) < 0.1_dp .and. abs(scpp_max - sgc_max) < 0.1_dp, nfail)
  end subroutine test_trace_and_tracking

  subroutine test_loss_propagation(nfail)
    ! A z(1) > 1 must short-circuit the production wrapper to ierr/=0 (loss), as
    ! times_lost/confined_fraction expect.
    integer, intent(inout) :: nfail
    type(tracer_t) :: edge
    real(dp) :: zedge(5)
    integer :: ierr

    zedge = [0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.3_dp]
    call init_sympl(edge%si, edge%f, zedge, dtaumin, dtaumin, relerr, integmode)
    call init_cpp(edge%cpp, edge%f, zedge, dtaumin)
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
