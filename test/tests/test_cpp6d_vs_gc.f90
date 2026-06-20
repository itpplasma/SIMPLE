program test_cpp6d_vs_gc
  ! Genuine 6D canonical CPP (orbit_model=ORBIT_CPP6D) wired into the production
  ! alpha-loss pipeline through native Boozer coordinates on the real reactor-scale
  ! equilibrium test_data/wout.nc (a QA stellarator, rho* ~ 1/200), validated
  ! against the production guiding center.
  !
  ! Acceptance gates:
  !   (a) METRIC CONSISTENCY -- on the production Boozer chart h_i g^ij h_j =
  !       |h|^2 = 1.
  !   (b) The 6D canonical-midpoint scheme conserves energy and holds mu fixed over
  !       a short resolved trace (the symplectic / fixed-mu signature).
  !   (c) The 6D->GC reduction stays on a bounded flux band overlapping the GC band
  !       (the chart-independent s label), i.e. the 6D orbit follows the GC surface.
  !   (d) The s>=1 loss propagates through the production wrapper to ierr/=0.
  !
  ! The wire keeps the SIMPLE GC normalization (mass=1, qc=sqrt(2)/ro0,
  ! dt=dtaumin/sqrt(2)). The Boozer field-spline h_i g^ij h_j floor is tested at
  ! the same 2e-4 level as test_boozer_field_metric. The native Boozer CPP pusher
  ! is run at npoiper2=1024 here: still a production microstep, but small enough
  ! for the first-derivative Newton solve to complete this long trapped trace.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use parmot_mod, only: ro0
  use simple, only: init_sympl, init_cpp, init_params, tracer_t, &
    orbit_timestep_cpp_canonical
  use simple_main, only: init_field
  use orbit_symplectic, only: orbit_timestep_sympl
  use orbit_cpp_canonical, only: cpp_canon_energy
  use boozer_field_metric, only: boozer_field_metric_eval
  use params, only: field_input, coord_input, integmode, relerr, dtaumin, orbit_coord
  use velo_mod, only: isw_field_type
  use magfie_sub, only: BOOZER
  use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
  use boozer_sub, only: get_boozer_coordinates

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
  orbit_coord = 1
  integmode = 1
  relerr = 1.0d-13
  call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, integmode)
  use_B_r = .true.
  use_del_tp_B = .true.
  call get_boozer_coordinates
  call init_params(norb, 2, 4, 3.5e6_dp, 1024, 1, 1.0d-13)
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
    ! h_i g^ij h_j must be 1: h is the covariant unit field and g^ij raises it.
    real(dp), intent(in) :: z0(5)
    integer, intent(inout) :: nfail
    real(dp) :: u(3), g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: hgh, hcon(3)
    integer :: i, j

    u = [z0(1), z0(2), z0(3)]
    call boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
         Bctr, Bcov, Bmod, dBmod, hcov)

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
    print '(A,F18.15)', '  h_i g^ij h_j (must be ~1) = ', hgh
    print '(A,ES12.4)', '  |B| (Gauss)               = ', Bmod
    call check('Boozer metric consistent (|h_i g^ij h_j - 1| < 2e-4)', &
        abs(hgh - 1.0_dp) < 2.0e-4_dp, nfail)
  end subroutine test_metric_consistency

  subroutine test_trace_and_tracking(norb, z0, nfail)
    ! Drive the production GC (orbit_timestep_sympl on the BOOZER chart) and the
    ! genuine 6D CPP through the production wrapper at the bare GC
    ! macrostep -- no sub-cycling -- from the SAME (s,theta,phi,v,lambda) start.
    ! Both must stay confined and conserve their invariants; the 6D s band must
    ! overlap the GC band. Both paths use Boozer angles here.
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
      if (ierr /= 0) then
        print '(A,I0,A,I0,A,ES12.4)', '  CPP6D stopped at step ', it, &
          ' ierr=', ierr, ' s=', zcpp(1)
        cpp_lost = .true.
        exit
      end if
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
    ! they must overlap and neither may eject. The 6D banana is a FULL orbit with a
    ! finite Larmor radius, so its turning point can sit further out than the
    ! zero-width GC banana tip; the bound is "not lost to the edge" (s < 1), and the
    ! "tracks GC band (edges within 0.1)" check below enforces that the excess
    ! stays within the FLR tolerance. (The single-source metric tracks the GC
    ! banana tip less tightly than the old dual-source metric, whose dg was NOT the
    ! derivative of its g -- that inconsistency made an analytic Jacobian diverge;
    ! a self-consistent dg is required for the smooth no-FD-ejection Jacobian.)
    call check('GC stays confined (0.05 < s < 0.95)', &
        sgc_min > 0.05_dp .and. sgc_max < 0.95_dp, nfail)
    call check('CPP6D stays confined (not lost: 0.05 < s < 1.0)', &
        scpp_min > 0.05_dp .and. scpp_max < 1.0_dp, nfail)
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
