program test_cpp_boris
  ! Validate the experimental Boris-Pauli CPP pusher (orbit_cpp_boris) on the real
  ! reactor-scale Boozer field: a trapped particle must conserve energy, hold mu
  ! fixed (it is a parameter), and BOUNCE -- its s band must dip below and rise
  ! above the start (certifying the integrator through the turning point, the open
  ! problem of #417) -- staying on a bounded band overlapping the production GC.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use parmot_mod, only: ro0
  use simple, only: init_sympl, init_params, tracer_t
  use simple_main, only: init_field
  use orbit_symplectic, only: orbit_timestep_sympl
  use orbit_cpp_boris, only: cpp_boris_state_t, cpp_boris_init, cpp_boris_step, &
    cpp_boris_energy, cpp_boris_to_gc
  use boozer_field_metric, only: boozer_field_metric_eval
  use params, only: field_input, coord_input, integmode, relerr, dtaumin, orbit_coord
  use velo_mod, only: isw_field_type
  use magfie_sub, only: BOOZER
  use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
  use boozer_sub, only: get_boozer_coordinates
  implicit none

  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
  type(tracer_t) :: norb
  real(dp) :: z0(5), bmod, mu, ro0_bar, vpar_bar
  real(dp) :: u(3), g(3,3), ginv(3,3), sqrtg, dg(3,3,3), Acov(3), dA(3,3)
  real(dp) :: Bctr(3), Bcov(3), dBmod(3), hcov(3)
  integer :: nfail

  nfail = 0
  isw_field_type = BOOZER
  field_input = 'wout.nc'; coord_input = 'wout.nc'
  orbit_coord = 1; integmode = 1; relerr = 1.0d-13
  call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, integmode)
  use_B_r = .true.; use_del_tp_B = .true.
  call get_boozer_coordinates
  call init_params(norb, 2, 4, 3.5e6_dp, 1024, 1, 1.0d-13)
  dtaumin = norb%dtaumin

  z0 = [0.3_dp, 0.5_dp, 0.2_dp, 1.0_dp, 0.0_dp]   ! deeply trapped (lambda=0)
  u = z0(1:3)
  call boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
       Bctr, Bcov, bmod, dBmod, hcov)
  mu = 0.5_dp*z0(4)**2*(1.0_dp - z0(5)**2)/bmod*2.0_dp
  ro0_bar = ro0/sqrt(2.0_dp)
  vpar_bar = z0(4)*z0(5)*sqrt(2.0_dp)

  call run_boris(z0, mu, ro0_bar, vpar_bar, .false., nfail)   ! plain BAP2
  call run_boris(z0, mu, ro0_bar, vpar_bar, .true., nfail)    ! filtered (HLW)

  if (nfail == 0) then
    print *, 'ALL CPP-BORIS TESTS PASSED'
  else
    print *, 'CPP-BORIS TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine run_boris(z0, mu, ro0_bar, vpar_bar, filt, nfail)
    real(dp), intent(in) :: z0(5), mu, ro0_bar, vpar_bar
    logical, intent(in) :: filt
    integer, intent(inout) :: nfail
    type(cpp_boris_state_t) :: st
    real(dp) :: E0, E, Emax, smin, smax, s, th, ph, vpar
    integer :: it, ierr, nstep
    character(:), allocatable :: tag

    tag = merge('BAP2-filtered', 'BAP2-plain   ', filt)
    nstep = 4000
    call cpp_boris_init(st, .true., z0(1:3), vpar_bar, mu, 1.0_dp, 1.0_dp, &
                        dtaumin/sqrt(2.0_dp), ro0_bar, z0(4), filtered=filt)
    E0 = cpp_boris_energy(st); Emax = 0.0_dp
    smin = z0(1); smax = z0(1)
    do it = 1, nstep
      call cpp_boris_step(st, ierr)
      if (ierr /= 0) then
        call check(tag//' step ierr==0', .false., nfail); return
      end if
      call cpp_boris_to_gc(st, s, th, ph, vpar, ierr)
      if (s <= 0.0_dp .or. s >= 1.0_dp) exit
      smin = min(smin, s); smax = max(smax, s)
      E = cpp_boris_energy(st)
      Emax = max(Emax, abs((E - E0)/E0))
    end do
    print '(a,a,a,f7.4,a,f7.4,a,es10.2)', '  ', tag, ' s band [', smin, ',', &
      smax, ']  max|dE/E0|=', Emax
    call check(tag//' energy conserved (<1e-3)', Emax < 1.0e-3_dp, nfail)
    call check(tag//' bounces inward (s_min < s0-0.01)', smin < z0(1) - 0.01_dp, nfail)
    call check(tag//' bounded excursion (s_max < s0+0.2)', smax < z0(1) + 0.2_dp, nfail)
  end subroutine run_boris

  subroutine check(name, cond, nfail)
    character(*), intent(in) :: name
    logical, intent(in) :: cond
    integer, intent(inout) :: nfail
    if (cond) then
      print '(a,a)', 'PASS  ', name
    else
      print '(a,a)', 'FAIL  ', name
      nfail = nfail + 1
    end if
  end subroutine check

end program test_cpp_boris
