program test_cpp_jacobian_fd
  ! Finite-difference self-check of the full-metric CPP Newton Jacobian on the
  ! COORD_BOOZER path. The analytic jacobian() must match a central finite
  ! difference of the public residual() to high relative accuracy -- this proves
  ! the analytic second-derivative terms (d2g, d2A, d2Bmod) in
  ! jacobian_vmec_analytic are correct.
  !
  ! Setup mirrors test_cpp6d_vs_gc: production BOOZER chart on the real
  ! reactor-scale equilibrium test_data/wout.nc, a trapped-class IC, GC-sized dt.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use parmot_mod, only: ro0
  use simple, only: init_sympl, init_cpp, init_params, tracer_t
  use simple_main, only: init_field
  use orbit_cpp_canonical, only: cpp_canon_state_t, residual, jacobian
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

  ! Trapped-class IC in flux coords (s, theta, phi, v/v0, lambda); s=0.5 mid-radius.
  z0 = [0.5_dp, 0.5_dp, 0.2_dp, 1.0_dp, 0.3_dp]

  call test_jacobian_fd(norb, z0, nfail)

  if (nfail == 0) then
    print *, 'ALL CPP-JACOBIAN-FD TESTS PASSED'
  else
    print *, 'CPP-JACOBIAN-FD TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine test_jacobian_fd(norb, z0, nfail)
    type(tracer_t), intent(inout) :: norb
    real(dp), intent(in) :: z0(5)
    integer, intent(inout) :: nfail
    type(tracer_t) :: cpp
    type(cpp_canon_state_t) :: st
    real(dp) :: zcpp(5), zold(6), z(6)
    real(dp) :: jac(6,6), jfd(6,6)
    real(dp) :: fplus(6), fminus(6), zp(6), zm(6)
    real(dp) :: h, denom, reldiff, maxrel, scale
    integer :: m, i, worst_i, worst_j

    zcpp = z0
    call init_sympl(cpp%si, cpp%f, zcpp, dtaumin, dtaumin, relerr, integmode)
    call init_cpp(cpp%cpp, cpp%f, zcpp, dtaumin)
    st = cpp%cpp

    ! Evaluate the Jacobian at a partially-advanced iterate (z != zold) so the
    ! velocity vmid is nonzero and every d2 term is exercised. zold is the start
    ! state; z is one small Newton-scale displacement away.
    zold = st%z
    z = zold
    z(1) = z(1) + 1.0e-4_dp
    z(2) = z(2) + 2.0e-4_dp
    z(3) = z(3) + 3.0e-4_dp
    z(4) = z(4) + 1.0e-3_dp
    z(5) = z(5) - 2.0e-3_dp
    z(6) = z(6) + 1.5e-3_dp

    call jacobian(st, zold, z, jac)

    ! Central FD of residual w.r.t. each z component.
    h = 1.0e-6_dp
    do m = 1, 6
      zp = z; zm = z
      zp(m) = zp(m) + h
      zm(m) = zm(m) - h
      call residual(st, zold, zp, fplus)
      call residual(st, zold, zm, fminus)
      do i = 1, 6
        jfd(i,m) = (fplus(i) - fminus(i))/(2.0_dp*h)
      end do
    end do

    ! Relative agreement, scaled by the column magnitude so small entries are not
    ! penalized against round-off.
    maxrel = 0.0_dp
    worst_i = 0; worst_j = 0
    do m = 1, 6
      scale = 0.0_dp
      do i = 1, 6
        scale = max(scale, abs(jfd(i,m)), abs(jac(i,m)))
      end do
      if (scale <= 0.0_dp) cycle
      do i = 1, 6
        denom = max(abs(jfd(i,m)), 1.0e-3_dp*scale)
        reldiff = abs(jac(i,m) - jfd(i,m))/denom
        if (reldiff > maxrel) then
          maxrel = reldiff; worst_i = i; worst_j = m
        end if
      end do
    end do

    print '(A,ES12.4)', '  max relative |jac - jfd| = ', maxrel
    print '(A,I0,A,I0,A,ES12.4,A,ES12.4)', '  worst component (', worst_i, ',', &
        worst_j, '): jac=', jac(worst_i,worst_j), ' jfd=', jfd(worst_i,worst_j)

    call check('analytic Jacobian matches central FD (< 1e-5 relative)', &
        maxrel < 1.0e-5_dp, nfail)
  end subroutine test_jacobian_fd

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

end program test_cpp_jacobian_fd
