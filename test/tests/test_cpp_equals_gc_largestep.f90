program test_cpp_equals_gc_largestep
  ! CPP (flux-canonical, mu fixed) must reproduce the GC symplectic trajectory
  ! at GC-sized steps. The CPP Gauss residual is the GC degenerate-Lagrangian
  ! Euler-Lagrange system specialized to fixed mu, so on the same canonical
  ! chart, same stage, same dt the two integrators advance the 4D state
  ! z=(r,theta,phi,pphi) identically up to the solver tolerance.
  !
  ! Cheap real field: BOOZER chart on the QA wout. Cross-check (strongest
  ! available): max over the orbit of |z_CPP - z_GC| < 1e-10.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: twopi
  use simple_main, only: init_field
  use simple, only: tracer_t, init_sympl, init_params
  use params, only: isw_field_type, field_input, coord_input, integmode
  use new_vmec_stuff_mod, only: rmajor
  use magfie_sub, only: BOOZER
  use field_can_mod, only: field_can_t, evaluate
  use orbit_symplectic_base, only: symplectic_integrator_t, GAUSS1, GAUSS2
  use orbit_symplectic, only: orbit_timestep_sympl
  use orbit_cpp, only: orbit_timestep_cpp, orbit_cpp_init, cpp_stages_from_mode

  implicit none

  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
  integer :: nfail
  type(tracer_t) :: norb

  nfail = 0

  isw_field_type = BOOZER
  field_input = 'wout.nc'
  coord_input = 'wout.nc'
  call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, GAUSS1)
  call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0e-12_dp)
  if (.not. associated(evaluate)) then
    print *, 'evaluate pointer not associated for BOOZER field'
    error stop 1
  end if

  call run_compare(norb, GAUSS1, 'GAUSS1', nfail)
  call run_compare(norb, GAUSS2, 'GAUSS2', nfail)

  if (nfail == 0) then
    print *, 'ALL CPP==GC LARGE-STEP TESTS PASSED'
  else
    print *, 'CPP==GC TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine run_compare(norb, mode, tag, nfail)
    type(tracer_t), intent(inout) :: norb
    integer, intent(in) :: mode
    character(*), intent(in) :: tag
    integer, intent(inout) :: nfail

    type(symplectic_integrator_t) :: si_gc, si_cpp
    type(field_can_t) :: f_gc, f_cpp
    real(dp) :: z0(5), dtau, dtaumin, relerr, rbig, maxdiff, d
    integer :: it, ierr_gc, ierr_cpp, s, nstep

    rbig = rmajor*1.0e2_dp
    ! GC-sized step: 256 points per torus (production-class resolution).
    dtau    = twopi*rbig/256.0_dp
    dtaumin = dtau
    relerr  = 1.0e-13_dp
    nstep   = 300
    s = cpp_stages_from_mode(mode)

    ! Trapped-particle-class IC.
    z0 = [0.4_dp, 0.7_dp, 0.1_dp, 1.0_dp, 0.2_dp]

    integmode = mode
    norb%relerr = relerr

    call init_sympl(si_gc, f_gc, z0, dtau, dtaumin, relerr, mode)
    call init_sympl(si_cpp, f_cpp, z0, dtau, dtaumin, relerr, mode)
    call orbit_cpp_init(si_cpp, f_cpp)

    maxdiff = 0.0_dp
    do it = 1, nstep
      call orbit_timestep_sympl(si_gc, f_gc, ierr_gc)
      call orbit_timestep_cpp(si_cpp, f_cpp, s, ierr_cpp)
      if (ierr_gc /= 0 .or. ierr_cpp /= 0) then
        print '(A,A,A,I0,A,I0)', '  ', tag, ': early exit ierr_gc=', ierr_gc, &
            ' ierr_cpp=', ierr_cpp
        exit
      end if
      d = maxval(abs(si_cpp%z - si_gc%z))
      maxdiff = max(maxdiff, d)
    end do

    print '(A,A,A,ES12.4)', '  ', tag, ': max |z_CPP - z_GC| = ', maxdiff
    call check(tag//': CPP matches GC to Newton tol', maxdiff < 1.0e-10_dp, nfail)
  end subroutine run_compare

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

end program test_cpp_equals_gc_largestep
