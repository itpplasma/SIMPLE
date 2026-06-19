program test_cpp_canonical_device
  ! Verify the COORD_TOK 6D canonical step runs on the OpenACC device and matches
  ! the host result to round-off. cpp_canon_step_tok and its whole kernel chain
  ! (residual_tok -> eval_block_tok/eval_field_correct_test/dLdq/raise/residual_blk,
  ! jacobian_analytic -> grad_jacobian_tok, rk_solve) are !$acc routine seq with
  ! fixed-size state and integer dispatch. One particle per gang/vector lane, all
  ! three models (CP, CPP_SYM, CPP_VAR).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use orbit_cpp_canonical, only: cpp_canon_state_t, cpp_canon_init, &
       cpp_canon_step, cpp_canon_step_tok, &
       MODEL_CP, MODEL_CPP_SYM, MODEL_CPP_VAR, COORD_TOK
  implicit none

  ! nstep is short on purpose. The device kernel reproduces the host step bit for
  ! bit, but the dt=800 variational orbit (MODEL_CPP_VAR) is Lyapunov-unstable, so
  ! the last-bit difference between host (x86 FMA) and device (GPU FMA) ordering
  ! amplifies after ~5 steps. Five steps validate device==host step-for-step
  ! agreement across all three models without entering that chaotic regime.
  integer, parameter :: npt = 256, nstep = 5
  real(dp), parameter :: mu = 1.0e-5_dp, mass = 1.0_dp, charge = 1.0_dp
  integer :: nfail
  nfail = 0

  call run_model(MODEL_CP, 1.0_dp, 'CP', nfail)
  call run_model(MODEL_CPP_SYM, 80.0_dp, 'CPP_SYM', nfail)
  call run_model(MODEL_CPP_VAR, 800.0_dp, 'CPP_VAR', nfail)

  if (nfail == 0) then
    print *, 'ALL CPP CANONICAL DEVICE TESTS PASSED'
  else
    print *, 'CPP CANONICAL DEVICE TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine run_model(model, dt, name, nfail)
    integer, intent(in) :: model
    real(dp), intent(in) :: dt
    character(*), intent(in) :: name
    integer, intent(inout) :: nfail
    type(cpp_canon_state_t) :: sh(npt), sd(npt)
    real(dp) :: x0(3), vperp0, maxdiff, dd
    integer :: ip, it, ierr, k

    ! Spread of initial radii/angles so the lanes are not all identical.
    do ip = 1, npt
      x0 = [0.08_dp + 0.20_dp*real(modulo(ip*7, 1000), dp)/1000.0_dp, &
            1.5_dp + 0.5_dp*real(modulo(ip*13, 1000), dp)/1000.0_dp, &
            0.0_dp]
      vperp0 = merge(1.0e-3_dp, 0.0_dp, model == MODEL_CP)
      call cpp_canon_init(sh(ip), model, COORD_TOK, x0, 0.0_dp, vperp0, mu, &
                          mass, charge, dt)
      sd(ip) = sh(ip)
    end do

    ! Host reference.
    do it = 1, nstep
      do ip = 1, npt
        call cpp_canon_step(sh(ip), ierr)
      end do
    end do

    ! Device: one particle per lane, the same number of steps.
    !$acc parallel loop gang vector copy(sd) private(ierr)
    do ip = 1, npt
      do it = 1, nstep
        call cpp_canon_step_tok(sd(ip), ierr)
      end do
    end do

    maxdiff = 0.0_dp
    do ip = 1, npt
      do k = 1, 6
        dd = abs(sh(ip)%z(k) - sd(ip)%z(k)); maxdiff = max(maxdiff, dd)
      end do
    end do

    print '(A,A,A,I0,A,I0,A,ES12.4)', '  ', name, ': ', npt, ' lanes x ', &
         nstep, ' steps, max|host-device| = ', maxdiff
    if (maxdiff <= 1.0e-12_dp) then
      print '(A,A)', 'PASS  ', name//' device == host'
    else
      print '(A,A)', 'FAIL  ', name//' device != host'
      nfail = nfail + 1
    end if
  end subroutine run_model

end program test_cpp_canonical_device
