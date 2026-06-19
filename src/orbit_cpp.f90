module orbit_cpp
  ! Classical Pauli Particle (CPP) pusher, flux-canonical realization
  ! (Xiao & Qin, CPC 265 (2021) 107981). The guiding-center motion is the slow
  ! manifold of the Pauli particle; a structure-preserving (variational /
  ! Gauss-collocation) integrator stays on it and reproduces GC at GC-sized
  ! (bounce-scale) steps.
  !
  ! State and field machinery are shared with the GC integrator verbatim:
  !   z(4) = (r, theta, phi, p_phi) in symplectic_integrator_t,
  !   field_can_t carries mu (fixed parameter), ro0, and the canonical field
  !   quantities Ath, Aph, hth, hph, Bmod with 1st and 2nd derivatives.
  !
  ! The discrete scheme is the degenerate-Lagrangian Euler-Lagrange system
  ! (implicit Gauss collocation) on field_can_t with mu held fixed. Because the
  ! GC canonical equations ARE the slow-manifold equations of the Pauli
  ! particle, the CPP Gauss residual coincides with the GC Gauss residual at
  ! fixed mu; test_cpp_equals_gc_largestep verifies the trajectories agree to
  ! Newton tolerance.
  !
  ! GPU portability: residual, Jacobian, the dense LU solve, and the Newton
  ! shell are pure, fixed-size, !$acc routine seq-able. No procedure pointers,
  ! no class() polymorphism in the per-step loop; the only runtime indirection
  ! is the field-evaluation pointer in field_can_mod (shared with GC and resolved
  ! once at init). The stage count s (<= S_CPP_MAX) parameterizes GAUSS1..4.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: pi
  use field_can_mod, only: field_can_t, get_derivatives2, eval_field => evaluate
  use orbit_symplectic_base, only: symplectic_integrator_t, coeff_rk_gauss, &
    GAUSS1, GAUSS2, GAUSS3, GAUSS4, extrap_field
  use orbit_rk_core, only: rk_gauss_newton, MODEL_CPP
  use diag_counters, only: count_event, EVT_RK_GAUSS_MAXIT, EVT_R_NEGATIVE

  implicit none
  private

  integer, parameter, public :: S_CPP_MAX = 4

  public :: orbit_cpp_init, orbit_timestep_cpp, cpp_stages_from_mode

contains

  ! Map an integmode-style mode (GAUSS1..GAUSS4) to a Gauss stage count.
  pure function cpp_stages_from_mode(mode) result(s)
    integer, intent(in) :: mode
    integer :: s
    select case (mode)
    case (GAUSS1)
      s = 1
    case (GAUSS2)
      s = 2
    case (GAUSS3)
      s = 3
    case (GAUSS4)
      s = 4
    case default
      s = 1
    end select
  end function cpp_stages_from_mode

  ! Initialize a CPP step from an already-prepared (si, f). The caller fills
  ! si%z, si%dt, si%ntau, si%atol, si%rtol and f%mu, f%ro0 (e.g. via
  ! simple%init_sympl, which sets mu from the initial pitch identically to GC).
  ! This seeds f%pth at the current state so the first residual is consistent.
  subroutine orbit_cpp_init(si, f)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f

    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives2(f, si%z(4))
  end subroutine orbit_cpp_init

  ! One CPP macro-step (si%ntau micro-steps of si%dt). Advances si%z and f
  ! exactly as orbit_timestep_sympl_rk_gauss does for GC, on the CPP residual.
  recursive subroutine orbit_timestep_cpp(si, f, s, ierr)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(in) :: s
    integer, intent(out) :: ierr

    integer, parameter :: maxit = 32
    real(dp) :: x(4*s), xlast(4*s)
    real(dp) :: a(s,s), b(s), c(s), Hprime(s), dz(4)
    type(field_can_t) :: fs(s)
    integer :: k, l, ktau
    logical :: hit_maxit

    do k = 1, s
      fs(k) = f
    end do

    ierr = 0
    ktau = 0
    do while (ktau < si%ntau)
      si%pthold = f%pth

      do k = 1, s
        x((4*k-3):(4*k)) = si%z
      end do

      call rk_gauss_newton(MODEL_CPP, si, fs, s, x, si%atol, si%rtol, maxit, &
                           xlast, hit_maxit)
      if (hit_maxit) call count_event(EVT_RK_GAUSS_MAXIT)

      if (x(1) > 1.0d0) then
        ierr = 1
        return
      end if

      if (x(1) < 0.0d0) then
        x(1) = -x(1)
        x(2) = x(2) + pi
        call count_event(EVT_R_NEGATIVE)
        if (x(1) > 1.0d0) then
          ierr = 1
          return
        end if
      end if

      call coeff_rk_gauss(s, a, b, c)

      if (extrap_field) then
        do k = 1, s
          dz(1) = x(4*k-3) - xlast(4*k-3)
          dz(2) = x(4*k-2) - xlast(4*k-2)
          dz(3) = x(4*k-1) - xlast(4*k-1)
          dz(4) = x(4*k)   - xlast(4*k)

          fs(k)%pth = fs(k)%pth + fs(k)%dpth(1)*dz(1) + fs(k)%dpth(2)*dz(2) &
                    + fs(k)%dpth(3)*dz(3) + fs(k)%dpth(4)*dz(4)
          fs(k)%dpth(1) = fs(k)%dpth(1) + fs(k)%d2pth(1)*dz(1) + fs(k)%d2pth(2)*dz(2) &
                        + fs(k)%d2pth(3)*dz(3) + fs(k)%d2pth(7)*dz(4)
          fs(k)%dpth(2) = fs(k)%dpth(2) + fs(k)%d2pth(2)*dz(1) + fs(k)%d2pth(4)*dz(2) &
                        + fs(k)%d2pth(5)*dz(3) + fs(k)%d2pth(8)*dz(4)
          fs(k)%dpth(3) = fs(k)%dpth(3) + fs(k)%d2pth(2)*dz(1) + fs(k)%d2pth(5)*dz(2) &
                        + fs(k)%d2pth(6)*dz(3) + fs(k)%d2pth(9)*dz(4)
          fs(k)%dH(1) = fs(k)%dH(1) + fs(k)%d2H(1)*dz(1) + fs(k)%d2H(2)*dz(2) &
                      + fs(k)%d2H(3)*dz(3) + fs(k)%d2H(7)*dz(4)
          fs(k)%dH(2) = fs(k)%dH(2) + fs(k)%d2H(2)*dz(1) + fs(k)%d2H(4)*dz(2) &
                      + fs(k)%d2H(5)*dz(3) + fs(k)%d2H(8)*dz(4)
          fs(k)%dH(3) = fs(k)%dH(3) + fs(k)%d2H(3)*dz(1) + fs(k)%d2H(5)*dz(2) &
                      + fs(k)%d2H(6)*dz(3) + fs(k)%d2H(9)*dz(4)
          fs(k)%vpar = fs(k)%vpar + fs(k)%dvpar(1)*dz(1) + fs(k)%dvpar(2)*dz(2) &
                     + fs(k)%dvpar(3)*dz(3)
          fs(k)%hth = fs(k)%hth + fs(k)%dhth(1)*dz(1) + fs(k)%dhth(2)*dz(2) &
                    + fs(k)%dhth(3)*dz(3)
          fs(k)%hph = fs(k)%hph + fs(k)%dhph(1)*dz(1) + fs(k)%dhph(2)*dz(2) &
                    + fs(k)%dhph(3)*dz(3)
        end do
      else
        do k = 1, s
          call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 2)
          call get_derivatives2(fs(k), x(4*k))
        end do
      end if

      do k = 1, s
        Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
      end do

      f = fs(s)
      f%pth = si%pthold
      si%z(1) = x(4*s-3)

      do l = 1, s
        f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
        si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
        si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph
        si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
      end do

      ktau = ktau + 1
    end do
  end subroutine orbit_timestep_cpp

end module orbit_cpp
