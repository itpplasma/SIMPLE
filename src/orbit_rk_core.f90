module orbit_rk_core
  ! Device-callable shared core for the implicit Gauss-collocation step used by
  ! the guiding-center (GC) integrator and the Classical Pauli Particle (CPP)
  ! pusher. The CPP guiding-center motion is the slow manifold of the Pauli
  ! particle, so its degenerate-Lagrangian Gauss residual on field_can_t (mu
  ! fixed) coincides with the GC Gauss residual; both map to one concrete
  ! kernel here, selected by an integer model code. No procedure pointers and
  ! no class() polymorphism in the hot path: dispatch is select case to
  ! inlinable !$acc routine seq kernels, the way orbit_symplectic_euler1
  ! already shares the symplectic-Euler residual/Jacobian/Newton between the
  ! CPU integrator and the OpenACC GPU kernel.
  !
  ! The arithmetic is lifted byte-identically from the former f_rk_gauss /
  ! jac_rk_gauss bodies in orbit_symplectic and the matching orbit_cpp bodies.
  ! coeff_rk_gauss stays in orbit_symplectic_base (plain data) and is called
  ! inside the kernels. The Newton shell here uses the device LU solve rk_solve;
  ! the GC CPU path keeps dgesv in orbit_symplectic for byte-identical results.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: twopi
  use field_can_mod, only: field_can_t, get_derivatives2, eval_field => evaluate
  use orbit_symplectic_base, only: symplectic_integrator_t, coeff_rk_gauss
  use linalg_lu_device, only: rk_solve

  implicit none
  private

  ! Integer model codes for select-case dispatch (no proc pointers / class()).
  integer, parameter, public :: MODEL_GC = 0, MODEL_CPP = 1, MODEL_FO_SYMPL = 2

  public :: rk_gauss_residual, rk_gauss_jacobian, rk_gauss_newton, rk_solve

contains

  ! Model-dispatched Gauss residual. Both MODEL_GC and MODEL_CPP evaluate the
  ! field-canonical degenerate-Lagrangian residual (the CPP slow manifold equals
  ! the GC equations at fixed mu), so they share one inlinable kernel.
  subroutine rk_gauss_residual(model, si, fs, s, x, fvec)
    !$acc routine seq
    integer, intent(in) :: model
    type(symplectic_integrator_t), intent(inout) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(inout) :: fs(:)
    real(dp), intent(in) :: x(4*s)
    real(dp), intent(out) :: fvec(4*s)

    select case (model)
    case (MODEL_GC, MODEL_CPP)
      call gauss_canfield_residual(si, fs, s, x, fvec)
    case default
      fvec = 0d0
    end select
  end subroutine rk_gauss_residual

  ! Model-dispatched Gauss Jacobian, analytic from the 2nd derivatives of
  ! field_can_t. Shared kernel for MODEL_GC and MODEL_CPP.
  subroutine rk_gauss_jacobian(model, si, fs, s, jac)
    !$acc routine seq
    integer, intent(in) :: model
    type(symplectic_integrator_t), intent(in) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(in) :: fs(:)
    real(dp), intent(out) :: jac(4*s, 4*s)

    select case (model)
    case (MODEL_GC, MODEL_CPP)
      call gauss_canfield_jacobian(si, fs, s, jac)
    case default
      jac = 0d0
    end select
  end subroutine rk_gauss_jacobian

  ! Concrete residual kernel: layout x(4*k-3:4*k) = (r,theta,phi,pphi) per stage.
  ! Bodies lifted verbatim from f_rk_gauss / f_rk_gauss_cpp.
  subroutine gauss_canfield_residual(si, fs, s, x, fvec)
    !$acc routine seq
    type(symplectic_integrator_t), intent(inout) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(inout) :: fs(:)
    real(dp), intent(in) :: x(4*s)
    real(dp), intent(out) :: fvec(4*s)

    real(dp) :: a(s,s), b(s), c(s), Hprime(s)
    integer :: k, l

    call coeff_rk_gauss(s, a, b, c)

    do k = 1, s
      call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 2)
      call get_derivatives2(fs(k), x(4*k))
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
    end do

    do k = 1, s
      fvec(4*k-3) = fs(k)%pth - si%pthold
      fvec(4*k-2) = x(4*k-2)  - si%z(2)
      fvec(4*k-1) = x(4*k-1)  - si%z(3)
      fvec(4*k)   = x(4*k)    - si%z(4)
      do l = 1, s
        fvec(4*k-3) = fvec(4*k-3) + si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
        fvec(4*k-2) = fvec(4*k-2) - si%dt*a(k,l)*Hprime(l)                                      ! thdot
        fvec(4*k-1) = fvec(4*k-1) - si%dt*a(k,l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
        fvec(4*k)   = fvec(4*k)   + si%dt*a(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
      end do
    end do
  end subroutine gauss_canfield_residual

  ! Concrete Jacobian kernel. Bodies lifted verbatim from jac_rk_gauss /
  ! jac_rk_gauss_cpp.
  subroutine gauss_canfield_jacobian(si, fs, s, jac)
    !$acc routine seq
    type(symplectic_integrator_t), intent(in) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(in) :: fs(:)
    real(dp), intent(out) :: jac(4*s, 4*s)

    real(dp) :: a(s,s), b(s), c(s), Hprime(s), dHprime(4*s)
    integer :: k, l, m

    call coeff_rk_gauss(s, a, b, c)

    do k = 1, s
      m = 4*k
      Hprime(k)    = fs(k)%dH(1)/fs(k)%dpth(1)
      dHprime(m-3) = (fs(k)%d2H(1) - Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)  ! d/dr
      dHprime(m-2) = (fs(k)%d2H(2) - Hprime(k)*fs(k)%d2pth(2))/fs(k)%dpth(1)  ! d/dth
      dHprime(m-1) = (fs(k)%d2H(3) - Hprime(k)*fs(k)%d2pth(3))/fs(k)%dpth(1)  ! d/dph
      dHprime(m)   = (fs(k)%d2H(7) - Hprime(k)*fs(k)%d2pth(7))/fs(k)%dpth(1)  ! d/dpph
    end do

    jac = 0d0

    do k = 1, s
      m = 4*k
      jac(m-3, m-3) = fs(k)%dpth(1)
      jac(m-3, m-2) = fs(k)%dpth(2)
      jac(m-3, m-1) = fs(k)%dpth(3)
      jac(m-3, m)   = fs(k)%dpth(4)
      jac(m-2, m-2) = 1d0
      jac(m-1, m-1) = 1d0
      jac(m, m) = 1d0
    end do

    do l = 1, s
      do k = 1, s
          m = 4*k
          jac(m-3, 4*l-3) = jac(m-3, 4*l-3) & ! d/dr
            + si%dt*a(k,l)*(fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-3))
          jac(m-3, 4*l-2) = jac(m-3, 4*l-2) & ! d/dth
            + si%dt*a(k,l)*(fs(l)%d2H(4) - fs(l)%d2pth(4)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-2))
          jac(m-3, 4*l-1) = jac(m-3, 4*l-1) & ! d/dph
            + si%dt*a(k,l)*(fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-1))
          jac(m-3, 4*l) = jac(m-3, 4*l) & ! d/dpph
            + si%dt*a(k,l)*(fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l))

        jac(m-2, 4*l-3)   = jac(m-2, 4*l-3)   - si%dt*a(k,l)*dHprime(4*l-3)   ! d/dr
        jac(m-2, 4*l-2) = jac(m-2, 4*l-2) - si%dt*a(k,l)*dHprime(4*l-2) ! d/dth
        jac(m-2, 4*l-1) = jac(m-2, 4*l-1) - si%dt*a(k,l)*dHprime(4*l-1) ! d/dph
        jac(m-2, 4*l) = jac(m-2, 4*l) - si%dt*a(k,l)*dHprime(4*l) ! d/dpph

        jac(m-1, 4*l-3) = jac(m-1, 4*l-3) & ! d/dr
          - si%dt*a(k,l)*(-fs(l)%dhph(1)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
            + (fs(l)%dvpar(1) - dHprime(4*l-3)*fs(l)%hth - Hprime(l)*fs(l)%dhth(1))/fs(l)%hph)
        jac(m-1, 4*l-2) = jac(m-1, 4*l-2) & ! d/dth
          - si%dt*a(k,l)*(-fs(l)%dhph(2)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
            + (fs(l)%dvpar(2) - dHprime(4*l-2)*fs(l)%hth - Hprime(l)*fs(l)%dhth(2))/fs(l)%hph)
        jac(m-1, 4*l-1) = jac(m-1, 4*l-1) & ! d/dph
          - si%dt*a(k,l)*(-fs(l)%dhph(3)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
            + (fs(l)%dvpar(3) - dHprime(4*l-1)*fs(l)%hth - Hprime(l)*fs(l)%dhth(3))/fs(l)%hph)
        jac(m-1, 4*l) = jac(m-1, 4*l) & ! d/dpph
          - si%dt*a(k,l)*((fs(l)%dvpar(4) - dHprime(4*l)*fs(l)%hth)/fs(l)%hph)

        jac(m, 4*l-3) = jac(m, 4*l-3) & ! d/dr
          + si%dt*a(k,l)*(fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-3))
        jac(m, 4*l-2) = jac(m, 4*l-2) & ! d/dth
          + si%dt*a(k,l)*(fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-2))
        jac(m, 4*l-1) = jac(m, 4*l-1) & ! d/dph
          + si%dt*a(k,l)*(fs(l)%d2H(6) - fs(l)%d2pth(6)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-1))
        jac(m, 4*l) = jac(m, 4*l) & ! d/dpph
          + si%dt*a(k,l)*(fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l))
      end do
    end do
  end subroutine gauss_canfield_jacobian

  ! Device-callable Newton iteration for the Gauss step. Mirrors the
  ! newton_rk_gauss control flow (atol/rtol/tolref, boundary guards, maxit) with
  ! the device LU solver rk_solve in place of dgesv. No event counters here:
  ! hit_maxit is returned so the CPU caller can record EVT_RK_GAUSS_MAXIT,
  ! exactly as the GPU newton1 path omits fort.6601.
  subroutine rk_gauss_newton(model, si, fs, s, x, atol, rtol, maxit, xlast, &
                             hit_maxit)
    !$acc routine seq
    integer, intent(in) :: model
    type(symplectic_integrator_t), intent(inout) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(inout) :: fs(:)
    real(dp), intent(inout) :: x(4*s)
    real(dp), intent(in) :: atol, rtol
    integer, intent(in) :: maxit
    real(dp), intent(out) :: xlast(4*s)
    logical, intent(out) :: hit_maxit

    real(dp) :: fvec(4*s), fjac(4*s, 4*s)
    real(dp) :: xabs(4*s), tolref(4*s), fabs(4*s)
    integer :: kit, ks, info
    logical :: conv

    hit_maxit = .false.
    do kit = 1, maxit

      ! Check if radius left the boundary
      do ks = 1, s
        if (x(4*ks-3) > 1d0) return
        ! Transient guard for intermediate iterates; the converged-negative
        ! case is handled by the caller via a chart switch (#370).
        if (x(4*ks-3) < 0.0d0) x(4*ks-3) = 0.01d0
      end do

      call rk_gauss_residual(model, si, fs, s, x, fvec)
      call rk_gauss_jacobian(model, si, fs, s, fjac)
      fabs = abs(fvec)
      xlast = x
      call rk_solve(4*s, fjac, fvec, info)
      if (info /= 0) return
      ! after solution: fvec = (xold-xnew)_Newton
      x = x - fvec
      xabs = abs(x - xlast)

      do ks = 1, s
        tolref(4*ks-3) = 1d0
        tolref(4*ks-2) = twopi
        tolref(4*ks-1) = twopi
        tolref(4*ks)   = abs(xlast(4*ks))
      end do

      conv = .true.
      do ks = 1, 4*s
        if (fabs(ks) >= atol) conv = .false.
      end do
      if (conv) return
      conv = .true.
      do ks = 1, 4*s
        if (xabs(ks) >= rtol*tolref(ks)) conv = .false.
      end do
      if (conv) return
    end do
    hit_maxit = .true.
  end subroutine rk_gauss_newton

end module orbit_rk_core
