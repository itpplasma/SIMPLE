module orbit_cpp
  ! Classical Pauli Particle (CPP) pusher, flux-canonical realization
  ! (Xiao & Qin, CPC 265 (2021) 107981). The guiding-center motion is the slow
  ! manifold of the Pauli particle; a structure-preserving (variational /
  ! Gauss-collocation) integrator stays on it and reproduces GC at GC-sized
  ! (bounce-scale) steps.
  !
  ! TAUTOLOGICAL BY DESIGN. This residual is the GC degenerate-Lagrangian
  ! Euler-Lagrange system specialized to fixed mu. Because the GC canonical
  ! equations ARE the slow-manifold equations of the Pauli particle in these
  ! flux-canonical coordinates, this residual is byte-identical to the GC Gauss
  ! residual. So "CPP == GC" here is an IDENTITY, not a physics cross-check: the
  ! accompanying tests (test_cpp_equals_gc_largestep, test_cpp_invariants) are
  ! refactor / code-motion correctness ORACLES on the shared symplectic core,
  ! verifying that this device-portable Newton/LU realization reproduces the GC
  ! trajectory it is built from. The genuine, non-tautological CPP -- a full 6D
  ! particle that carries real gyration and is a DIFFERENT method from GC -- is
  ! orbit_cpp_pauli; its banana matches GC only to O(rho*).
  !
  ! State and field machinery are shared with the GC integrator verbatim:
  !   z(4) = (r, theta, phi, p_phi) in symplectic_integrator_t,
  !   field_can_t carries mu (fixed parameter), ro0, and the canonical field
  !   quantities Ath, Aph, hth, hph, Bmod with 1st and 2nd derivatives.
  !
  ! The discrete scheme is the degenerate-Lagrangian Euler-Lagrange system
  ! (implicit Gauss collocation) on field_can_t with mu held fixed.
  !
  ! GPU portability: residual, Jacobian, the dense LU solve, and the Newton
  ! shell are pure, fixed-size, !$acc routine seq-able. No procedure pointers,
  ! no class() polymorphism in the per-step loop; the only runtime indirection
  ! is the field-evaluation pointer in field_can_mod (shared with GC and resolved
  ! once at init). The stage count s (<= S_CPP_MAX) parameterizes GAUSS1..4.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: pi, twopi
  use field_can_mod, only: field_can_t, get_derivatives2, eval_field => evaluate
  use orbit_symplectic_base, only: symplectic_integrator_t, coeff_rk_gauss, &
    GAUSS1, GAUSS2, GAUSS3, GAUSS4, extrap_field
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

  ! Degenerate-Lagrangian Euler-Lagrange residual on field_can_t, mu fixed.
  ! Layout x(4*k-3:4*k) = (r,theta,phi,pphi) for stage k. fvec same layout.
  subroutine f_rk_gauss_cpp(si, fs, s, x, fvec)
    !$acc routine seq
    type(symplectic_integrator_t), intent(in) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(inout) :: fs(s)
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
      fvec(4*k-2) = x(4*k-2) - si%z(2)
      fvec(4*k-1) = x(4*k-1) - si%z(3)
      fvec(4*k)   = x(4*k)   - si%z(4)
      do l = 1, s
        fvec(4*k-3) = fvec(4*k-3) + si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
        fvec(4*k-2) = fvec(4*k-2) - si%dt*a(k,l)*Hprime(l)
        fvec(4*k-1) = fvec(4*k-1) - si%dt*a(k,l)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph
        fvec(4*k)   = fvec(4*k)   + si%dt*a(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
      end do
    end do
  end subroutine f_rk_gauss_cpp

  ! Analytic Jacobian of f_rk_gauss_cpp using the 2nd derivatives of field_can_t.
  subroutine jac_rk_gauss_cpp(si, fs, s, jac)
    !$acc routine seq
    type(symplectic_integrator_t), intent(in) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(in) :: fs(s)
    real(dp), intent(out) :: jac(4*s, 4*s)

    real(dp) :: a(s,s), b(s), c(s), Hprime(s), dHprime(4*s)
    integer :: k, l, m

    call coeff_rk_gauss(s, a, b, c)

    do k = 1, s
      m = 4*k
      Hprime(k)    = fs(k)%dH(1)/fs(k)%dpth(1)
      dHprime(m-3) = (fs(k)%d2H(1) - Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)
      dHprime(m-2) = (fs(k)%d2H(2) - Hprime(k)*fs(k)%d2pth(2))/fs(k)%dpth(1)
      dHprime(m-1) = (fs(k)%d2H(3) - Hprime(k)*fs(k)%d2pth(3))/fs(k)%dpth(1)
      dHprime(m)   = (fs(k)%d2H(7) - Hprime(k)*fs(k)%d2pth(7))/fs(k)%dpth(1)
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
        jac(m-3, 4*l-3) = jac(m-3, 4*l-3) &
          + si%dt*a(k,l)*(fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-3))
        jac(m-3, 4*l-2) = jac(m-3, 4*l-2) &
          + si%dt*a(k,l)*(fs(l)%d2H(4) - fs(l)%d2pth(4)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-2))
        jac(m-3, 4*l-1) = jac(m-3, 4*l-1) &
          + si%dt*a(k,l)*(fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-1))
        jac(m-3, 4*l) = jac(m-3, 4*l) &
          + si%dt*a(k,l)*(fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l))

        jac(m-2, 4*l-3) = jac(m-2, 4*l-3) - si%dt*a(k,l)*dHprime(4*l-3)
        jac(m-2, 4*l-2) = jac(m-2, 4*l-2) - si%dt*a(k,l)*dHprime(4*l-2)
        jac(m-2, 4*l-1) = jac(m-2, 4*l-1) - si%dt*a(k,l)*dHprime(4*l-1)
        jac(m-2, 4*l)   = jac(m-2, 4*l)   - si%dt*a(k,l)*dHprime(4*l)

        jac(m-1, 4*l-3) = jac(m-1, 4*l-3) &
          - si%dt*a(k,l)*(-fs(l)%dhph(1)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
            + (fs(l)%dvpar(1) - dHprime(4*l-3)*fs(l)%hth - Hprime(l)*fs(l)%dhth(1))/fs(l)%hph)
        jac(m-1, 4*l-2) = jac(m-1, 4*l-2) &
          - si%dt*a(k,l)*(-fs(l)%dhph(2)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
            + (fs(l)%dvpar(2) - dHprime(4*l-2)*fs(l)%hth - Hprime(l)*fs(l)%dhth(2))/fs(l)%hph)
        jac(m-1, 4*l-1) = jac(m-1, 4*l-1) &
          - si%dt*a(k,l)*(-fs(l)%dhph(3)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
            + (fs(l)%dvpar(3) - dHprime(4*l-1)*fs(l)%hth - Hprime(l)*fs(l)%dhth(3))/fs(l)%hph)
        jac(m-1, 4*l) = jac(m-1, 4*l) &
          - si%dt*a(k,l)*((fs(l)%dvpar(4) - dHprime(4*l)*fs(l)%hth)/fs(l)%hph)

        jac(m, 4*l-3) = jac(m, 4*l-3) &
          + si%dt*a(k,l)*(fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-3))
        jac(m, 4*l-2) = jac(m, 4*l-2) &
          + si%dt*a(k,l)*(fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-2))
        jac(m, 4*l-1) = jac(m, 4*l-1) &
          + si%dt*a(k,l)*(fs(l)%d2H(6) - fs(l)%d2pth(6)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-1))
        jac(m, 4*l) = jac(m, 4*l) &
          + si%dt*a(k,l)*(fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l))
      end do
    end do
  end subroutine jac_rk_gauss_cpp

  ! Dense LU solve A x = rhs with partial pivoting, in place on rhs. Device
  ! portable replacement for dgesv (n <= 4*S_CPP_MAX). info=0 on success.
  pure subroutine lu_solve(n, A, rhs, info)
    !$acc routine seq
    integer, intent(in) :: n
    real(dp), intent(inout) :: A(n,n), rhs(n)
    integer, intent(out) :: info
    integer :: i, j, k, ipiv
    real(dp) :: piv, amax, factor, tmp

    info = 0
    do k = 1, n
      ipiv = k
      amax = abs(A(k,k))
      do i = k+1, n
        if (abs(A(i,k)) > amax) then
          amax = abs(A(i,k))
          ipiv = i
        end if
      end do
      if (amax == 0d0) then
        info = k
        return
      end if
      if (ipiv /= k) then
        do j = 1, n
          tmp = A(k,j); A(k,j) = A(ipiv,j); A(ipiv,j) = tmp
        end do
        tmp = rhs(k); rhs(k) = rhs(ipiv); rhs(ipiv) = tmp
      end if
      piv = A(k,k)
      do i = k+1, n
        factor = A(i,k)/piv
        A(i,k) = factor
        do j = k+1, n
          A(i,j) = A(i,j) - factor*A(k,j)
        end do
        rhs(i) = rhs(i) - factor*rhs(k)
      end do
    end do

    do i = n, 1, -1
      tmp = rhs(i)
      do j = i+1, n
        tmp = tmp - A(i,j)*rhs(j)
      end do
      rhs(i) = tmp/A(i,i)
    end do
  end subroutine lu_solve

  ! Newton iteration for the CPP Gauss step. Mirrors newton_rk_gauss control
  ! flow (atol/rtol/tolref, boundary guards, maxit) with the device LU solver.
  subroutine newton_cpp(si, fs, s, x, atol, rtol, maxit, xlast)
    !$acc routine seq
    type(symplectic_integrator_t), intent(inout) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(inout) :: fs(s)
    real(dp), intent(inout) :: x(4*s)
    real(dp), intent(in) :: atol, rtol
    integer, intent(in) :: maxit
    real(dp), intent(out) :: xlast(4*s)

    real(dp) :: fvec(4*s), fjac(4*s, 4*s)
    real(dp) :: xabs(4*s), tolref(4*s), fabs(4*s)
    integer :: kit, ks, info
    logical :: conv

    do kit = 1, maxit
      do ks = 1, s
        if (x(4*ks-3) > 1d0) return
        if (x(4*ks-3) < 0.0d0) x(4*ks-3) = 0.01d0
      end do

      call f_rk_gauss_cpp(si, fs, s, x, fvec)
      call jac_rk_gauss_cpp(si, fs, s, fjac)
      fabs = abs(fvec)
      xlast = x
      call lu_solve(4*s, fjac, fvec, info)
      if (info /= 0) return
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
    call count_event(EVT_RK_GAUSS_MAXIT)
  end subroutine newton_cpp

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

      call newton_cpp(si, fs, s, x, si%atol, si%rtol, maxit, xlast)

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
