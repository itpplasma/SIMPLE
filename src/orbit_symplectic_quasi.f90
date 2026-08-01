module orbit_symplectic_quasi

use util, only: pi
use field_can_mod, only: eval_field => evaluate, field_can_t, get_derivatives, &
  get_derivatives2
use orbit_symplectic_base, only: symplectic_integrator_t, multistage_integrator_t, &
  orbit_timestep_quasi_i, coeff_rk_gauss, coeff_rk_lobatto, f_rk_lobatto, sympl_rmax
use fortnum_multiroot, only: multiroot_hybrids
use fortnum_status, only: fortnum_status_t
use diag_counters, only: count_event, EVT_R_NEGATIVE

implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)
save

logical, parameter :: exact_steps = .False.
integer, parameter :: S_MAX_GAUSS = 3

type(field_can_t) :: f
type(field_can_t) :: fs(S_MAX_GAUSS)
type(symplectic_integrator_t) :: si
  !$omp threadprivate(f, si)

procedure(orbit_timestep_quasi_i), pointer :: orbit_timestep_quasi => null()

! Step size carried across macro-steps by the adaptive explicit integrators.
! Without this each macro-step restarts the step-size search from a default
! guess and spends most of its evaluations re-discovering a step it already
! knew, which inflates the evaluation count by more than an order of magnitude
! and would make any cost comparison meaningless.
real(dp) :: explicit_h_carry = 0d0
  !$omp threadprivate(explicit_h_carry)

contains

  !
  ! Wrapper routines for ODEPACK
  !

subroutine f_exact_quasi(x, fvec, ctx)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: fvec(:)
  class(*), intent(in), optional :: ctx

  type(field_can_t) :: f2

  f2 = f

  call eval_field(f2, x(1), si%z(2), si%z(3), 0)
  call get_derivatives(f2, si%z(4))

  fvec(1) = f%pth - f2%pth

  f%H = f2%H

end subroutine f_exact_quasi


subroutine f_euler1_quasi(x, fvec, ctx)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: fvec(:)
  class(*), intent(in), optional :: ctx

  call eval_field(f, x(1), si%z(2), si%z(3), 0)
  call get_derivatives(f, x(2))

  fvec(1) = f%dpth(1)*(f%pth - si%pthold) &
    + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
  fvec(2) = f%dpth(1)*(x(2) - si%z(4)) &
    + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))

end subroutine f_euler1_quasi


subroutine f_euler2_quasi(x, fvec, ctx)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: fvec(:)
  class(*), intent(in), optional :: ctx

  call eval_field(f, x(1), x(2), x(3), 0)
  call get_derivatives(f, si%z(4))

  fvec(1) = f%pth - si%pthold
  fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
  fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) &
    - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)

end subroutine f_euler2_quasi


subroutine f_midpoint_quasi(x, fvec, ctx)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: fvec(:)
    class(*), intent(in), optional :: ctx

    real(dp) :: dpthmid, pthdotbar

    call eval_field(f, x(5), 0.5*(x(2) + si%z(2)), 0.5*(x(3) + si%z(3)), 0)
    call get_derivatives(f, 0.5*(x(4) + si%z(4)))

    fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
    fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) &
      - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)
    fvec(4) = f%dpth(1)*(x(4) - si%z(4)) &
      + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
    fvec(5) = f%dpth(1)*(f%pth - si%pthold) &
      + si%dt/2.0d0*(f%dpth(1)*f%dH(2)-f%dpth(2)*f%dH(1))

    ! save evaluation from midpoint
    dpthmid = f%dpth(1)
    pthdotbar = f%dpth(1)*f%dH(2) - f%dpth(2)*f%dH(1)

    ! evaluate at endpoint
    call eval_field(f, x(1), x(2), x(3), 0)
    call get_derivatives(f, x(4))
    fvec(1) = dpthmid*(f%pth - si%pthold) + si%dt*pthdotbar

end subroutine f_midpoint_quasi


subroutine f_rk_gauss_quasi(x, fvec, ctx)

  real(dp), intent(in) :: x(:)  ! = (rend, thend, phend, pphend, rmid)
  real(dp), intent(out) :: fvec(:)
  class(*), intent(in), optional :: ctx

  integer :: n
  real(dp) :: a(size(x)/4,size(x)/4), b(size(x)/4), c(size(x)/4), Hprime(size(x)/4)
  integer :: k,l  ! counters

  n = size(x)
  call coeff_rk_gauss(n/4, a, b, c)  ! TODO: move this to preprocessing

  ! evaluate stages
  do k = 1, n/4
    call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 0)
    call get_derivatives(fs(k), x(4*k))
    Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
  end do

  do k = 1, n/4
    fvec(4*k-3) = fs(k)%pth - si%pthold
    fvec(4*k-2) = x(4*k-2)  - si%z(2)
    fvec(4*k-1) = x(4*k-1)  - si%z(3)
    fvec(4*k)   = x(4*k)    - si%z(4)
    do l = 1, n/4
      fvec(4*k-3) = fvec(4*k-3) &
        + si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))       ! pthdot
      fvec(4*k-2) = fvec(4*k-2) - si%dt*a(k,l)*Hprime(l)             ! thdot
      fvec(4*k-1) = fvec(4*k-1) &
        - si%dt*a(k,l)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
      fvec(4*k)   = fvec(4*k) &
        + si%dt*a(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))       ! pphdot
    end do
  end do

end subroutine f_rk_gauss_quasi

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! Lobatto (IIIA)-(IIIB) Runge-Kutta method, s stages (n=4*s-2 variables)
  !
subroutine f_rk_lobatto_quasi(x, fvec, ctx)
  !
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: fvec(:)
  class(*), intent(in), optional :: ctx

  call f_rk_lobatto(si, fs, (size(x)+2)/4, x, fvec, 0)

  end subroutine f_rk_lobatto_quasi


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine orbit_timestep_multi_quasi(mi, ierr)
  !
    type(multistage_integrator_t), intent(inout) :: mi

    integer, intent(out) :: ierr

    integer :: kstage

    si = mi%stages(1)
    call orbit_timestep_quasi(ierr)

    do kstage = 2, 2*mi%s
      mi%stages(kstage)%z = si%z
      mi%stages(kstage)%pthold = f%pth
      si = mi%stages(kstage)
      call orbit_timestep_quasi(ierr)
    end do
    mi%stages(1)%z = si%z
    mi%stages(1)%pthold = f%pth

  end subroutine orbit_timestep_multi_quasi

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine timestep_midpoint_quasi(ierr)
  !
  integer, intent(out) :: ierr

  integer, parameter :: n = 5

  real(dp), dimension(n) :: x, x0
  type(fortnum_status_t) :: status
  integer :: ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x0(1:4) = si%z
    x0(5) = si%z(1)

    call multiroot_hybrids(f_midpoint_quasi, n, x0, x, status, &
      xtol=si%rtol, ftol=si%rtol)

    if (x(1) > sympl_rmax) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      ! Axis crossing: continue on the opposite ray, (r, theta) ->
      ! (|r|, theta + pi), instead of the former teleport to r = 0.01 (#370).
      call count_event(EVT_R_NEGATIVE)
      x(1) = -x(1)
      x(2) = x(2) + pi
      if (x(1) > sympl_rmax) then
        ierr = 1
        return
      end if
    end if

    si%z = x(1:4)

    ktau = ktau+1
  enddo

end subroutine timestep_midpoint_quasi


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine timestep_expl_impl_euler_quasi(ierr)
  !
  integer, intent(out) :: ierr

  integer, parameter :: n = 2

  real(dp), dimension(n) :: x, x0
  real(dp), dimension(1) :: xe0, xe
  type(fortnum_status_t) :: status
  integer :: ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x0(1)=si%z(1)
    x0(2)=si%z(4)

    call multiroot_hybrids(f_euler1_quasi, n, x0, x, status, &
      xtol=si%rtol, ftol=si%rtol)

    if (x(1) > sympl_rmax) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      ! Axis crossing: continue on the opposite ray, (r, theta) ->
      ! (|r|, theta + pi), instead of the former teleport to r = 0.01 (#370).
      call count_event(EVT_R_NEGATIVE)
      x(1) = -x(1)
      x(2) = x(2) + pi
      if (x(1) > sympl_rmax) then
        ierr = 1
        return
      end if
    end if

    si%z(1) = x(1)
    si%z(4) = x(2)

    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives(f, si%z(4))

    si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
    si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph

    if (exact_steps) then
      xe0(1) = si%z(1)
      call multiroot_hybrids(f_exact_quasi, 1, xe0, xe, status, &
        xtol=si%rtol, ftol=si%rtol)
      si%z(1) = xe(1)
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    end if

    ktau = ktau+1
  enddo

end subroutine timestep_expl_impl_euler_quasi



  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine timestep_impl_expl_euler_quasi(ierr)
  !
  integer, intent(out) :: ierr

  integer, parameter :: n = 3

  real(dp), dimension(n) :: x, x0
  real(dp), dimension(1) :: xe0, xe
  type(fortnum_status_t) :: status
  integer :: ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x0 = si%z(1:3)

    call multiroot_hybrids(f_euler2_quasi, n, x0, x, status, &
      xtol=si%rtol, ftol=si%rtol)

    if (x(1) > sympl_rmax) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      ! Axis crossing: continue on the opposite ray, (r, theta) ->
      ! (|r|, theta + pi), instead of the former teleport to r = 0.01 (#370).
      call count_event(EVT_R_NEGATIVE)
      x(1) = -x(1)
      x(2) = x(2) + pi
      if (x(1) > sympl_rmax) then
        ierr = 1
        return
      end if
    end if

    si%z(1:3) = x

    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives(f, si%z(4))

    f%pth = si%pthold - si%dt*(f%dH(2) - f%dH(1)*f%dpth(2)/f%dpth(1))
    si%z(4) = si%z(4) - si%dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

    if (exact_steps) then
      xe0(1) = si%z(1)
      call multiroot_hybrids(f_exact_quasi, 1, xe0, xe, status, &
        xtol=si%rtol, ftol=si%rtol)
      si%z(1) = xe(1)
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    end if

    ktau = ktau+1
  enddo

end subroutine timestep_impl_expl_euler_quasi

subroutine timestep_rk_gauss_quasi(s, ierr)
  !
  integer, intent(out) :: ierr


  integer, intent(in) :: s
  real(dp), dimension(4*s) :: x, x0
  type(fortnum_status_t) :: status
  integer :: ktau, k, l

  real(dp) :: a(s,s), b(s), c(s), Hprime(s)

  do k = 1,s
    fs(k) = f
  end do

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    do k = 1,s
      x0((4*k-3):(4*k)) = si%z
    end do

    call multiroot_hybrids(f_rk_gauss_quasi, 4*s, x0, x, status, &
      xtol=si%rtol, ftol=si%rtol)

    if (x(1) > sympl_rmax) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      ! Axis crossing: continue on the opposite ray, (r, theta) ->
      ! (|r|, theta + pi), instead of the former teleport to r = 0.01 (#370).
      call count_event(EVT_R_NEGATIVE)
      x(1) = -x(1)
      x(2) = x(2) + pi
      if (x(1) > sympl_rmax) then
        ierr = 1
        return
      end if
    end if

    call coeff_rk_gauss(s, a, b, c)

    f = fs(s)
    f%pth = si%pthold
    si%z(1) = x(4*s-3)

    do l = 1, s
      Hprime(l) = fs(l)%dH(1)/fs(l)%dpth(1)
      f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
      si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
      si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar-Hprime(l)*fs(l)%hth)/fs(l)%hph
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
    end do

    ktau = ktau+1
  enddo

end subroutine timestep_rk_gauss_quasi


subroutine timestep_rk_lobatto_quasi(s, ierr)
  !
  integer, intent(out) :: ierr


  integer, intent(in) :: s
  real(dp), dimension(4*s-2) :: x, x0
  type(fortnum_status_t) :: status
  integer :: ktau, k, l

  real(dp) :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s)

  do k = 1,s
    fs(k) = f
  end do

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x0(1) = si%z(1)
    x0(2) = si%z(4)
    do k = 2,s
      x0((4*k-3-2):(4*k-2)) = si%z
    end do

    !
    ! q1 = q0
    ! q2 = q0 + h*(a21*f(w1) + a22*f(z2) + a23*f(z3)
    ! q3 = q0 + h*(a31*f(w1) + a32*f(z2) + a33*f(z3))
    !
    ! p1(w1) = p0 + h*(o11*f(z1) + o12*f(z2))
    ! p2(z2) = p0 + h*(o21*f(z1) + o22*f(z2))
    ! p3(z3) = p0 + h*(o31*f(z1) + o32*f(z2))
    !
    ! 0 0 0 0 0 0
    ! x x x x x x
    ! x x x x x x
    ! x x x x 0 0
    ! x x x x 0 0
    ! x x x x 0 0
    !
    ! effectively
    !
    ! x x x x x
    ! x x x x x
    ! x x x 0 0
    ! x x x 0 0
    ! x x x 0 0

    call multiroot_hybrids(f_rk_lobatto_quasi, 4*s-2, x0, x, status, &
      xtol=si%rtol, ftol=si%rtol)

    call coeff_rk_lobatto(s, a, ahat, b, c)

    call eval_field(fs(1), x(1), si%z(2), si%z(3), 0)
    call get_derivatives(fs(1), x(2))
    Hprime(1) = fs(1)%dH(1)/fs(1)%dpth(1)

    do k = 2, s
      call eval_field(fs(k), x(4*k-3-2), x(4*k-2-2), x(4*k-1-2), 0)
      call get_derivatives(fs(k), x(4*k-2))
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
    end do

    f = fs(s)
    f%pth   = si%pthold
    si%z(1) = x(4*s-3-2)

    do l = 1, s
      f%pth   = f%pth   - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
      si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
      si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar-Hprime(l)*fs(l)%hth)/fs(l)%hph
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
    end do

    ktau = ktau+1
  enddo

end subroutine timestep_rk_lobatto_quasi


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine f_ode(tau, z, zdot)
  !

  real(dp), intent(in)  :: tau
  real(dp), intent(in)  :: z(:)
  real(dp), intent(out) :: zdot(:)
  real(dp) :: Hprime

  call eval_field(f, z(1), z(2), z(3), 0)
  call get_derivatives(f, z(4))

  Hprime = f%dH(1)/f%dpth(1)

  zdot(1) = -(f%dH(2) - f%hth/f%hph*f%dH(3))/f%dpth(1)
  zdot(2) = Hprime
  zdot(3) = (f%vpar-Hprime*f%hth)/f%hph
  zdot(4) = -(f%dH(3) - Hprime*f%dpth(3))

end subroutine f_ode

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
! Index into the packed symmetric second-derivative arrays d2H, d2pth, d2vpar.
!
! field_can_base documents the order as
!   1:(r,r) 2:(r,th) 3:(r,ph) 4:(th,th) 5:(th,ph) 6:(ph,ph)
!   7:(pph,r) 8:(pph,th) 9:(pph,ph) 10:(pph,pph)
! over the variables (r, th, ph, pph) = (1, 2, 3, 4).
pure function d2idx(i, j) result(k)
  integer, intent(in) :: i, j
  integer :: k
  integer :: lo, hi

  lo = min(i, j)
  hi = max(i, j)

  if (hi == 4) then
    k = 6 + lo            ! (r,pph)=7, (th,pph)=8, (ph,pph)=9, (pph,pph)=10
  else if (lo == 1) then
    k = hi                ! (r,r)=1, (r,th)=2, (r,ph)=3
  else if (lo == 2) then
    k = hi + 2            ! (th,th)=4, (th,ph)=5
  else
    k = 6                 ! (ph,ph)
  end if
end function d2idx

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
! Second time derivative of the canonical guiding-centre state,
!
!   G(z) = F'(z) F(z),
!
! which is what a two-derivative Runge-Kutta method integrates in place of F.
!
! The derivation and its verification live in the study repo at
! derivation/tdrk_guiding_center_G.wl: Mathematica differentiates f_ode's F
! symbolically and compares against this closed form. The symbolic residual is
! exactly zero and a numeric spot-check with nonlinear, fully coupled test
! functions agrees to 3.6e-15.
!
! The cost argument rests on which derivatives appear. That check is in the same
! file and reports second derivatives of H and pth, first of vpar, hth and hph,
! and no third derivatives anywhere -- so G needs exactly one mode_secders = 2
! evaluation and no extra F evaluations. Note also that only d2H indices 1 to 9
! are touched: F depends on H solely through dH(1:3), so differentiating it can
! never raise index 10, which get_derivatives2 does not fill.
!
! hth and hph are geometric and carry no pph dependence, hence the explicit zero
! in the fourth slot of their gradients rather than an out-of-bounds read of a
! three-element array.
subroutine g_ode(tau, z, g)
  real(dp), intent(in)  :: tau
  real(dp), intent(in)  :: z(:)
  real(dp), intent(out) :: g(:)

  real(dp) :: Hprime, dHprime(4), dhth4(4), dhph4(4)
  real(dp) :: zdot(4), dF(4, 4)
  integer  :: j

  call eval_field(f, z(1), z(2), z(3), 2)
  call get_derivatives2(f, z(4))

  Hprime = f%dH(1)/f%dpth(1)

  dhth4(1:3) = f%dhth
  dhth4(4)   = 0d0
  dhph4(1:3) = f%dhph
  dhph4(4)   = 0d0

  zdot(1) = -(f%dH(2) - f%hth/f%hph*f%dH(3))/f%dpth(1)
  zdot(2) = Hprime
  zdot(3) = (f%vpar - Hprime*f%hth)/f%hph
  zdot(4) = -(f%dH(3) - Hprime*f%dpth(3))

  ! Hprime is a quotient of two state-dependent quantities, so both numerator
  ! and denominator contribute. This is the term a hand derivation usually drops.
  do j = 1, 4
    dHprime(j) = (f%d2H(d2idx(1, j))*f%dpth(1) &
                - f%dH(1)*f%d2pth(d2idx(1, j)))/f%dpth(1)**2
  end do

  do j = 1, 4
    dF(1, j) = -((f%d2H(d2idx(2, j)) &
               - (dhth4(j)/f%hph - f%hth*dhph4(j)/f%hph**2)*f%dH(3) &
               - f%hth/f%hph*f%d2H(d2idx(3, j)))*f%dpth(1) &
               - (f%dH(2) - f%hth/f%hph*f%dH(3))*f%d2pth(d2idx(1, j))) &
               /f%dpth(1)**2

    dF(2, j) = dHprime(j)

    dF(3, j) = ((f%dvpar(j) - dHprime(j)*f%hth - Hprime*dhth4(j))*f%hph &
               - (f%vpar - Hprime*f%hth)*dhph4(j))/f%hph**2

    dF(4, j) = -(f%d2H(d2idx(3, j)) - dHprime(j)*f%dpth(3) &
               - Hprime*f%d2pth(d2idx(3, j)))
  end do

  g(1:4) = matmul(dF, zdot)
end subroutine g_ode

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
! Adapter from g_ode to fortnum's ode_rhs2_t.
subroutine g_ode_fortnum(t, y, g, ctx)
  real(dp), intent(in)  :: t
  real(dp), intent(in)  :: y(:)
  real(dp), intent(out) :: g(:)
  class(*), intent(in), optional :: ctx

  call g_ode(t, y, g)
end subroutine g_ode_fortnum

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine orbit_timestep_rk45(ierr)
  !
  use odeint_allroutines_sub, only : odeint_allroutines
  integer, intent(out) :: ierr
  integer :: ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    call odeint_allroutines(si%z, 4, ktau*si%dt, (ktau+1)*si%dt, si%rtol, f_ode)
    ktau = ktau+1
  end do
end subroutine orbit_timestep_rk45

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
! Adapter from f_ode to fortnum's ode_rhs_t, which carries an optional context
! argument that this right-hand side does not need.
subroutine f_ode_fortnum(t, y, dydt, ctx)
  real(dp), intent(in)  :: t
  real(dp), intent(in)  :: y(:)
  real(dp), intent(out) :: dydt(:)
  class(*), intent(in), optional :: ctx

  call f_ode(t, y, dydt)
end subroutine f_ode_fortnum

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
! Gauss-Radau (IAS15-style, 15th order) on the same canonical Hamiltonian the
! symplectic schemes use.
subroutine orbit_timestep_radau15(ierr)
  !
  use fortnum_ode_gauss_radau, only : ode_integrate_radau
  use fortnum_ode, only : ode_problem_t, ode_solution_t
  use fortnum_status, only : fortnum_status_t, FORTNUM_OK
  integer, intent(out) :: ierr
  integer :: ktau, nlast
  type(ode_problem_t) :: problem
  type(ode_solution_t) :: solution
  type(fortnum_status_t) :: status

  ierr = 0
  ktau = 0
  problem%rhs => f_ode_fortnum
  problem%rtol = si%rtol
  problem%atol = si%atol
  allocate(problem%y0(4))
  do while(ktau .lt. si%ntau)
    problem%t0 = ktau*si%dt
    problem%t1 = (ktau+1)*si%dt
    problem%y0 = si%z(1:4)
    problem%h0 = explicit_h_carry
    call ode_integrate_radau(problem, solution, status)
    if (status%code /= FORTNUM_OK) then
      ierr = 1
      return
    end if
    nlast = size(solution%t)
    si%z(1:4) = solution%y(:, nlast)
    if (nlast > 1) explicit_h_carry = solution%h(nlast)
    ktau = ktau+1
  end do
  call sync_field_to_state
end subroutine orbit_timestep_radau15

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
! Re-evaluate the field at the converged end-of-step state.
!
! Only z(1:4) are integrated; the caller reconstructs the parallel velocity as
! z(5) = f%vpar/(pabs*sqrt(2)), reading f%vpar off the module-global field_can_t
! that f_ode fills as a side effect. After a multi-stage step that side effect
! belongs to the last interior stage node, not to the end of the step, so z(5)
! carries an error the tolerance cannot remove: measured on the sweep,
! components 1-4 of the final state converged to 5e-9 while z(5) stalled at
! 3.3e-4 and got no better from rtol 1e-8 to 1e-10.
!
! One evaluation per macro-step, against the hundreds a step already costs.
subroutine sync_field_to_state
  real(dp) :: zdot_discard(4)

  call f_ode(0d0, si%z(1:4), zdot_discard)
end subroutine sync_field_to_state

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
! Gragg-Bulirsch-Stoer extrapolation on the same canonical Hamiltonian.
subroutine orbit_timestep_gbs16(ierr)
  !
  use fortnum_ode_extrapolation, only : ode_integrate_gbs
  use fortnum_ode, only : ode_problem_t, ode_solution_t
  use fortnum_status, only : fortnum_status_t, FORTNUM_OK
  integer, intent(out) :: ierr
  integer :: ktau, nlast
  type(ode_problem_t) :: problem
  type(ode_solution_t) :: solution
  type(fortnum_status_t) :: status

  ierr = 0
  ktau = 0
  problem%rhs => f_ode_fortnum
  problem%rtol = si%rtol
  problem%atol = si%atol
  allocate(problem%y0(4))
  do while(ktau .lt. si%ntau)
    problem%t0 = ktau*si%dt
    problem%t1 = (ktau+1)*si%dt
    problem%y0 = si%z(1:4)
    problem%h0 = explicit_h_carry
    call ode_integrate_gbs(problem, solution, status)
    if (status%code /= FORTNUM_OK) then
      ierr = 1
      return
    end if
    nlast = size(solution%t)
    si%z(1:4) = solution%y(:, nlast)
    if (nlast > 1) explicit_h_carry = solution%h(nlast)
    ktau = ktau+1
  end do
  call sync_field_to_state
end subroutine orbit_timestep_gbs16

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
! Two-derivative Runge-Kutta on the canonical Hamiltonian: the Runge-Kutta-
! Nystrom transplant.
!
! RKN integrates y'' = f(t,y) and does not apply to the guiding-centre system as
! written, which is first order. Differentiating once gives zdd = F'(z)F(z) =
! G(z), which is in special second-order form because zd is itself determined by
! z, and re-synchronising the velocity as F(z) at the start of every step turns
! the scheme into a special two-derivative Runge-Kutta method (Chan & Tsai,
! Numer. Algorithms 53 (2010) 171). Its order conditions coincide with the
! Nystrom ones, which is what lets an RKN tableau be used directly here.
!
! Because the velocity is recomputed rather than propagated, only the position
! order conditions are needed, so order 4 costs two stages instead of three.
!
! Fixed step, like the symplectic schemes: fortnum's TDRK has no embedded error
! estimate yet, so the accuracy parameter is the step count rather than a
! tolerance. npoiper2 therefore sets resolution here exactly as it does for the
! symplectic integrators, which also makes the two directly comparable.
subroutine orbit_timestep_tdrk24(ierr)
  !
  use fortnum_ode_tdrk, only : tdrk_integrate_fixed
  use fortnum_status, only : fortnum_status_t, FORTNUM_OK
  integer, intent(out) :: ierr
  integer :: ktau, nfev_f, nfev_g
  real(dp) :: yend(4)
  type(fortnum_status_t) :: status

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    call tdrk_integrate_fixed(f_ode_fortnum, g_ode_fortnum, &
      ktau*si%dt, (ktau+1)*si%dt, si%z(1:4), 1, 4, .false., yend, &
      nfev_f, nfev_g, status)
    if (status%code /= FORTNUM_OK) then
      ierr = 1
      return
    end if
    si%z(1:4) = yend
    ktau = ktau+1
  end do
  call sync_field_to_state
end subroutine orbit_timestep_tdrk24

end module orbit_symplectic_quasi
