module orbit_symplectic

use util, only: pi, twopi
use field_can_mod, only: FieldCan, eval_field, get_val, get_derivatives, get_derivatives2

implicit none
save

public

integer, parameter :: NLAG_MAX = 2
integer, parameter :: NBUF_MAX = 16*NLAG_MAX

type :: SymplecticIntegrator
  integer :: nlag          ! Lagrange polynomial order
  integer :: nbuf          ! values to store back
  logical :: extrap_field  ! do extrapolation after final iteration

  double precision :: atol
  double precision :: rtol

! Current phase-space coordinates z and old pth
  double precision, dimension(4) :: z  ! z = (r, th, ph, pphi)
  double precision :: pthold

  ! Buffer for Lagrange polynomial interpolation
  integer :: kbuf
  integer :: kt
  integer :: k
  integer :: bufind(0:NLAG_MAX)
  double precision, dimension(8, 0:NBUF_MAX) :: zbuf
  double precision, dimension(0:0, NLAG_MAX+1) :: coef

  ! Timestep and variables from z0
  integer :: ntau
  double precision :: dt
  double precision :: pabs

  ! Integrator mode
  integer :: mode  ! 1 = euler1, 2 = euler2
end type SymplecticIntegrator

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Composition method with 2s internal stages according to Hairer, 2002 V.3.1
!
integer, parameter :: S_MAX = 32
type :: MultistageIntegrator
  integer :: s
  double precision :: alpha(S_MAX), beta(S_MAX)
  type(SymplecticIntegrator) stages(2*S_MAX)
end type MultistageIntegrator

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init(si, f, z, dt, ntau, rtol_init, mode_init, nlag)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  double precision, intent(in) :: z(:)
  double precision, intent(in) :: dt
  integer, intent(in) :: ntau
  double precision, intent(in) :: rtol_init
  integer, intent(in) :: mode_init
  integer, intent(in) :: nlag ! Lagrangian polynomials

  integer :: k

  si%mode = mode_init
  si%atol = 1d-15
  si%rtol = rtol_init

  si%kbuf = 0
  si%kt = 0
  si%k = 0

  si%bufind = 0

  si%ntau = ntau
  si%dt = dt

  si%z = z
  si%nlag = nlag
  si%nbuf = 4*nlag

  si%extrap_field = .True.

  call eval_field(f, z(1), z(2), z(3), 0)
  call get_val(f, si%z(4)) ! for pth
  si%pthold = f%pth
  if(si%nlag > 0) call plag_coeff(si%nlag+1, 0, 1d0, 1d0*(/(k,k=-si%nlag,0)/), si%coef)
end subroutine orbit_sympl_init


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_sympl_euler1(si, f, n, x, fvec, iflag)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(f, x(1), si%z(2), si%z(3), 2)
  call get_derivatives2(f, x(2))

  fvec(1) = f%dpth(1)*(f%pth - si%pthold) + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
  fvec(2) = f%dpth(1)*(x(2) - si%z(4))  + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))

end subroutine f_sympl_euler1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_sympl_euler1(si, f, x, jac)
!
  type(SymplecticIntegrator), intent(in) :: si
  type(FieldCan), intent(inout) :: f

  double precision, intent(in)  :: x(2)
  double precision, intent(out) :: jac(2, 2)

  jac(1,1) = f%d2pth(1)*(f%pth - si%pthold) + f%dpth(1)**2 &
    + si%dt*(f%d2H(2)*f%dpth(1) + f%dH(2)*f%d2pth(1) - f%d2H(1)*f%dpth(2) - f%dH(1)*f%d2pth(2))
  jac(1,2) = f%d2pth(7)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(4) &
    + si%dt*(f%d2H(8)*f%dpth(1) + f%dH(2)*f%d2pth(7) - f%d2H(7)*f%dpth(2) - f%dH(1)*f%d2pth(8))
  jac(2,1) = f%d2pth(1)*(x(2) - si%z(4)) &
    + si%dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
  jac(2,2) = f%d2pth(7)*(x(2) - si%z(4)) + f%dpth(1) &
    + si%dt*(f%d2H(9)*f%dpth(1) + f%dH(3)*f%d2pth(7) - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9))

end subroutine jac_sympl_euler1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_sympl_euler2(si, f, n, x, fvec, iflag)
!
  type(SymplecticIntegrator), intent(in) :: si
  type(FieldCan), intent(inout) :: f
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(f, x(1), x(2), x(3), 2)
  call get_derivatives2(f, si%z(4))

  fvec(1) = f%pth - si%pthold
  fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
  fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)

end subroutine f_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_sympl_euler2(si, f, x, jac)
!
  type(SymplecticIntegrator), intent(in) :: si
  type(FieldCan), intent(inout) :: f
  double precision, intent(in)  :: x(3)
  double precision, intent(out) :: jac(3, 3)

  jac(1,:) = f%dpth(1:3)
  jac(2,:) = f%d2pth(1:3)*(x(2) - si%z(2)) - si%dt*f%d2H(1:3)
  jac(2,2) = jac(2,2) + f%dpth(1)
  jac(3,:) = (f%d2pth(1:3)*f%hph + f%dpth(1)*f%dhph)*(x(3) - si%z(3)) &
    - si%dt*(f%d2pth(1:3)*f%vpar + f%dpth(1)*f%dvpar(1:3) - f%d2H(1:3)*f%hth - f%dH(1)*f%dhth)
  jac(3,3) = jac(3,3) + f%dpth(1)*f%hph

end subroutine jac_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_midpoint_part1(si, f, n, x, fvec)
  !
    type(SymplecticIntegrator), intent(in) :: si
    type(FieldCan), intent(inout) :: f
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    double precision, intent(out) :: fvec(n)

    ! evaluate at midpoint
    call eval_field(f, x(5), 0.5d0*(x(2) + si%z(2)), 0.5d0*(x(3) + si%z(3)), 2)
    call get_derivatives2(f, 0.5d0*(x(4) + si%z(4)))

    fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
    fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)
    fvec(4) = f%dpth(1)*(x(4) - si%z(4)) + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
    fvec(5) = f%dpth(1)*(f%pth - si%pthold) + 0.5d0*si%dt*(f%dpth(1)*f%dH(2)-f%dpth(2)*f%dH(1))

  end subroutine f_midpoint_part1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_midpoint_part2(si, f, n, x, fvec)
  !
    type(SymplecticIntegrator), intent(in) :: si
    type(FieldCan), intent(inout) :: f
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    double precision, intent(out) :: fvec(n)

    double precision :: dpthmid, pthdotbar

    ! save evaluation from midpoint
    dpthmid = f%dpth(1)
    pthdotbar = f%dpth(1)*f%dH(2) - f%dpth(2)*f%dH(1)

    ! evaluate at endpoint
    call eval_field(f, x(1), x(2), x(3), 2)
    call get_derivatives2(f, x(4))
    fvec(1) = dpthmid*(f%pth - si%pthold) + si%dt*pthdotbar

  end subroutine f_midpoint_part2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_midpoint_part1(si, f, x, jac)
  !
    type(SymplecticIntegrator), intent(in) :: si
    type(FieldCan), intent(inout) :: f
    double precision, intent(in)  :: x(5)
    double precision, intent(out) :: jac(5, 5)

    jac(2,1) = 0d0
    jac(2,5) = f%d2pth(1)*(x(2) - si%z(2)) - si%dt*f%d2H(1)
    jac(2,2:3) = 0.5d0*(f%d2pth(2:3)*(x(2) - si%z(2)) - si%dt*f%d2H(2:3))
    jac(2,2) = jac(2,2) + f%dpth(1)
    jac(2,4) = 0.5d0*(f%d2pth(7)*(x(2) - si%z(2)) - si%dt*f%d2H(7))

    jac(3,1) = 0d0
    jac(3,5) = (f%d2pth(1)*f%hph + f%dpth(1)*f%dhph(1))*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(1)*f%vpar + f%dpth(1)*f%dvpar(1) &
      - f%d2H(1)*f%hth - f%dH(1)*f%dhth(1))
    jac(3,2:3) = 0.5d0*((f%d2pth(2:3)*f%hph + f%dpth(1)*f%dhph(2:3))*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(2:3)*f%vpar + f%dpth(1)*f%dvpar(2:3) &
      - f%d2H(2:3)*f%hth - f%dH(1)*f%dhth(2:3)))
    jac(3,3) = jac(3,3) + f%dpth(1)*f%hph
    jac(3,4) = 0.5d0*(f%d2pth(7)*f%hph*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(7)*f%vpar + f%dpth(1)*f%dvpar(4) - f%d2H(7)*f%hth))

    jac(4,1) = 0d0
    jac(4,5) = f%d2pth(1)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) &
      - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
    jac(4,2) = 0.5d0*(f%d2pth(2)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(5)*f%dpth(1) + f%dH(3)*f%d2pth(2) &
      - f%d2H(2)*f%dpth(3) - f%dH(1)*f%d2pth(5)))
    jac(4,3) = 0.5d0*(f%d2pth(3)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(6)*f%dpth(1) + f%dH(3)*f%d2pth(3) &
      - f%d2H(3)*f%dpth(3) - f%dH(1)*f%d2pth(6)))
    jac(4,4) = 0.5d0*(f%d2pth(7)*(x(4) - si%z(4)) &
      + si%dt*(f%dH(3)*f%d2pth(7) + f%d2H(9)*f%dpth(1)&
      - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9)))
    jac(4,4) = jac(4,4) + f%dpth(1)

    jac(5,1) = 0d0
    jac(5,5) = f%d2pth(1)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(1)&
      + si%dt/2.0d0*(f%d2pth(1)*f%dH(2) + f%dpth(1)*f%d2H(2) &
      - f%d2pth(2)*f%dH(1) - f%dpth(2)*f%d2H(1))
    jac(5,2) = 0.5d0*(f%d2pth(2)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(2) &
      + si%dt/2.0d0*(f%d2pth(2)*f%dH(2) + f%dpth(1)*f%d2H(4) &
      - f%d2pth(4)*f%dH(1) - f%dpth(2)*f%d2H(2)))
    jac(5,3) = 0.5d0*(f%d2pth(3)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(3) &
      + si%dt/2.0d0*(f%d2pth(3)*f%dH(2) + f%dpth(1)*f%d2H(5) &
      - f%d2pth(5)*f%dH(1) - f%dpth(2)*f%d2H(3)))
    jac(5,4) = 0.5d0*(f%d2pth(7)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(4) &
      + si%dt/2.0d0*(f%d2pth(7)*f%dH(2) + f%dpth(1)*f%d2H(8) &
      - f%d2pth(8)*f%dH(1) - f%dpth(2)*f%d2H(7)))

end subroutine jac_midpoint_part1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_midpoint_part2(si, f, fmid, x, jac)
  !
    type(SymplecticIntegrator), intent(in) :: si
    type(FieldCan), intent(inout) :: f
    type(FieldCan), intent(inout) :: fmid
    double precision, intent(in)  :: x(5)
    double precision, intent(out) :: jac(5, 5)

    ! fmid%dpth(1)*(f%pth - si%pthold) + si%dt*(fmid%dpth(1)*fmid%dH(2)-fmid%dpth(2)*fmid%dH(1))

    jac(1,1) = fmid%dpth(1)*f%dpth(1)
    jac(1,2) = 0.5d0*(fmid%d2pth(2)*(f%pth - si%pthold) &
      + si%dt*(fmid%d2pth(2)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(4) &
      - fmid%dpth(2)*fmid%d2H(2) - fmid%d2pth(4)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(2)
    jac(1,3) = 0.5d0*(fmid%d2pth(3)*(f%pth - si%pthold) &
      + si%dt*(fmid%d2pth(3)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(5) &
      - fmid%dpth(2)*fmid%d2H(3) - fmid%d2pth(5)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(3)
    jac(1,4) = 0.5d0*(fmid%d2pth(7)*(f%pth - si%pthold) &
        + si%dt*(fmid%d2pth(7)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(8) &
        - fmid%dpth(2)*fmid%d2H(7) - fmid%d2pth(8)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(4)
    jac(1,5) = fmid%d2pth(1)*(f%pth - si%pthold) &
        + si%dt*(fmid%d2pth(1)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(2) &
        - fmid%dpth(2)*fmid%d2H(1) - fmid%d2pth(2)*fmid%dH(1))

end subroutine jac_midpoint_part2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine newton1(si, f, x, maxit, xlast)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  integer, parameter :: n = 2

  double precision, intent(inout) :: x(n)
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(n)

  double precision :: fvec(n), fjac(n,n), ijac(n,n)
  double precision :: tolref(n)
  integer :: kit

  tolref(1) = 1d0
  tolref(2) = dabs(x(2))

  do kit = 1, maxit
    if(x(1) > 1d0) return
    if(x(1) < 0d0) x(1) = 0.01d0

    call f_sympl_euler1(si, f, n, x, fvec, 1)
    call jac_sympl_euler1(si, f, x, fjac)
    ijac(1,1) = 1d0/(fjac(1,1) - fjac(1,2)*fjac(2,1)/fjac(2,2))
    ijac(1,2) = -1d0/(fjac(1,1)*fjac(2,2)/fjac(1,2) - fjac(2,1))
    ijac(2,1) = -1d0/(fjac(1,1)*fjac(2,2)/fjac(2,1) - fjac(1,2))
    ijac(2,2) = 1d0/(fjac(2,2) - fjac(1,2)*fjac(2,1)/fjac(1,1))
    xlast = x
    x = x - matmul(ijac, fvec)

    ! Don't take too small values in pphi as tolerance reference
    tolref(2) = max(dabs(x(2)), tolref(2))
    tolref(2) = max(dabs(x(2)), tolref(2))

    if (all(dabs(fvec) < si%atol)) return
    if (all(dabs(x-xlast) < si%rtol*tolref)) return
  enddo
  print *, 'newton1: maximum iterations reached: ', maxit
  write(6601,*) x(1), x(2)
  write(6601,*) x-xlast
  write(6601,*) fvec
  write(6601,*) ''
  write(6601,*) fjac(1,1), fjac(1,2)
  write(6601,*) fjac(2,1), fjac(2,2)
  write(6601,*) ''
  write(6601,*) ijac(1,1), ijac(1,2)
  write(6601,*) ijac(2,1), ijac(2,2)
  write(6601,*) ''
  write(6601,*) si%z(2), si%z(3)
  write(6601,*) ''
  write(6601,*) ''
end subroutine

subroutine newton2(si, f, x, atol, rtol, maxit, xlast)
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f

  integer, parameter :: n = 3
  integer :: kit

  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(n)

  double precision :: fvec(n), fjac(n,n), jinv(n,n)
  integer :: pivot(n), info

  double precision :: xabs(n), tolref(n), fabs(n)
  double precision :: det

  do kit = 1, maxit
    if(x(1) > 1.0) return
    if(x(1) < 0.0) x(1) = 0.01
    call f_sympl_euler2(si, f, n, x, fvec, 1)
    fabs = dabs(fvec)
    call jac_sympl_euler2(si, f, x, fjac)
    xlast = x

    det = fjac(1,1)*fjac(2,2)*fjac(3,3) + fjac(1,2)*fjac(2,3)*fjac(3,1) + fjac(1,3)*fjac(2,1)*fjac(3,2) &
        - fjac(1,3)*fjac(2,2)*fjac(3,1) - fjac(1,1)*fjac(2,3)*fjac(3,2) - fjac(1,2)*fjac(2,1)*fjac(3,3)

    jinv(1,1) = fjac(2,2)*fjac(3,3) - fjac(2,3)*fjac(3,2)
    jinv(1,2) = fjac(1,3)*fjac(3,2) - fjac(1,2)*fjac(3,3)
    jinv(1,3) = fjac(1,2)*fjac(2,3) - fjac(1,3)*fjac(2,2)

    jinv(2,1) = fjac(2,3)*fjac(3,1) - fjac(2,1)*fjac(3,3)
    jinv(2,2) = fjac(1,1)*fjac(3,3) - fjac(1,3)*fjac(3,1)
    jinv(2,3) = fjac(1,3)*fjac(2,1) - fjac(1,1)*fjac(2,3)

    jinv(3,1) = fjac(2,1)*fjac(3,2) - fjac(3,1)*fjac(2,2)
    jinv(3,2) = fjac(1,2)*fjac(3,1) - fjac(3,2)*fjac(1,1)
    jinv(3,3) = fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1)

    x = x - matmul(jinv, fvec)/det

    !call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    ! after solution: fvec = (xold-xnew)_Newton
    !x = x - fvec

    xabs = dabs(x-xlast)

    tolref(1) = 1d0
    tolref(2) = twopi
    tolref(3) = twopi

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton2: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6602,*) x(1), x(2), x(3), si%z(4), xabs
end subroutine

subroutine newton_midpoint(si, f, x, atol, rtol, maxit, xlast)
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f

  type(FieldCan) :: fmid

  integer, parameter :: n = 5
  integer :: kit

  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(n)

  double precision :: fvec(n), fjac(n,n)
  integer :: pivot(n), info

  double precision :: xabs(n), tolref(n), fabs(n)

  do kit = 1, maxit
    if(x(1) > 1.0) return
    if(x(1) < 0.0) x(1) = 0.01
    if(x(5) < 0.0) x(5) = 0.01
    call f_midpoint_part1(si, f, n, x, fvec)
    call jac_midpoint_part1(si, f, x, fjac)
    fmid = f
    call f_midpoint_part2(si, f, n, x, fvec)
    call jac_midpoint_part2(si, f, fmid, x, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)

    tolref = dabs(xlast)
    tolref(2) = twopi
    tolref(3) = twopi

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton_midpoint: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6603,*) x(1), x(2), x(3), x(4), x(5), xabs, fvec
end subroutine

subroutine coeff_rk_gauss(n, a, b, c)
  integer, intent(in) :: n
  double precision, intent(inout) :: a(n,n), b(n), c(n)

  if (n == 1) then
    a(1,1) = 0.5d0
    b(1) = 1.0d0
    c(1) = 0.5d0
  elseif (n == 2) then
    a(1,1) =  0.25d0
    a(1,2) = -0.038675134594812d0
    a(2,1) =  0.538675134594812d0
    a(2,2) =  0.25d0

    b(1) = 0.5d0
    b(2) = 0.5d0

    c(1) = 0.211324865405187d0
    c(2) = 0.788675134594812d0
  elseif (n == 3) then
    a(1,1) =  0.1388888888888889d0
    a(1,2) = -0.03597666752493894d0
    a(1,3) =  0.009789444015308318d0
    a(2,1) =  0.3002631949808646d0
    a(2,2) =  0.2222222222222222d0
    a(2,3) = -0.022485417203086805d0
    a(3,1) = 0.26798833376246944d0
    a(3,2) = 0.48042111196938336d0
    a(3,3) = 0.1388888888888889d0

    b(1) = 0.2777777777777778d0
    b(2) = 0.4444444444444444d0
    b(3) = 0.2777777777777778d0

    c(1) = 0.1127016653792583d0
    c(2) = 0.5d0
    c(3) = 0.8872983346207417d0
  elseif (n == 4) then  ! with help of coefficients from GeometricIntegrators.jl of Michael Kraus
    a(1,1) = 0.086963711284363462428182d0
    a(1,2) = -0.026604180084998794303397d0
    a(1,3) = 0.012627462689404725035280d0
    a(1,4) = -0.003555149685795683332096d0

    a(2,1) = 0.188118117499868064967927d0
    a(2,2) = 0.163036288715636523694030d0
    a(2,3) = -0.027880428602470894855481d0
    a(2,4) = 0.006735500594538155853808d0

    a(3,1) = 0.167191921974188778543535d0
    a(3,2) = 0.353953006033743966529670d0
    a(3,3) = 0.163036288715636523694030d0
    a(3,4) = -0.014190694931141143581010d0

    a(4,1) = 0.177482572254522602550608d0
    a(4,2) = 0.313445114741868369190314d0
    a(4,3) = 0.352676757516271865977586d0
    a(4,4) = 0.086963711284363462428182d0

    b(1) = 0.173927422568726924856364d0
    b(2) = 0.326072577431273047388061d0
    b(3) = 0.326072577431273047388061d0
    b(4) = 0.173927422568726924856364d0

    c(1) = 0.069431844202973713731097d0
    c(2) = 0.330009478207571871344328d0
    c(3) = 0.669990521792428128655672d0
    c(4) = 0.930568155797026341780054d0
  else
    ! not implemented
    a = 0d0
    b = 0d0
    c = 0d0
  endif
end subroutine coeff_rk_gauss


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Gauss-Legendre Runge-Kutta method with s internal stages (n=4*s variables)
!
subroutine f_rk_gauss(si, fs, s, x, fvec)
  !
  type(SymplecticIntegrator), intent(in) :: si
  integer, intent(in) :: s
  type(FieldCan), intent(inout) :: fs(:)
  double precision, intent(in) :: x(4*s)  ! = (rend, thend, phend, pphend)
  double precision, intent(out) :: fvec(4*s)

  double precision :: a(s,s), b(s), c(s), Hprime(s)
  integer :: k,l  ! counters

  call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

  ! evaluate stages
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

  end subroutine f_rk_gauss



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_rk_gauss(si, fs, s, jac)
!
  type(SymplecticIntegrator), intent(in) :: si
  integer, intent(in) :: s
  type(FieldCan), intent(in) :: fs(:)
  double precision, intent(out) :: jac(4*s, 4*s)

  double precision :: a(s,s), b(s), c(s), Hprime(s), dHprime(4*s)
  integer :: k,l,m  ! counters

  call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

  ! evaluate stages
  do k = 1, s
    m=4*k
    Hprime(k)      = fs(k)%dH(1)/fs(k)%dpth(1)
    dHprime(m-3) = (fs(k)%d2H(1) - Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)  ! d/dr
    dHprime(m-2) = (fs(k)%d2H(2) - Hprime(k)*fs(k)%d2pth(2))/fs(k)%dpth(1)  ! d/dth
    dHprime(m-1) = (fs(k)%d2H(3) - Hprime(k)*fs(k)%d2pth(3))/fs(k)%dpth(1)  ! d/dph
    dHprime(m)   = (fs(k)%d2H(7) - Hprime(k)*fs(k)%d2pth(7))/fs(k)%dpth(1)  ! d/dpph
  end do

  jac = 0d0

  do k = 1, s
    m=4*k
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
        m=4*k
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

end subroutine jac_rk_gauss

subroutine newton_rk_gauss(si, fs, s, x, atol, rtol, maxit, xlast)
  type(SymplecticIntegrator), intent(inout) :: si

  integer, intent(in) :: s
  type(FieldCan), intent(inout) :: fs(:)
  integer :: kit, ks

  double precision, intent(inout) :: x(4*s)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(4*s)

  double precision :: fvec(4*s), fjac(4*s, 4*s)
  integer :: pivot(4*s), info

  double precision :: xabs(4*s), tolref(4*s), fabs(4*s)

  do kit = 1, maxit

    ! Check if radius left the boundary
    do ks = 1, s
      if (x(4*ks-3) > 1d0) return
      if (x(4*ks-3) < 0.0) x(4*ks-3) = 0.01d0
    end do

    call f_rk_gauss(si, fs, s, x, fvec)
    call jac_rk_gauss(si, fs, s, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(4*s, 1, fjac, 4*s, pivot, fvec, 4*s, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)

    do ks = 1, s
      tolref(4*ks-3) = 1d0
      tolref(4*ks-2) = twopi
      tolref(4*ks-1) = twopi
      tolref(4*ks) = dabs(xlast(4*ks))
    end do

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton_rk_gauss: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6603,*) x, xabs, fvec
end subroutine newton_rk_gauss


subroutine fixpoint_rk_gauss(si, fs, s, x, atol, rtol, maxit, xlast)
  ! TODO: this doesn't work well yet
  type(SymplecticIntegrator), intent(inout) :: si

  integer, intent(in) :: s
  type(FieldCan), intent(inout) :: fs(:)
  integer :: kit, ks

  double precision, intent(inout) :: x(4*s)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(4*s)

  double precision :: fvec(4*s)

  double precision :: xabs(4*s), tolref(4*s), fabs(4*s)

  double precision :: a(s,s), b(s), c(s), Hprime(s), dHprimedr(s)
  integer :: k, l

  double precision :: pthnew, dpthnewdr, damp

  call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

  do kit = 1, 4096

    ! Check if radius left the boundary
    do ks = 1, s
      if (x(4*ks-3) > 1d0) return
      if (x(4*ks-3) < 0.0) x(4*ks-3) = 0.01d0
    end do

    call f_rk_gauss(si, fs, s, x, fvec)
    fabs = dabs(fvec)
    xlast = x

    do k = 1, s
      call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 2)
      call get_derivatives2(fs(k), x(4*k))
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
      dHprimedr(k) = (fs(k)%d2H(1) - Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)
    end do

    k = 1
    l = 1

    pthnew = si%pthold - si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
    dpthnewdr = -si%dt*a(k,l)*( &
    fs(l)%d2H(2) - Hprime(l)*fs(l)%d2pth(2) - dHprimedr(l)*fs(l)%dpth(2))
    !print *, pthnew, dpthnewdr, pthnew/dpthnewdr

    damp = 0.99d0

    x(1) = x(1) - (1d0-damp)*pthnew/dpthnewdr
    x(2) = damp*x(2) + (1d0-damp)*(si%z(2) + si%dt*a(k,l)*si%dt*a(k,l)*Hprime(l))
    x(3) = damp*x(3) + (1d0-damp)*(si%z(3) + si%dt*a(k,l)*(fs(l)%vpar &
                                    - Hprime(l)*fs(l)%hth)/fs(l)%hph)
    x(4) = damp*x(4) + (1d0-damp)*(si%z(4) - si%dt*a(k,l)*(fs(l)%dH(3) &
                                    - Hprime(l)*fs(l)%dpth(3)))

    print *, x

    xabs = dabs(x - xlast)

    do ks = 1, s
      xabs(4*ks-2) = modulo(xabs(4*ks-2), pi)
      xabs(4*ks-1) = modulo(xabs(4*ks-1), pi)
    end do

    do ks = 1, s
      tolref(4*ks-3) = 1d0
      tolref(4*ks-2) = twopi
      tolref(4*ks-1) = twopi
      tolref(4*ks) = dabs(x(4*ks))
    end do

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'fixpoint_rk_gauss: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6603,*) x, xabs, fvec
end subroutine fixpoint_rk_gauss


subroutine coeff_rk_lobatto(n, a, ahat, b, c)
  integer, intent(in) :: n
  double precision, intent(inout) :: a(n,n), ahat(n,n), b(n), c(n)

  if (n == 3) then
    a(1,1) =  0d0
    a(1,2) =  0d0
    a(1,3) =  0d0

    a(2,1) =  0.20833333333333334d0
    a(2,2) =  0.33333333333333333d0
    a(2,3) = -0.041666666666666664d0

    a(3,1) =  0.16666666666666667d0
    a(3,2) =  0.66666666666666667d0
    a(3,3) =  0.16666666666666667d0

    ahat(1,1) =  0.16666666666666667d0
    ahat(1,2) = -0.16666666666666667d0
    ahat(1,3) =  0d0

    ahat(2,1) =  0.16666666666666667d0
    ahat(2,2) =  0.33333333333333333d0
    ahat(2,3) =  0d0

    ahat(3,1) =  0.16666666666666667d0
    ahat(3,2) =  0.83333333333333333d0
    ahat(3,3) =  0d0

    b(1) =  0.16666666666666667d0
    b(2) =  0.66666666666666667d0
    b(3) =  0.16666666666666667d0

    c(1) = 0d0
    c(2) = 0.5d0
    c(3) = 1.0d0

  else
    ! not implemented
    a = 0d0
    b = 0d0
    c = 0d0
  endif
end subroutine coeff_rk_lobatto

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Lobatto (IIIA)-(IIIB) Runge-Kutta method with s internal stages (n=4*s variables)
!


subroutine f_rk_lobatto(si, fs, s, x, fvec, jactype)
  !
  type(SymplecticIntegrator), intent(in) :: si
  integer, intent(in) :: s
  type(FieldCan), intent(inout) :: fs(:)
  double precision, intent(in) :: x(4*s)  ! = (rend, thend, phend, pphend)
  double precision, intent(out) :: fvec(4*s)
  integer, intent(in) :: jactype  ! 0 = no second derivatives, 2 = second derivatives

  double precision :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s)
  integer :: k,l  ! counters

  call coeff_rk_lobatto(s, a, ahat, b, c)

  call eval_field(fs(1), x(1), si%z(2), si%z(3), jactype)
  call get_derivatives(fs(1), x(2))

  do k = 2, s
    call eval_field(fs(k), x(4*k-3-2), x(4*k-2-2), x(4*k-1-2), jactype)
    call get_derivatives(fs(k), x(4*k-2))
  end do

  Hprime = fs%dH(1)/fs%dpth(1)

  fvec(1) = fs(1)%pth - si%pthold
  fvec(2) = x(2) - si%z(4)

  do l = 1, s
      fvec(1) = fvec(1) + si%dt*ahat(1,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
      fvec(2) = fvec(2) + si%dt*ahat(1,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
  end do

  do k = 2, s
    fvec(4*k-3-2) = fs(k)%pth - si%pthold
    fvec(4*k-2-2) = x(4*k-2-2)  - si%z(2)
    fvec(4*k-1-2) = x(4*k-1-2)  - si%z(3)
    fvec(4*k-2)   = x(4*k-2)    - si%z(4)
  end do

  do l = 1, s
    do k = 2, s
      fvec(4*k-3-2) = fvec(4*k-3-2) + si%dt*ahat(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
      fvec(4*k-2-2) = fvec(4*k-2-2) - si%dt*a(k,l)*Hprime(l)                                         ! thdot
      fvec(4*k-1-2) = fvec(4*k-1-2) - si%dt*a(k,l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph     ! phdot
      fvec(4*k-2)   = fvec(4*k-2)   + si%dt*ahat(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
    end do
  end do

end subroutine f_rk_lobatto


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_rk_lobatto(si, fs, s, jac)
!
  type(SymplecticIntegrator), intent(in) :: si
  integer, intent(in) :: s
  type(FieldCan), intent(in) :: fs(:)
  double precision, intent(out) :: jac(4*s-2, 4*s-2)

  double precision :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s), dHprime(4*s-2)
  integer :: k,l,m,n  ! counters

  call coeff_rk_lobatto(s, a, ahat, b, c)
  jac = 0d0

  Hprime = fs%dH(1)/fs%dpth(1)
  dHprime(1) = (fs(1)%d2H(1)-Hprime(1)*fs(1)%d2pth(1))/fs(1)%dpth(1)  ! d/dr
  dHprime(2) = (fs(1)%d2H(7)-Hprime(1)*fs(1)%d2pth(7))/fs(1)%dpth(1)  ! d/dpph
  do k = 2, s
    m = 4*k-2
    dHprime(m-3)=(fs(k)%d2H(1)-Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)  ! d/dr
    dHprime(m-2)=(fs(k)%d2H(2)-Hprime(k)*fs(k)%d2pth(2))/fs(k)%dpth(1)  ! d/dth
    dHprime(m-1)=(fs(k)%d2H(3)-Hprime(k)*fs(k)%d2pth(3))/fs(k)%dpth(1)  ! d/dph
    dHprime(m)  =(fs(k)%d2H(7)-Hprime(k)*fs(k)%d2pth(7))/fs(k)%dpth(1)  ! d/dpph
  end do

  ! evaluate first stage with only r and pph as unknowns
  jac(1, 1) = fs(1)%dpth(1)
  jac(1, 2) = fs(1)%dpth(4)
  jac(2, 2) = 1d0

  jac(1, 1) = jac(1, 1) + si%dt*ahat(1,1) * &
    (fs(1)%d2H(2) - fs(1)%d2pth(2)*Hprime(1) - fs(1)%dpth(2)*dHprime(1))
  jac(1, 2) = jac(1, 2) + si%dt*ahat(1,1)* &
    (fs(1)%d2H(8) - fs(1)%d2pth(8)*Hprime(1) - fs(1)%dpth(2)*dHprime(2))
  jac(2, 1) = jac(2, 1) + si%dt*ahat(1,1) * &
    (fs(1)%d2H(3) - fs(1)%d2pth(3)*Hprime(1) - fs(1)%dpth(3)*dHprime(1))
  jac(2, 2) = jac(2, 2) + si%dt*ahat(1,1) * &
    (fs(1)%d2H(9) - fs(1)%d2pth(9)*Hprime(1) - fs(1)%dpth(3)*dHprime(2))

  do l = 2, s
    n = 4*l-2
    jac(1, n-3) = jac(1, n-3) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-3))
    jac(1, n-2) = jac(1, n-2) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(4) - fs(l)%d2pth(4)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-2))
    jac(1, n-1) = jac(1, n-1) + si%dt*ahat(1,l)* &
      (fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-1))
    jac(1, n) = jac(1, n) + si%dt*ahat(1,l)* &
      (fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(n))

    jac(2, n-3) = jac(2, n-3) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-3))
    jac(2, n-2) = jac(2, n-2) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-2))
    jac(2, n-1) = jac(2, n-1) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(6) - fs(l)%d2pth(6)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-1))
    jac(2, n) = jac(2, n) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(n))
  end do

  ! evaluate remaining stages with r, th, ph, pph as unknowns

  do k = 2, s
    m = 4*k-2
    jac(m-3, m-3) = fs(k)%dpth(1)
    jac(m-3, m-2) = fs(k)%dpth(2)
    jac(m-3, m-1) = fs(k)%dpth(3)
    jac(m-3, m)   = fs(k)%dpth(4)
    jac(m-2, m-2) = 1d0
    jac(m-1, m-1) = 1d0
    jac(m, m) = 1d0
  end do

  l = 1
  do k = 2, s
    m = 4*k-2
    jac(m-3, 1) = jac(m-3, 1) + si%dt*ahat(k,l) * &              ! d/dr
      (fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(1))
    jac(m-3, 2) = jac(m-3, 2) + si%dt*ahat(k,l)* &                   ! d/dpph
      (fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(2))

    jac(m-2, 1) = jac(m-2, 1) - si%dt*a(k,l)*dHprime(1)  ! d/dr
    jac(m-2, 2) = jac(m-2, 2) - si%dt*a(k,l)*dHprime(2)  ! d/dpph

    jac(m-1, 1) = jac(m-1, 1) - si%dt*a(k,l) * & ! d/dr
      (-fs(l)%dhph(1)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
      +(fs(l)%dvpar(1)-dHprime(1)*fs(l)%hth-Hprime(l)*fs(l)%dhth(1))/fs(l)%hph)
    jac(m-1, 2) = jac(m-1, 2) - si%dt*a(k,l) * & ! d/dpph
      (fs(l)%dvpar(4)-dHprime(2)*fs(l)%hth)/fs(l)%hph

    jac(m, 1) = jac(m, 1) + si%dt*ahat(k,l) * & ! d/dr
      (fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(1))
    jac(m, 2) = jac(m, 2) + si%dt*ahat(k,l) * & ! d/dpph
      (fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(2))
  end do

  do l = 2, s
    do k = 2, s
      m = 4*k-2
      n = 4*l-2
      jac(m-3, n-3) = jac(m-3, n-3) + si%dt*ahat(k,l) * &              ! d/dr
        (fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-3))
      jac(m-3, n-2) = jac(m-3, n-2) + si%dt*ahat(k,l) * &              ! d/dth
        (fs(l)%d2H(4) - fs(l)%d2pth(4)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-2))
      jac(m-3, n-1) = jac(m-3, n-1) + si%dt*ahat(k,l)* &               ! d/dph
        (fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-1))
      jac(m-3, n) = jac(m-3, n) + si%dt*ahat(k,l)* &                   ! d/dpph
        (fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(n))

      jac(m-2, n-3) = jac(m-2, n-3) - si%dt*a(k,l)*dHprime(n-3)
      jac(m-2, n-2) = jac(m-2, n-2) - si%dt*a(k,l)*dHprime(n-2)
      jac(m-2, n-1) = jac(m-2, n-1) - si%dt*a(k,l)*dHprime(n-1)
      jac(m-2, n)   = jac(m-2, n)   - si%dt*a(k,l)*dHprime(n)

      jac(m-1, n-3) = jac(m-1, n-3) - si%dt*a(k,l) * & ! d/dr
        (-fs(l)%dhph(1)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
        +(fs(l)%dvpar(1)-dHprime(n-3)*fs(l)%hth-Hprime(l)*fs(l)%dhth(1))/fs(l)%hph)
      jac(m-1, n-2) = jac(m-1, n-2) - si%dt*a(k,l) * & ! d/dth
        (-fs(l)%dhph(2)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
        +(fs(l)%dvpar(2)-dHprime(n-2)*fs(l)%hth-Hprime(l)*fs(l)%dhth(2))/fs(l)%hph)
      jac(m-1, n-1) = jac(m-1, n-1) - si%dt*a(k,l) * & ! d/dph
        (-fs(l)%dhph(3)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
        +(fs(l)%dvpar(3)-dHprime(n-1)*fs(l)%hth-Hprime(l)*fs(l)%dhth(3))/fs(l)%hph)
      jac(m-1, n) = jac(m-1, n) - si%dt*a(k,l) * & ! d/dpph
        (fs(l)%dvpar(4)-dHprime(n)*fs(l)%hth)/fs(l)%hph

      jac(m, n-3) = jac(m, n-3) + si%dt*ahat(k,l) * & ! d/dr
        (fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-3))
      jac(m, n-2) = jac(m, n-2) + si%dt*ahat(k,l) * & ! d/dth
        (fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-2))
      jac(m, n-1) = jac(m, n-1) + si%dt*ahat(k,l) * & ! d/dph
        (fs(l)%d2H(6) - fs(l)%d2pth(6)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-1))
      jac(m, n) = jac(m, n) + si%dt*ahat(k,l) * & ! d/dpph
        (fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(n))
    end do
  end do
end subroutine jac_rk_lobatto


subroutine newton_rk_lobatto(si, fs, s, x, atol, rtol, maxit, xlast)
  type(SymplecticIntegrator), intent(inout) :: si

  integer, intent(in) :: s
  type(FieldCan), intent(inout) :: fs(:)
  integer :: kit, ks

  double precision, intent(inout) :: x(4*s-2)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(4*s-2)

  double precision :: fvec(4*s-2), fjac(4*s-2, 4*s-2)
  integer :: pivot(4*s-2), info

  double precision :: xabs(4*s-2), tolref(4*s-2), fabs(4*s-2)

  do kit = 1, maxit

    ! Check if radius left the boundary
    if (x(1) > 1d0) return
    if (x(1) < 0.0) x(1) = 0.01d0
    do ks = 2, s
      if (x(4*ks-2-3) > 1d0) return
      if (x(4*ks-2-3) < 0.0) x(4*ks-3) = 0.01d0
    end do

    call f_rk_lobatto(si, fs, s, x, fvec, 2)
    call jac_rk_lobatto(si, fs, s, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(4*s-2, 1, fjac, 4*s-2, pivot, fvec, 4*s-2, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)

    tolref(1) = 1d0
    tolref(2) = dabs(xlast(2))
    do ks = 2, s
      tolref(4*ks-2-3) = 1d0
      tolref(4*ks-2-2) = twopi
      tolref(4*ks-2-1) = twopi
      tolref(4*ks-2) = dabs(xlast(4*ks-2))
    end do

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton_rk_lobatto: maximum iterations reached: ', &
           maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6603,*) x, xabs, fvec
end subroutine newton_rk_lobatto


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl(si, f, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f

  integer, intent(out) :: ierr

  select case (si%mode)
   case (1)
      call orbit_timestep_sympl_euler1(si, f, ierr)
   case (2)
      call orbit_timestep_sympl_euler2(si, f, ierr)
   case (3)
      call orbit_timestep_sympl_midpoint(si, f, ierr)
   case (4)
      call orbit_timestep_sympl_rk_gauss(si, f, 1, ierr)
   case (5)
      call orbit_timestep_sympl_rk_gauss(si, f, 2, ierr)
   case (6)
      call orbit_timestep_sympl_rk_gauss(si, f, 3, ierr)
   case (7)
      call orbit_timestep_sympl_rk_gauss(si, f, 4, ierr)
   case (15)
      call orbit_timestep_sympl_rk_lobatto(si, f, 3, ierr)
   case default
      print *, 'invalid mode for orbit_timestep_sympl: ', si%mode
      stop
  end select

end subroutine orbit_timestep_sympl


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_multi(mi, f, ierr)
!
  type(MultistageIntegrator), intent(inout) :: mi
  type(FieldCan), intent(inout) :: f

  integer, intent(out) :: ierr

  integer :: kstage

  call orbit_timestep_sympl(mi%stages(1), f, ierr)

  do kstage = 2, 2*mi%s
    mi%stages(kstage)%z = mi%stages(kstage-1)%z
    mi%stages(kstage)%pthold = f%pth
    call orbit_timestep_sympl(mi%stages(kstage), f, ierr)
  end do

  mi%stages(1)%z = mi%stages(2*mi%s)%z
  mi%stages(1)%pthold = f%pth

end subroutine orbit_timestep_sympl_multi


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
!
  type(MultistageIntegrator), intent(inout) :: mi
  type(FieldCan), intent(inout) :: f

  double precision, intent(in) :: z(:)
  double precision, intent(in) :: dtau
  integer, intent(in) :: ntau
  double precision, intent(in) :: rtol_init

  double precision, intent(in) :: alpha(:), beta(:)

  integer :: ks

  mi%s = size(alpha)

  do ks = 1, mi%s
    call orbit_sympl_init(mi%stages(2*ks-1), f, z, &
      alpha(ks)*dtau, ntau, rtol_init, 1, 1)
    call orbit_sympl_init(mi%stages(2*ks), f, z, &
      beta(ks)*dtau, ntau, rtol_init, 2, 1)
  end do
end subroutine orbit_sympl_init_multi


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init_verlet(mi, f, z, dtau, ntau, rtol_init)
!
  type(MultistageIntegrator), intent(inout) :: mi
  type(FieldCan), intent(inout) :: f

  double precision, intent(in) :: z(:)
  double precision, intent(in) :: dtau
  integer, intent(in) :: ntau
  double precision, intent(in) :: rtol_init

  double precision :: alpha(1), beta(1)

  alpha(1) = 0.5d0
  beta(1)  = 0.5d0

  call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
end subroutine orbit_sympl_init_verlet


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init_order4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=3
  !
  !
    type(MultistageIntegrator), intent(inout) :: mi
    type(FieldCan), intent(inout) :: f

    double precision, intent(in) :: z(:)
    double precision, intent(in) :: dtau
    integer, intent(in) :: ntau
    double precision, intent(in) :: rtol_init

    double precision :: alpha(5), beta(5)

    alpha(1) = 1d0/(2d0*(2d0 - 2d0**(1d0/3d0)))
    alpha(2) = 2d0**(1d0/3d0)/(2d0*(2d0 - 2d0**(1d0/3d0)))
    alpha(3) = alpha(1)

    beta = alpha

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_order4


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init_mclachlan4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=5 by McLachlan (1995)
  !
  !
    type(MultistageIntegrator), intent(inout) :: mi
    type(FieldCan), intent(inout) :: f

    double precision, intent(in) :: z(:)
    double precision, intent(in) :: dtau
    integer, intent(in) :: ntau
    double precision, intent(in) :: rtol_init

    double precision :: alpha(5), beta(5)

    alpha(1) = (146d0 + 5d0*dsqrt(19d0))/540d0
    alpha(2) = (-2d0 + 10d0*dsqrt(19d0))/135d0
    alpha(3) = 1d0/5d0
    alpha(4) = (-23d0 + 20d0*dsqrt(19d0))/270d0
    alpha(5) = (14d0 - dsqrt(19d0))/108d0

    beta(5)  = alpha(1)
    beta(4) = alpha(2)
    beta(3) = alpha(3)
    beta(2) = alpha(4)
    beta(1) = alpha(5)

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_mclachlan4


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init_blanes4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=6 by Blanes&Moan (2002)
  ! with coefficients in the form of Hairer (2002)
  !
  !
    type(MultistageIntegrator), intent(inout) :: mi
    type(FieldCan), intent(inout) :: f

    double precision, intent(in) :: z(:)
    double precision, intent(in) :: dtau
    integer, intent(in) :: ntau
    double precision, intent(in) :: rtol_init

    double precision :: alpha(6), beta(6)

    alpha(1) = 0.16231455076687d0
    alpha(2) = 0.37087741497958d0
    alpha(3) = 0.059762097006575d0
    alpha(4) = 0.40993371990193d0
    alpha(5) = 0.23399525073150d0
    alpha(6) = 0.082984406417405d0

    beta(6) = alpha(1)
    beta(5) = alpha(2)
    beta(4) = alpha(3)
    beta(3) = alpha(4)
    beta(2) = alpha(5)
    beta(1) = alpha(6)

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_blanes4


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init_kahan6(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 6 with s=9 by Kahan&Li (1995)
  ! with coefficients in the form of Hairer (2002)
  !
  !
    type(MultistageIntegrator), intent(inout) :: mi
    type(FieldCan), intent(inout) :: f

    double precision, intent(in) :: z(:)
    double precision, intent(in) :: dtau
    integer, intent(in) :: ntau
    double precision, intent(in) :: rtol_init

    double precision :: gam(9)

    gam(1) = 0.39216144400731413927925056d0
    gam(2) = 0.33259913678935943859974864d0
    gam(3) = -0.70624617255763935980996482d0
    gam(4) = 0.08221359629355080023149045d0
    gam(5) = 0.79854399093482996339895035d0
    gam(6) = gam(4)
    gam(7) = gam(3)
    gam(8) = gam(2)
    gam(9) = gam(1)

    call orbit_sympl_init_multi(mi,f,z,dtau,ntau,rtol_init,gam/2d0,gam/2d0)
  end subroutine orbit_sympl_init_kahan6


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine orbit_sympl_init_kahan8(mi, f, z, dtau, ntau, rtol_init)
    !
    ! Composition method of order 8 with s=17 by Kahan&Li (1995)
    ! with coefficients in the form of Hairer (2002)
    !
    !
      type(MultistageIntegrator), intent(inout) :: mi
      type(FieldCan), intent(inout) :: f

      double precision, intent(in) :: z(:)
      double precision, intent(in) :: dtau
      integer, intent(in) :: ntau
      double precision, intent(in) :: rtol_init

      double precision :: gam(17)

      gam(1) =  0.13020248308889008087881763d0
      gam(2) =  0.56116298177510838456196441d0
      gam(3) = -0.38947496264484728640807860d0
      gam(4) =  0.15884190655515560089621075d0
      gam(5) = -0.39590389413323757733623154d0
      gam(6) =  0.18453964097831570709183254d0
      gam(7) =  0.25837438768632204729397911d0
      gam(8) =  0.29501172360931029887096624d0
      gam(9) = -0.60550853383003451169892108d0
      gam(10) = gam(8)
      gam(11) = gam(7)
      gam(12) = gam(6)
      gam(13) = gam(5)
      gam(14) = gam(4)
      gam(15) = gam(3)
      gam(16) = gam(2)
      gam(17) = gam(1)

      call orbit_sympl_init_multi(mi,f,z,dtau,ntau,rtol_init,gam/2d0,gam/2d0)
    end subroutine orbit_sympl_init_kahan8

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_euler1(si, f, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  integer, intent(out) :: ierr

  integer, parameter :: n = 2
  integer, parameter :: maxit = 32

  double precision, dimension(n) :: x, xlast
  integer :: k, ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    ! Initial guess with Lagrange extrapolation
    if (si%nlag>0) then
      do k=0, si%nlag
        si%bufind(k) = si%kbuf-si%nlag+k
        if (si%bufind(k)<1) si%bufind(k) = si%bufind(k) + si%nbuf
      end do
    end if

    if (si%nlag>0 .and. si%kt>si%nlag) then
      x(1)=sum(si%zbuf(1,si%bufind)*si%coef(0,:))
      x(2)=sum(si%zbuf(4,si%bufind)*si%coef(0,:))
    else
      x(1)=si%z(1)
      x(2)=si%z(4)
    end if

    ! correct if Lagrange messed up
    if (x(1) < 0.0 .or. x(1) > 1.0) then
      x(1) = si%z(1)
      x(2) = si%z(4)
    end if

    call newton1(si, f, x, maxit, xlast)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    si%z(1) = x(1)
    si%z(4) = x(2)

    if (si%extrap_field) then
      f%pth = f%pth + f%dpth(1)*(x(1)-xlast(1))  + f%dpth(4)*(x(2)-xlast(2))
      f%dH(1) = f%dH(1) + f%d2H(1)*(x(1)-xlast(1)) + f%d2H(7)*(x(2)-xlast(2))
      f%dpth(1)=f%dpth(1)+f%d2pth(1)*(x(1)-xlast(1))+f%d2pth(7)*(x(2)-xlast(2))
      f%vpar = f%vpar + f%dvpar(1)*(x(1)-xlast(1)) + f%dvpar(4)*(x(2)-xlast(2))
      f%hth = f%hth + f%dhth(1)*(x(1)-xlast(1))
      f%hph = f%hph + f%dhph(1)*(x(1)-xlast(1))
    else
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    endif

    si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
    si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph

    if (si%nbuf > 0) then
      si%kbuf = mod(si%kt, si%nbuf) + 1
      si%zbuf(1:4,si%kbuf) = si%z
    endif

    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_euler1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_euler2(si, f, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: n = 3
  integer, parameter :: maxit = 32

  double precision, dimension(n) :: x, xlast, dz
  integer :: k, ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    ! Initial guess with Lagrange extrapolation
    if (si%nlag>0) then
      do k=0, si%nlag
        si%bufind(k) = si%kbuf-si%nlag+k
        if (si%bufind(k)<1) si%bufind(k) = si%bufind(k) + si%nbuf
      end do
    endif

    if (si%nlag>0 .and. si%kt>si%nlag) then
      x(1)=sum(si%zbuf(1,si%bufind)*si%coef(0,:))
      x(2)=sum(si%zbuf(2,si%bufind)*si%coef(0,:))
      x(3)=sum(si%zbuf(3,si%bufind)*si%coef(0,:))
    else
      x = si%z(1:3)
    end if

    call newton2(si, f, x, si%atol, si%rtol, maxit, xlast)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    si%z(1:3) = x

    if (si%extrap_field) then
      dz(1) = x(1)-xlast(1)
      dz(2) = x(2)-xlast(2)
      dz(3) = x(3)-xlast(3)

      f%dH(1) = f%dH(1) + f%d2H(1)*dz(1) + f%d2H(2)*dz(2) + f%d2H(3)*dz(3)
      f%dH(2) = f%dH(2) + f%d2H(2)*dz(1) + f%d2H(4)*dz(2) + f%d2H(5)*dz(3)
      f%dH(3) = f%dH(3) + f%d2H(3)*dz(1) + f%d2H(5)*dz(2) + f%d2H(6)*dz(3)

      f%dpth(1) = f%dpth(1)+f%d2pth(1)*dz(1)+f%d2pth(2)*dz(2)+f%d2pth(3)*dz(3)
      f%dpth(2) = f%dpth(2)+f%d2pth(2)*dz(1)+f%d2pth(4)*dz(2)+f%d2pth(5)*dz(3)
      f%dpth(3) = f%dpth(3)+f%d2pth(3)*dz(1)+f%d2pth(5)*dz(2)+f%d2pth(6)*dz(3)
    else
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    endif

    f%pth = si%pthold - si%dt*(f%dH(2) - f%dH(1)*f%dpth(2)/f%dpth(1))
    si%z(4) = si%z(4) - si%dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

    if (si%nbuf > 0) then
      si%kbuf = mod(si%kt, si%nbuf) + 1
      si%zbuf(1:4,si%kbuf) = si%z
    endif

    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_midpoint(si, f, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: n = 5
  integer, parameter :: maxit = 32

  double precision, dimension(n) :: x, xlast
  integer :: k, ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    ! Initial guess with Lagrange extrapolation
    if (si%nlag>0) then
      do k=0, si%nlag
        si%bufind(k) = si%kbuf-si%nlag+k
        if (si%bufind(k)<1) si%bufind(k) = si%bufind(k) + si%nbuf
      end do
    endif

    if (si%nlag>0 .and. si%kt>si%nlag) then
      x(1)=sum(si%zbuf(1,si%bufind)*si%coef(0,:))
      x(2)=sum(si%zbuf(2,si%bufind)*si%coef(0,:))
      x(3)=sum(si%zbuf(3,si%bufind)*si%coef(0,:))
      x(4)=sum(si%zbuf(4,si%bufind)*si%coef(0,:))
      x(5)=sum(si%zbuf(5,si%bufind)*si%coef(0,:))
    else
      x(1:4) = si%z
      x(5) = si%z(1)
    end if

    call newton_midpoint(si, f, x, si%atol, si%rtol, maxit, xlast)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    si%z = x(1:4)

    if (si%extrap_field) then
      f%pth = f%pth + f%dpth(1)*(x(1)-xlast(1) + x(5) - xlast(5)) &  ! d/dr
                    + f%dpth(2)*(x(2)-xlast(2)) &  ! d/dth
                    + f%dpth(3)*(x(3)-xlast(3)) &  ! d/dph
                    + f%dpth(4)*(x(4)-xlast(4))    ! d/dpph
    else
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_val(f, si%z(4))
    endif

    if (si%nbuf > 0) then
      si%kbuf = mod(si%kt, si%nbuf) + 1
      si%zbuf(1:4,si%kbuf) = si%z
      si%zbuf(5,si%kbuf) = x(5)  ! midpoint radius
    endif

    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_midpoint


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_rk_gauss(si, f, s, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: maxit = 32

  integer, intent(in) :: s
  double precision, dimension(4*s) :: x, xlast
  integer :: k, l, ktau

  type(FieldCan) :: fs(s)
  double precision :: a(s,s), b(s), c(s), Hprime(s), dz(4*s)

  do k = 1,s
    fs(k) = f
    fs(k) = f
  end do

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    ! Initial guess with Lagrange extrapolation
    if (si%nlag>0) then
      do k=0, si%nlag
        si%bufind(k) = si%kbuf-si%nlag+k
        if (si%bufind(k)<1) si%bufind(k) = si%bufind(k) + si%nbuf/(4*s)
      end do
    end if

    if (si%nlag>0 .and. si%kt>si%nlag) then
      do k = 1,s
        x(1)=sum(si%zbuf(1,si%bufind+(4*k-3))*si%coef(0,:))
        x(1)=sum(si%zbuf(2,si%bufind+(4*k-2))*si%coef(0,:))
        x(1)=sum(si%zbuf(3,si%bufind+(4*k-1))*si%coef(0,:))
        x(2)=sum(si%zbuf(4,si%bufind+(4*k))*si%coef(0,:))
      end do
    else
      do k = 1,s
        x((4*k-3):(4*k)) = si%z
      end do
    end if

    call newton_rk_gauss(si, fs, s, x, si%atol, si%rtol, maxit, xlast)
    !optionally try fixed point iterations, doesn't work yet
    !call fixpoint_rk_gauss(si, fs, s, x, si%atol, si%rtol, maxit, xlast)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

    if (si%extrap_field) then
      do k = 1, s
        dz(1) = x(4*k-3)-xlast(4*k-3)
        dz(2) = x(4*k-2)-xlast(4*k-2)
        dz(3) = x(4*k-1)-xlast(4*k-1)
        dz(4) = x(4*k)-xlast(4*k)

        fs(k)%pth = fs(k)%pth + fs(k)%dpth(1)*dz(1) &  ! d/dr
                      + fs(k)%dpth(2)*dz(2) &          ! d/dth
                      + fs(k)%dpth(3)*dz(3) &          ! d/dph
                      + fs(k)%dpth(4)*dz(4)            ! d/dpph

        fs(k)%dpth(1) = fs(k)%dpth(1) + fs(k)%d2pth(1)*dz(1) &
                      + fs(k)%d2pth(2)*dz(2) &
                      + fs(k)%d2pth(3)*dz(3) &
                      + fs(k)%d2pth(7)*dz(4)

        fs(k)%dpth(2) = fs(k)%dpth(2) + fs(k)%d2pth(2)*dz(1)&
                      + fs(k)%d2pth(4)*dz(2) &
                      + fs(k)%d2pth(5)*dz(3) &
                      + fs(k)%d2pth(8)*dz(4)

        fs(k)%dpth(3) = fs(k)%dpth(3) + fs(k)%d2pth(2)*dz(1) &
                      + fs(k)%d2pth(5)*dz(2) &
                      + fs(k)%d2pth(6)*dz(3) &
                      + fs(k)%d2pth(9)*dz(4)

        fs(k)%dH(1) = fs(k)%dH(1) + fs(k)%d2H(1)*dz(1) &
                      + fs(k)%d2H(2)*dz(2) &
                      + fs(k)%d2H(3)*dz(3) &
                      + fs(k)%d2H(7)*dz(4)

        fs(k)%dH(2) = fs(k)%dH(2) + fs(k)%d2H(2)*dz(1) &
                      + fs(k)%d2H(4)*dz(2) &
                      + fs(k)%d2H(5)*dz(3) &
                      + fs(k)%d2H(8)*dz(4)

        fs(k)%dH(3) = fs(k)%dH(3) + fs(k)%d2H(3)*dz(1) &
                      + fs(k)%d2H(5)*dz(2) &
                      + fs(k)%d2H(6)*dz(3) &
                      + fs(k)%d2H(9)*dz(4)

        fs(k)%vpar = fs(k)%vpar + fs(k)%dvpar(1)*dz(1) &  ! d/dr
                   + fs(k)%dvpar(2)*dz(2) &
                   + fs(k)%dvpar(3)*dz(3)

        fs(k)%hth = fs(k)%hth + fs(k)%dhth(1)*dz(1) &  ! d/dr
                      + fs(k)%dhth(2)*dz(2) &
                      + fs(k)%dhth(3)*dz(3)

        fs(k)%hph = fs(k)%hph + fs(k)%dhph(1)*dz(1) &  ! d/dr
                      + fs(k)%dhph(2)*dz(2) &
                      + fs(k)%dhph(3)*dz(3)
      end do
    else
      do k = 1, s
        call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 0)
        call get_derivatives(fs(k), x(4*k))
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
      si%z(3) = si%z(3)+si%dt*b(l)*(fs(l)%vpar-Hprime(l)*fs(l)%hth)/fs(l)%hph
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
    end do

    if (si%nbuf > 0) then
      si%kbuf = mod(si%kt, si%nbuf) + 1
      do k = 1, s
        si%zbuf(1, si%kbuf+(4*k-3)) = si%z(1)
        si%zbuf(2, si%kbuf+(4*k-2)) = si%z(2)
        si%zbuf(3, si%kbuf+(4*k-1)) = si%z(3)
        si%zbuf(4, si%kbuf+(4*k)) = si%z(4)
      end do
    endif

    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_rk_gauss


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_rk_lobatto(si, f, s, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  integer, intent(out) :: ierr

  integer, parameter :: maxit = 32

  integer, intent(in) :: s
  double precision, dimension(4*s-2) :: x, xlast
  integer :: ktau, k, l

  type(FieldCan) :: fs(s)
  double precision :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s)

  do k = 1,s
    fs(k) = f
  end do

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x(1) = si%z(1)
    x(2) = si%z(4)
    do k = 2,s
      x((4*k-3-2):(4*k-2)) = si%z
    end do

    call newton_rk_lobatto(si, fs, s, x, si%atol, si%rtol, maxit, xlast)

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

    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_rk_lobatto


subroutine debug_root(si, f, x0)
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), pointer :: f

  double precision :: x0(2), x(2)
  integer :: k, l, iflag
  integer, parameter :: n = 100
  double precision, parameter :: eps = 1d-15

  double precision :: fvec(2)

  do k = -n,n
    do l = -n,n
      x = x0 + l*eps/n*(/x0(1),0d0/) + k*eps/n*(/0d0,x0(2)/)
      call f_sympl_euler1(si, f, 2, x, fvec, iflag)
      write(5001,*) x(1), x(2), fvec
    end do
  end do

end subroutine debug_root

end module orbit_symplectic

module orbit_symplectic_global
  use field_can_mod, only: FieldCan
  use orbit_symplectic, only: SymplecticIntegrator, MultistageIntegrator, &
    orbit_sympl_init, orbit_sympl_init_multi, orbit_timestep_sympl, &
    orbit_timestep_sympl_multi

  implicit none
  save

  type(SymplecticIntegrator) :: si
  type(MultistageIntegrator) :: mi
  type(FieldCan) :: f
  !$omp threadprivate(si, f)

contains
  ! subroutine init(z, dt, ntau, rtol_init, mode_init, nlag)
  !   double precision, dimension(:), intent(in) :: z
  !   double precision, intent(in) :: dt
  !   integer, intent(in) :: ntau
  !   double precision, intent(in) :: rtol_init
  !   integer, intent(in) :: mode_init
  !   integer, intent(in) :: nlag

  !   call orbit_sympl_init(si, f, z, dt, ntau, rtol_init, mode_init, nlag)
  ! end subroutine init

  subroutine init_multi(z, dtau, ntau, rtol_init, alpha, beta)
    double precision, dimension(:), intent(in) :: z
    double precision, intent(in) :: dtau
    integer, intent(in) :: ntau
    double precision, intent(in) :: rtol_init
    double precision, dimension(:), intent(in) :: alpha
    double precision, dimension(:), intent(in) :: beta

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine init_multi

  ! subroutine timestep(z, ierr)
  !   double precision, intent(out) :: z(:)
  !   integer, intent(out) :: ierr

  !   call orbit_timestep_sympl(si, f, ierr)
  !   z = si%z
  ! end subroutine timestep

  subroutine timesteps(n, zs, ierr)
    integer, intent(in) :: n
    double precision, intent(out) :: zs(:,:)
    integer, intent(out) :: ierr

    integer :: k

    do k = 1, n
      call orbit_timestep_sympl(si, f, ierr)
      zs(:, k) = si%z
    end do
  end subroutine timesteps


end module orbit_symplectic_global
