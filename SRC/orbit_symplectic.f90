module orbit_symplectic

use common, only: pi
use field_can_mod, only: FieldCan, eval_field, get_val, get_derivatives, get_derivatives2

implicit none
save

public

integer, parameter :: NLAG_MAX = 2
integer, parameter :: NBUF_MAX = 4*NLAG_MAX

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
integer, parameter :: S_MAX = 64
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
  integer, intent(in) :: mode_init ! 1 = expl.-impl. Euler, 2 = impl.-expl. Euler
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
subroutine f_midpoint_part1(si, f, n, x, fvec, iflag)
  !
    type(SymplecticIntegrator), intent(in) :: si
    type(FieldCan), intent(inout) :: f
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    double precision, intent(out) :: fvec(n)
    integer, intent(in) :: iflag
    
    double precision :: thmid, phmid


    ! evaluate at midpoint
    call eval_field(f, x(5), 0.5*(x(2) + si%z(2)), 0.5*(x(3) + si%z(3)), 2)
    call get_derivatives2(f, 0.5*(x(4) + si%z(4)))

    fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
    fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)
    fvec(4) = f%dpth(1)*(x(4) - si%z(4)) + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
    fvec(5) = f%dpth(1)*(f%pth - si%pthold) + si%dt/2.0d0*(f%dpth(1)*f%dH(2)-f%dpth(2)*f%dH(1))
  
  end subroutine f_midpoint_part1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_midpoint_part2(si, f, n, x, fvec, iflag)
  !
    type(SymplecticIntegrator), intent(in) :: si
    type(FieldCan), intent(inout) :: f
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    double precision, intent(out) :: fvec(n)
    integer, intent(in) :: iflag

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
  integer :: kit

  do kit = 1, maxit
    if(x(1) > 1.0) return
    if(x(1) < 0.0) x(1) = 0.2
    call f_sympl_euler1(si, f, n, x, fvec, 1)
    call jac_sympl_euler1(si, f, x, fjac)
    ijac(1,1) = fjac(2,2)
    ijac(1,2) = -fjac(1,2)
    ijac(2,1) = -fjac(2,1)
    ijac(2,2) = fjac(1,1)
    ijac = ijac/(fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1))
    xlast = x
    x = x - matmul(ijac, fvec)
    if (all(dabs(fvec) < si%atol)) return
    if (all(dabs(x-xlast) < si%rtol*dabs(x))) return
  enddo
  print *, 'newton1: maximum iterations reached: ', maxit
  write(6601,*) x(1), si%z(2), si%z(3), x(2), x-xlast, fvec
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

  double precision :: fvec(n), fjac(n,n)
  integer :: pivot(n), info
  
  double precision :: xabs(n), tolref(n), fabs(n)

  do kit = 1, maxit
    if(x(1) > 1.0) return
    if(x(1) < 0.0) x(1) = 0.01
    call f_sympl_euler2(si, f, n, x, fvec, 1)
    fabs = dabs(fvec)
    call jac_sympl_euler2(si, f, x, fjac)
    xlast = x
    call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x-xlast)
    xabs(2) = modulo(xabs(2), pi)
    xabs(3) = modulo(xabs(3), pi)
    
    tolref = dabs(x)
    tolref(2) = 2d0*pi
    tolref(3) = 2d0*pi
    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton2: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6602,*) x(1), x(2), x(3), si%z(4), xabs, fvec
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
    call f_midpoint_part1(si, f, n, x, fvec, 1)
    call jac_midpoint_part1(si, f, x, fjac)
    fmid = f
    call f_midpoint_part2(si, f, n, x, fvec, 1)
    call jac_midpoint_part2(si, f, fmid, x, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)
    xabs(2) = modulo(xabs(2), pi)
    xabs(3) = modulo(xabs(3), pi)
    
    tolref = dabs(x)
    tolref(2) = 2d0*pi
    tolref(3) = 2d0*pi
    
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
  else
    ! not implemented
    a = 0d0
    b = 0d0
    c = 0d0
  endif
end subroutine coeff_rk_gauss


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_rk_gauss(si, fs, n, x, fvec, iflag)
  !
    type(SymplecticIntegrator), intent(in) :: si
    integer, intent(in) :: n
    type(FieldCan), intent(inout) :: fs(n)
    double precision, intent(in) :: x(4*n)  ! = (rend, thend, phend, pphend)
    double precision, intent(out) :: fvec(4*n)
    integer, intent(in) :: iflag

    double precision :: a(n,n), b(n), c(n), Hprime(n)
    integer :: k,l  ! counters

    call coeff_rk_gauss(n, a, b, c)  ! TODO: move this to preprocessing

    ! evaluate stages
    do k = 1, n
      call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 2)
      call get_derivatives2(fs(k), x(4*k))
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
    end do

    do k = 1, n
      fvec(4*k-3) = fs(k)%pth - si%pthold
      fvec(4*k-2) = x(4*k-2)  - si%z(2)
      fvec(4*k-1) = x(4*k-1)  - si%z(3)
      fvec(4*k)   = x(4*k)    - si%z(4)
      do l = 1, n
        fvec(4*k-3) = fvec(4*k-3) + si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
        fvec(4*k-2) = fvec(4*k-2) - si%dt*a(k,l)*Hprime(l)                                      ! thdot
        fvec(4*k-1) = fvec(4*k-1) - si%dt*a(k,l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
        fvec(4*k)   = fvec(4*k)   + si%dt*a(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
      end do
    end do
  
  end subroutine f_rk_gauss



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine jac_rk_gauss(si, fs, n, x, jac)
    !
      type(SymplecticIntegrator), intent(in) :: si
      integer, intent(in) :: n
      type(FieldCan), intent(in) :: fs(n)
      double precision, intent(in)  :: x(4*n)
      double precision, intent(out) :: jac(4*n, 4*n)

      double precision :: a(n,n), b(n), c(n), Hprime(n), dHprime(4*n)
      integer :: k,l  ! counters
  
      call coeff_rk_gauss(n, a, b, c)  ! TODO: move this to preprocessing
    
      ! evaluate stages
      do k = 1,n
        Hprime(k)      = fs(k)%dH(1)/fs(k)%dpth(1)
        dHprime(4*k-3) = (fs(k)%d2H(1) - Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)  ! d/dr
        dHprime(4*k-2) = (fs(k)%d2H(2) - Hprime(k)*fs(k)%d2pth(2))/fs(k)%dpth(1)  ! d/dth 
        dHprime(4*k-1) = (fs(k)%d2H(3) - Hprime(k)*fs(k)%d2pth(3))/fs(k)%dpth(1)  ! d/dph
        dHprime(4*k)   = (fs(k)%d2H(7) - Hprime(k)*fs(k)%d2pth(7))/fs(k)%dpth(1)  ! d/dpph 
      end do

      jac = 0d0
        
      do k = 1, n
        jac(4*k-3, 4*k-3) = fs(k)%dpth(1)
        jac(4*k-3, 4*k-2) = fs(k)%dpth(2)
        jac(4*k-3, 4*k-1) = fs(k)%dpth(3)
        jac(4*k-3, 4*k)   = fs(k)%dpth(4)
        jac(4*k-2, 4*k-2) = 1d0
        jac(4*k-1, 4*k-1) = 1d0
        jac(4*k, 4*k) = 1d0

        do l = 1, n
           jac(4*k-3, 4*l-3) = jac(4*k-3, 4*l-3) & ! d/dr
             + si%dt*a(k,l)*(fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-3))
           jac(4*k-3, 4*l-2) = jac(4*k-3, 4*l-2) & ! d/dth
             + si%dt*a(k,l)*(fs(l)%d2H(4) - fs(l)%d2pth(4)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-2))
           jac(4*k-3, 4*l-1) = jac(4*k-3, 4*l-1) & ! d/dph
             + si%dt*a(k,l)*(fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-1))
           jac(4*k-3, 4*l) = jac(4*k-3, 4*l) & ! d/dpph
             + si%dt*a(k,l)*(fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l))
          
          jac(4*k-2, 4*l-3)   = jac(4*k-2, 4*l-3)   - si%dt*a(k,l)*dHprime(4*l-3)   ! d/dr
          jac(4*k-2, 4*l-2) = jac(4*k-2, 4*l-2) - si%dt*a(k,l)*dHprime(4*l-2) ! d/dth
          jac(4*k-2, 4*l-1) = jac(4*k-2, 4*l-1) - si%dt*a(k,l)*dHprime(4*l-1) ! d/dph
          jac(4*k-2, 4*l) = jac(4*k-2, 4*l) - si%dt*a(k,l)*dHprime(4*l) ! d/dpph

          jac(4*k-1, 4*l-3) = jac(4*k-1, 4*l-3) & ! d/dr
            - si%dt*a(k,l)*(-fs(l)%dhph(1)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
              + (fs(l)%dvpar(1) - dHprime(4*l-3)*fs(l)%hth - Hprime(l)*fs(l)%dhth(1))/fs(l)%hph)
          jac(4*k-1, 4*l-2) = jac(4*k-1, 4*l-2) & ! d/dth
            - si%dt*a(k,l)*(-fs(l)%dhph(2)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
              + (fs(l)%dvpar(2) - dHprime(4*l-2)*fs(l)%hth - Hprime(l)*fs(l)%dhth(2))/fs(l)%hph)
          jac(4*k-1, 4*l-1) = jac(4*k-1, 4*l-1) & ! d/dph
            - si%dt*a(k,l)*(-fs(l)%dhph(3)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
              + (fs(l)%dvpar(3) - dHprime(4*l-1)*fs(l)%hth - Hprime(l)*fs(l)%dhth(3))/fs(l)%hph)
          jac(4*k-1, 4*l) = jac(4*k-1, 4*l) & ! d/dpph
            - si%dt*a(k,l)*((fs(l)%dvpar(4) - dHprime(4*l)*fs(l)%hth)/fs(l)%hph)

          jac(4*k, 4*l-3) = jac(4*k, 4*l-3) & ! d/dr
            + si%dt*a(k,l)*(fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-3))
          jac(4*k, 4*l-2) = jac(4*k, 4*l-2) & ! d/dth
            + si%dt*a(k,l)*(fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-2))
          jac(4*k, 4*l-1) = jac(4*k, 4*l-1) & ! d/dph
            + si%dt*a(k,l)*(fs(l)%d2H(6) - fs(l)%d2pth(6)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-1))
          jac(4*k, 4*l) = jac(4*k, 4*l) & ! d/dpph
            + si%dt*a(k,l)*(fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l)) 
        end do
  end do
    
  end subroutine jac_rk_gauss

  subroutine newton_rk_gauss(si, fs, x, atol, rtol, maxit, xlast)
    type(SymplecticIntegrator), intent(inout) :: si
  
    integer, parameter :: n = 2
    type(FieldCan), intent(inout) :: fs(n)
    integer :: kit
  
    double precision, intent(inout) :: x(4*n)
    double precision, intent(in) :: atol, rtol
    integer, intent(in) :: maxit
    double precision, intent(out) :: xlast(4*n)
  
    double precision :: fvec(4*n), fjac(4*n,4*n)
    integer :: pivot(4*n), info
    
    double precision :: xabs(4*n), tolref(4*n), fabs(4*n)

    do kit = 1, maxit
      if(x(1) > 1.0) return
      if(x(1) < 0.0) x(1) = 0.01
      if(x(5) < 0.0) x(5) = 0.01
      call f_rk_gauss(si, fs, n, x, fvec, 1)
      call jac_rk_gauss(si, fs, n, x, fjac)
      fabs = dabs(fvec)
      xlast = x
      call dgesv(4*n, 1, fjac, 4*n, pivot, fvec, 4*n, info)
      ! after solution: fvec = (xold-xnew)_Newton
      x = x - fvec
      xabs = dabs(x - xlast)
      xabs(2) = modulo(xabs(2), pi)
      xabs(3) = modulo(xabs(3), pi)
      xabs(6) = modulo(xabs(6), pi)
      xabs(7) = modulo(xabs(7), pi)
      
      tolref = dabs(x)
      tolref(2) = 2d0*pi
      tolref(3) = 2d0*pi
      tolref(6) = 2d0*pi
      tolref(7) = 2d0*pi
      
      if (all(fabs < atol)) return
      if (all(xabs < rtol*tolref)) return
    enddo
    print *, 'newton_rk_gauss: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
    write(6603,*) x, xabs, fvec
  end subroutine newton_rk_gauss

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl(si, f, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f

  integer, intent(out) :: ierr

  ! si%z(1:3) = z0(1:3)  ! r, th, ph
  ! si%z(4) = f%vpar*f%hph + f%Aph/f%ro0 ! pphi
  ! print *, si%z(1:3)-z0(1:3)
  ! print *, si%z(4)-(f%vpar*f%hph + f%Aph/f%ro0)

  select case (si%mode)
   case (1)
      call orbit_timestep_sympl_euler1(si, f, ierr)
   case (2)
      call orbit_timestep_sympl_euler2(si, f, ierr)
   case (3)
      call orbit_timestep_sympl_midpoint(si, f, ierr)
   case (4)
      call orbit_timestep_sympl_rk_gauss(si, f, ierr)
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
    mi%stages(kstage)%pthold = mi%stages(kstage-1)%pthold
    call orbit_timestep_sympl(mi%stages(kstage), f, ierr)
  end do
  mi%stages(1)%z = mi%stages(2*mi%s)%z
  mi%stages(1)%pthold = mi%stages(2*mi%s)%pthold

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
      alpha(ks)*dtau, ntau, rtol_init, 1, 0)
    call orbit_sympl_init(mi%stages(2*ks), f, z, &
      beta(ks)*dtau, ntau, rtol_init, 2, 0)
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
subroutine orbit_timestep_sympl_euler1(si, f, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: n = 2
  integer, parameter :: maxit = 256

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
      f%dpth(1) = f%dpth(1) + f%d2pth(1)*(x(1)-xlast(1)) + f%d2pth(7)*(x(2)-xlast(2))
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
  integer, parameter :: maxit = 256

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
      f%dH(1) = f%dH(1) + f%d2H(1)*(x(1)-xlast(1)) &  ! d/dr
                        + f%d2H(2)*(x(2)-xlast(2)) &  ! d/dth
                        + f%d2H(3)*(x(3)-xlast(3))    ! d/dph

      f%dH(2) = f%dH(2) + f%d2H(2)*(x(1)-xlast(1)) &
                        + f%d2H(4)*(x(2)-xlast(2)) & 
                        + f%d2H(5)*(x(3)-xlast(3))

      f%dH(3) = f%dH(3) + f%d2H(3)*(x(1)-xlast(1)) &
                        + f%d2H(5)*(x(2)-xlast(2)) &
                        + f%d2H(6)*(x(3)-xlast(3))

      f%dpth(1) = f%dpth(1) + f%d2pth(1)*(x(1)-xlast(1)) &
                            + f%d2pth(2)*(x(2)-xlast(2)) &
                            + f%d2pth(3)*(x(3)-xlast(3))

      f%dpth(2) = f%dpth(2) + f%d2pth(2)*(x(1)-xlast(1)) &
                            + f%d2pth(4)*(x(2)-xlast(2)) &
                            + f%d2pth(5)*(x(3)-xlast(3))

      f%dpth(3) = f%dpth(3) + f%d2pth(3)*(x(1)-xlast(1)) &
                            + f%d2pth(5)*(x(2)-xlast(2)) &
                            + f%d2pth(6)*(x(3)-xlast(3))
      
      f%vpar = f%vpar + f%dvpar(1)*(x(1)-xlast(1)) &
                      + f%dvpar(2)*(x(2)-xlast(2)) &
                      + f%dvpar(3)*(x(3)-xlast(3))

      f%hth = f%hth + f%dhth(1)*(x(1)-xlast(1)) &
                    + f%dhth(2)*(x(2)-xlast(2)) &
                    + f%dhth(3)*(x(3)-xlast(3))

      f%hph = f%hph + f%dhph(1)*(x(1)-xlast(1)) &
                    + f%dhph(2)*(x(2)-xlast(2)) &
                    + f%dhph(3)*(x(3)-xlast(3))
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
  integer, parameter :: maxit = 256

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
subroutine orbit_timestep_sympl_rk_gauss(si, f, ierr)
  !
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
  
    integer, intent(out) :: ierr
  
    integer, parameter :: n = 2
    integer, parameter :: maxit = 256
  
    double precision, dimension(4*n) :: x, xlast
    integer :: k, l, ktau

    type(FieldCan) :: fs(n)
    double precision :: a(n,n), b(n), c(n), Hprime(n)

    fs(1) = f
    fs(2) = f
  
    ierr = 0
    ktau = 0
    do while(ktau .lt. si%ntau)
      si%pthold = f%pth
      
      x(1:4) = si%z
      x(5:8) = si%z
  
      call newton_rk_gauss(si, fs, x, si%atol, si%rtol, maxit, xlast)
  
      if (x(1) > 1.0) then
        ierr = 1
        return
      end if
  
      if (x(1) < 0.0) then
        print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
        x(1) = 0.01
      end if
  
      call coeff_rk_gauss(n, a, b, c)  ! TODO: move this to preprocessing

      f = fs(2)
      f%pth = si%pthold
      si%z(1) = x(5)

      do l = 1, n
        Hprime(l) = fs(l)%dH(1)/fs(l)%dpth(1)
        f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))            ! pthdot
        si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)                                      ! thdot
        si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
        si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
      end do
    
      si%kt = si%kt+1
      ktau = ktau+1
    enddo
  
  end subroutine orbit_timestep_sympl_rk_gauss

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
