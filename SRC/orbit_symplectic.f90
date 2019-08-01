module orbit_symplectic

use parmot_mod, only : ro0_parmot => ro0
use field_can_mod, only: FieldCan, eval_field, get_val, get_derivatives, get_derivatives2

implicit none
save

integer, parameter :: nlag = 3 ! order
integer, parameter :: nbuf = 16 ! values to store back
logical, parameter :: extrap_field = .true.

public

type :: SymplecticIntegrator
  double precision :: atol = 1d-15
  double precision :: rtol

! Current phase-space coordinates z and old pth
  double precision, dimension(4) :: z  ! z = (r, th, ph, pphi)
  double precision :: pthold

  ! Buffer for Lagrange polynomial interpolation
  integer :: kbuf
  integer :: kt
  integer :: k
  integer :: bufind(0:nlag)
  double precision, dimension(4, nbuf) :: zbuf
  double precision, dimension(0:0, nlag+1) :: coef

  ! Timestep and variables from z0
  integer :: ntau
  double precision :: dt
  double precision :: pabs

  ! For initial conditions
  double precision :: z0init(5)

  ! Integrator mode
  integer :: mode  ! 1 = euler1, 2 = euler2, 3 = verlet

  ! Data for evaluating field, Hamiltonian and canonical poloidal momentum
  type(FieldCan), pointer :: f
end type SymplecticIntegrator

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init(si, z0, dtau, dtaumin, rtol_init, mode_init)
!
  type(SymplecticIntegrator), intent(inout) :: si
  double precision, intent(in) :: z0(5)
  double precision, intent(in) :: dtau, dtaumin
  double precision, intent(in) :: rtol_init
  integer, intent(in) :: mode_init ! 1 = euler1, 2 = euler2, 3 = verlet

  integer, parameter :: n = 2
  integer :: k

  si%mode = mode_init
  si%rtol = rtol_init

  si%kbuf = 0
  si%kt = 0
  si%k = 0

  if(min(dabs(mod(dtau, dtaumin)), dabs(mod(dtau, dtaumin)-dtaumin)) > 1d-9*dtaumin) then
    stop 'orbit_sympl_init - error: dtau/dtaumin not integer'
  endif

  si%ntau = nint(dtau/dtaumin)
  si%dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  if (si%mode==3) si%dt = si%dt/2d0 ! Verlet out of two Euler steps

  allocate(si%f)
  call eval_field(si%f, z0(1), z0(2), z0(3), 0)

  si%pabs = z0(4)

  si%f%mu = .5d0*si%pabs**2*(1.d0-z0(5)**2)/si%f%Bmod*2d0 ! mu by factor 2 from other modules
  si%f%ro0 = ro0_parmot/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
  si%f%vpar = si%pabs*z0(5)*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

  ! initialize canonical variables
  si%z(1:3) = z0(1:3)  ! r, th, ph
  si%z(4) = si%f%vpar*si%f%hph + si%f%Aph/si%f%ro0 ! pphi
  call get_val(si%f, si%z(4)) ! for pth

  call plag_coeff(nlag+1, 0, 1d0, 1d0*(/(k,k=-nlag,0)/), si%coef)
end subroutine orbit_sympl_init


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_sympl_euler1(si, n, x, fvec, iflag)
!
  type(SymplecticIntegrator), intent(inout) :: si
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  type(FieldCan), pointer :: f

  f => si%f

  call eval_field(f, x(1), si%z(2), si%z(3), 2)
  call get_derivatives2(f, x(2))

  fvec(1) = f%dpth(1)*(f%pth - si%pthold) + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
  fvec(2) = f%dpth(1)*(x(2) - si%z(4))  + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))

end subroutine f_sympl_euler1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_sympl_euler1(si, x, jac)
!
  type(SymplecticIntegrator), intent(in) :: si
  double precision, intent(in)  :: x(2)
  double precision, intent(out) :: jac(2, 2)

  type(FieldCan), pointer :: f

  f => si%f

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
subroutine f_sympl_euler2(si, n, x, fvec, iflag)
!
  type(SymplecticIntegrator), intent(in) :: si
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  type(FieldCan), pointer :: f

  f => si%f

  call eval_field(si%f, x(1), x(2), x(3), 2)
  call get_derivatives2(f, si%z(4))

  fvec(1) = f%pth - si%pthold
  fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
  fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)

end subroutine f_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_sympl_euler2(si, x, jac)
!
  type(SymplecticIntegrator), intent(in) :: si
  double precision, intent(in)  :: x(3)
  double precision, intent(out) :: jac(3, 3)

  type(FieldCan), pointer :: f

  f => si%f

  jac(1,:) = f%dpth(1:3)
  jac(2,:) = f%d2pth(1:3)*(x(2) - si%z(2)) - si%dt*f%d2H(1:3)
  jac(2,2) = jac(2,2) + f%dpth(1)
  jac(3,:) = (f%d2pth(1:3)*f%hph + f%dpth(1)*f%dhph)*(x(3) - si%z(3)) &
    - si%dt*(f%d2pth(1:3)*f%vpar + f%dpth(1)*f%dvpar(1:3) - f%d2H(1)*f%hth - f%dH(1)*f%dhth)
  jac(3,3) = jac(3,3) + f%dpth(1)*f%hph

end subroutine jac_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine newton1(si, x, maxit, xlast)
!
  type(SymplecticIntegrator), intent(inout) :: si
  integer, parameter :: n = 2

  double precision, intent(inout) :: x(n)
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(n)

  double precision :: fvec(n), fjac(n,n), ijac(n,n)
  integer :: kit

  do kit = 1, maxit
    if(x(1) > 1.0) return
    if(x(1) < 0.0) x(1) = 0.2
    call f_sympl_euler1(si, n, x, fvec, 1)
    call jac_sympl_euler1(si, x, fjac)
    ijac(1,1) = fjac(2,2)
    ijac(1,2) = -fjac(1,2)
    ijac(2,1) = -fjac(2,1)
    ijac(2,2) = fjac(1,1)
    ijac = ijac/(fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1))
    xlast = x
    x = x - matmul(ijac, fvec)
    if (all(dabs(fvec) < si%atol) &
      .or. all(dabs(x-xlast) < si%rtol*dabs(x))) return
  enddo
  print *, 'newton1: maximum iterations reached: ', maxit
  write(6601,*) x(1), si%z(2), si%z(3), x(2), x-xlast, fvec
end subroutine

subroutine newton2(si, x, atol, rtol, maxit, xlast)
  type(SymplecticIntegrator), intent(inout) :: si
  integer, parameter :: n = 3
  integer :: kit

  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(n)

  double precision :: fvec(n), fjac(n,n)
  integer :: pivot(n), info

  do kit = 1, maxit
    if(x(1) > 1.0) return
    if(x(1) < 0.0) x(1) = 0.01
    call f_sympl_euler2(si, n, x, fvec, 1)
    call jac_sympl_euler2(si, x, fjac)
    call dgesv(n, 1, fjac, n, pivot, fvec, 3, info)
    xlast = x
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    if (all(dabs(fvec) < atol) &
      .or. all(dabs(x-xlast) < rtol*dabs(x))) return
  enddo
  print *, 'newton2: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6602,*) x(1), x(2), x(3), si%z(4), x-xlast, fvec
end subroutine


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl(si, z0, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0

  select case (si%mode)
   case (1)
      call orbit_timestep_sympl_euler1(si, z0, ierr)
   case (2)
      call orbit_timestep_sympl_euler2(si, z0, ierr)
   case (3)
      call orbit_timestep_sympl_verlet(si, z0, ierr)
   case default
      stop 'invalid mode for orbit_timestep_sympl'
  end select

end subroutine orbit_timestep_sympl


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_euler1(si, z0, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0

  integer, parameter :: n = 2
  integer, parameter :: maxit = 256

  double precision, dimension(n) :: x, xlast
  integer :: k, ktau

  type(FieldCan), pointer :: f

  f => si%f
  ierr = 0

  if (z0(1) < 0.0 .or. z0(1) > 1.0) then
    ierr = 1
    return
  end if

  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    ! Initial guess with Lagrange extrapolation
    do k=0, nlag
      si%bufind(k) = si%kbuf-nlag+k
      if (si%bufind(k)<1) si%bufind(k) = si%bufind(k) + nbuf
    end do

    if (nlag>0 .and. si%kt>nlag) then
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

    call newton1(si, x, maxit, xlast)

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

    if (extrap_field) then
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

    si%kbuf = mod(si%kt, nbuf) + 1
    si%zbuf(:,si%kbuf) = si%z

    si%kt = si%kt+1
    ktau = ktau+1
  enddo
  z0(1:3) = si%z(1:3)
  !z0(4) = H
  z0(5) = f%vpar/(si%pabs*dsqrt(2d0))  ! alambda

end subroutine orbit_timestep_sympl_euler1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_euler2(si, z0, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0

  integer, parameter :: n = 3
  integer, parameter :: maxit = 256

  double precision, dimension(n) :: x, xlast
  integer :: ktau

  type(FieldCan), pointer :: f

  f => si%f

stop 'need to implement behavior on axis and for lost particles'

  ! ierr = 0
  ! ktau = 0
  ! do while(ktau .lt. si%ntau)
  !   pthold = f%pth

  !   ! Initial guess with Lagrange extrapolation
  !   do k=0, nlag
  !     bufind(k) = kbuf-nlag+k
  !     if (bufind(k)<1) bufind(k) = bufind(k) + nbuf
  !   end do

  !   if (nlag>0 .and. kt>nlag) then
  !     x(1)=sum(zbuf(1,bufind)*coef(0,:))
  !     x(2)=sum(zbuf(2,bufind)*coef(0,:))
  !     x(3)=sum(zbuf(3,bufind)*coef(0,:))
  !   else
  !     x = z(1:3)
  !   end if

  !   call newton2(x, atol, rtol, maxit, xlast)

  !   z(1:3) = x

  !   call eval_field(f, z(1), z(2), z(3), 0)
  !   call get_derivatives(f, z(4))

  !   f%pth = pthold - dt*(f%dH(2) - f%dH(1)*f%dpth(2)/f%dpth(1))
  !   z(4) = z(4) - dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

  !   kbuf = mod(kt, nbuf) + 1
  !   zbuf(:,kbuf) = z

  !   kt = kt+1
  !   ktau = ktau+1
  ! enddo
  ! z0(1:3) = z(1:3)
  ! !z0(4) = H
  ! z0(5) = f%vpar/(pabs*dsqrt(2d0))  ! alambda

end subroutine orbit_timestep_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_verlet(si, z0, ierr)
!
  type(SymplecticIntegrator), intent(inout) :: si
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0

  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 2
  integer, parameter :: maxit = 256

  double precision, dimension(n1) :: x1, xlast1
  double precision, dimension(n2) :: x2, xlast2
  integer :: ktau

stop 'need to implement behavior on axis and for lost particles'

  ! ierr = 0
  ! ktau = 0
  ! do while(ktau .lt. ntau)
  !   pthold = f%pth

  !   ! Verlet part one: impl/expl Euler

  !   ! Initial guess with Lagrange extrapolation

  !   if (nlag>0 .and. kt>nlag) then
  !     do k=0, nlag
  !       bufind(k) = kbuf-2*nlag+2*k-1
  !       if (bufind(k)<1) bufind(k) = bufind(k) + nbuf
  !     end do
  ! !    print *, bufind
  !     x1(1)=sum(zbuf(1,bufind)*coef(0,:))
  !     x1(2)=sum(zbuf(2,bufind)*coef(0,:))
  !     x1(3)=sum(zbuf(3,bufind)*coef(0,:))
  !   else
  !     x1 = z(1:3)
  !   end if

  !   ! correct if Lagrange messed up
  !   if (x1(1) < 0.0 .or. x1(1) > 1.0) then
  !     x1 = z(1:3)
  !   end if

  !   call newton2(x1, atol, rtol, maxit, xlast1)

  !   if (x1(1) < 0.0 .or. x1(1) > 1.0) then
  !     ierr = 1
  !     return
  !   end if

  !   z(1:3) = x1

  !   call eval_field(f, z(1), z(2), z(3), 0)
  !   call get_derivatives(f, z(4))

  !   pthold = pthold - dt*(f%dH(2) - f%dH(1)*f%dpth(2)/f%dpth(1))
  !   z(4) = z(4) - dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

  !   kbuf = mod(2*kt, nbuf) + 1
  !   zbuf(:,kbuf) = z

  !   ! Verlet part two: expl/impl Euler

  !   ! Initial guess with Lagrange extrapolation

  !   if (nlag>0 .and. kt>nlag) then
  !     do k=0, nlag
  !       bufind(k) = kbuf-2*nlag+2*k
  !       if (bufind(k)<1) bufind(k) = bufind(k) + nbuf
  !     end do
  !     x2(1)=sum(zbuf(1,bufind)*coef(0,:))
  !     x2(2)=sum(zbuf(4,bufind)*coef(0,:))
  !   else
  !     x2(1)=z(1)
  !     x2(2)=z(4)
  !   end if

  !   ! correct if Lagrange messed up
  !   if (x2(1) < 0.0 .or. x2(1) > 1.0) then
  !     x2(1) = z(1)
  !     x2(2) = z(4)
  !   end if

  !   call newton1(x2, atol, rtol, maxit, xlast2)

  !   if (x2(1) < 0.0 .or. x2(1) > 1.0) then
  !     ierr = 1
  !     return
  !   end if

  !   z(1) = x2(1)
  !   z(4) = x2(2)

  !   call eval_field(f, z(1), z(2), z(3), 0)
  !   call get_derivatives(f, z(4))
  !   z(2) = z(2) + dt*f%dH(1)/f%dpth(1)
  !   z(3) = z(3) + dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph

  !   kbuf = mod(2*kt+1, nbuf) + 1
  !   zbuf(:,kbuf) = z

  !   kt = kt+1
  !   ktau = ktau+1
  ! enddo
  ! z0(1:3) = z(1:3)
  ! !z0(4) = H
  ! z0(5) = f%vpar/(pabs*dsqrt(2d0))  ! alambda

end subroutine orbit_timestep_sympl_verlet

subroutine debug_root(si, x0)
  type(SymplecticIntegrator), intent(inout) :: si
  double precision :: x0(2), x(2)
  integer :: k, l, iflag
  integer, parameter :: n = 100
  double precision, parameter :: eps = 1d-15

  double precision :: fvec(2)

  do k = -n,n
    do l = -n,n
      x = x0 + l*eps/n*(/x0(1),0d0/) + k*eps/n*(/0d0,x0(2)/)
      call f_sympl_euler1(si, 2, x, fvec, iflag)
      write(5001,*) x(1), x(2), fvec
    end do
  end do

end subroutine debug_root

end module orbit_symplectic
