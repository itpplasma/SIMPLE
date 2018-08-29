module orbit_symplectic

use parmot_mod, only : ro0_parmot => ro0
use field_can_mod, only: field_can, d_field_can, d2_field_can, eval_field, &
  get_val, get_derivatives, H, pth, vpar, dvpar, dH, dpth, f, df, d2f, ro0, mu

implicit none
save


double precision, dimension(4) :: z  ! z = (r, th, ph, pphi)
double precision :: pthold

integer, parameter :: nbuf = 16 ! values to store back
integer :: kbuf = 0
integer :: kt = 0
double precision, dimension(4, nbuf) :: zbuf

double precision :: dt

double precision :: coala
double precision :: derphi(3)
double precision :: alambd, pabs

interface orbit_timestep_sympl
  module procedure orbit_timestep_sympl_verlet
end interface

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init(z0)
!
  double precision, intent(in) :: z0(5)

  call eval_field(z0(1), z0(2), z0(3), 0)

  pabs=z0(4)
  alambd=z0(5)

  mu = .5d0*pabs**2*(1.d0-alambd**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
  ro0 = ro0_parmot/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
  vpar = pabs*alambd*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

  z(1:3) = z0(1:3)  ! r, th, ph
  z(4) = vpar*f%Bph/f%Bmod + f%Aph/ro0 ! pphi
  pth = vpar*f%Bth/f%Bmod + f%Ath/ro0 ! pth

end subroutine orbit_sympl_init


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_sympl_euler1(n, x, fvec, iflag)
!
  integer, parameter :: mode = 2

  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(x(1), z(2), z(3), 0)
  call get_derivatives(x(2))

  fvec(1) = dpth(1)*(pth - pthold) + dt*(dH(2)*dpth(1) - dH(1)*dpth(2))
  fvec(2) = dpth(1)*(x(2) - z(4))  + dt*(dH(3)*dpth(1) - dH(1)*dpth(3))

  !print *, fvec

end subroutine f_sympl_euler1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_sympl_euler2(n, x, fvec, iflag)
!
  integer, parameter :: mode = 2

  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(out) :: iflag

  call eval_field(x(1), x(2), x(3), 0)
  call get_derivatives(z(4))

  fvec(1) = pth - pthold
  fvec(2) = dpth(1)*(x(2) - z(2)) - dt*dH(1)
  fvec(3) = dpth(1)*(x(3) - z(3)) - dt*(dpth(1)*vpar*f%Bmod - dH(1)*f%Bth)/f%Bph

end subroutine f_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_euler1(z0, dtau, dtaumin, ierr)
!
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  integer, parameter :: n = 2
  double precision :: tol

  double precision, dimension(n) :: fvec, x

  double precision :: tau2

  ! for Lagrange interpolation
  integer, parameter :: nlag = 2 ! order
  integer :: bufind(0:nlag), k
  double precision, dimension(0:0, nlag+1) :: coef

  call plag_coeff(nlag+1, 0, 1d0, 1d0*(/(k,k=-nlag,0)/), coef)
  
  ierr = 0

  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0d0
  do while(tau2.lt.dtau)
    pthold = pth

    ! Initial guess with Lagrange extrapolation
    do k=0, nlag
      bufind(k) = kbuf-nlag+k
      if (bufind(k)<1) bufind(k) = bufind(k) + nbuf 
    end do
    
    if (nlag>0 .and. kt>nlag) then
      x(1)=sum(zbuf(1,bufind)*coef(0,:))
      x(2)=sum(zbuf(4,bufind)*coef(0,:))
    else
      x(1)=z(1)
      x(2)=z(4)
    end if

    !print *, z(1), x(1) 
    !print *, z(4), x(2)

    tol = 1d-8
    call hybrd1 (f_sympl_euler1, n, x, fvec, tol, ierr)
    if(ierr > 1) stop 'error: root finding'
    
    z(1) = x(1)
    z(4) = x(2)

    call eval_field(z(1), z(2), z(3), 0)
    call get_derivatives(z(4))
    z(2) = z(2) + dt*dH(1)/dpth(1)
    z(3) = z(3) + dt*(vpar*f%Bmod - dH(1)/dpth(1)*f%Bth)/f%Bph

    kbuf = mod(kt, nbuf) + 1
    zbuf(:,kbuf) = z
    kt = kt+1

    tau2 = tau2 + dtaumin
  enddo
  z0(1:3) = z(1:3)
  z0(4) = H
  z0(5) = vpar/sqrt(H) ! alambda

end subroutine orbit_timestep_sympl_euler1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_euler2(z0, dtau, dtaumin, ierr)
!
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  integer, parameter :: n = 3
  double precision :: tol

  double precision, dimension(n) :: x, fvec

  double precision :: tau2

  ! for Lagrange interpolation
  integer, parameter :: nlag = 2 ! order
  integer :: bufind(0:nlag), k
  double precision, dimension(0:0, nlag+1) :: coef

  call plag_coeff(nlag+1, 0, 1d0, 1d0*(/(k,k=-nlag,0)/), coef)
  
  ierr = 0

  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0d0
  do while(tau2.lt.dtau)
    pthold = pth

    ! Initial guess with Lagrange extrapolation
    do k=0, nlag
      bufind(k) = kbuf-nlag+k
      if (bufind(k)<1) bufind(k) = bufind(k) + nbuf 
    end do
    
    if (nlag>0 .and. kt>nlag) then
      x(1)=sum(zbuf(1,bufind)*coef(0,:))
      x(2)=sum(zbuf(2,bufind)*coef(0,:))
      x(3)=sum(zbuf(3,bufind)*coef(0,:))
    else
      x = z(1:3)
    end if
    
    tol = 1d-8
    
    call hybrd1 (f_sympl_euler2, n, x, fvec, tol, ierr)
    if(ierr > 1) stop 'error: root finding'
    z(1:3) = x

    call eval_field(z(1), z(2), z(3), 0)
    call get_derivatives(z(4))

    pth = pthold - dt*(dH(2) - dH(1)*dpth(2)/dpth(1))
    z(4) = z(4) - dt*(dH(3) - dH(1)*dpth(3)/dpth(1))

    kbuf = mod(kt, nbuf) + 1
    zbuf(:,kbuf) = z
    kt = kt+1

    tau2 = tau2 + dtaumin
  enddo
  z0(1:3) = z(1:3)
  z0(4) = H
  z0(5) = vpar/sqrt(H) ! alambda

end subroutine orbit_timestep_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_verlet(z0, dtau, dtaumin, ierr)
!
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 2
  double precision :: tol

  double precision, dimension(n1) :: x1, fvec1
  double precision, dimension(n2) :: x2, fvec2

  double precision :: tau2
  
  ! for Lagrange interpolation
  integer, parameter :: nlag = 2 ! order
  integer :: bufind(0:nlag), k
  double precision, dimension(0:0, nlag+1) :: coef

  call plag_coeff(nlag+1, 0, 1d0, 1d0*(/(k,k=-nlag,0)/), coef)

  ierr = 0

  dt = 0.5d0*dtaumin/dsqrt(2d0)
  tau2 = 0d0
  do while(tau2.lt.dtau)
    pthold = pth

    ! Verlet part one: impl/expl Euler

    ! Initial guess with Lagrange extrapolation
    
    if (nlag>0 .and. kt>nlag) then
      do k=0, nlag
        bufind(k) = kbuf-2*nlag+2*k-1
        if (bufind(k)<1) bufind(k) = bufind(k) + nbuf 
      end do
  !    print *, bufind
      x1(1)=sum(zbuf(1,bufind)*coef(0,:))
      x1(2)=sum(zbuf(2,bufind)*coef(0,:))
      x1(3)=sum(zbuf(3,bufind)*coef(0,:))
    else
      x1 = z(1:3)
    end if

    tol = 1d-8
    call hybrd1 (f_sympl_euler2, n1, x1, fvec1, tol, ierr)
    if(ierr > 1) stop 'error in root finding'
    z(1:3) = x1

    call eval_field(z(1), z(2), z(3), 0)
    call get_derivatives(z(4))

    pthold = pthold - dt*(dH(2) - dH(1)*dpth(2)/dpth(1))
    z(4) = z(4) - dt*(dH(3) - dH(1)*dpth(3)/dpth(1))

    kbuf = mod(2*kt, nbuf) + 1
    zbuf(:,kbuf) = z

    ! Verlet part two: expl/impl Euler

    ! Initial guess with Lagrange extrapolation
    
    if (nlag>0 .and. kt>nlag) then
      do k=0, nlag
        bufind(k) = kbuf-2*nlag+2*k
        if (bufind(k)<1) bufind(k) = bufind(k) + nbuf 
      end do
      x2(1)=sum(zbuf(1,bufind)*coef(0,:))
      x2(2)=sum(zbuf(4,bufind)*coef(0,:))
!      print *, bufind
!      print *, x2, z(1), z(4)
    else
      x2(1)=z(1)
      x2(2)=z(4)
    end if

    tol = 1d-8
    call hybrd1 (f_sympl_euler1, n2, x2, fvec2, tol, ierr)
    if(ierr > 1) stop 'error in root finding'
    z(1) = x2(1)
    z(4) = x2(2)
!    print *, z(1), z(4)

    call eval_field(z(1), z(2), z(3), 0)
    call get_derivatives(z(4))
    z(2) = z(2) + dt*dH(1)/dpth(1)
    z(3) = z(3) + dt*(vpar*f%Bmod - dH(1)/dpth(1)*f%Bth)/f%Bph

    kbuf = mod(2*kt+1, nbuf) + 1
    zbuf(:,kbuf) = z
    kt = kt+1
    
    tau2 = tau2 + dtaumin
  enddo
  z0(1:3) = z(1:3)
  z0(4) = H
  z0(5) = vpar/(pabs*dsqrt(2d0))  ! alambda

end subroutine orbit_timestep_sympl_verlet



subroutine debug_root(x0)
  double precision :: x0(2), x(2)
  integer :: k, l, iflag
  integer, parameter :: n = 100
  double precision, parameter :: eps = 1e-15

  double precision :: fvec(2)

  do k = -n,n
    do l = -n,n
      x = x0 + l*eps/n*(/x0(1),0d0/) + k*eps/n*(/0d0,x0(2)/) 
      call f_sympl_euler1(2, x, fvec, iflag)
      write(5001,*) x(1), x(2), fvec
    end do
  end do

end subroutine debug_root

end module orbit_symplectic


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE plag_coeff(npoi,nder,x,xp,coef)
  !
  ! npoi - number of points (determines the order of Lagrange
  ! polynomial
  ! which is equal npoi-1)
  ! nder - number of derivatives computed 0 - function only, 1 - first
  ! derivative
  ! x - actual point where function and derivatives are evaluated
  ! xp(npoi) - array of points where function is known
  ! coef(0:nder,npoi) - weights for computation of function and
  ! derivatives,
  ! f=sum(fun(1:npoi)*coef(0,1:npoi) gives the function value
  ! df=sum(fun(1:npoi)*coef(1,1:npoi) gives the derivative value value
  !
  implicit none
  !
  INTEGER, INTENT(in)                                :: npoi,nder
  double precision, INTENT(in)                          :: x
  double precision, DIMENSION(npoi), INTENT(in)         :: xp
  double precision, DIMENSION(0:nder,npoi), INTENT(out) :: coef
  double precision, DIMENSION(:), ALLOCATABLE           :: dummy
  !
  INTEGER                                            :: i,k,j
  double precision                                      :: fac
  !
  DO i=1,npoi
      coef(0,i)=1.d0
      DO k=1,npoi
        IF(k.EQ.i) CYCLE
        coef(0,i)=coef(0,i)*(x-xp(k))/(xp(i)-xp(k))
      ENDDO
  ENDDO
  !
  IF(nder.EQ.0) RETURN
  !
  ALLOCATE(dummy(npoi))
  !
  DO i=1,npoi
      dummy=1.d0
      dummy(i)=0.d0
      DO k=1,npoi
        IF(k.EQ.i) CYCLE
        fac=(x-xp(k))/(xp(i)-xp(k))
        DO j=1,npoi
            IF(j.EQ.k) THEN
              dummy(j)=dummy(j)/(xp(i)-xp(k))
            ELSE
              dummy(j)=dummy(j)*fac
            ENDIF
        ENDDO
      ENDDO
      coef(1,i)=SUM(dummy)
  ENDDO
  !
  DEALLOCATE(dummy)
  !
  RETURN
END SUBROUTINE plag_coeff

