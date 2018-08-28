module orbit_symplectic

use parmot_mod, only : ro0_parmot => ro0
use field_can_mod, only: field_can, d_field_can, d2_field_can, eval_field, &
  get_val, get_derivatives, H, pth, vpar, dvpar, dH, dpth, f, df, d2f, ro0, mu

implicit none
save

double precision, dimension(2) :: qold, wold
double precision, dimension(2) :: q, w  ! q = (th, ph), w = (r, pphi)
double precision :: pthold

double precision :: dt

double precision :: coala
double precision :: derphi(3)
double precision :: alambd, pabs

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init(z)
  double precision, intent(in) :: z(5)

  call eval_field(z(1), z(2), z(3), 0)

  pabs=z(4)
  alambd=z(5)

  mu = .5d0*pabs**2*(1.d0-alambd**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
  ro0 = ro0_parmot/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
  vpar = pabs*alambd*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

  q = z(2:3)  ! theta, varphi
  w(1) = z(1) ! r
  w(2) = vpar*f%Bph/f%Bmod + f%Aph/ro0 ! p_ph

  pth = vpar*f%Bth/f%Bmod + f%Ath/ro0 ! p_th
end subroutine orbit_sympl_init

subroutine f_sympl_euler_2(n, x, fvec, iflag)
  integer, parameter :: mode = 2

  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(out) :: iflag

  call eval_field(x(1), x(2), x(3), 0)
  call get_derivatives(w(2))

  fvec(1) = pth - pthold
  fvec(2) = dpth(1)*(x(2) - qold(1)) - dt*dH(1)
  fvec(3) = dpth(1)*(x(3) - qold(2)) - dt*(dpth(1)*vpar*f%Bmod - dH(1)*f%Bth)/f%Bph

end subroutine f_sympl_euler_2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl(z, dtau, dtaumin, ierr)
!
!
  implicit none
  
  integer, parameter :: ndim = 5
  integer, parameter :: n = 3
  double precision :: tol

  integer, intent(out) :: ierr
  double precision, dimension(ndim), intent(inout) :: z
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  double precision, dimension(n) :: x, fvec

  double precision :: tau2
  
  ierr = 0

  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0d0
  do while(tau2.lt.dtau)
    qold = q
    wold = w
    pthold = pth
    ! TODO: initial guess with Lagrange
    tol = 1d-10
    x(1) = w(1)
    x(2:3) = q
    call hybrd1 (f_sympl_euler_2, n, x, fvec, tol, ierr)
    if(ierr > 1) stop 'error in root finding'
    call eval_field(x(1), x(2), x(3), 0)
    call get_derivatives(w(2))

    pth = pthold - dt*(dH(2) - dH(1)*dpth(2)/dpth(1))
    w(2) = wold(2) - dt*(dH(3) - dH(1)*dpth(3)/dpth(1))

    w(1) = x(1)
    q = x(2:3)
    tau2 = tau2 + dtaumin
  enddo
  z(1:3) = x
  z(4) = H
  z(5) = vpar/sqrt(H) ! alambda
end

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

