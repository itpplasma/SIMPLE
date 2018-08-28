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

