#ifndef _OPENMP
  function omp_get_thread_num()
    integer :: omp_get_thread_num
    omp_get_thread_num = 0
  end function omp_get_thread_num
#endif

module common

implicit none

double precision, parameter  :: pi=3.14159265358979d0
double precision, parameter  :: c=2.9979d10
double precision, parameter  :: e_charge=4.8032d-10
double precision, parameter  :: e_mass=9.1094d-28
double precision, parameter  :: p_mass=1.6726d-24
double precision, parameter  :: ev=1.6022d-12

contains

! From: http://fortranwiki.org/fortran/show/newunit
! This is a simple function to search for an available unit.
! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
! The UNIT value is returned by the function, and also by the optional
! argument. This allows the function to be used directly in an OPEN
! statement, and optionally save the result in a local variable.
! If no units are available, -1 is returned.
integer function newunit(unit)
  integer, intent(out), optional :: unit
! local
  integer, parameter :: LUN_MIN=10, LUN_MAX=1000
  logical :: opened
  integer :: lun
! begin
  newunit=-1
  do lun=LUN_MIN,LUN_MAX
    inquire(unit=lun,opened=opened)
    if (.not. opened) then
      newunit=lun
      exit
    end if
  end do
  if (present(unit)) unit=newunit
end function newunit  

end module common


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