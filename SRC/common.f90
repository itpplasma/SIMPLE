#ifndef _OPENMP
  function omp_get_thread_num()
    integer :: omp_get_thread_num
    omp_get_thread_num = 0
  end function omp_get_thread_num
#endif

module common

implicit none

double precision, parameter  :: pi=3.14159265358979d0
double precision, parameter  :: twopi=6.28318530717958d0
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

