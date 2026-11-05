#ifndef _OPENMP
  function omp_get_thread_num()
      integer :: omp_get_thread_num
      omp_get_thread_num = 0
  end function omp_get_thread_num
#endif

  module util

      implicit none

      double precision, parameter  :: pi = 3.141592653589793d0
      double precision, parameter  :: twopi = 6.283185307179586d0
! CODATA 2018 in Gaussian units; c and e are exact in the 2019 SI.
      double precision, parameter  :: c = 2.99792458d10
      double precision, parameter  :: e_charge = 4.803204712570263d-10
      double precision, parameter  :: e_mass = 9.1093837015d-28
      double precision, parameter  :: p_mass = 1.67262192369d-24
      double precision, parameter  :: ev = 1.602176634d-12
      double precision, parameter  :: hbar = 1.054571817d-27
      double precision, parameter  :: sqrt2 = dsqrt(2.0d0)

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
          integer, parameter :: LUN_MIN = 10, LUN_MAX = 1000
          logical :: opened
          integer :: lun
! begin
          newunit = -1
          do lun = LUN_MIN, LUN_MAX
              inquire (unit=lun, opened=opened)
              if (.not. opened) then
                  newunit = lun
                  exit
              end if
          end do
          if (present(unit)) unit = newunit
      end function newunit

  end module util
