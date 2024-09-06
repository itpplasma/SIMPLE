program test_field_can_albert

use, intrinsic :: iso_fortran_env, only: dp => real64

use simple, only: Tracer, init_field
use magfie_sub, only: ALBERT
use velo_mod, only: isw_field_type
use field, only: VmecField
use field_can_albert, only: init_albert

implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

type(Tracer) :: norb
class(VmecField), allocatable :: magfie

isw_field_type = ALBERT
magfie = VmecField()

print *, 'init_field'
call init_field(norb, 'wout.nc', 5, 5, 3, 0)

end program test_field_can_albert
