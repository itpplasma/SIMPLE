program test_field_can_albert

use, intrinsic :: iso_fortran_env, only: dp => real64

use simple, only: tracer_t
use simple_main, only: init_field
use magfie_sub, only: ALBERT
use velo_mod, only: isw_field_type
use field, only: vmec_field_t
use field_can_albert, only: init_albert

implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

type(tracer_t) :: norb
class(vmec_field_t), allocatable :: magfie

isw_field_type = ALBERT
magfie = vmec_field_t()

print *, 'init_field'
call init_field(norb, 'wout.nc', 5, 5, 3, 0)

end program test_field_can_albert
