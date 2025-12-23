program test_field_can_albert

use, intrinsic :: iso_fortran_env, only: dp => real64

use simple, only: tracer_t
use simple_main, only: init_field
use magfie_sub, only: ALBERT
use velo_mod, only: isw_field_type
use field, only: vmec_field_t, create_vmec_field
use field_can_albert, only: init_albert

implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

type(tracer_t) :: norb
type(vmec_field_t) :: magfie

isw_field_type = ALBERT
call create_vmec_field(magfie)

print *, 'init_field'
call init_field(norb, 'wout.nc', 5, 5, 3, 0)

end program test_field_can_albert
