program test_magfie_coils

use, intrinsic :: iso_fortran_env, only: dp => real64
use simple, only : init_vmec
use magfie_sub, only : VMEC
use velo_mod, only: isw_field_type
use simple_magfie, only: VmecField, CoilsField, create_coils_field

implicit none

class(VmecField), allocatable :: vmec_field
class(CoilsField), allocatable :: coils_field
real(dp) :: dummy, x(3), Acov(3), hcov(3), Bmod

isw_field_type = VMEC

call init_vmec('wout.nc', 5, 5, 5, dummy)
allocate(vmec_field)
coils_field = create_coils_field('coils.5C')

x = [0.3d0, 0.2d0, 0.1d0]

print *, 'x = ', x

call vmec_field%evaluate(x, Acov, hcov, Bmod)
print *, 'vmec_field%evaluate'
print *, 'A = ', Acov
print *, 'h = ', hcov
print *, 'B = ', Bmod

call coils_field%evaluate_direct(x, Acov, hcov, Bmod)
print *, 'coils_field%evaluate_direct'
print *, 'A = ', Acov
print *, 'h = ', hcov
print *, 'B = ', Bmod

call coils_field%evaluate(x, Acov, hcov, Bmod)
print *, 'coils_field%evaluate'
print *, 'A = ', Acov
print *, 'h = ', hcov
print *, 'B = ', Bmod


end program test_magfie_coils
