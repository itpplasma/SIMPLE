program test_field_can_transforms
use, intrinsic :: iso_fortran_env, only: dp => real64

use magfie_sub, only : TEST, CANFLUX, BOOZER, MEISS, ALBERT
use field_can_mod, only: init_field_can, can_to_ref, ref_to_can, twopi
use simple, only: init_vmec

implicit none

real(dp), parameter :: TOL = 1d-13

real(dp) :: fper
real(dp) :: xstart(3), xref(3), xcan(3)

call init_vmec('wout.nc', 5, 5, 5, fper)
call init_field_can(BOOZER)

xstart = [0.3, 10.4, 20.5]
xref = xstart

call ref_to_can(xref, xcan)
call can_to_ref(xcan, xref)

if (any(abs(mod(xstart - xref, twopi)) > TOL)) then
    error stop "test_field_can_transforms failed"
end if

end program test_field_can_transforms
