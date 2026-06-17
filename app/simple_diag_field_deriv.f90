program diag_field_deriv
!> Dump the Boozer field and its radial derivatives at fixed (s, theta_B, phi_B)
!> points, using whichever field the config selects (internal get_boozer_coordinates
!> or a loaded booz_xform chartmap). The symplectic step consumes the radial
!> derivatives (dAph, dhth, dhph, dBmod and the second derivatives); comparing them
!> between the two eval paths at identical points isolates whether the
!> internal-vs-chartmap symplectic loss difference comes from the field
!> derivatives the two spline constructions produce (#398). The field VALUES are
!> known to agree to <0.15%; this checks the derivatives.
!>
!> Usage: ./diag_field_deriv.x [config_file] [out_file]

use, intrinsic :: iso_fortran_env, only: dp => real64
use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, integmode
use simple, only: tracer_t
use simple_main, only: init_field
use field_can_mod, only: field_can_t, eval_field => evaluate

implicit none

character(256) :: config_file, out_file, arg
type(tracer_t) :: norb
type(field_can_t) :: f
integer :: i, unit, nargs
real(dp) :: s, theta, phi, lo, hi
integer, parameter :: ns = 60

config_file = 'simple.in'
out_file = 'field_deriv.dat'
nargs = command_argument_count()
if (nargs >= 1) call get_command_argument(1, config_file)
if (nargs >= 2) call get_command_argument(2, out_file)

call read_config(config_file)
call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)

theta = 0.3_dp
phi = 0.2_dp
lo = log(1.0e-6_dp); hi = log(0.3_dp)

open(newunit=unit, file=trim(out_file), status='replace')
write(unit, '(A)') '# s  Ath  Aph  dAph_ds  hth  dhth_ds  hph  dhph_ds  '// &
    'Bmod  dBmod_ds  d2Bmod_ds2  d2hth_ds2  d2hph_ds2'
do i = 1, ns
    s = exp(lo + (hi - lo)*real(i - 1, dp)/real(ns - 1, dp))
    call eval_field(f, s, theta, phi, 2)   ! mode_secders=2: all second derivatives
    write(unit, '(13ES16.8)') s, f%Ath, f%Aph, f%dAph(1), f%hth, f%dhth(1), &
        f%hph, f%dhph(1), f%Bmod, f%dBmod(1), f%d2Bmod(1), f%d2hth(1), f%d2hph(1)
end do
close(unit)
print '(A,A)', 'field+derivatives written: ', trim(out_file)
print '(A,I0)', 'field type (isw_field_type via config), integmode=', integmode

end program diag_field_deriv
