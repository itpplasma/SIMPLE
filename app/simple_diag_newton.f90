program diag_newton_main
!> Diagnostic application for newton_midpoint function analysis
!> 
!> Reads configuration from simple.in (default) or specified file,
!> validates that integration mode is set to midpoint, initializes the field,
!> and provides detailed Newton iteration diagnostics for trapped orbit

use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
    integmode, params_init, isw_field_type, dtaumin, dtau, ntau, relerr
use simple, only: Tracer, init_sympl
use simple_main, only: init_field
use timing, only: init_timer, print_phase_time
use diag_newton, only: integrate_orbit_with_newton_debug
use orbit_symplectic_base, only: SymplecticIntegrator, MIDPOINT
use field_can_mod, only: FieldCan, get_val, eval_field => evaluate

implicit none

! Define dp kind parameter
integer, parameter :: dp = kind(1.0d0)

character(256) :: config_file
type(Tracer) :: norb
type(SymplecticIntegrator) :: si
type(FieldCan) :: field_can
integer, parameter :: NUM_DEBUG_STEPS = 10
real(dp), dimension(5) :: z

! Initialize timing
call init_timer()

! Read configuration file name from command line arguments
if (command_argument_count() == 0) then
    config_file = 'simple.in'
else
    call get_command_argument(1, config_file)
end if
call print_phase_time('Command line parsing completed')

! Initialize the system following simple.x pattern
call read_config(config_file)
call print_phase_time('Configuration reading completed')

! Validate that integration mode is set to midpoint
if (integmode /= MIDPOINT) then
    write(*, '(A)') 'ERROR: Newton midpoint diagnostic requires midpoint integration!'
    write(*, '(A,I0)') 'Current integmode = ', integmode
    write(*, '(A)') 'Please set integmode = 3 in simple.in for midpoint rule'
    write(*, '(A)') 'Integration mode mapping:'
    write(*, '(A)') '  0: RK45'
    write(*, '(A)') '  1: Explicit-Implicit Euler'
    write(*, '(A)') '  2: Implicit-Explicit Euler'
    write(*, '(A)') '  3: Midpoint (required for this diagnostic)'
    write(*, '(A)') '  4: GAUSS1'
    write(*, '(A)') '  5: GAUSS2'
    write(*, '(A)') '  6: GAUSS3'
    write(*, '(A)') '  7: GAUSS4'
    write(*, '(A)') '  15: LOBATTO3'
    error stop 'Incorrect integration mode for Newton midpoint diagnostic'
end if
call print_phase_time('Midpoint integration mode validation completed')

! Validate coordinate system
if (isw_field_type /= 4) then
    write(*, '(A)') 'WARNING: This diagnostic was designed for Albert coordinates (isw_field_type = 4)'
    write(*, '(A,I0)') 'Current isw_field_type = ', isw_field_type
    write(*, '(A)') 'The Newton iterations will still be analyzed but results may differ'
end if

! Use the complete field initialization from simple_main
call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
call print_phase_time('Complete field initialization completed')

call params_init
call print_phase_time('Parameter initialization completed')

print *, "Newton Midpoint Diagnostic Program"
print *, "=================================="
print *
print *, "This program provides detailed analysis of Newton iteration convergence"
print *, "during orbit integration using the midpoint rule."
print *
print *, "Focus: Trapped orbit at s=0.6, th=0, ph=0, vpar=0 using Midpoint method"
print '(A,I0,A)', "Integration: ", NUM_DEBUG_STEPS, " time steps"
print *

! Initialize trapped orbit coordinates exactly like simple_main.f90
z(1) = 0.6_dp      ! s = 0.6
z(2) = 0.0_dp      ! theta = 0
z(3) = 0.0_dp      ! phi = 0
z(4) = 0.0_dp      ! pphi = 0 (trapped orbit with vpar=0)  
z(5) = 0.0_dp      ! lambda = 0 (vpar=0)

! Initialize symplectic integrator EXACTLY like simple_main.f90
if (integmode > 0) call init_sympl(si, field_can, z, dtaumin, dtaumin, relerr, integmode)

! Transfer coordinates to integrator
si%z = z(1:4)

print '(A,ES12.5)', 'dtau (large time step): ', dtau
print '(A,ES12.5)', 'dtaumin (integration time step): ', si%dt
print '(A,I0)', 'ntau (substeps per dtau): ', si%ntau
print '(A,ES12.5)', 'Absolute tolerance: ', si%atol
print '(A,ES12.5)', 'Relative tolerance: ', si%rtol
print *

! Initialize field at starting position
call eval_field(field_can, si%z(1), si%z(2), si%z(3), 0)
call get_val(field_can, si%z(4))

print *, "Field initialization completed"
print *

! Perform diagnostic integration
call integrate_orbit_with_newton_debug(si, field_can, NUM_DEBUG_STEPS)

call print_phase_time('Newton midpoint diagnostic analysis completed')

print *
print *, "Newton Midpoint Diagnostic completed successfully!"
print *, "Analysis shows detailed convergence behavior of Newton iterations"
print *, "during symplectic orbit integration using the midpoint rule."

end program diag_newton_main