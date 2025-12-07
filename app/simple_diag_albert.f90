program diag_albert_main
!> Diagnostic application for Albert canonical coordinate analysis
!> 
!> Reads configuration from simple.in (default) or specified file,
!> validates that coordinates are set to Albert, initializes the field,
!> and generates contour plots of Aph_of_xc and Bmod_of_xc

use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
    integmode, params_init, isw_field_type
use simple, only: tracer_t
use simple_main, only: init_field
use timing, only: init_timer, print_phase_time
use diag_albert, only: plot_albert_contours

implicit none

character(256) :: config_file
type(tracer_t) :: norb
integer, parameter :: ALBERT_FIELD_TYPE = 4

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

! Validate that coordinates are set to Albert
if (isw_field_type /= ALBERT_FIELD_TYPE) then
    write(*, '(A)') 'ERROR: Albert coordinate diagnostic requires Albert coordinates!'
    write(*, '(A,I0)') 'Current isw_field_type = ', isw_field_type
    write(*, '(A)') 'Please set isw_field_type = 4 in simple.in for Albert coordinates'
    write(*, '(A)') 'Field type mapping:'
    write(*, '(A)') '  -1: Testing'
    write(*, '(A)') '   0: Canonical (flux)'
    write(*, '(A)') '   1: VMEC'
    write(*, '(A)') '   2: Boozer'
    write(*, '(A)') '   3: Meiss'
    write(*, '(A)') '   4: Albert'
    error stop 'Incorrect coordinate system for Albert diagnostic'
end if
call print_phase_time('Albert coordinate validation completed')

! Use the complete field initialization from simple_main
! This properly handles both netcdffile and field_input parameters
call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
call print_phase_time('Complete field initialization completed')

call params_init
call print_phase_time('Parameter initialization completed')

print *, "Generating Albert coordinate diagnostic plots..."
print *, "This will create line plots of Aph_of_xc and Bmod_of_xc"
print *, "vs theta and phi at three radial slices (inner, middle, outer)..."

! Generate Albert coordinate diagnostic plots
call plot_albert_contours()

call print_phase_time('Albert diagnostic analysis completed')

end program diag_albert_main