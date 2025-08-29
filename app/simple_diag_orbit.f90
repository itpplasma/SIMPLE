program diag_orbit_main
!> Diagnostic application for orbit trajectory analysis
!> 
!> Reads configuration from simple.in (default) or specified file,
!> validates the integration setup, initializes the field exactly like simple.x,
!> and provides detailed orbit trajectory plotting for the Nth particle

use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
    integmode, params_init, isw_field_type, dtaumin, relerr, ntestpart, ntimstep, ntau, &
    zstart, startmode, grid_density, special_ants_file, reuse_batch, num_surf, sbeg
use simple, only: Tracer, init_sympl
use simple_main, only: init_field
use magfie_sub, only: init_magfie, VMEC
use samplers, only: init_starting_surf, sample, START_FILE
use timing, only: init_timer, print_phase_time
use diag_orbit, only: integrate_orbit_with_trajectory_debug
use orbit_symplectic_base, only: SymplecticIntegrator
use field_can_mod, only: FieldCan, get_val, eval_field => evaluate

implicit none

! Define dp kind parameter
integer, parameter :: dp = kind(1.0d0)

character(256) :: config_file, particle_arg
type(Tracer) :: norb
type(SymplecticIntegrator) :: si
type(FieldCan) :: field_can
integer :: particle_number
real(dp), dimension(5) :: z

! Initialize timing
call init_timer()

! Read configuration file name from command line arguments
if (command_argument_count() == 0) then
    config_file = 'simple.in'
    particle_number = 1  ! Default to first particle
elseif (command_argument_count() == 1) then
    call get_command_argument(1, particle_arg)
    read(particle_arg, *) particle_number
    config_file = 'simple.in'
elseif (command_argument_count() == 2) then
    call get_command_argument(1, config_file)
    call get_command_argument(2, particle_arg)
    read(particle_arg, *) particle_number
else
    print *, 'Usage: ./diag_orbit.x [config_file] [particle_number]'
    print *, '  or:  ./diag_orbit.x [particle_number]'
    print *, 'Example: ./diag_orbit.x simple.in 2'
    print *, '         ./diag_orbit.x 3'
    stop
end if
call print_phase_time('Command line parsing completed')

! Initialize the system following simple.x pattern
call read_config(config_file)
call print_phase_time('Configuration reading completed')

! Validate particle number against ntestpart
if (particle_number < 1 .or. particle_number > ntestpart) then
    print *, 'ERROR: Invalid particle number!'
    print '(A,I0)', 'Requested particle: ', particle_number
    print '(A,I0)', 'Available particles (ntestpart): ', ntestpart
    print *, 'Please adjust particle number or ntestpart in config file.'
    error stop 'Invalid particle number for orbit trajectory diagnostic'
end if
call print_phase_time('Particle number validation completed')

! Use the complete field initialization from simple_main
call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
call print_phase_time('Complete field initialization completed')

call params_init
call print_phase_time('Parameter initialization completed')

! Initialize VMEC magnetic field (required for sampling)
call init_magfie(VMEC)
call print_phase_time('VMEC magnetic field initialization completed')

! Initialize starting surfaces (required before sampling)
call init_starting_surf
call print_phase_time('Starting surface initialization completed')

! Perform particle sampling exactly like simple_main.f90
if (1 == startmode) then
    if ((0d0 < grid_density) .and. (1d0 > grid_density)) then
        call sample(zstart, grid_density)
    else
        call sample(zstart)
    endif
elseif (2 == startmode) then
    call sample(zstart, START_FILE)
elseif (3 == startmode) then
    call sample(special_ants_file)
elseif (4 == startmode) then
    call sample(zstart, reuse_batch)
elseif (5 == startmode) then
    if (0 == num_surf) then
        call sample(zstart, 0.0d0, 1.0d0)
    elseif (1 == num_surf) then
        call sample(zstart, 0.0d0, sbeg(1))
    elseif (2 == num_surf) then
        call sample(zstart, sbeg(1), sbeg(num_surf))
    else
        print *, 'Invalid surface range for volume sample defined.'
        error stop 'Invalid surface range for volume sample'
    endif
else
    print *, 'Invalid startmode: ', startmode
    error stop 'Invalid startmode'
endif
call print_phase_time('Particle sampling completed')

print *, "Orbit Trajectory Diagnostic Program"
print *, "==================================="
print *
print *, "This program provides detailed orbit trajectory plotting"
print *, "for individual particles using real physics integration."
print *
print '(A,A)', "Configuration file: ", trim(config_file)
print '(A,I0,A,I0)', "Selected particle: ", particle_number, " out of ", ntestpart
print '(A,I0)', "Field type (isw_field_type): ", isw_field_type
print '(A,I0)', "Integration mode: ", integmode
print '(A,I0,A,I0,A,I0,A)', "Integration: ", ntimstep, " macrosteps × ", ntau, " substeps = ", ntimstep*ntau, " total timesteps"
print *

! Initialize default starting conditions (will be overridden by sampling)
z(1) = 0.6_dp      ! s = 0.6
z(2) = 0.0_dp      ! theta = 0
z(3) = 0.0_dp      ! phi = 0
z(4) = 0.0_dp      ! pphi = 0 (trapped orbit with vpar=0)  
z(5) = 0.0_dp      ! lambda = 0 (vpar=0)

! Initialize symplectic integrator EXACTLY like simple_main.f90
if (integmode > 0) call init_sympl(si, field_can, z, dtaumin, dtaumin, relerr, integmode)

print '(A,ES12.5)', 'dtaumin (integration time step): ', si%dt
print '(A,I0)', 'ntau (substeps per dtau): ', si%ntau
print '(A,ES12.5)', 'Absolute tolerance: ', si%atol
print '(A,ES12.5)', 'Relative tolerance: ', si%rtol
print *

! Perform orbit trajectory diagnostic integration
call integrate_orbit_with_trajectory_debug(si, field_can, particle_number)

call print_phase_time('Orbit trajectory diagnostic analysis completed')

print *
print *, "Orbit Trajectory Diagnostic completed successfully!"
print *, "Generated detailed trajectory plots showing real particle"
print *, "motion in the magnetic field using symplectic integration."

end program diag_orbit_main