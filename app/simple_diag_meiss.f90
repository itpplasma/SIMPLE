program diag_meiss_main
!> Diagnostic application for field_can_meiss analysis
!> 
!> Reads configuration from simple.in (default) or specified file,
!> initializes the field, and generates diagnostic plots for canonical coordinates

use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, integmode, params_init, isw_field_type
use simple, only: Tracer, init_vmec
use timing, only: init_timer, print_phase_time
use diag_meiss, only: plot_rh_can_vs_rc
use field_can_mod, only: field_can_from_id
use field, only: field_from_file, vmec_field_t
use field_can_meiss, only: init_transformation_arrays

implicit none

character(256) :: config_file
type(Tracer) :: norb

! Initialize timing
call init_timer()

! Read configuration file name from command line arguments
if (command_argument_count() == 0) then
    config_file = 'simple.in'
else
    call get_command_argument(1, config_file)
end if
call print_phase_time('Command line parsing completed')

! Initialize the system following simple.x pattern BUT stop before init_field_can
call read_config(config_file)
call print_phase_time('Configuration reading completed')

! Call init_vmec directly (like init_field does) but skip init_field_can
call init_vmec(netcdffile, ns_s, ns_tp, multharm, norb%fper)
call print_phase_time('VMEC initialization completed')

norb%integmode = integmode
call print_phase_time('Integration mode set')

! Initialize field_can system up to the point before expensive computation
if (norb%integmode >= 0) then
    call field_can_from_id(isw_field_type, vmec_field_t())
    call print_phase_time('Field canonical setup completed')
    
    ! Initialize only the transformation arrays (without expensive computation)
    call init_transformation_arrays()
    call print_phase_time('Transformation arrays initialized')
end if

call params_init
call print_phase_time('Parameter initialization completed')

print *, "Generating diagnostic plots..."

! Generate plots for different grid indices
print *, "Creating plot for i_th=1, i_phi=1..."
call plot_rh_can_vs_rc(1, 1, "diag_meiss_1_1.pdf")

print *, "Creating plot for i_th=2, i_phi=2..."  
call plot_rh_can_vs_rc(1, 2, "diag_meiss_1_2.pdf")

print *, "Creating plot for i_th=2, i_phi=2..."  
call plot_rh_can_vs_rc(2, 1, "diag_meiss_2_1.pdf")

print *, "Creating plot for i_th=2, i_phi=2..."
call plot_rh_can_vs_rc(2, 2, "diag_meiss_2_2.pdf")

print *, "Diagnostic plots completed successfully!"
print *, "Generated files: diag_meiss.pdf, diag_meiss_1_1.pdf, diag_meiss_2_2.pdf"

call print_phase_time('Diagnostic analysis completed')

end program diag_meiss_main