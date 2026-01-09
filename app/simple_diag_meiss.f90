program diag_meiss_main
!> Diagnostic application for field_can_meiss analysis
!>
!> Reads configuration from simple.in (default) or specified file,
!> initializes the field (including chartmap if specified), and generates
!> diagnostic plots for canonical coordinates

use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, integmode, &
                  params_init, isw_field_type, coord_input, field_input
use simple, only: tracer_t, init_vmec
use timing, only: init_timer, print_phase_time
use diag_meiss, only: plot_rh_can_vs_rc
use field_can_mod, only: field_can_from_id
use field_base, only: magnetic_field_t
use field, only: field_from_file, vmec_field_t, create_vmec_field
use field_splined, only: splined_field_t
use field_can_meiss, only: init_transformation_arrays
use reference_coordinates, only: init_reference_coordinates

implicit none

character(256) :: config_file
type(tracer_t) :: norb
type(vmec_field_t) :: vmec_field
class(magnetic_field_t), allocatable :: field_temp

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

! Initialize VMEC
call init_vmec(netcdffile, ns_s, ns_tp, multharm, norb%fper)
call print_phase_time('VMEC initialization completed')

! Initialize reference coordinates (chartmap) if specified
call init_reference_coordinates(coord_input)
call print_phase_time('Reference coordinate system initialization completed')

norb%integmode = integmode

! Initialize field_can system
if (norb%integmode >= 0) then
    if (trim(coord_input) /= '') then
        ! Chartmap case: load field from coils and use splined field
        call field_from_file(field_input, field_temp)
        call print_phase_time('Field from file loading completed')

        select type (fld => field_temp)
        type is (splined_field_t)
            call field_can_from_id(isw_field_type, fld)
        class default
            print *, "Warning: expected splined_field_t for chartmap case"
            call create_vmec_field(vmec_field)
            call field_can_from_id(isw_field_type, vmec_field)
        end select
    else
        ! Pure VMEC case
        call create_vmec_field(vmec_field)
        call field_can_from_id(isw_field_type, vmec_field)
    end if
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