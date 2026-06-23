program export_boozer_chartmap_app
!> Export SIMPLE's VMEC->Boozer transform as a Boozer chartmap NetCDF.
!>
!> Reads field configuration from simple.in (netcdffile, ns_s, ns_tp,
!> multharm, axis healing), builds the internal Boozer field exactly as a
!> tracing run does, and writes the chartmap in the current schema (rho and
!> s coordinates, A_phi on the s abscissa, B_theta/B_phi on rho, Bmod on the
!> rho/theta/zeta grid). Lets SIMPLE's own field be compared against an
!> external booz_xform chartmap on equal footing.
!>
!> Usage: export_boozer_chartmap.x <out.nc> [config.in]   (default simple.in)

use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
    integmode, params_init, isw_field_type
use simple, only: tracer_t
use simple_main, only: init_field
use boozer_chartmap, only: export_boozer_chartmap
use timing, only: init_timer, print_phase_time

implicit none

integer, parameter :: BOOZER_FIELD_TYPE = 2
character(1024) :: out_file
character(256) :: config_file
type(tracer_t) :: norb

call init_timer()

if (command_argument_count() < 1) then
    write(*, '(A)') 'Usage: export_boozer_chartmap.x <out.nc> [config.in]'
    error stop 'missing output file argument'
end if
call get_command_argument(1, out_file)
if (command_argument_count() >= 2) then
    call get_command_argument(2, config_file)
else
    config_file = 'simple.in'
end if

call read_config(config_file)
call print_phase_time('Configuration reading completed')

if (isw_field_type /= BOOZER_FIELD_TYPE) then
    write(*, '(A,I0)') 'ERROR: export_boozer_chartmap needs isw_field_type = 2 &
        &(Boozer); got ', isw_field_type
    error stop 'incorrect field type for Boozer chartmap export'
end if

call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
call params_init
call print_phase_time('Field initialization completed')

call export_boozer_chartmap(out_file)
call print_phase_time('Chartmap export completed')
write(*, '(A)') 'Written '//trim(out_file)

end program export_boozer_chartmap_app
