program test_tokamak_alpha_diagnostic

    use, intrinsic :: iso_fortran_env, only : dp => real64
    use pyplot_module, only : pyplot, pyplot_wp
    use geoflux_coordinates, only : geoflux_to_cyl, init_geoflux_coordinates

    implicit none

    character(len=256), parameter :: template_cfg = '../../../examples/tokamak_alpha_confinement/simple.in'
    character(len=256), parameter :: diag_cfg = 'tokamak_alpha_diag.in'
    character(len=256), parameter :: simple_exec = '../../simple.x'
    character(len=256), parameter :: orbit_plot_png = 'tokamak_orbits_RZ.png'
    character(len=256), parameter :: orbit_plot_py = 'tokamak_orbits_RZ.py'

    real(dp), parameter :: tol_trace = 1.0d-8
    real(dp), parameter :: loss_fraction = 1.0d-3

    integer :: exit_code
    integer :: lost_particle, confined_particle
    real(dp) :: trace_time, loss_tolerance
    type(pyplot) :: plt
    real(dp), allocatable :: time_lost(:), R_lost(:), Z_lost(:)
    real(dp), allocatable :: time_conf(:), R_conf(:), Z_conf(:)

    call prepare_diagnostic_config(template_cfg, diag_cfg)
    call copy_file('../../../examples/tokamak_alpha_confinement/tokamak.in', 'tokamak.in')
    call copy_file('../../../examples/tokamak/EQDSK_I.geqdsk', 'EQDSK_I.geqdsk')

    call execute_command_line(trim(simple_exec)//' '//trim(diag_cfg), exitstat=exit_code)
    if (exit_code /= 0) then
        write(*,*) 'simple.x failed with exit code ', exit_code
        error stop 'Tokamak diagnostic example execution failed'
    end if

    call read_trace_time(diag_cfg, trace_time)
    loss_tolerance = max(trace_time * loss_fraction, tol_trace)

    call select_particles('times_lost.dat', trace_time, loss_tolerance, &
        lost_particle, confined_particle)

    if (lost_particle < 0) error stop 'No lost particle found for diagnostic plot'
    if (confined_particle < 0) error stop 'No confined or long-lived particle found for diagnostic plot'

    call init_geoflux_coordinates('EQDSK_I.geqdsk')

    call load_orbit(lost_particle, time_lost, R_lost, Z_lost)
    call load_orbit(confined_particle, time_conf, R_conf, Z_conf)

    call plt%initialize(grid=.true., xlabel='R [cm]', ylabel='Z [cm]', &
        title='Tokamak Poloidal Orbits', legend=.true., axis_equal=.true.)
    call plt%add_plot(real(R_conf, pyplot_wp), real(Z_conf, pyplot_wp), 'Confined', '-', linewidth=2)
    call plt%add_plot(real(R_lost, pyplot_wp), real(Z_lost, pyplot_wp), 'Lost', '--', linewidth=2)
    call plt%savefig(trim(orbit_plot_png), pyfile=trim(orbit_plot_py))
    call plt%destroy()

contains

    subroutine prepare_diagnostic_config(template_path, output_path)
        character(*), intent(in) :: template_path
        character(*), intent(in) :: output_path

        character(len=512) :: line, trimmed
        integer :: in_unit, out_unit, ios
        logical :: inserted

        inserted = .false.
        open(newunit=in_unit, file=trim(template_path), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'Unable to read tokamak template configuration'
        open(newunit=out_unit, file=trim(output_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) error stop 'Unable to create diagnostic configuration'

        do
            read(in_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            trimmed = adjustl(line)
            if (index(trimmed, 'output_orbits_macrostep') > 0) then
                write(out_unit, '(A)') 'output_orbits_macrostep = .True.'
                inserted = .true.
            else if (.not. inserted .and. trim(trimmed) == '/') then
                write(out_unit, '(A)') 'output_orbits_macrostep = .True.'
                write(out_unit, '(A)') trim(line)
                inserted = .true.
            else
                write(out_unit, '(A)') trim(line)
            end if
        end do

        if (.not. inserted) then
            write(out_unit, '(A)') 'output_orbits_macrostep = .True.'
            write(out_unit, '(A)') '/'
        end if

        close(in_unit)
        close(out_unit)
    end subroutine prepare_diagnostic_config


    subroutine read_trace_time(config_path, value)
        character(*), intent(in) :: config_path
        real(dp), intent(out) :: value

        character(len=512) :: line, working
        integer :: unit_cfg, ios_cfg, comment_pos, eq_pos, pos

        value = -1.0_dp
        open(newunit=unit_cfg, file=trim(config_path), status='old', action='read', iostat=ios_cfg)
        if (ios_cfg /= 0) error stop 'Unable to read diagnostic configuration'

        do
            read(unit_cfg, '(A)', iostat=ios_cfg) line
            if (ios_cfg /= 0) exit
            working = adjustl(line)
            comment_pos = index(working, '!')
            if (comment_pos > 0) working = working(:comment_pos-1)
            pos = index(working, 'trace_time')
            if (pos > 0) then
                eq_pos = index(working, '=')
                if (eq_pos > 0 .and. eq_pos < len_trim(working)) then
                    read(working(eq_pos+1:), *, iostat=ios_cfg) value
                    if (ios_cfg == 0) exit
                end if
            end if
        end do

        close(unit_cfg)
        if (value <= 0.0_dp) error stop 'trace_time not found in diagnostic configuration'
    end subroutine read_trace_time


    subroutine select_particles(filename, trace_time, loss_tolerance, lost_idx, confined_idx)
        character(*), intent(in) :: filename
        real(dp), intent(in) :: trace_time, loss_tolerance
        integer, intent(out) :: lost_idx, confined_idx

        integer :: unit, ios, particle_index
        real(dp) :: t_lost, trap_par_val, s_initial, perp_inv_val
        real(dp) :: zend(5)
        real(dp) :: best_survival_time
        integer :: best_survival_idx

        lost_idx = -1
        confined_idx = -1
        best_survival_time = -1.0_dp
        best_survival_idx = -1

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'times_lost.dat not produced for diagnostic'

        do
            read(unit, *, iostat=ios) particle_index, t_lost, trap_par_val, &
                s_initial, perp_inv_val, zend
            if (ios /= 0) exit
            if (t_lost > loss_tolerance .and. t_lost < trace_time - loss_tolerance) then
                if (lost_idx < 0) lost_idx = particle_index
            else if (abs(t_lost - trace_time) <= loss_tolerance) then
                if (confined_idx < 0) confined_idx = particle_index
            end if
            if (t_lost > best_survival_time) then
                best_survival_time = t_lost
                best_survival_idx = particle_index
            end if
            if (lost_idx >= 0 .and. confined_idx >= 0) exit
        end do

        close(unit)
        if (confined_idx < 0) confined_idx = best_survival_idx
    end subroutine select_particles


    subroutine load_orbit(particle_index, time, R, Z)
        integer, intent(in) :: particle_index
        real(dp), allocatable, intent(out) :: time(:), R(:), Z(:)

        character(len=32) :: filename
        integer :: unit, ios, count, i
        real(dp) :: t, s, theta, phi, aux1, aux2
        real(dp) :: coords(3), cyl(3)

        write(filename, '(A,I0)') 'fort.', 90000 + particle_index
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'Orbit diagnostic file missing for selected particle'

        count = 0
        do
            read(unit, *, iostat=ios) t, s, theta, phi, aux1, aux2
            if (ios /= 0) exit
            count = count + 1
        end do
        if (count <= 1) error stop 'Insufficient orbit samples for diagnostic plot'

        rewind(unit)
        allocate(time(count), R(count), Z(count))
        do i = 1, count
            read(unit, *, iostat=ios) time(i), s, theta, phi, aux1, aux2
            if (ios /= 0) exit
            coords = [s, theta, phi]
            call geoflux_to_cyl(coords, cyl)
            R(i) = cyl(1)
            Z(i) = cyl(3)
        end do
        close(unit)
    end subroutine load_orbit


    subroutine copy_file(source_path, dest_path)
        character(*), intent(in) :: source_path
        character(*), intent(in) :: dest_path

        integer :: in_unit, out_unit, ios
        character(len=512) :: line

        open(newunit=in_unit, file=trim(source_path), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'Unable to read tokamak parameter file for diagnostic'
        open(newunit=out_unit, file=trim(dest_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) error stop 'Unable to create tokamak parameter file copy'

        do
            read(in_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            write(out_unit, '(A)') trim(line)
        end do

        close(in_unit)
        close(out_unit)
    end subroutine copy_file

end program test_tokamak_alpha_diagnostic
