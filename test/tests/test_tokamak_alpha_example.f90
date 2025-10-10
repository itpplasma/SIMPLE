program test_tokamak_alpha_example

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    character(len=256) :: simple_exec
    character(len=256) :: example_cfg
    character(len=512) :: command
    integer :: exit_code
    integer :: unit, ios, particle_index
    real(dp) :: t_lost, trap_par_val, s_initial, perp_inv_val
    real(dp) :: zend(5)
    real(dp) :: trace_time, loss_tolerance

    simple_exec = '../../simple.x'
    example_cfg = '../../../examples/tokamak_alpha_confinement/simple.in'
    command = trim(simple_exec)//' '//trim(example_cfg)
    call execute_command_line(command, exitstat=exit_code)
    if (exit_code /= 0) then
        write(*,*) 'simple.x failed with exit code ', exit_code
        error stop 'Tokamak example execution failed'
    end if

    trace_time = read_trace_time(example_cfg)
    loss_tolerance = max(trace_time*1.0e-3_dp, 1.0e-8_dp)

    open(newunit=unit, file='times_lost.dat', status='old', action='read', &
        iostat=ios)
    if (ios /= 0) then
        error stop 'times_lost.dat not produced by example run'
    end if

    do
        read(unit, *, iostat=ios) particle_index, t_lost, trap_par_val, &
            s_initial, perp_inv_val, zend
        if (ios /= 0) exit
        if (t_lost > loss_tolerance .and. t_lost < trace_time - loss_tolerance) then
            write(*,*) 'Particle ', particle_index, ' lost at t = ', t_lost
            error stop 'Alpha confinement example should have no losses'
        end if
    end do

    close(unit, status='delete')

    call cleanup_output('confined_fraction.dat')
    call cleanup_output('avg_inverse_t_lost.dat')

contains

    function read_trace_time(config_path) result(value)
        character(*), intent(in) :: config_path
        real(dp) :: value
        character(len=512) :: line
        character(len=512) :: working
        integer :: unit_cfg, ios_cfg, comment_pos, eq_pos, pos

        value = -1.0_dp
        open(newunit=unit_cfg, file=trim(config_path), status='old', action='read', &
            iostat=ios_cfg)
        if (ios_cfg /= 0) then
            error stop 'Unable to read example configuration'
        end if

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

        if (value <= 0.0_dp) then
            error stop 'trace_time not found in configuration'
        end if
    end function read_trace_time

    subroutine cleanup_output(filename)
        character(*), intent(in) :: filename
        integer :: unit_local, status_local

        open(newunit=unit_local, file=filename, status='old', &
            action='read', iostat=status_local)
        if (status_local == 0) then
            close(unit_local, status='delete')
        end if
    end subroutine cleanup_output

end program test_tokamak_alpha_example
