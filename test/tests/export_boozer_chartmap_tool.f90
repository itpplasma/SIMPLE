program export_boozer_chartmap_tool
    !> Tool to export a Boozer chartmap from VMEC and convert start.dat
    !> from VMEC reference coordinates to Boozer/chartmap reference coordinates.
    !>
    !> Usage: export_boozer_chartmap_tool.x <wout.nc> <chartmap.nc> <start_vmec.dat> <start_boozer.dat>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: netcdffile, multharm, ns_A, ns_s, ns_tp, nper
    use parmot_mod, only: rmu
    use velo_mod, only: isw_field_type
    use boozer_coordinates_mod, only: use_B_r
    use boozer_sub, only: get_boozer_coordinates, vmec_to_boozer, &
        export_boozer_chartmap
    use spline_vmec_sub, only: spline_vmec_data
    use vmecin_sub, only: stevvo

    implicit none

    real(dp), parameter :: twopi = 2.0_dp * 3.14159265358979_dp
    character(len=1024) :: wout_file, chartmap_file, start_vmec, start_boozer
    character(len=1024) :: line
    integer :: nargs, ipart, npart, ios
    real(dp) :: s, theta_v, phi_v, v, lam, theta_b, phi_b
    real(dp) :: RT0, R0i, cbfi, bz0i, bf0, fper, rho
    integer :: L1i

    nargs = command_argument_count()
    if (nargs /= 4) then
        print *, 'Usage: export_boozer_chartmap_tool.x <wout.nc> <chartmap.nc>', &
                 ' <start_vmec.dat> <start_boozer.dat>'
        error stop
    end if

    call get_command_argument(1, wout_file)
    call get_command_argument(2, chartmap_file)
    call get_command_argument(3, start_vmec)
    call get_command_argument(4, start_boozer)

    ! Initialize VMEC
    netcdffile = wout_file
    multharm = 3
    ns_A = 5
    ns_s = 5
    ns_tp = 5
    isw_field_type = 2
    rmu = 1.0e8_dp

    call spline_vmec_data
    call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0)
    fper = twopi / real(L1i, dp)

    ! Compute Boozer coordinates
    use_B_r = .false.
    call get_boozer_coordinates

    ! Export chartmap
    call export_boozer_chartmap(chartmap_file)

    ! Count particles in start_vmec
    npart = 0
    open(unit=10, file=trim(start_vmec), status='old', iostat=ios)
    if (ios /= 0) then
        print *, 'Cannot open ', trim(start_vmec)
        error stop
    end if
    do
        read(10, *, iostat=ios)
        if (ios /= 0) exit
        npart = npart + 1
    end do
    close(10)

    print *, 'Converting', npart, ' particles from VMEC to Boozer coords'

    ! Convert start.dat coordinates
    open(unit=10, file=trim(start_vmec), status='old')
    open(unit=11, file=trim(start_boozer), status='replace', recl=1024)

    do ipart = 1, npart
        read(10, *) s, theta_v, phi_v, v, lam

        ! Transform VMEC angles to Boozer angles
        call vmec_to_boozer(s, theta_v, phi_v, theta_b, phi_b)

        ! In chartmap reference coords: x(1) = rho = sqrt(s)
        rho = sqrt(max(s, 0.0_dp))

        write(11, *) rho, theta_b, phi_b, v, lam
    end do

    close(10)
    close(11)

    print *, 'Written ', trim(start_boozer)

end program export_boozer_chartmap_tool
