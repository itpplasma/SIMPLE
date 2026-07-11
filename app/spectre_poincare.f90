program spectre_poincare
    !> Symplectic field-line Poincare sections and rotational-transform profiles
    !> for a SPECTRE/SPEC equilibrium, one file per volume.
    !>
    !> Reads a namelist (default spectre_poincare.in), loads the SPECTRE HDF5
    !> equilibrium, seeds field lines per volume at theta = 0.2 with s spread
    !> across the volume, and integrates each line with the semi-implicit
    !> symplectic Euler map in field_line_spectre. Field periodicity puts every
    !> section exactly on zeta0 + k*2*pi/Nfp, so no interpolation is needed.
    !>
    !> Output per volume: poincare_vol<N>.dat (seed, section, s, theta, R, Z)
    !> and iota_vol<N>.dat (seed, fitted rotational transform).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spectre_reader, only: spectre_data_t, load_spectre, free_spectre
    use field_line_spectre, only: field_line_psi, field_line_step, spectre_rz

    implicit none

    real(dp), parameter :: pi = 3.141592653589793238462643383279_dp
    real(dp), parameter :: seed_theta = 0.2_dp
    real(dp), parameter :: seed_s_min = -0.85_dp
    real(dp), parameter :: seed_s_max = 0.85_dp
    real(dp), parameter :: s_terminate = 1.0_dp - 1.0e-10_dp
    integer, parameter :: max_volumes = 128

    character(len=1024) :: spectre_file
    real(dp) :: zeta0
    integer :: nturns, nsteps_per_period, nseeds_per_volume
    integer :: volumes(max_volumes)
    namelist /poincare/ spectre_file, zeta0, nturns, &
        nsteps_per_period, nseeds_per_volume, volumes

    type(spectre_data_t) :: data
    character(len=1024) :: config_file
    integer :: ierr, iv, lvol, n_terminated
    integer, allocatable :: vol_list(:)
    real(dp) :: h

    spectre_file = 'spectre.h5'
    zeta0 = 0.0_dp
    nturns = 200
    nsteps_per_period = 256
    nseeds_per_volume = 4
    volumes = -1

    config_file = 'spectre_poincare.in'
    if (command_argument_count() >= 1) call get_command_argument(1, config_file)

    call read_input(config_file)
    call load_spectre(trim(spectre_file), data, ierr)
    if (ierr /= 0) error stop 'spectre_poincare: load_spectre failed'

    call select_volumes(data%Mvol, volumes, vol_list)
    h = (2.0_dp*pi/real(data%Nfp, dp))/real(nsteps_per_period, dp)

    write (*, '(A,I0,A,I0,A,I0)') 'spectre_poincare: Nfp = ', data%Nfp, &
        ', Mvol = ', data%Mvol, ', tracing volumes: ', size(vol_list)

    n_terminated = 0
    do iv = 1, size(vol_list)
        lvol = vol_list(iv)
        call trace_volume(data, lvol, h, n_terminated)
    end do

    write (*, '(A,I0,A)') 'spectre_poincare: ', n_terminated, &
        ' seed(s) terminated at |s| > 1 - 1e-10'

    call free_spectre(data)

contains

    subroutine read_input(fname)
        character(len=*), intent(in) :: fname

        integer :: unit, ios
        logical :: exists

        inquire (file=fname, exist=exists)
        if (.not. exists) then
            write (*, '(A)') 'spectre_poincare: input file not found: '//trim(fname)
            error stop 'spectre_poincare: missing input file'
        end if

        open (newunit=unit, file=fname, status='old', action='read')
        read (unit, nml=poincare, iostat=ios)
        close (unit)
        if (ios /= 0) error stop 'spectre_poincare: failed to read namelist'
    end subroutine read_input

    subroutine select_volumes(mvol, requested, list)
        integer, intent(in) :: mvol
        integer, intent(in) :: requested(:)
        integer, allocatable, intent(out) :: list(:)

        integer :: i, k, count

        count = 0
        do i = 1, size(requested)
            if (requested(i) >= 1) count = count + 1
        end do

        if (count == 0) then
            allocate (list(mvol))
            do i = 1, mvol
                list(i) = i
            end do
            return
        end if

        allocate (list(count))
        k = 0
        do i = 1, size(requested)
            if (requested(i) >= 1) then
                if (requested(i) > mvol) then
                    error stop 'spectre_poincare: requested volume exceeds Mvol'
                end if
                k = k + 1
                list(k) = requested(i)
            end if
        end do
    end subroutine select_volumes

    pure function seed_s(iseed, nseeds) result(s)
        integer, intent(in) :: iseed, nseeds
        real(dp) :: s

        if (nseeds <= 1) then
            s = 0.5_dp*(seed_s_min + seed_s_max)
        else
            s = seed_s_min + (seed_s_max - seed_s_min) &
                *real(iseed - 1, dp)/real(nseeds - 1, dp)
        end if
    end function seed_s

    subroutine trace_volume(data, lvol, h, n_terminated)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: h
        integer, intent(inout) :: n_terminated

        integer :: nsec, iseed, ksec, istep, uposec, uiota
        real(dp) :: s, theta, zeta, psi, R, Z, iota
        real(dp), allocatable :: zeta_sec(:), theta_sec(:)
        character(len=64) :: poi_name, iota_name
        logical :: terminated

        nsec = nturns*data%Nfp
        allocate (zeta_sec(0:nsec), theta_sec(0:nsec))

        write (poi_name, '(A,I0,A)') 'poincare_vol', lvol, '.dat'
        write (iota_name, '(A,I0,A)') 'iota_vol', lvol, '.dat'
        open (newunit=uposec, file=trim(poi_name), status='replace', action='write')
        open (newunit=uiota, file=trim(iota_name), status='replace', action='write')
        write (uposec, '(A)') '# seed  section  s  theta  R  Z'
        write (uiota, '(A)') '# seed  iota'

        do iseed = 1, nseeds_per_volume
            s = seed_s(iseed, nseeds_per_volume)
            theta = seed_theta
            zeta = zeta0
            psi = field_line_psi(data, lvol, s, theta, zeta)
            terminated = .false.

            call spectre_rz(data, lvol, s, theta, zeta0, R, Z)
            call write_section(uposec, iseed, 0, s, theta, R, Z)
            zeta_sec(0) = zeta
            theta_sec(0) = theta

            do ksec = 1, nsec
                do istep = 1, nsteps_per_period
                    call field_line_step(data, lvol, h, s, theta, zeta, psi)
                    if (abs(s) > s_terminate) then
                        terminated = .true.
                        exit
                    end if
                end do
                if (terminated) exit
                call spectre_rz(data, lvol, s, theta, zeta0, R, Z)
                call write_section(uposec, iseed, ksec, s, theta, R, Z)
                zeta_sec(ksec) = zeta
                theta_sec(ksec) = theta
            end do

            if (terminated) then
                n_terminated = n_terminated + 1
                write (uiota, '(I6,1x,A)') iseed, 'terminated'
            else
                iota = fit_iota(zeta_sec, theta_sec)
                write (uiota, '(I6,1x,ES23.15)') iseed, iota
            end if
        end do

        close (uposec)
        close (uiota)
        deallocate (zeta_sec, theta_sec)
    end subroutine trace_volume

    subroutine write_section(unit, iseed, ksec, s, theta, R, Z)
        integer, intent(in) :: unit, iseed, ksec
        real(dp), intent(in) :: s, theta, R, Z

        write (unit, '(I6,1x,I8,4(1x,ES23.15))') iseed, ksec, s, theta, R, Z
    end subroutine write_section

    pure function fit_iota(zeta_sec, theta_sec) result(iota)
        real(dp), intent(in) :: zeta_sec(:), theta_sec(:)
        real(dp) :: iota

        integer :: n
        real(dp) :: zbar, tbar, szz, szt

        n = size(zeta_sec)
        zbar = sum(zeta_sec)/real(n, dp)
        tbar = sum(theta_sec)/real(n, dp)
        szz = sum((zeta_sec - zbar)**2)
        szt = sum((zeta_sec - zbar)*(theta_sec - tbar))
        iota = szt/szz
    end function fit_iota

end program spectre_poincare
