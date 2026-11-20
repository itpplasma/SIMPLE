program test_nbi_birth_chartmap
    !> Neutral beam deposition on a chartmap of a VMEC equilibrium equals the
    !> deposition in VMEC coordinates. The chart radius is rho = sqrt(s), so
    !> profiles read at the chart radius instead of s change the shine-through.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_api, only: simple_reset_config, simple_init, simple_sample
    use params, only: field_input, coord_input, startmode, isw_field_type, integmode, &
                      nper, npoiper, ntimstep, deterministic, n_d, n_e, &
                      particle_energy_eV, nbi_file, ntestpart, contr_pp, nbi_torque_bins
    use simple_profiles, only: ni1_coef, ni2_coef, te_coef
    use profile_radius, only: profile_s
    use timing, only: init_timer
    implicit none
    real(dp) :: ionized(2), shine(2), smax
    real(dp), allocatable :: z(:, :)
    character(1024) :: wout, chartmap
    integer :: i

    call get_command_argument(1, wout)
    call get_command_argument(2, chartmap)
    call init_timer
    open (newunit=i, file='beamlets.dat', status='replace')
    ! Tangential pencil through the core, 0.5 MeV deuterons
    write (i, '(a)') 'T 2 5.0e5 1.0e6 1.0 0.0 0.0 10.2265706505403 -8.0 0.0 0.0 1.0 0.0'
    close (i)

    call births(trim(wout), 1, ionized(1), shine(1), smax)
    call births(trim(chartmap), 5, ionized(2), shine(2), smax)
    print '(a,2es14.6)', 'shine-through VMEC, chartmap [1/s] ', shine
    print '(a,2es14.6)', 'ionized rate VMEC, chartmap [1/s]  ', ionized
    print '(a,f9.5)', 'largest birth s on the chartmap     ', smax
    if (abs(shine(2)/shine(1) - 1.0_dp) > 1.0e-3_dp) &
        error stop 'chartmap shine-through differs from VMEC'
    if (abs(ionized(2)/ionized(1) - 1.0_dp) > 1.0e-3_dp) &
        error stop 'chartmap ionized rate differs from VMEC'
    if (smax >= 1.0_dp) error stop 'birth outside the LCFS'
    print '(a)', 'PASS'

contains

    subroutine births(coords, field_type, ionized, shine, smax)
        character(*), intent(in) :: coords
        integer, intent(in) :: field_type
        real(dp), intent(out) :: ionized, shine, smax
        integer :: k

        call simple_reset_config
        field_input = trim(wout)
        coord_input = coords
        isw_field_type = field_type
        integmode = 0
        startmode = 8
        deterministic = .true.
        nper = 1
        npoiper = 4
        ntimstep = 2
        n_d = 2
        n_e = 1
        particle_energy_eV = 5.0e5_dp
        ntestpart = 200
        contr_pp = -1.0e10_dp
        nbi_file = 'beamlets.dat'
        nbi_torque_bins = 0
        ! Steep profiles in s: deposition depends on the radius they are read at.
        te_coef(1:2) = [1.0_dp, -0.9_dp]
        ni1_coef(1:2) = [1.0_dp, -0.95_dp]
        ni2_coef(1:2) = [1.0_dp, -0.95_dp]
        call simple_init
        if (allocated(z)) deallocate (z)
        allocate (z(5, ntestpart))
        call simple_sample(ntestpart, z)
        smax = 0.0_dp
        do k = 1, ntestpart
            smax = max(smax, profile_s(z(1:3, k)))
        end do
        call read_rates(ionized, shine)
    end subroutine births

    subroutine read_rates(ionized, shine)
        real(dp), intent(out) :: ionized, shine
        character(256) :: line, key
        real(dp) :: value
        integer :: unit, ios

        ionized = -1.0_dp
        shine = -1.0_dp
        open (newunit=unit, file='nbi_source.dat', status='old')
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:1) == '#') cycle
            read (line, *, iostat=ios) key, value
            if (ios /= 0) cycle
            if (key == 'ionized_rate') ionized = value
            if (key == 'shine_rate') shine = value
        end do
        close (unit)
        if (min(ionized, shine) <= 0.0_dp) error stop 'rates missing in nbi_source.dat'
    end subroutine read_rates

end program test_nbi_birth_chartmap
