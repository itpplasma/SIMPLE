program test_boozer_chartmap_roundtrip
    !> Roundtrip validation: VMEC -> Boozer -> export chartmap -> reimport ->
    !> compare field values and symplectic orbits between direct and file paths.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: nper, rmajor, netcdffile, multharm, ns_A, ns_s, ns_tp
    use parmot_mod, only: rmu, ro0
    use velo_mod, only: isw_field_type
    use boozer_coordinates_mod, only: use_B_r
    use boozer_sub, only: splint_boozer_coord, get_boozer_coordinates, &
        vmec_to_boozer, export_boozer_chartmap, load_boozer_from_chartmap, &
        reset_boozer_batch_splines
    use spline_vmec_sub, only: spline_vmec_data
    use vmecin_sub, only: stevvo
    use field_can_mod, only: field_can_from_name, field_can_init, &
        eval_field => evaluate, field_can_t, get_val
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl, &
        symplectic_integrator_t

    implicit none

    real(dp), parameter :: pi = 3.14159265358979_dp
    real(dp), parameter :: twopi = 2.0_dp * pi

    ! Field comparison
    integer, parameter :: n_test = 50
    real(dp) :: s_test(n_test), th_test(n_test), ph_test(n_test)
    real(dp) :: Bmod_ref(n_test), Bmod_new(n_test)

    ! Orbit comparison
    integer, parameter :: n_orbit = 500
    real(dp) :: orbit_direct(n_orbit, 5), orbit_chartmap(n_orbit, 5)

    ! Work variables
    real(dp) :: fper, dtau, RT0
    real(dp) :: A_theta, A_phi_val, dA_theta_dr, dA_phi_dr
    real(dp) :: d2A_phi_dr2, d3A_phi_dr3
    real(dp) :: Bth, dBth, d2Bth, Bph, dBph, d2Bph
    real(dp) :: Bmod, dBmod(3), d2Bmod(6), Br, dBr(3), d2Br(6)
    real(dp) :: phi_period, z(4), z0(5), vartheta, varphi
    real(dp) :: rel_err, max_err_bmod, max_err_orbit
    integer :: i, ierr, nfail, n_steps_done
    character(len=256) :: chartmap_file
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f

    nfail = 0
    chartmap_file = 'roundtrip_test.nc'

    print *, 'Starting roundtrip test...'

    ! =========================================================
    ! Step 1: Initialize VMEC and compute Boozer coordinates
    ! =========================================================
    isw_field_type = 2
    rmu = 1.0e8_dp
    netcdffile = 'wout.nc'
    multharm = 5
    ns_A = 5
    ns_s = 5
    ns_tp = 5

    call spline_vmec_data
    block
        integer :: L1i
        real(dp) :: R0i, cbfi, bz0i, bf0
        call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0)
        fper = twopi / real(L1i, dp)
    end block
    phi_period = fper
    RT0 = rmajor

    use_B_r = .false.
    call get_boozer_coordinates
    call field_can_from_name("boozer")

    ! ro0 for orbit integration: rlarm * bmod00
    ro0 = 2.23e-2_dp * 2.81e5_dp

    print *, '=== Step 1: VMEC + Boozer initialized ==='
    print *, '  nper=', nper, ' R0=', RT0, ' fper=', fper

    ! =========================================================
    ! Step 2: Evaluate fields at test points (direct path)
    ! =========================================================
    call init_test_points(phi_period)

    do i = 1, n_test
        call splint_boozer_coord(s_test(i), th_test(i), ph_test(i), 0, &
            A_theta, A_phi_val, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, &
            d3A_phi_dr3, Bth, dBth, d2Bth, Bph, dBph, d2Bph, &
            Bmod, dBmod, d2Bmod, Br, dBr, d2Br)
        Bmod_ref(i) = Bmod
    end do

    print *, '=== Step 2: Direct field evaluation done ==='

    ! =========================================================
    ! Step 3: Export to Boozer chartmap NetCDF
    ! =========================================================
    call export_boozer_chartmap(chartmap_file)
    print *, '=== Step 3: Exported chartmap ==='

    ! =========================================================
    ! Step 4: Reimport from chartmap and evaluate fields
    ! =========================================================
    call reset_boozer_batch_splines
    call load_boozer_from_chartmap(chartmap_file)

    do i = 1, n_test
        call splint_boozer_coord(s_test(i), th_test(i), ph_test(i), 0, &
            A_theta, A_phi_val, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, &
            d3A_phi_dr3, Bth, dBth, d2Bth, Bph, dBph, d2Bph, &
            Bmod, dBmod, d2Bmod, Br, dBr, d2Br)
        Bmod_new(i) = Bmod
    end do

    print *, '=== Step 4: Chartmap field evaluation done ==='

    ! =========================================================
    ! Step 5: Compare fields
    ! =========================================================
    max_err_bmod = 0.0_dp
    open(unit=20, file='/tmp/boozer_field_comparison.dat', status='replace')
    write(20, '(a)') '# s  theta  phi  Bmod_ref  Bmod_chartmap  rel_err'

    do i = 1, n_test
        if (abs(Bmod_ref(i)) > 0.0_dp) then
            rel_err = abs(Bmod_new(i) - Bmod_ref(i)) / abs(Bmod_ref(i))
        else
            rel_err = abs(Bmod_new(i))
        end if
        max_err_bmod = max(max_err_bmod, rel_err)
        write(20, '(6es18.10)') s_test(i), th_test(i), ph_test(i), &
            Bmod_ref(i), Bmod_new(i), rel_err
    end do
    close(20)

    print *, '  max relative error Bmod:', max_err_bmod
    if (max_err_bmod > 2.0e-3_dp) then
        print *, 'FAIL: Bmod roundtrip error too large:', max_err_bmod
        nfail = nfail + 1
    else
        print *, 'PASS: Bmod roundtrip error < 2e-3'
    end if

    ! =========================================================
    ! Step 6: Run symplectic orbit (direct path)
    ! =========================================================
    call reset_boozer_batch_splines
    call get_boozer_coordinates
    call field_can_from_name("boozer")

    call vmec_to_boozer(0.35_dp, 0.33_dp, 0.97_dp, vartheta, varphi)

    dtau = fper * RT0 / 400.0_dp

    z0(1) = 0.35_dp
    z0(2) = vartheta
    z0(3) = varphi
    z0(4) = 1.0_dp
    z0(5) = 0.1_dp

    call run_symplectic_orbit(z0, dtau, n_orbit, orbit_direct, n_steps_done)

    open(unit=21, file='/tmp/orbit_direct.dat', status='replace')
    write(21, '(a)') '# time  s  theta  phi  pphi'
    do i = 1, n_steps_done
        write(21, '(5es18.10)') orbit_direct(i, :)
    end do
    close(21)

    print *, '=== Step 6: Direct orbit done, steps=', n_steps_done, ' ==='

    ! =========================================================
    ! Step 7: Run symplectic orbit (chartmap path)
    ! =========================================================
    call reset_boozer_batch_splines
    call load_boozer_from_chartmap(chartmap_file)
    call field_can_from_name("boozer")

    call run_symplectic_orbit(z0, dtau, n_orbit, orbit_chartmap, n_steps_done)

    open(unit=22, file='/tmp/orbit_chartmap.dat', status='replace')
    write(22, '(a)') '# time  s  theta  phi  pphi'
    do i = 1, n_steps_done
        write(22, '(5es18.10)') orbit_chartmap(i, :)
    end do
    close(22)

    print *, '=== Step 7: Chartmap orbit done, steps=', n_steps_done, ' ==='

    ! =========================================================
    ! Step 8: Compare orbits
    ! =========================================================
    max_err_orbit = 0.0_dp
    do i = 1, n_steps_done
        rel_err = abs(orbit_direct(i, 2) - orbit_chartmap(i, 2))
        max_err_orbit = max(max_err_orbit, rel_err)
    end do

    print *, '  max |s_direct - s_chartmap|:', max_err_orbit
    if (max_err_orbit > 1.0e-6_dp) then
        print *, 'FAIL: orbit roundtrip error too large:', max_err_orbit
        nfail = nfail + 1
    else
        print *, 'PASS: orbit roundtrip error < 1e-6'
    end if

    ! =========================================================
    ! Summary
    ! =========================================================
    print *, ''
    if (nfail == 0) then
        print *, 'All roundtrip tests passed.'
    else
        print *, nfail, ' roundtrip tests failed.'
        error stop 'Roundtrip test failures'
    end if

contains

    subroutine init_test_points(phi_per)
        real(dp), intent(in) :: phi_per
        integer :: j
        real(dp) :: frac

        do j = 1, n_test
            frac = real(j - 1, dp) / real(n_test - 1, dp)
            s_test(j) = 0.1_dp + 0.7_dp * frac
            th_test(j) = twopi * frac
            ph_test(j) = phi_per * mod(real(j * 7, dp), real(n_test, dp)) &
                          / real(n_test, dp)
        end do
    end subroutine init_test_points

    subroutine run_symplectic_orbit(z0_in, dt, nsteps, orbit_out, steps_done)
        real(dp), intent(in) :: z0_in(5), dt
        integer, intent(in) :: nsteps
        real(dp), intent(out) :: orbit_out(nsteps, 5)
        integer, intent(out) :: steps_done

        type(symplectic_integrator_t) :: si_loc
        type(field_can_t) :: f_loc
        real(dp) :: z_loc(4), pphi
        integer :: j, ierr_loc

        ! Initialize field_can_t
        call eval_field(f_loc, z0_in(1), z0_in(2), z0_in(3), 0)

        f_loc%mu = 0.5_dp * z0_in(4)**2 * (1.0_dp - z0_in(5)**2) / &
                   f_loc%Bmod * 2.0_dp
        f_loc%ro0 = ro0 / sqrt(2.0_dp)
        f_loc%vpar = z0_in(4) * z0_in(5) * sqrt(2.0_dp)

        z_loc(1:3) = z0_in(1:3)
        pphi = f_loc%vpar * f_loc%hph + f_loc%Aph / f_loc%ro0
        z_loc(4) = pphi

        ! Midpoint integrator (mode=3), single step per call
        call orbit_sympl_init(si_loc, f_loc, z_loc, &
                               dt / sqrt(2.0_dp), 1, 1.0e-10_dp, 3)

        steps_done = 0
        do j = 1, nsteps
            call orbit_timestep_sympl(si_loc, f_loc, ierr_loc)
            if (ierr_loc /= 0) exit
            steps_done = j
            orbit_out(j, 1) = dt * real(j, dp)
            orbit_out(j, 2) = si_loc%z(1)   ! s
            orbit_out(j, 3) = si_loc%z(2)   ! theta
            orbit_out(j, 4) = si_loc%z(3)   ! phi
            orbit_out(j, 5) = si_loc%z(4)   ! pphi
        end do
    end subroutine run_symplectic_orbit

end program test_boozer_chartmap_roundtrip
