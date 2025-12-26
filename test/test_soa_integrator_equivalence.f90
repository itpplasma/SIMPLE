program test_soa_integrator_equivalence
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_symplectic_soa, only: trace_orbit_soa_omp1
    use orbit_symplectic, only: orbit_timestep_sympl
    use orbit_symplectic_base, only: symplectic_integrator_t, EXPL_IMPL_EULER
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_derivatives
    use boozer_sub, only: get_boozer_coordinates, vmec_to_boozer
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
        params_init, ntau, dtaumin
    use simple_main, only: init_field
    use simple, only: tracer_t, init_sympl
    use parmot_mod, only: ro0

    implicit none

    real(dp), parameter :: RTOL_FIRST_STEPS = 1.0d-12
    type(tracer_t) :: norb
    character(len=256) :: config_file
    integer :: i, npts, nerrors, ntimstep, ierr_single, kt

    real(dp), allocatable :: zstart(:,:), z_final_soa(:,:), z_final_old(:,:)
    real(dp), allocatable :: times_lost_soa(:), times_lost_old(:)
    integer, allocatable :: ierr_soa(:), ierr_old(:)
    real(dp) :: z(5), atol, rtol_newton, theta_B, phi_B
    integer :: maxit, ktau, it
    real(dp) :: max_reldiff_r, max_reldiff_th, max_reldiff_ph

    print *, '=========================================='
    print *, 'Testing SoA vs Old integrator equivalence'
    print *, '=========================================='
    print *, 'Verifying trajectories match to floating point precision'
    print *, ''

    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)
    call params_init
    call get_boozer_coordinates()

    npts = 1
    ntimstep = 5
    maxit = 32
    atol = 1.0d-15
    rtol_newton = 1.0d-10

    allocate(zstart(5, npts), z_final_soa(5, npts), z_final_old(5, npts))
    allocate(times_lost_soa(npts), times_lost_old(npts))
    allocate(ierr_soa(npts), ierr_old(npts))

    zstart(1, 1) = 0.3d0
    zstart(2, 1) = 0.1d0
    zstart(3, 1) = 0.2d0
    zstart(4, 1) = 1.0d0
    zstart(5, 1) = 0.5d0

    print *, 'Test parameters:'
    print *, '  npts = ', npts
    print *, '  ntimstep = ', ntimstep
    print *, '  ntau = ', ntau
    print *, '  dtaumin = ', dtaumin
    print *, '  ro0 = ', ro0
    print *, ''

    do i = 1, npts
        z(1) = zstart(1, i)
        call vmec_to_boozer(zstart(1, i), mod(zstart(2, i), 6.283185307179586d0), &
            mod(zstart(3, i), 6.283185307179586d0), theta_B, phi_B)
        z(2) = mod(theta_B, 6.283185307179586d0)
        z(3) = mod(phi_B, 6.283185307179586d0)
        z(4) = zstart(4, i)
        z(5) = zstart(5, i)

        call init_sympl(norb%si, norb%f, z, dtaumin, dtaumin, rtol_newton, EXPL_IMPL_EULER)

        print *, 'OLD initial state:'
        print *, '  z = ', norb%si%z
        print *, '  mu = ', norb%f%mu
        print *, '  vpar = ', norb%f%vpar

        ierr_single = 0
        times_lost_old(i) = 0.0d0
        kt = 0
        do it = 1, ntimstep
            do ktau = 1, ntau
                call orbit_timestep_sympl(norb%si, norb%f, ierr_single)
                if (ierr_single /= 0) exit
                kt = kt + 1
            end do
            if (ierr_single /= 0) then
                times_lost_old(i) = real(kt, dp) * dtaumin
                exit
            end if
        end do
        ierr_old(i) = ierr_single

        z_final_old(1, i) = norb%si%z(1)
        z_final_old(2, i) = norb%si%z(2)
        z_final_old(3, i) = norb%si%z(3)
        z_final_old(4, i) = dsqrt(norb%f%mu * norb%f%Bmod + 0.5d0 * norb%f%vpar**2)
        z_final_old(5, i) = norb%f%vpar / (z_final_old(4, i) * dsqrt(2.0d0))
    end do

    print *, 'OLD after ', ntimstep * ntau, ' substeps:'
    print *, '  z_final = ', z_final_old(:, 1)
    print *, ''

    call trace_orbit_soa_omp1(npts, zstart, ntimstep, ntau, dtaumin, ro0, &
        atol, rtol_newton, maxit, z_final_soa, times_lost_soa, ierr_soa)

    print *, 'SoA after ', ntimstep * ntau, ' substeps:'
    print *, '  z_final = ', z_final_soa(:, 1)
    print *, ''

    nerrors = 0
    max_reldiff_r = 0.0d0
    max_reldiff_th = 0.0d0
    max_reldiff_ph = 0.0d0

    do i = 1, npts
        max_reldiff_r = max(max_reldiff_r, reldiff(z_final_soa(1, i), z_final_old(1, i)))
        max_reldiff_th = max(max_reldiff_th, reldiff(z_final_soa(2, i), z_final_old(2, i)))
        max_reldiff_ph = max(max_reldiff_ph, reldiff(z_final_soa(3, i), z_final_old(3, i)))

        if (reldiff(z_final_soa(1, i), z_final_old(1, i)) > RTOL_FIRST_STEPS) then
            print *, 'Mismatch z_r at i=', i
            print *, '  OLD: ', z_final_old(1, i)
            print *, '  SoA: ', z_final_soa(1, i)
            print *, '  reldiff: ', reldiff(z_final_soa(1, i), z_final_old(1, i))
            nerrors = nerrors + 1
        end if
        if (reldiff(z_final_soa(2, i), z_final_old(2, i)) > RTOL_FIRST_STEPS) then
            print *, 'Mismatch z_th at i=', i
            print *, '  OLD: ', z_final_old(2, i)
            print *, '  SoA: ', z_final_soa(2, i)
            print *, '  reldiff: ', reldiff(z_final_soa(2, i), z_final_old(2, i))
            nerrors = nerrors + 1
        end if
        if (reldiff(z_final_soa(3, i), z_final_old(3, i)) > RTOL_FIRST_STEPS) then
            print *, 'Mismatch z_ph at i=', i
            print *, '  OLD: ', z_final_old(3, i)
            print *, '  SoA: ', z_final_soa(3, i)
            print *, '  reldiff: ', reldiff(z_final_soa(3, i), z_final_old(3, i))
            nerrors = nerrors + 1
        end if
    end do

    print *, 'Results after ', ntimstep, ' macrosteps (', ntimstep * ntau, ' substeps):'
    print *, '  max relative diff z_r:  ', max_reldiff_r
    print *, '  max relative diff z_th: ', max_reldiff_th
    print *, '  max relative diff z_ph: ', max_reldiff_ph
    print *, ''

    if (nerrors > 0) then
        print *, 'FAILED: ', nerrors, ' errors found'
        print *, 'Trajectories diverge beyond floating point tolerance'
        stop 1
    end if

    print *, 'PASSED: SoA integrator matches old integrator to floating point precision'
    print *, '        (relative difference < ', RTOL_FIRST_STEPS, ')'

contains

    pure function reldiff(a, b) result(rd)
        real(dp), intent(in) :: a, b
        real(dp) :: rd, denom
        denom = max(abs(a), abs(b), 1.0d-30)
        rd = abs(a - b) / denom
    end function reldiff

end program test_soa_integrator_equivalence
