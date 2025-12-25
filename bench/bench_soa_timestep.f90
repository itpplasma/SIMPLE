program bench_soa_timestep
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_symplectic_soa, only: trace_orbit_soa
    use orbit_symplectic, only: orbit_timestep_sympl
    use orbit_symplectic_base, only: symplectic_integrator_t
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_derivatives
    use field_can_boozer, only: evaluate_boozer
    use boozer_sub, only: get_boozer_coordinates, vmec_to_boozer
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm
    use simple_main, only: init_field
    use simple, only: tracer_t, init_sympl
    use parmot_mod, only: ro0

    implicit none

    type(tracer_t) :: norb
    type(field_can_t) :: f
    type(symplectic_integrator_t) :: si
    character(len=256) :: config_file
    integer :: i, npts, nerrors, ntimstep, ntau, ierr_single, kt

    real(dp), allocatable :: zstart(:, :), z_final(:, :), z_final_ref(:, :)
    real(dp), allocatable :: times_lost(:), times_lost_ref(:)
    integer, allocatable :: ierr(:), ierr_ref(:)
    real(dp) :: z(5), dt, atol, rtol_newton, theta_B, phi_B
    integer :: maxit, ktau, it, rep, nreps
    real(dp) :: t_soa, t_single
    integer :: count_start, count_end, count_rate

    print *, '=========================================='
    print *, 'Benchmarking SoA vs single-point tracing'
    print *, '=========================================='

    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)
    call get_boozer_coordinates()

    npts = 32
    ntimstep = 10
    ntau = 5
    dt = 1.0d-5
    maxit = 32
    atol = 1.0d-15
    rtol_newton = 1.0d-10
    nreps = 10

    print *, 'npts = ', npts
    print *, 'ntimstep = ', ntimstep
    print *, 'ntau = ', ntau
    print *, 'nreps = ', nreps
    print *, 'Total timesteps per particle per rep = ', ntimstep*ntau

    allocate (zstart(5, npts), z_final(5, npts), z_final_ref(5, npts))
    allocate (times_lost(npts), times_lost_ref(npts))
    allocate (ierr(npts), ierr_ref(npts))

    do i = 1, npts
        zstart(1, i) = 0.3d0 + 0.3d0*real(i - 1, dp)/real(npts - 1, dp)
        zstart(2, i) = 6.28318530718d0*real(i, dp)/real(npts, dp)
        zstart(3, i) = 3.14159265359d0*real(i, dp)/real(npts, dp)
        zstart(4, i) = 1.0d0
        zstart(5, i) = 0.3d0 + 0.4d0*real(i - 1, dp)/real(npts - 1, dp)
    end do

    print *, ''
    print *, 'Benchmarking single-point loop version...'
    call system_clock(count_start, count_rate)
    do rep = 1, nreps
        do i = 1, npts
            z(1) = zstart(1, i)
            call vmec_to_boozer(zstart(1, i), mod(zstart(2, i), 6.283185307179586d0), &
                                mod(zstart(3, i), 6.283185307179586d0), theta_B, phi_B)
            z(2) = mod(theta_B, 6.283185307179586d0)
            z(3) = mod(phi_B, 6.283185307179586d0)
            z(4) = zstart(4, i)
            z(5) = zstart(5, i)

            call init_sympl(norb%si, norb%f, z, dt, dt, rtol_newton, 1)

            ierr_single = 0
            times_lost_ref(i) = 0.0d0
            kt = 0
            do it = 1, ntimstep
                do ktau = 1, ntau
                    call orbit_timestep_sympl(norb%si, norb%f, ierr_single)
                    if (ierr_single /= 0) exit
                    kt = kt + 1
                end do
                if (ierr_single /= 0) then
                    times_lost_ref(i) = real(kt, dp)*dt
                    exit
                end if
            end do
            ierr_ref(i) = ierr_single

            z_final_ref(1, i) = norb%si%z(1)
            z_final_ref(2, i) = norb%si%z(2)
            z_final_ref(3, i) = norb%si%z(3)
            z_final_ref(4, i) = dsqrt(norb%f%mu*norb%f%Bmod + 0.5d0*norb%f%vpar**2)
            z_final_ref(5, i) = norb%f%vpar/(z_final_ref(4, i)*dsqrt(2.0d0))
        end do
    end do
    call system_clock(count_end)
    t_single = real(count_end - count_start, dp)/real(count_rate, dp)

    print *, ''
    print *, 'Benchmarking SoA version (trace_orbit_soa)...'
    call system_clock(count_start, count_rate)
    do rep = 1, nreps
        call trace_orbit_soa(npts, zstart, ntimstep, ntau, dt, ro0, &
                             atol, rtol_newton, maxit, z_final, times_lost, ierr)
    end do
    call system_clock(count_end)
    t_soa = real(count_end - count_start, dp)/real(count_rate, dp)

    print *, ''
    print *, '=========================================='
    print *, 'RESULTS'
    print *, '=========================================='
    print *, 'SoA time:         ', t_soa, ' s'
    print *, 'Single-point time:', t_single, ' s'
    if (t_soa > 0.0d0) then
        print *, 'Speedup:          ', t_single/t_soa, 'x'
    end if
    print *, ''
    print *, 'Time per particle-trace:'
    print *, '  SoA:         ', t_soa/real(nreps*npts, dp)*1.0d6, ' us'
    print *, '  Single-point:', t_single/real(nreps*npts, dp)*1.0d6, ' us'

end program bench_soa_timestep
