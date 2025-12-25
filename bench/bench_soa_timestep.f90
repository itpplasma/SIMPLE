program bench_soa_timestep
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_symplectic_soa, only: orbit_timestep_euler1_soa
    use orbit_symplectic, only: orbit_timestep_sympl_expl_impl_euler
    use orbit_symplectic_base, only: symplectic_integrator_t
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_derivatives
    use field_can_boozer, only: evaluate_boozer
    use boozer_sub, only: get_boozer_coordinates
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm
    use simple_main, only: init_field
    use simple, only: tracer_t, init_sympl

    implicit none

    type(tracer_t) :: norb
    type(field_can_t) :: f
    type(symplectic_integrator_t) :: si
    character(len=256) :: config_file
    integer :: i, npts, ntau, ierr_single, rep, nreps

    real(dp), allocatable :: z_r(:), z_th(:), z_ph(:), z_pphi(:)
    real(dp), allocatable :: z_r0(:), z_th0(:), z_ph0(:), z_pphi0(:)
    real(dp), allocatable :: mu(:)
    logical, allocatable :: escaped(:)
    integer, allocatable :: ierr(:)
    real(dp) :: z(5), dt, ro0, mu_scalar, atol, rtol_newton
    integer :: maxit
    real(dp) :: t_soa, t_single
    integer :: count_start, count_end, count_rate

    print *, '=========================================='
    print *, 'Benchmarking SoA vs single-point timestep'
    print *, '=========================================='

    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)
    call get_boozer_coordinates()

    npts = 100
    ntau = 5
    dt = 1.0d-5
    maxit = 32
    atol = 1.0d-15
    rtol_newton = 1.0d-10
    nreps = 20

    print *, 'npts = ', npts
    print *, 'ntau = ', ntau
    print *, 'nreps = ', nreps

    allocate(z_r(npts), z_th(npts), z_ph(npts), z_pphi(npts))
    allocate(z_r0(npts), z_th0(npts), z_ph0(npts), z_pphi0(npts))
    allocate(escaped(npts), ierr(npts), mu(npts))

    do i = 1, npts
        z_r0(i) = 0.3d0 + 0.4d0 * real(i-1, dp) / real(npts-1, dp)
        z_th0(i) = 6.28318530718d0 * real(i, dp) / real(npts, dp)
        z_ph0(i) = 3.14159265359d0 * real(i, dp) / real(npts, dp)
    end do

    call eval_field(f, z_r0(1), z_th0(1), z_ph0(1), 0)
    ro0 = f%ro0
    mu_scalar = f%mu
    mu = mu_scalar

    do i = 1, npts
        call eval_field(f, z_r0(i), z_th0(i), z_ph0(i), 0)
        z_pphi0(i) = f%Aph / ro0 + 0.1d0 * f%hph
    end do

    print *, ''
    print *, 'Warming up...'
    z_r = z_r0; z_th = z_th0; z_ph = z_ph0; z_pphi = z_pphi0
    call orbit_timestep_euler1_soa(npts, dt, ntau, ro0, mu, atol, rtol_newton, maxit, &
        z_r, z_th, z_ph, z_pphi, escaped, ierr)

    print *, ''
    print *, 'Benchmarking SoA version...'
    call system_clock(count_start, count_rate)
    do rep = 1, nreps
        z_r = z_r0; z_th = z_th0; z_ph = z_ph0; z_pphi = z_pphi0
        call orbit_timestep_euler1_soa(npts, dt, ntau, ro0, mu, atol, rtol_newton, maxit, &
            z_r, z_th, z_ph, z_pphi, escaped, ierr)
    end do
    call system_clock(count_end)
    t_soa = real(count_end - count_start, dp) / real(count_rate, dp)

    print *, ''
    print *, 'Benchmarking single-point loop version...'
    call system_clock(count_start, count_rate)
    do rep = 1, nreps
        do i = 1, npts
            z(1) = z_r0(i)
            z(2) = z_th0(i)
            z(3) = z_ph0(i)
            z(4) = 1.0d0
            z(5) = 0.5d0

            call init_sympl(si, f, z, dt, dt, rtol_newton, 1)
            call orbit_timestep_sympl_expl_impl_euler(si, f, ierr_single)
        end do
    end do
    call system_clock(count_end)
    t_single = real(count_end - count_start, dp) / real(count_rate, dp)

    print *, ''
    print *, '=========================================='
    print *, 'RESULTS'
    print *, '=========================================='
    print *, 'SoA time:         ', t_soa, ' s'
    print *, 'Single-point time:', t_single, ' s'
    print *, 'Speedup:          ', t_single / t_soa, 'x'
    print *, ''
    print *, 'Time per particle-timestep:'
    print *, '  SoA:         ', t_soa / real(nreps * npts * ntau, dp) * 1.0d6, ' us'
    print *, '  Single-point:', t_single / real(nreps * npts * ntau, dp) * 1.0d6, ' us'

end program bench_soa_timestep
