program test_soa_timestep
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_symplectic_soa, only: orbit_timestep_euler1_soa
    use orbit_symplectic, only: orbit_timestep_sympl_expl_impl_euler
    use orbit_symplectic_base, only: symplectic_integrator_t, EXPL_IMPL_EULER
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_derivatives
    use field_can_boozer, only: evaluate_boozer
    use boozer_sub, only: get_boozer_coordinates
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm
    use simple_main, only: init_field
    use simple, only: tracer_t

    implicit none

    real(dp), parameter :: RTOL = 1.0d-10
    type(tracer_t) :: norb
    type(field_can_t) :: f
    type(symplectic_integrator_t) :: si
    character(len=256) :: config_file
    integer :: i, npts, nerrors, ntau, ierr_single

    real(dp), allocatable :: z_r(:), z_th(:), z_ph(:), z_pphi(:)
    real(dp), allocatable :: z_r_ref(:), z_th_ref(:), z_ph_ref(:), z_pphi_ref(:)
    real(dp), allocatable :: mu(:)
    logical, allocatable :: escaped(:)
    integer, allocatable :: ierr(:)
    real(dp) :: dt, ro0, mu_scalar, atol, rtol_newton
    integer :: maxit

    print *, '=========================================='
    print *, 'Testing orbit_timestep_euler1_soa'
    print *, '=========================================='

    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)
    call get_boozer_coordinates()

    npts = 10
    ntau = 5
    dt = 1.0d-5
    maxit = 32
    atol = 1.0d-15
    rtol_newton = 1.0d-10

    allocate(z_r(npts), z_th(npts), z_ph(npts), z_pphi(npts))
    allocate(z_r_ref(npts), z_th_ref(npts), z_ph_ref(npts), z_pphi_ref(npts))
    allocate(escaped(npts), ierr(npts))
    allocate(mu(npts))

    do i = 1, npts
        z_r(i) = 0.3d0 + 0.4d0 * real(i-1, dp) / real(npts-1, dp)
        z_th(i) = 6.28318530718d0 * real(i, dp) / real(npts, dp)
        z_ph(i) = 3.14159265359d0 * real(i, dp) / real(npts, dp)
    end do

    call eval_field(f, z_r(1), z_th(1), z_ph(1), 0)
    ro0 = f%ro0
    mu_scalar = f%mu
    mu = mu_scalar

    do i = 1, npts
        call eval_field(f, z_r(i), z_th(i), z_ph(i), 0)
        z_pphi(i) = f%Aph / ro0 + 0.1d0 * f%hph
    end do

    z_r_ref = z_r
    z_th_ref = z_th
    z_ph_ref = z_ph
    z_pphi_ref = z_pphi

    do i = 1, npts
        si%ntau = ntau
        si%dt = dt
        si%z(1) = z_r_ref(i)
        si%z(2) = z_th_ref(i)
        si%z(3) = z_ph_ref(i)
        si%z(4) = z_pphi_ref(i)
        si%atol = atol
        si%rtol = rtol_newton

        call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
        call get_derivatives(f, si%z(4))

        call orbit_timestep_sympl_expl_impl_euler(si, f, ierr_single)

        z_r_ref(i) = si%z(1)
        z_th_ref(i) = si%z(2)
        z_ph_ref(i) = si%z(3)
        z_pphi_ref(i) = si%z(4)
    end do

    call orbit_timestep_euler1_soa(npts, dt, ntau, ro0, mu, atol, rtol_newton, maxit, &
        z_r, z_th, z_ph, z_pphi, escaped, ierr)

    nerrors = 0
    do i = 1, npts
        if (reldiff(z_r(i), z_r_ref(i)) > RTOL) then
            print *, 'Mismatch z_r at i=', i, z_r(i), z_r_ref(i), reldiff(z_r(i), z_r_ref(i))
            nerrors = nerrors + 1
        end if
        if (reldiff(z_th(i), z_th_ref(i)) > RTOL) then
            print *, 'Mismatch z_th at i=', i, z_th(i), z_th_ref(i), reldiff(z_th(i), z_th_ref(i))
            nerrors = nerrors + 1
        end if
        if (reldiff(z_ph(i), z_ph_ref(i)) > RTOL) then
            print *, 'Mismatch z_ph at i=', i, z_ph(i), z_ph_ref(i), reldiff(z_ph(i), z_ph_ref(i))
            nerrors = nerrors + 1
        end if
        if (reldiff(z_pphi(i), z_pphi_ref(i)) > RTOL) then
            print *, 'Mismatch z_pphi at i=', i, z_pphi(i), z_pphi_ref(i), reldiff(z_pphi(i), z_pphi_ref(i))
            nerrors = nerrors + 1
        end if
    end do

    if (nerrors > 0) then
        print *, 'FAILED: ', nerrors, ' errors found'
        stop 1
    end if

    print *, 'PASSED: orbit_timestep_euler1_soa matches single-point version'

contains

    pure function reldiff(a, b) result(rd)
        real(dp), intent(in) :: a, b
        real(dp) :: rd, denom
        denom = max(abs(a), abs(b), 1.0d-30)
        rd = abs(a - b) / denom
    end function reldiff

end program test_soa_timestep
