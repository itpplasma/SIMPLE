program test_soa_newton
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_symplectic_soa, only: newton1_soa
    use orbit_symplectic, only: newton1, f_sympl_euler1, jac_sympl_euler1
    use orbit_symplectic_base, only: symplectic_integrator_t
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_derivatives2
    use field_can_boozer, only: evaluate_boozer
    use boozer_sub, only: get_boozer_coordinates
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm
    use simple_main, only: init_field
    use simple, only: tracer_t
    use vector_potentail_mod, only: torflux

    implicit none

    real(dp), parameter :: RTOL = 1.0d-10
    type(tracer_t) :: norb
    type(field_can_t) :: f
    type(symplectic_integrator_t) :: si
    character(len=256) :: config_file
    integer :: i, npts, nerrors

    real(dp), allocatable :: z_r(:), z_th(:), z_ph(:), z_pphi(:), pthold(:)
    real(dp), allocatable :: x_r(:), x_pphi(:), x_r_ref(:), x_pphi_ref(:)
    real(dp), allocatable :: xlast(:,:)
    logical, allocatable :: converged(:)
    real(dp) :: x(2), xlast_single(2)
    real(dp) :: dt, ro0, mu, atol, rtol_newton
    integer :: maxit

    print *, '=========================================='
    print *, 'Testing newton1_soa vs newton1'
    print *, '=========================================='

    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)
    call get_boozer_coordinates()

    npts = 20
    dt = 1.0d-5
    maxit = 32
    atol = 1.0d-15
    rtol_newton = 1.0d-10

    allocate(z_r(npts), z_th(npts), z_ph(npts), z_pphi(npts), pthold(npts))
    allocate(x_r(npts), x_pphi(npts), x_r_ref(npts), x_pphi_ref(npts))
    allocate(xlast(2, npts))
    allocate(converged(npts))

    do i = 1, npts
        z_r(i) = 0.2d0 + 0.6d0 * real(i-1, dp) / real(npts-1, dp)
        z_th(i) = 6.28318530718d0 * real(i, dp) / real(npts, dp)
        z_ph(i) = 3.14159265359d0 * real(i, dp) / real(npts, dp)
    end do

    call eval_field(f, z_r(1), z_th(1), z_ph(1), 0)
    ro0 = f%ro0
    mu = f%mu

    do i = 1, npts
        call eval_field(f, z_r(i), z_th(i), z_ph(i), 0)
        z_pphi(i) = f%Aph / ro0 + 0.1d0 * f%hph
        call eval_field(f, z_r(i), z_th(i), z_ph(i), 2)
        call get_derivatives2(f, z_pphi(i))
        pthold(i) = f%pth
    end do

    x_r = z_r
    x_pphi = z_pphi

    do i = 1, npts
        si%dt = dt
        si%z(1) = z_r(i)
        si%z(2) = z_th(i)
        si%z(3) = z_ph(i)
        si%z(4) = z_pphi(i)
        si%pthold = pthold(i)
        si%atol = atol
        si%rtol = rtol_newton

        call eval_field(f, z_r(i), z_th(i), z_ph(i), 2)
        call get_derivatives2(f, z_pphi(i))

        x(1) = z_r(i)
        x(2) = z_pphi(i)
        call newton1(si, f, x, maxit, xlast_single)
        x_r_ref(i) = x(1)
        x_pphi_ref(i) = x(2)
    end do

    x_r = z_r
    x_pphi = z_pphi
    call newton1_soa(npts, dt, ro0, mu, atol, rtol_newton, maxit, &
        z_r, z_th, z_ph, z_pphi, pthold, x_r, x_pphi, converged)

    nerrors = 0
    do i = 1, npts
        if (reldiff(x_r(i), x_r_ref(i)) > RTOL) then
            print *, 'Mismatch x_r at i=', i, x_r(i), x_r_ref(i), reldiff(x_r(i), x_r_ref(i))
            nerrors = nerrors + 1
        end if
        if (reldiff(x_pphi(i), x_pphi_ref(i)) > RTOL) then
            print *, 'Mismatch x_pphi at i=', i, x_pphi(i), x_pphi_ref(i), reldiff(x_pphi(i), x_pphi_ref(i))
            nerrors = nerrors + 1
        end if
    end do

    if (nerrors > 0) then
        print *, 'FAILED: ', nerrors, ' errors found'
        stop 1
    end if

    print *, 'PASSED: newton1_soa matches single-point newton1'

contains

    pure function reldiff(a, b) result(rd)
        real(dp), intent(in) :: a, b
        real(dp) :: rd, denom
        denom = max(abs(a), abs(b), 1.0d-30)
        rd = abs(a - b) / denom
    end function reldiff

end program test_soa_newton
