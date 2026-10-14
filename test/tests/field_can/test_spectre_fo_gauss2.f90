program test_spectre_fo_gauss2
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use field_spectre, only: create_spectre_field
    use magfie_sub, only: spectre_field
    use spectre_fo_gauss2, only: gauss2_prepare, gauss2_geometry, gauss2_potential, &
        gauss2_encode, gauss2_decode, gauss2_hamiltonian, gauss2_step
    implicit none
    real(dp), parameter :: WIDTH = 0.001_dp, RO0 = 2000.0_dp
    integer, parameter :: NSTEP = 160000
    real(dp) :: q(3), x(3), e(3, 3), metric(3, 3), dg(3, 3, 3), d2g(3, 3, 3, 3)
    real(dp) :: a(3), da(3, 3), d2a(3, 3, 3), b(3), hcov(3), bm, sqgb(3)
    real(dp) :: z(6), z1(6), zp(6), zm(6), energy, ep, em, grad(6), hess(6, 6)
    real(dp) :: gp(6), gm(6), hh(6, 6), tangent(6, 6), js(6, 6), defect(6, 6)
    real(dp) :: dt, total_time, first, last, shell, h, length, vmax(3), scale(6)
    real(dp) :: velocity(3), pz0, pzspread, derivative_error, delta
    real(dp) :: periodic_x(3), periodic_v(3), reference_x(3), reference_v(3)
    integer :: ierr, i, j, k, crossings, previous
    logical :: ok
    character(1024) :: filename
    character(32) :: mode

    call get_command_argument(1, filename)
    call get_command_argument(2, mode)
    allocate(spectre_field)
    call create_spectre_field(spectre_field, trim(filename), ierr)
    if (ierr /= 0) error stop 'field load'
    if (trim(mode) == 'incompatible') then
        do i = 1, spectre_field%data%mn
            if (spectre_field%data%im(i) > 0) exit
        end do
        if (i > spectre_field%data%mn) error stop 'missing test mode'
        ! A nonclosed tangential jump changes the interface normal field.
        spectre_field%data%vol(2)%Aze(0, i) = 0.01_dp
    end if
    call gauss2_prepare
    if (trim(mode) == 'incompatible') stop 0
    ! Outside the layer the embedding and curl agree with the existing field.
    do k = 1, spectre_field%data%Mvol
        q = [real(k, dp) - 0.3_dp, 0.7_dp, 0.4_dp]
        call gauss2_geometry(q, WIDTH, x, e, metric, dg, d2g)
        call spectre_field%coords%evaluate_cart(q, b)
        if (norm2(x - 100.0_dp*b) > 1.0e-10_dp) error stop 'embedding oracle'
        call gauss2_potential(q, WIDTH, RO0, a, da, d2a)
        call spectre_field%evaluate(q, a, hcov, bm, sqgb)
        b = [da(3, 2) - da(2, 3), da(1, 3) - da(3, 1), da(2, 1) - da(1, 2)]
        if (norm2(b*RO0/1.0e8_dp - sqgb) > 1.0e-12_dp) error stop 'curl oracle'
    end do
    q = [1.0_dp + WIDTH/8.0_dp, 0.7_dp, 0.4_dp]
    call gauss2_geometry(q, WIDTH, x, e, metric, dg, d2g)
    ! Transversal velocity ensures repeated gyro crossings of the layer.
    velocity = e(:, 1)/norm2(e(:, 1))
    z(1:3) = q
    call gauss2_encode(x, velocity, WIDTH, RO0, z, ok)
    if (.not. ok) error stop 'encode'
    call gauss2_hamiltonian(z, WIDTH, RO0, energy, grad, hess)
    scale = [1.0_dp, 1.0_dp, 1.0_dp, 1000.0_dp, 1000.0_dp, 1000.0_dp]
    derivative_error = 0.0_dp
    do j = 1, 6
        delta = 1.0e-7_dp*scale(j)
        if (j == 1) delta = 1.0e-8_dp
        zp = z
        zm = z
        zp(j) = zp(j) + delta
        zm(j) = zm(j) - delta
        call gauss2_hamiltonian(zp, WIDTH, RO0, ep, gp, hh)
        call gauss2_hamiltonian(zm, WIDTH, RO0, em, gm, hh)
        if (abs((ep-em)/(2.0_dp*delta) - grad(j))*scale(j) > &
            1.0e-5_dp*max(1.0_dp, abs(grad(j))*scale(j))) &
            error stop 'Hamiltonian gradient oracle'
        derivative_error = max(derivative_error, &
            maxval(abs((gp - gm)/(2.0_dp*delta) - hess(:, j))* &
            scale*scale(j))/max(1.0_dp, maxval(abs(hess(:, j))*scale*scale(j))))
    end do
    print *, 'relative Hessian oracle error: ', derivative_error
    if (derivative_error > 1.0e-5_dp) error stop 'analytic Hessian'
    h = 0.015_dp
    length = 0.02_dp
    call gauss2_step(z, WIDTH, RO0, energy, length, 0.025_dp, h, z1, dt, ok, tangent)
    if (.not. ok) error stop 'layer step'
    js = 0.0_dp
    do i = 1, 3
        js(i, i + 3) = 1.0_dp
        js(i + 3, i) = -1.0_dp
    end do
    ! A constant canonical normalization balances cm momenta and logical q.
    do j = 1, 6
        do i = 1, 6
            tangent(i, j) = tangent(i, j)*scale(j)/scale(i)
        end do
    end do
    defect = matmul(transpose(tangent), matmul(js, tangent)) - js
    print *, 'analytic tangent symplectic defect: ', maxval(abs(defect))
    if (maxval(abs(defect)) > 1.0e-10_dp) error stop 'symplecticity'
    ! The same physical map at large unwrapped angles must still converge.
    call gauss2_decode(z1, WIDTH, RO0, reference_x, reference_v)
    zp = z
    zp(2:3) = zp(2:3) + 1600.0_dp*atan(1.0_dp)
    call gauss2_step(zp, WIDTH, RO0, energy, length, 0.025_dp, h, zm, dt, ok)
    if (.not. ok) error stop 'unwrapped-angle Newton convergence'
    call gauss2_decode(zm, WIDTH, RO0, periodic_x, periodic_v)
    if (norm2(periodic_x - reference_x) > 1.0e-8_dp .or. &
        norm2(periodic_v - reference_v) > 1.0e-9_dp) &
        error stop 'discrete-map periodicity'

    if (trim(mode) == 'field') stop 0
    ! E, h, g and layer width remain frozen throughout this long FO excursion.
    pz0 = z(6)
    first = 0.0_dp
    last = 0.0_dp
    shell = 0.0_dp
    pzspread = 0.0_dp
    crossings = 0
    previous = merge(1, -1, z(1) > 1.0_dp)
    total_time = 0.0_dp
    do k = 1, NSTEP
        call gauss2_step(z, WIDTH, RO0, energy, length, 0.025_dp, h, z1, dt, ok)
        if (.not. ok) then
            print *, 'failed step and canonical state: ', k, z
            error stop 'long run Newton failure'
        end if
        z = z1
        total_time = total_time + dt
        call gauss2_decode(z, WIDTH, RO0, x, vmax)
        if (.not. all(ieee_is_finite(vmax))) error stop 'nonfinite velocity'
        delta = abs(0.5_dp*dot_product(vmax, vmax)/energy - 1.0_dp)
        shell = max(shell, delta)
        if (k <= NSTEP/4) first = max(first, delta)
        if (k > 3*NSTEP/4) last = max(last, delta)
        pzspread = max(pzspread, abs(z(6) - pz0)/max(abs(pz0), 1.0_dp))
        i = merge(1, -1, z(1) > 1.0_dp)
        if (i /= previous) crossings = crossings + 1
        previous = i
    end do
    print *, 'steps crossings physical time: ', NSTEP, crossings, total_time
    print *, 'energy first quarter / last quarter / max: ', first, last, shell
    print *, 'relative canonical p_zeta spread: ', pzspread
    if (crossings < 100) error stop 'insufficient crossings'
    if (shell > 5.0e-6_dp) error stop 'energy envelope'
    if (last > 1.25_dp*first + 1.0e-10_dp) error stop 'growing energy error'
    if (pzspread > 6.0e-8_dp) error stop 'canonical p_zeta'
end program test_spectre_fo_gauss2
