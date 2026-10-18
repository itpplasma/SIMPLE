program test_candidates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use candidate_field, only: model_t, local_t, field, geometry, potential, &
        side_of, cross, frame, connection, canonical_pphi
    use transforms, only: lift, inverse, gc_energy, kernel_available
    use candidate_orbits, only: options_t, orbit_t, initialize, advance, &
        particle, helix, plane_step
    use rk_reference, only: rk_advance
    use sheet_field, only: psi_value
    implicit none
    type(model_t) :: pm
    type(local_t) :: lf
    type(options_t) :: opt, unsupported
    type(orbit_t) :: orbit
    real(dp) :: x(3), b(3), jac(3, 3), num(3, 3), xp(3), xm(3), bp(3), bm(3)
    real(dp) :: y(6), yn(6), g(6), gn(6), err, speed, rr, psi, z(4), ht
    real(dp) :: out(6), v(3), xr(3), vr(3), xn(3), vn(3)
    real(dp) :: tangent(6, 6), ofo(6, 6), ogc(6, 6), diff(6, 6), scales(6)
    real(dp) :: rs(3), cr(3), bst(3), p(3), n(3), hh, defects(3), exponent
    real(dp) :: chunk_x(3), chunk_v(3), target_x(3), target_v(3)
    integer :: i, j, order, region, nc, level, events, total_events
    logical :: ok
    pm%torus%a_in = 0.25_dp
    pm%torus%a_out = 2.0_dp
    pm%delta = 0.1_dp
    x = [10.3_dp, 0.7_dp, 0.2_dp]
    call field(pm, x, 1, b, jac)
    do j = 1, 3
        xp = x
        xm = x
        xp(j) = xp(j) + 1.0e-5_dp
        xm(j) = xm(j) - 1.0e-5_dp
        call field(pm, xp, 1, bp)
        call field(pm, xm, 1, bm)
        num(:, j) = (bp - bm)/2.0e-5_dp
    end do
    call require(maxval(abs(num - jac)) < 2.0e-10_dp, 'analytic field Jacobian')
    ! Curl A from an independent finite-difference oracle.
    do j = 1, 3
        xp = x
        xm = x
        xp(j) = xp(j) + 1.0e-5_dp
        xm(j) = xm(j) - 1.0e-5_dp
        num(:, j) = (potential(pm, xp) - potential(pm, xm))/2.0e-5_dp
    end do
    v = [num(3, 2) - num(2, 3), num(1, 3) - num(3, 1), &
         num(2, 1) - num(1, 2)]
    call require(norm2(v - b) < 1.0e-8_dp, 'B=curl A numerical oracle')
    call require(abs(jac(1, 1) + jac(2, 2) + jac(3, 3)) < 1.0e-14_dp, &
        'divergence numerical oracle')
    call connection(pm, x, 1, rs, cr)
    do j = 1, 3
        xp = x
        xm = x
        xp(j) = xp(j) + 2.0e-5_dp
        xm(j) = xm(j) - 2.0e-5_dp
        call geometry(pm, xp, 1, lf)
        bp = lf%gauge + 0.5_dp*lf%tau*lf%b
        call geometry(pm, xm, 1, lf)
        bm = lf%gauge + 0.5_dp*lf%tau*lf%b
        num(:, j) = (bp - bm)/4.0e-5_dp
    end do
    v = [num(3, 2) - num(2, 3), num(1, 3) - num(3, 1), &
         num(2, 1) - num(1, 2)]
    call require(norm2(v - cr) < 2.0e-10_dp, 'generated curl Rstar oracle')
    do order = 0, 1
        g = [x, 0.017_dp, 0.00013_dp, 0.37_dp]
        call lift(pm, g, order, y, ok)
        call require(ok, 'lift success')
        call inverse(pm, y, order, gn, ok, err)
        call require(ok .and. err < 3.0e-11_dp, 'Newton inverse convergence')
        call lift(pm, gn, order, yn, ok)
        call require(maxval(abs(yn - y)) < 1.0e-11_dp, 'lift inverse round trip')
        call require(abs(sum(y(4:6)**2)/(2*gc_energy(pm, g)) - 1) < 2.0e-15_dp, &
            'switch energy exact')
    end do
    call require(.not. kernel_available(2), 'order 2 pending, no silent fallback')
    unsupported%method = 3
    unsupported%order = 2
    call initialize(pm, unsupported, y, orbit)
    call require(orbit%failed, 'unavailable order 2 rejects orbit initialization')

    ! Full 3-D reference versus an independently reduced Hamiltonian prototype.
    pm%delta = 0.0_dp
    speed = 0.04_dp
    y = [10.49_dp, 0.0_dp, 0.0_dp, 0.035_dp, 0.016_dp, &
         sqrt(speed**2 - 0.035_dp**2 - 0.016_dp**2)]
    rr = y(1)
    psi = psi_value(pm%torus, 1, [rr, y(3)])
    pm%torus%p_phi = rr*y(5) + psi
    z = [rr, y(3), y(4), y(6) - pm%torus%f_tor*log(rr)]
    region = 1
    ht = 0.03_dp
    nc = 0
    call rk_advance(pm%torus, region, z, 300.0_dp, 1.0e-12_dp, ht, nc)
    opt%tol = 1.0e-11_dp
    call initialize(pm, opt, y, orbit)
    do i = 1, 600
        call advance(pm, opt, orbit, 0.5_dp)
    end do
    call particle(pm, opt, orbit, yn, ok)
    out(1:4) = [hypot(yn(1), yn(2)), yn(3), &
        dot_product(yn(1:2), yn(4:5))/hypot(yn(1), yn(2)), &
        yn(6) - pm%torus%f_tor*log(hypot(yn(1), yn(2)))]
    call require(maxval(abs(out(1:4) - z)) < 2.0e-8_dp, '3D vs reduced FO')
    call require(orbit%crossings > 10, 'reference repeated sheet events exercised')
    do i = 1, 600
        call advance(pm, opt, orbit, -0.5_dp)
    end do
    call particle(pm, opt, orbit, yn, ok)
    call require(norm2(yn - y)/speed < 3.0e-7_dp, 'FO backward-time reversibility')

    ! A crossing and return concentrated near a radial turn; both the end
    ! and half-step remain inside. This exercises the extremum event bracket.
    y = [10.4999999_dp, 0.0_dp, 0.0_dp, 0.0001_dp, 0.008_dp, 0.018_dp]
    opt%tol = 1.0e-6_dp
    call initialize(pm, opt, y, orbit)
    call advance(pm, opt, orbit, 0.05_dp)
    call require(orbit%crossings == 2, 'shallow crossing-return pair detected')

    x = [0.0_dp, 0.0_dp, 0.12_dp]
    v = [0.4_dp, 0.7_dp, -0.3_dp]
    b = [0.4_dp, 0.7_dp, 0.2_dp]
    call helix(x, v, b, 4.37_dp, xn, vn)
    call helix(xn, vn, b, -4.37_dp, xr, vr)
    call require(norm2(xr - x) + norm2(vr - v) < 1.0e-14_dp, 'helix inverse')
    n = [0.0_dp, 0.0_dp, 1.0_dp]
    p = 0.0_dp
    bm = [1.0_dp, 0.0_dp, 0.0_dp]
    bp = [0.7_dp, 1.1_dp, 0.0_dp]
    call plane_step(p, n, bm, bp, x, v, 70.0_dp, xn, vn, nc, ok)
    call require(ok .and. nc > 10, 'analytic plane repeated crossings')
    call require(abs(sum(vn**2) - sum(v**2)) < 1.0e-13_dp, 'plane energy')
    call plane_step(p, n, bm, bp, xn, vn, -70.0_dp, xr, vr, nc, ok)
    call require(norm2(xr - x) + norm2(vr - v) < 2.0e-11_dp, 'plane event inverse')
    x = [0.1_dp, -0.2_dp, 0.0_dp]
    v = [0.4_dp, 0.25_dp, 0.6_dp]
    bp = [1.04_dp, 0.78_dp, 0.0_dp]
    call plane_step(p, n, bm, bp, x, v, 2000.0_dp, &
        target_x, target_v, total_events, ok)
    call require(ok .and. total_events > 128, 'complete-cycle compression')
    xr = x
    vr = v
    events = 0
    do i = 1, 8000
        call plane_step(p, n, bm, bp, xr, vr, 0.25_dp, chunk_x, chunk_v, nc, ok)
        call require_quiet(ok, 'uncompressed short planar events')
        xr = chunk_x
        vr = chunk_v
        events = events + nc
    end do
    call require(events == total_events, 'compressed/direct crossing count')
    call require(norm2(xr - target_x) + norm2(vr - target_v) < 2.0e-9_dp, &
        'compressed cycles vs direct short events')

    ! General-field pullback scaling is measured here, not asserted as a proof.
    pm%delta = 0.1_dp
    print '(a)', 'pullback: order speed max_scaled_twoform_defect'
    do order = 0, 1
        do level = 1, 3
            speed = 0.01_dp*2**(level - 1)
            g = [10.3_dp, 0.0_dp, 0.2_dp, 0.4_dp*speed, &
                 0.42_dp*speed**2, 0.37_dp]
            call lift(pm, g, order, y, ok)
            scales = [speed, speed, speed, speed, speed**2, 1.0_dp]
            do j = 1, 6
                hh = 1.0e-4_dp*scales(j)
                gn = g
                gn(j) = gn(j) + hh
                call lift(pm, gn, order, yn, ok)
                gn = g
                gn(j) = gn(j) - hh
                call lift(pm, gn, order, out, ok)
                tangent(:, j) = (yn - out)/(2*hh)
            end do
            call geometry(pm, y(1:3), side_of(pm, y(1:3)), lf)
            call magnetic_form(lf%bvec, ofo)
            call geometry(pm, g(1:3), side_of(pm, g(1:3)), lf)
            rs = lf%gauge
            cr = 0.0_dp
            if (order == 1) call connection(pm, g(1:3), side_of(pm, g(1:3)), rs, cr)
            bst = lf%bvec + g(4)*lf%curlb - g(5)*cr
            call magnetic_form(bst, ogc)
            ogc(1:3, 4:6) = 0.0_dp
            ogc(4:6, 1:3) = 0.0_dp
            ogc(1:3, 4) = -lf%b
            ogc(4, 1:3) = lf%b
            ogc(1:3, 5) = rs
            ogc(5, 1:3) = -rs
            ogc(5, 6) = -1.0_dp
            ogc(6, 5) = 1.0_dp
            diff = matmul(transpose(tangent), matmul(ofo, tangent)) - ogc
            do i = 1, 6
                do j = 1, 6
                    diff(i, j) = diff(i, j)*scales(i)*scales(j)/speed**2
                end do
            end do
            defects(level) = maxval(abs(diff))
            print '(i2,2es17.8)', order, speed, defects(level)
        end do
        exponent = log(defects(3)/defects(1))/log(4.0_dp)
        call require(abs(exponent - real(order + 1, dp)) < 0.08_dp, &
            'measured two-form order in non-axisymmetric field')
    end do
    print '(a)', 'PASS standalone behavioral checks'
contains
    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL '//label
            error stop 1
        end if
        print '(a)', 'PASS '//label
    end subroutine
    subroutine require_quiet(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL '//label
            error stop 1
        end if
    end subroutine
    subroutine magnetic_form(b, om)
        real(dp), intent(in) :: b(3)
        real(dp), intent(out) :: om(6, 6)
        integer :: i
        om = 0.0_dp
        om(1, 2) = b(3)
        om(1, 3) = -b(2)
        om(2, 3) = b(1)
        om(2, 1) = -b(3)
        om(3, 1) = b(2)
        om(3, 2) = -b(1)
        do i = 1, 3
            om(i, i + 3) = -1
            om(i + 3, i) = 1
        end do
    end subroutine
end program
