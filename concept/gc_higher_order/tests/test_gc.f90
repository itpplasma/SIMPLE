program test_gc
    use gc_api, only: dp, field_derivs_t, gc_transform, gc_inverse, gc_rhs, &
        gc_transform_jacobian, gc_diagnostics, GC_OK, GC_ORDER, GC_FIELD, GC_POLAR
    implicit none
    type(field_derivs_t) :: f, f2
    real(dp) :: e1(3), e2(3), localb(3), phi2
    real(dp) :: z(6), zz(6), xp(3), vp(3), xm(3), vm(3), jac(6, 6), fd(6, 6)
    real(dp) :: x0(3), v0(3), exactx(3), exactv(3), bh(3), perp(3), w(3)
    real(dp) :: xd(3), ud, md, pd, bs(3), bp, drift(3), coef(14), step, t, angle
    integer :: n, i, ierr, point
    f%b = [0.0_dp, 0.0_dp, 1.7_dp]
    f%epsilon = 0.06_dp
    z = [0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, 0.25_dp, 0.73_dp]
    bh = [0.0_dp, 0.0_dp, 1.0_dp]
    do n = 0, 1
        call gc_transform(n, f, z(1:3), z(4), z(5), z(6), x0, v0, ierr)
        call check(ierr == GC_OK, 'uniform transform')
        call gc_inverse(n, f, x0, v0, zz(1:3), zz(4), zz(5), zz(6), ierr)
        call check(ierr == GC_OK, 'uniform inverse status')
        call check(norm2(zz - z) < 2.0e-11_dp, 'uniform round trip')
        call gc_rhs(n, f, z(1:3), z(4), z(5), z(6), xd, ud, md, pd, bs, bp, drift, ierr)
        call check(ierr == GC_OK, 'uniform rhs')
        call check(norm2(xd - z(4)*bh) < 1.0e-14_dp, 'uniform parallel motion')
        call check(abs(pd + 1.7_dp/f%epsilon) < 1.0e-13_dp, 'uniform phase sign')
        t = 0.37_dp
        zz = z
        zz(1:3) = z(1:3) + t*xd
        zz(6) = z(6) + t*pd
        call gc_transform(n, f, zz(1:3), zz(4), zz(5), zz(6), xp, vp, ierr)
        perp = v0 - dot_product(v0, bh)*bh
        w = [-perp(2), perp(1), 0.0_dp]
        angle = 1.7_dp*t/f%epsilon
        exactv = z(4)*bh + cos(angle)*perp - sin(angle)*w
        exactx = x0 + t*z(4)*bh + f%epsilon/1.7_dp &
                 *(sin(angle)*perp + (cos(angle) - 1.0_dp)*w)
        call check(norm2(vp - exactv) < 1.0e-13_dp, 'uniform exact velocity')
        call check(norm2(xp - exactx) < 1.0e-13_dp, 'uniform exact position')
    end do
    ! Solenoidal linear field with varying direction and strength.
    f%b = [0.13_dp, 0.07_dp, 1.1_dp]
    f%db(:, 1) = [0.02_dp, 0.12_dp, 0.08_dp]
    f%db(:, 2) = [-0.04_dp, -0.01_dp, 0.03_dp]
    f%db(:, 3) = [0.02_dp, 0.01_dp, -0.01_dp]
    f%d2b(3, 1, 1) = 0.025_dp
    f%epsilon = 0.025_dp
    step = 1.0e-6_dp
    do n = 0, 1
        do point = 1, 7
            z = [0.03_dp*point, -0.02_dp, 0.01_dp, 0.4_dp, 0.3_dp, 0.2_dp*point]
            call gc_transform(n, f, z(1:3), z(4), z(5), z(6), x0, v0, ierr)
            call check(ierr == GC_OK, 'nonuniform transform')
            call gc_inverse(n, f, x0, v0, zz(1:3), zz(4), zz(5), zz(6), ierr)
            call check(ierr == GC_OK, 'nonuniform inverse')
            call check(norm2(zz - z) < 5.0e-10_dp, 'nonuniform round trip')
            call gc_transform_jacobian(n, f, z(1:3), z(4), z(5), z(6), jac, ierr)
            call check(ierr == GC_OK, 'Jacobian status')
            do i = 1, 6
                zz = z
                zz(i) = zz(i) + step
                call gc_transform(n, f, zz(1:3), zz(4), zz(5), zz(6), xp, vp, ierr)
                call check(ierr == GC_OK, 'positive finite difference')
                zz(i) = z(i) - step
                call gc_transform(n, f, zz(1:3), zz(4), zz(5), zz(6), xm, vm, ierr)
                call check(ierr == GC_OK, 'negative finite difference')
                fd(:, i) = ([xp, vp] - [xm, vm])/(2*step)
            end do
            call check(maxval(abs(jac - fd)) < 2.0e-8_dp, 'generated analytic Jacobian')
            call gc_rhs(n, f, z(1:3), z(4), z(5), z(6), &
                xd, ud, md, pd, bs, bp, drift, ierr)
            call check(ierr == GC_OK, 'nonuniform rhs')
        end do
    end do
    ! Independent chart change with the same physical perpendicular direction.
    f2 = f
    f2%gauge_axis = 2
    localb = f%b + matmul(f%db, z(1:3))
    localb(3) = localb(3) + 0.5_dp*f%d2b(3, 1, 1)*z(1)**2
    bh = localb/norm2(localb)
    call gc_transform(0, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    perp = vp - z(4)*bh
    e1 = [0.0_dp, 1.0_dp, 0.0_dp] - bh(2)*bh
    e1 = e1/norm2(e1)
    e2 = [bh(2)*e1(3) - bh(3)*e1(2), bh(3)*e1(1) - bh(1)*e1(3), &
          bh(1)*e1(2) - bh(2)*e1(1)]
    phi2 = atan2(dot_product(perp, e2), dot_product(perp, e1))
    call gc_transform(1, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    call check(ierr == GC_OK, 'original gauge')
    call gc_transform(1, f2, z(1:3), z(4), z(5), phi2, xm, vm, ierr)
    call check(ierr == GC_OK, 'changed gauge')
    call check(norm2(xp-xm)+norm2(vp-vm) < 2.0e-12_dp,'transform gyrogauge invariance')
    call gc_diagnostics(f, z(1:3), z(4), z(5), z(6), coef, ierr)
    call check(ierr == GC_OK, 'second-order diagnostic status')
    ! Independently derived linear-slab H2 and averaged G2 coefficients.
    f2 = field_derivs_t()
    f2%b(3) = 1.3_dp
    f2%db(3,1) = 0.4_dp
    f2%epsilon = 0.02_dp
    zz = [0.0_dp,0.0_dp,0.0_dp,0.35_dp,0.25_dp,0.7_dp]
    call gc_diagnostics(f2,zz(1:3),zz(4),zz(5),zz(6),coef,ierr)
    call check(ierr == GC_OK,'slab diagnostics')
    call check(abs(coef(10)+0.75_dp*zz(5)**2*0.4_dp**2/1.3_dp**2) &
        < 1.0e-14_dp,'independent slab normal-form H2')
    call check(abs(coef(11)-zz(5)**2*0.4_dp**2/(4*1.3_dp**3)) &
        < 1.0e-14_dp,'slab averaged G2_mu')
    f2%b(3) = 1.0_dp
    f2%db = 0.0_dp
    f2%db(2,1) = 0.4_dp
    call gc_rhs(1,f2,zz(1:3),zz(4),zz(5),zz(6),xd,ud,md,pd,bs,bp,drift,ierr)
    call check(ierr == GC_OK,'shear Banos status')
    call check(abs(pd+1.0_dp/f2%epsilon+0.2_dp*zz(4)) < 1.0e-12_dp, &
        'shear Banos phase coefficient')
    call gc_transform(2, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    call check(ierr == GC_ORDER, 'unimplemented T2 explicitly refused')
    call gc_inverse(2, f, x0, v0, zz(1:3), zz(4), zz(5), zz(6), ierr)
    call check(ierr == GC_ORDER, 'unimplemented inverse T2 refused')
    call gc_transform_jacobian(2, f, z(1:3), z(4), z(5), z(6), jac, ierr)
    call check(ierr == GC_ORDER, 'unimplemented Jacobian T2 refused')
    call gc_rhs(2, f, z(1:3), z(4), z(5), z(6), xd, ud, md, pd, bs, bp, drift, ierr)
    call check(ierr == GC_ORDER, 'unimplemented RHS2 refused')
    call gc_transform(1, f, z(1:3), z(4), 0.0_dp, z(6), xp, vp, ierr)
    call check(ierr == GC_POLAR, 'polar singularity refused')
    f%b = 0.0_dp
    f%db = 0.0_dp
    f%d2b = 0.0_dp
    call gc_transform(0, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    call check(ierr == GC_FIELD, 'zero field refused')
    f%b = [1.0_dp, 0.0_dp, 0.0_dp]
    call gc_transform(0, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    call check(ierr == GC_FIELD, 'singular gauge chart refused')
    print '(a)', 'PASS: uniform exactness, both inverses, generated Jacobians, refusals'
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_gc
