program test_gc
    use gc_api, only: dp, field_derivs_t, gc_transform, gc_inverse, gc_rhs, &
        gc_transform_jacobian, gc_diagnostics, GC_OK, GC_ORDER, GC_FIELD, GC_POLAR, GC_FIRST_PLUS
    implicit none
    type(field_derivs_t) :: f, f2
    real(dp) :: e1(3), e2(3), localb(3), phi2
    real(dp) :: z(6), zz(6), xp(3), vp(3), xm(3), vm(3), jac(6, 6), fd(6, 6)
    real(dp) :: x0(3), v0(3), exactx(3), exactv(3), bh(3), perp(3), w(3)
    real(dp) :: xd(3), ud, md, pd, bs(3), bp, drift(3), coef(14), step, t, angle
    real(dp) :: residual(2, 2), zdot(6), field_particle(3)
    integer :: n, i, ierr, point, selector, scale
    integer, parameter :: selectors(4) = [0, 1, GC_FIRST_PLUS, 2]
    f%b = [0.0_dp, 0.0_dp, 1.7_dp]
    f%epsilon = 0.06_dp
    z = [0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, 0.25_dp, 0.73_dp]
    bh = [0.0_dp, 0.0_dp, 1.0_dp]
    do selector = 1, size(selectors)
        n = selectors(selector)
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
    f%d3b(3, 1, 1, 1) = 0.08_dp
    f%epsilon = 0.025_dp
    step = 1.0e-6_dp
    do selector = 1, size(selectors)
        n = selectors(selector)
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
    ! Independent differential constraints of the original particle model.
    ! A third-order position correction contributes to velocity at eps^2.
    z = [0.12_dp,-0.02_dp,0.01_dp,0.4_dp,0.3_dp,0.7_dp]
    do scale = 1, 2
        f%epsilon = 0.04_dp/2**scale
        call gc_transform(2,f,z(1:3),z(4),z(5),z(6),xp,vp,ierr)
        call check(ierr == GC_OK,'T2 Lagrangian transform status')
        call gc_transform_jacobian(2,f,z(1:3),z(4),z(5),z(6),jac,ierr)
        call check(ierr == GC_OK,'T2 Lagrangian Jacobian status')
        call gc_rhs(2,f,z(1:3),z(4),z(5),z(6), &
            xd,ud,md,pd,bs,bp,drift,ierr)
        call check(ierr == GC_OK,'T2 Lagrangian RHS status')
        zdot = [xd,ud,md,pd]
        residual(1,scale) = norm2(matmul(jac(1:3,:),zdot)-vp)
        field_particle = polynomial_field(f,xp)
        residual(2,scale) = norm2(matmul(jac(4:6,:),zdot) &
            - cross_product(vp,field_particle)/f%epsilon)
    end do
    call check(residual(1,2) < 0.16_dp*residual(1,1), &
        'T2 position-velocity Lagrangian residual order three')
    call check(residual(2,2) < 0.30_dp*residual(2,1), &
        'T2 Lorentz residual order two')
    f%epsilon = 0.025_dp
    ! Independent chart change with the same physical perpendicular direction.
    f2 = f
    f2%gauge_axis = 2
    localb = f%b + matmul(f%db, z(1:3))
    localb(3) = localb(3) + 0.5_dp*f%d2b(3, 1, 1)*z(1)**2 &
                + f%d3b(3, 1, 1, 1)*z(1)**3/6
    bh = localb/norm2(localb)
    call gc_transform(0, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    perp = vp - z(4)*bh
    e1 = [0.0_dp, 1.0_dp, 0.0_dp] - bh(2)*bh
    e1 = e1/norm2(e1)
    e2 = [bh(2)*e1(3) - bh(3)*e1(2), bh(3)*e1(1) - bh(1)*e1(3), &
          bh(1)*e1(2) - bh(2)*e1(1)]
    phi2 = atan2(dot_product(perp, e2), dot_product(perp, e1))
    do selector = 2, size(selectors)
        n = selectors(selector)
        call gc_transform(n, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
        call check(ierr == GC_OK, 'original gauge')
        call gc_transform(n, f2, z(1:3), z(4), z(5), phi2, xm, vm, ierr)
        call check(ierr == GC_OK, 'changed gauge')
        call check(norm2(xp-xm)+norm2(vp-vm) < 2.0e-12_dp, &
            'transform gyrogauge invariance')
    end do
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
    call gc_rhs(GC_FIRST_PLUS,f2,zz(1:3),zz(4),zz(5),zz(6), &
        xd,ud,md,pd,bs,bp,drift,ierr)
    call check(ierr == GC_OK,'slab H2 RHS status')
    call check(abs(pd + 1.3_dp/f2%epsilon &
        - f2%epsilon*1.5_dp*zz(5)*0.4_dp**2/1.3_dp**2) < 1.0e-12_dp, &
        'independent slab second-order frequency')
    f2%b(3) = 1.0_dp
    f2%db = 0.0_dp
    f2%db(2,1) = 0.4_dp
    call gc_rhs(1,f2,zz(1:3),zz(4),zz(5),zz(6),xd,ud,md,pd,bs,bp,drift,ierr)
    call check(ierr == GC_OK,'shear Banos status')
    call check(abs(pd+1.0_dp/f2%epsilon+0.2_dp*zz(4)) < 1.0e-12_dp, &
        'shear Banos phase coefficient')
    call gc_transform(3, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    call check(ierr == GC_ORDER, 'unsupported T3 explicitly refused')
    call gc_inverse(3, f, x0, v0, zz(1:3), zz(4), zz(5), zz(6), ierr)
    call check(ierr == GC_ORDER, 'unsupported inverse T3 refused')
    call gc_transform_jacobian(3, f, z(1:3), z(4), z(5), z(6), jac, ierr)
    call check(ierr == GC_ORDER, 'unsupported Jacobian T3 refused')
    call gc_rhs(3, f, z(1:3), z(4), z(5), z(6), xd, ud, md, pd, bs, bp, drift, ierr)
    call check(ierr == GC_ORDER, 'unsupported RHS3 refused')
    call gc_transform(1, f, z(1:3), z(4), 0.0_dp, z(6), xp, vp, ierr)
    call check(ierr == GC_POLAR, 'polar singularity refused')
    f%b = 0.0_dp
    f%db = 0.0_dp
    f%d2b = 0.0_dp
    f%d3b = 0.0_dp
    call gc_transform(0, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    call check(ierr == GC_FIELD, 'zero field refused')
    f%b = [1.0_dp, 0.0_dp, 0.0_dp]
    call gc_transform(0, f, z(1:3), z(4), z(5), z(6), xp, vp, ierr)
    call check(ierr == GC_FIELD, 'singular gauge chart refused')
    print '(a)', 'PASS: all selectors, cubic-jet inverses/Jacobians, gyrogauge, slab phase'
contains
    function polynomial_field(jet,x) result(value)
        type(field_derivs_t), intent(in) :: jet
        real(dp), intent(in) :: x(3)
        real(dp) :: value(3), dx(3)
        integer :: j,k,l
        dx = x-jet%origin
        value = jet%b+matmul(jet%db,dx)
        do j=1,3
            do k=1,3
                value = value+jet%d2b(:,j,k)*dx(j)*dx(k)/2
                do l=1,3
                    value = value+jet%d3b(:,j,k,l)*dx(j)*dx(k)*dx(l)/6
                end do
            end do
        end do
    end function polynomial_field

    function cross_product(a,b) result(c)
        real(dp), intent(in) :: a(3),b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3), &
             a(1)*b(2)-a(2)*b(1)]
    end function cross_product

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_gc
