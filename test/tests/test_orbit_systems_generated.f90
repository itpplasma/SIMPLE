program test_orbit_systems_generated
    ! The generated residuals and Jacobians of the implicit symplectic schemes
    ! against the hand-derived formulas they replaced (SIMPLE af9c3346), on
    ! random field values. Every field_can_t entry is drawn independently, so
    ! a missing, extra or misindexed term changes the result at O(1).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: field_can_t
    use orbit_symplectic_base, only: coeff_rk_gauss, coeff_rk_lobatto
    use orbit_systems_generated
    implicit none

    real(dp), parameter :: tol = 1d-12
    integer :: errors, trial, s

    errors = 0
    call random_seed_fixed()
    do trial = 1, 20
        call check_euler1()
        call check_euler2()
        call check_midpoint()
        do s = 1, 4
            call check_gauss(s)
        end do
        call check_lobatto()
        call check_rates()
    end do

    if (errors /= 0) then
        print *, 'test_orbit_systems_generated: ', errors, ' mismatches'
        error stop 1
    end if
    print *, 'test_orbit_systems_generated: all checks passed'

contains

    subroutine random_seed_fixed()
        integer :: n
        integer, allocatable :: seed(:)
        call random_seed(size=n)
        allocate (seed(n))
        seed = 20261103
        call random_seed(put=seed)
    end subroutine random_seed_fixed

    real(dp) function uniform(lo, hi)
        real(dp), intent(in) :: lo, hi
        call random_number(uniform)
        uniform = lo + (hi - lo)*uniform
    end function uniform

    subroutine random_field(f)
        type(field_can_t), intent(out) :: f
        integer :: k

        f%H = uniform(-1d0, 1d0)
        f%pth = uniform(-1d0, 1d0)
        f%vpar = uniform(-1d0, 1d0)
        f%hth = uniform(-1d0, 1d0)
        f%hph = uniform(0.5d0, 1.5d0)
        do k = 1, 4
            f%dH(k) = uniform(-1d0, 1d0)
            f%dpth(k) = uniform(-1d0, 1d0)
            f%dvpar(k) = uniform(-1d0, 1d0)
        end do
        f%dpth(1) = uniform(0.5d0, 1.5d0)
        do k = 1, 3
            f%dhth(k) = uniform(-1d0, 1d0)
            f%dhph(k) = uniform(-1d0, 1d0)
        end do
        do k = 1, 10
            f%d2H(k) = uniform(-1d0, 1d0)
            f%d2pth(k) = uniform(-1d0, 1d0)
            f%d2vpar(k) = uniform(-1d0, 1d0)
        end do
    end subroutine random_field

    subroutine random_vector(v)
        real(dp), intent(out) :: v(:)
        integer :: k
        do k = 1, size(v)
            v(k) = uniform(-1d0, 1d0)
        end do
    end subroutine random_vector

    subroutine compare(name, new, old)
        character(*), intent(in) :: name
        real(dp), intent(in) :: new(:), old(:)
        real(dp) :: scale

        scale = max(1d0, maxval(abs(old)))
        if (maxval(abs(new - old)) > tol*scale) then
            errors = errors + 1
            print '(a,a,es10.2)', name, ': max |new - old| = ', maxval(abs(new - old))
        end if
    end subroutine compare

    subroutine check_euler1()
        type(field_can_t) :: f
        real(dp) :: dt, pthold, z(4), x(2), fvec(2), jac(2, 2), ref(2), rjac(2, 2)

        call random_field(f)
        dt = uniform(0.1d0, 1d0)
        pthold = uniform(-1d0, 1d0)
        call random_vector(z)
        call random_vector(x)
        call euler1_residual(dt, pthold, z, x, f, fvec)
        call euler1_jacobian(dt, pthold, z, x, f, jac)

        ref(1) = f%dpth(1)*(f%pth - pthold) + dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
        ref(2) = f%dpth(1)*(x(2) - z(4)) + dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
        rjac(1, 1) = f%d2pth(1)*(f%pth - pthold) + f%dpth(1)**2 &
                     + dt*(f%d2H(2)*f%dpth(1) + f%dH(2)*f%d2pth(1) &
                           - f%d2H(1)*f%dpth(2) - f%dH(1)*f%d2pth(2))
        rjac(1, 2) = f%d2pth(7)*(f%pth - pthold) + f%dpth(1)*f%dpth(4) &
                     + dt*(f%d2H(8)*f%dpth(1) + f%dH(2)*f%d2pth(7) &
                           - f%d2H(7)*f%dpth(2) - f%dH(1)*f%d2pth(8))
        rjac(2, 1) = f%d2pth(1)*(x(2) - z(4)) &
                     + dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) &
                           - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
        rjac(2, 2) = f%d2pth(7)*(x(2) - z(4)) + f%dpth(1) &
                     + dt*(f%d2H(9)*f%dpth(1) + f%dH(3)*f%d2pth(7) &
                           - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9))
        call compare('euler1 residual', fvec, ref)
        call compare('euler1 jacobian', reshape(jac, [4]), reshape(rjac, [4]))
    end subroutine check_euler1

    subroutine check_euler2()
        type(field_can_t) :: f
        real(dp) :: dt, pthold, z(4), x(3), fvec(3), jac(3, 3), ref(3), rjac(3, 3)

        call random_field(f)
        dt = uniform(0.1d0, 1d0)
        pthold = uniform(-1d0, 1d0)
        call random_vector(z)
        call random_vector(x)
        call euler2_residual(dt, pthold, z, x, f, fvec)
        call euler2_jacobian(dt, z, x, f, jac)

        ref(1) = f%pth - pthold
        ref(2) = f%dpth(1)*(x(2) - z(2)) - dt*f%dH(1)
        ref(3) = f%dpth(1)*f%hph*(x(3) - z(3)) - dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)
        rjac(1, :) = f%dpth(1:3)
        rjac(2, :) = f%d2pth(1:3)*(x(2) - z(2)) - dt*f%d2H(1:3)
        rjac(2, 2) = rjac(2, 2) + f%dpth(1)
        rjac(3, :) = (f%d2pth(1:3)*f%hph + f%dpth(1)*f%dhph)*(x(3) - z(3)) &
                     - dt*(f%d2pth(1:3)*f%vpar + f%dpth(1)*f%dvpar(1:3) &
                           - f%d2H(1:3)*f%hth - f%dH(1)*f%dhth)
        rjac(3, 3) = rjac(3, 3) + f%dpth(1)*f%hph
        call compare('euler2 residual', fvec, ref)
        call compare('euler2 jacobian', reshape(jac, [9]), reshape(rjac, [9]))
    end subroutine check_euler2

    subroutine check_midpoint()
        type(field_can_t) :: fm, fe
        real(dp) :: dt, pthold, z(4), x(5), fvec(5), jac(5, 5), ref(5), rjac(5, 5)

        call random_field(fm)
        call random_field(fe)
        dt = uniform(0.1d0, 1d0)
        pthold = uniform(-1d0, 1d0)
        call random_vector(z)
        call random_vector(x)
        call midpoint_residual(dt, pthold, z, x, fm, fe, fvec)
        call midpoint_jacobian(dt, pthold, z, x, fm, fe, jac)

  ref(1) = fm%dpth(1)*(fe%pth - pthold) + dt*(fm%dpth(1)*fm%dH(2) - fm%dpth(2)*fm%dH(1))
        ref(2) = fm%dpth(1)*(x(2) - z(2)) - dt*fm%dH(1)
    ref(3) = fm%dpth(1)*fm%hph*(x(3) - z(3)) - dt*(fm%dpth(1)*fm%vpar - fm%dH(1)*fm%hth)
      ref(4) = fm%dpth(1)*(x(4) - z(4)) + dt*(fm%dH(3)*fm%dpth(1) - fm%dH(1)*fm%dpth(3))
        ref(5) = fm%dpth(1)*(fm%pth - pthold) &
                 + 0.5d0*dt*(fm%dpth(1)*fm%dH(2) - fm%dpth(2)*fm%dH(1))
        call old_midpoint_jacobian(dt, pthold, z, x, fm, fe, rjac)
        call compare('midpoint residual', fvec, ref)
        call compare('midpoint jacobian', reshape(jac, [25]), reshape(rjac, [25]))
    end subroutine check_midpoint

    subroutine old_midpoint_jacobian(dt, pthold, z, x, f, fend, jac)
        real(dp), intent(in) :: dt, pthold, z(4), x(5)
        type(field_can_t), intent(in) :: f, fend
        real(dp), intent(out) :: jac(5, 5)

        jac(2, 1) = 0d0
        jac(2, 5) = f%d2pth(1)*(x(2) - z(2)) - dt*f%d2H(1)
        jac(2, 2:3) = 0.5d0*(f%d2pth(2:3)*(x(2) - z(2)) - dt*f%d2H(2:3))
        jac(2, 2) = jac(2, 2) + f%dpth(1)
        jac(2, 4) = 0.5d0*(f%d2pth(7)*(x(2) - z(2)) - dt*f%d2H(7))
        jac(3, 1) = 0d0
        jac(3, 5) = (f%d2pth(1)*f%hph + f%dpth(1)*f%dhph(1))*(x(3) - z(3)) &
                    - dt*(f%d2pth(1)*f%vpar + f%dpth(1)*f%dvpar(1) &
                          - f%d2H(1)*f%hth - f%dH(1)*f%dhth(1))
       jac(3, 2:3) = 0.5d0*((f%d2pth(2:3)*f%hph + f%dpth(1)*f%dhph(2:3))*(x(3) - z(3)) &
                             - dt*(f%d2pth(2:3)*f%vpar + f%dpth(1)*f%dvpar(2:3) &
                                   - f%d2H(2:3)*f%hth - f%dH(1)*f%dhth(2:3)))
        jac(3, 3) = jac(3, 3) + f%dpth(1)*f%hph
        jac(3, 4) = 0.5d0*(f%d2pth(7)*f%hph*(x(3) - z(3)) &
                       - dt*(f%d2pth(7)*f%vpar + f%dpth(1)*f%dvpar(4) - f%d2H(7)*f%hth))
        jac(4, 1) = 0d0
        jac(4, 5) = f%d2pth(1)*(x(4) - z(4)) + dt*(f%d2H(3)*f%dpth(1) &
                         + f%dH(3)*f%d2pth(1) - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
        jac(4, 2) = 0.5d0*(f%d2pth(2)*(x(4) - z(4)) + dt*(f%d2H(5)*f%dpth(1) &
                        + f%dH(3)*f%d2pth(2) - f%d2H(2)*f%dpth(3) - f%dH(1)*f%d2pth(5)))
        jac(4, 3) = 0.5d0*(f%d2pth(3)*(x(4) - z(4)) + dt*(f%d2H(6)*f%dpth(1) &
                        + f%dH(3)*f%d2pth(3) - f%d2H(3)*f%dpth(3) - f%dH(1)*f%d2pth(6)))
        jac(4, 4) = 0.5d0*(f%d2pth(7)*(x(4) - z(4)) + dt*(f%dH(3)*f%d2pth(7) &
            + f%d2H(9)*f%dpth(1) - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9))) + f%dpth(1)
        jac(5, 1) = 0d0
        jac(5, 5) = f%d2pth(1)*(f%pth - pthold) + f%dpth(1)*f%dpth(1) &
                    + 0.5d0*dt*(f%d2pth(1)*f%dH(2) + f%dpth(1)*f%d2H(2) &
                                - f%d2pth(2)*f%dH(1) - f%dpth(2)*f%d2H(1))
        jac(5, 2) = 0.5d0*(f%d2pth(2)*(f%pth - pthold) + f%dpth(1)*f%dpth(2) &
                           + 0.5d0*dt*(f%d2pth(2)*f%dH(2) + f%dpth(1)*f%d2H(4) &
                                       - f%d2pth(4)*f%dH(1) - f%dpth(2)*f%d2H(2)))
        jac(5, 3) = 0.5d0*(f%d2pth(3)*(f%pth - pthold) + f%dpth(1)*f%dpth(3) &
                           + 0.5d0*dt*(f%d2pth(3)*f%dH(2) + f%dpth(1)*f%d2H(5) &
                                       - f%d2pth(5)*f%dH(1) - f%dpth(2)*f%d2H(3)))
        jac(5, 4) = 0.5d0*(f%d2pth(7)*(f%pth - pthold) + f%dpth(1)*f%dpth(4) &
                           + 0.5d0*dt*(f%d2pth(7)*f%dH(2) + f%dpth(1)*f%d2H(8) &
                                       - f%d2pth(8)*f%dH(1) - f%dpth(2)*f%d2H(7)))
        jac(1, 1) = f%dpth(1)*fend%dpth(1)
        jac(1, 2) = 0.5d0*(f%d2pth(2)*(fend%pth - pthold) + dt*(f%d2pth(2)*f%dH(2) &
                      + f%dpth(1)*f%d2H(4) - f%dpth(2)*f%d2H(2) - f%d2pth(4)*f%dH(1))) &
                    + f%dpth(1)*fend%dpth(2)
        jac(1, 3) = 0.5d0*(f%d2pth(3)*(fend%pth - pthold) + dt*(f%d2pth(3)*f%dH(2) &
                      + f%dpth(1)*f%d2H(5) - f%dpth(2)*f%d2H(3) - f%d2pth(5)*f%dH(1))) &
                    + f%dpth(1)*fend%dpth(3)
        jac(1, 4) = 0.5d0*(f%d2pth(7)*(fend%pth - pthold) + dt*(f%d2pth(7)*f%dH(2) &
                      + f%dpth(1)*f%d2H(8) - f%dpth(2)*f%d2H(7) - f%d2pth(8)*f%dH(1))) &
                    + f%dpth(1)*fend%dpth(4)
        jac(1, 5) = f%d2pth(1)*(fend%pth - pthold) + dt*(f%d2pth(1)*f%dH(2) &
                         + f%dpth(1)*f%d2H(2) - f%dpth(2)*f%d2H(1) - f%d2pth(2)*f%dH(1))
    end subroutine old_midpoint_jacobian

    !> Contribution of stage l to the rates and their gradients with respect
    !> to the stage state (r, theta, phi, p_phi), as in the former
    !> gauss_canfield_jacobian and jac_rk_lobatto.
    subroutine stage_rates(f, hp, pthd, thd, phd, pphd, dpthd, dthd, dphd, dpphd)
        type(field_can_t), intent(in) :: f
        real(dp), intent(out) :: hp, pthd, thd, phd, pphd
        real(dp), intent(out) :: dpthd(4), dthd(4), dphd(4), dpphd(4)
        real(dp) :: dhp(4)
        integer, parameter :: i_r(4) = [1, 2, 3, 7], i_th(4) = [2, 4, 5, 8], &
                              i_ph(4) = [3, 5, 6, 9]
        integer :: i

        hp = f%dH(1)/f%dpth(1)
        do i = 1, 4
            dhp(i) = (f%d2H(i_r(i)) - hp*f%d2pth(i_r(i)))/f%dpth(1)
        end do
        pthd = f%dH(2) - hp*f%dpth(2)
        thd = hp
        phd = (f%vpar - hp*f%hth)/f%hph
        pphd = f%dH(3) - hp*f%dpth(3)
        do i = 1, 4
            dpthd(i) = f%d2H(i_th(i)) - f%d2pth(i_th(i))*hp - f%dpth(2)*dhp(i)
            dthd(i) = dhp(i)
            dpphd(i) = f%d2H(i_ph(i)) - f%d2pth(i_ph(i))*hp - f%dpth(3)*dhp(i)
        end do
        do i = 1, 3
            dphd(i) = -f%dhph(i)*(f%vpar - hp*f%hth)/f%hph**2 &
                      + (f%dvpar(i) - dhp(i)*f%hth - hp*f%dhth(i))/f%hph
        end do
        dphd(4) = (f%dvpar(4) - dhp(4)*f%hth)/f%hph
    end subroutine stage_rates

    subroutine check_gauss(s)
        integer, intent(in) :: s
        type(field_can_t) :: fs(s)
        real(dp) :: dt, pthold, z(4), x(4*s), fvec(4*s), jac(4*s, 4*s)
        real(dp) :: ref(4*s), rjac(4*s, 4*s), a(s, s), b(s), c(s)
        real(dp) :: hp(s), pthd(s), thd(s), phd(s), pphd(s)
        real(dp) :: dpthd(4, s), dthd(4, s), dphd(4, s), dpphd(4, s)
        integer :: k, l, m, n
        character(8) :: name

        do k = 1, s
            call random_field(fs(k))
            call stage_rates(fs(k), hp(k), pthd(k), thd(k), phd(k), pphd(k), &
                             dpthd(:, k), dthd(:, k), dphd(:, k), dpphd(:, k))
        end do
        dt = uniform(0.1d0, 1d0)
        pthold = uniform(-1d0, 1d0)
        call random_vector(z)
        call random_vector(x)
        call coeff_rk_gauss(s, a, b, c)
        select case (s)
        case (1)
            call gauss1_residual(dt, pthold, z, a, x, fs, fvec)
            call gauss1_jacobian(dt, a, fs, jac)
        case (2)
            call gauss2_residual(dt, pthold, z, a, x, fs, fvec)
            call gauss2_jacobian(dt, a, fs, jac)
        case (3)
            call gauss3_residual(dt, pthold, z, a, x, fs, fvec)
            call gauss3_jacobian(dt, a, fs, jac)
        case (4)
            call gauss4_residual(dt, pthold, z, a, x, fs, fvec)
            call gauss4_jacobian(dt, a, fs, jac)
        end select

        rjac = 0d0
        do k = 1, s
            m = 4*k
            ref(m - 3) = fs(k)%pth - pthold
            ref(m - 2) = x(m - 2) - z(2)
            ref(m - 1) = x(m - 1) - z(3)
            ref(m) = x(m) - z(4)
            rjac(m - 3, m - 3:m) = fs(k)%dpth(1:4)
            rjac(m - 2, m - 2) = 1d0
            rjac(m - 1, m - 1) = 1d0
            rjac(m, m) = 1d0
            do l = 1, s
                n = 4*l
                ref(m - 3) = ref(m - 3) + dt*a(k, l)*pthd(l)
                ref(m - 2) = ref(m - 2) - dt*a(k, l)*thd(l)
                ref(m - 1) = ref(m - 1) - dt*a(k, l)*phd(l)
                ref(m) = ref(m) + dt*a(k, l)*pphd(l)
                rjac(m - 3, n - 3:n) = rjac(m - 3, n - 3:n) + dt*a(k, l)*dpthd(:, l)
                rjac(m - 2, n - 3:n) = rjac(m - 2, n - 3:n) - dt*a(k, l)*dthd(:, l)
                rjac(m - 1, n - 3:n) = rjac(m - 1, n - 3:n) - dt*a(k, l)*dphd(:, l)
                rjac(m, n - 3:n) = rjac(m, n - 3:n) + dt*a(k, l)*dpphd(:, l)
            end do
        end do
        write (name, '(a,i0)') 'gauss', s
        call compare(trim(name)//' residual', fvec, ref)
        call compare(trim(name)//' jacobian', reshape(jac, [16*s*s]), &
                     reshape(rjac, [16*s*s]))
    end subroutine check_gauss

    subroutine check_lobatto()
        integer, parameter :: s = 3, n = 4*s - 2
        type(field_can_t) :: fs(s)
        real(dp) :: dt, pthold, z(4), x(n), fvec(n), jac(n, n), ref(n), rjac(n, n)
        real(dp) :: a(s, s), ahat(s, s), b(s), c(s)
        real(dp) :: hp(s), pthd(s), thd(s), phd(s), pphd(s)
        real(dp) :: dpthd(4, s), dthd(4, s), dphd(4, s), dpphd(4, s)
        integer :: cols(4, s), k, l, m, i

        do k = 1, s
            call random_field(fs(k))
            call stage_rates(fs(k), hp(k), pthd(k), thd(k), phd(k), pphd(k), &
                             dpthd(:, k), dthd(:, k), dphd(:, k), dpphd(:, k))
        end do
        dt = uniform(0.1d0, 1d0)
        pthold = uniform(-1d0, 1d0)
        call random_vector(z)
        call random_vector(x)
        call coeff_rk_lobatto(s, a, ahat, b, c)
        call lobatto3_residual(dt, pthold, z, a, ahat, x, fs, fvec)
        call lobatto3_jacobian(dt, a, ahat, fs, jac)

        ! Unknown column of each stage state component; the first stage has
        ! the angles of z (column 0: not an unknown).
        cols(:, 1) = [1, 0, 0, 2]
        do k = 2, s
            cols(:, k) = [4*k - 5, 4*k - 4, 4*k - 3, 4*k - 2]
        end do

        rjac = 0d0
        ref(1) = fs(1)%pth - pthold
        ref(2) = x(2) - z(4)
        rjac(1, 1) = fs(1)%dpth(1)
        rjac(1, 2) = fs(1)%dpth(4)
        rjac(2, 2) = 1d0
        do k = 2, s
            m = 4*k - 2
            ref(m - 3) = fs(k)%pth - pthold
            ref(m - 2) = x(m - 2) - z(2)
            ref(m - 1) = x(m - 1) - z(3)
            ref(m) = x(m) - z(4)
            rjac(m - 3, m - 3:m) = fs(k)%dpth(1:4)
            rjac(m - 2, m - 2) = 1d0
            rjac(m - 1, m - 1) = 1d0
            rjac(m, m) = 1d0
        end do
        do l = 1, s
            ref(1) = ref(1) + dt*ahat(1, l)*pthd(l)
            ref(2) = ref(2) + dt*ahat(1, l)*pphd(l)
            do k = 2, s
                m = 4*k - 2
                ref(m - 3) = ref(m - 3) + dt*ahat(k, l)*pthd(l)
                ref(m - 2) = ref(m - 2) - dt*a(k, l)*thd(l)
                ref(m - 1) = ref(m - 1) - dt*a(k, l)*phd(l)
                ref(m) = ref(m) + dt*ahat(k, l)*pphd(l)
            end do
            do i = 1, 4
                if (cols(i, l) == 0) cycle
                rjac(1, cols(i, l)) = rjac(1, cols(i, l)) + dt*ahat(1, l)*dpthd(i, l)
                rjac(2, cols(i, l)) = rjac(2, cols(i, l)) + dt*ahat(1, l)*dpphd(i, l)
                do k = 2, s
                    m = 4*k - 2
                    rjac(m - 3, cols(i, l)) = rjac(m - 3, cols(i, l)) &
                                              + dt*ahat(k, l)*dpthd(i, l)
                    rjac(m - 2, cols(i, l)) = rjac(m - 2, cols(i, l)) &
                                              - dt*a(k, l)*dthd(i, l)
                    rjac(m - 1, cols(i, l)) = rjac(m - 1, cols(i, l)) &
                                              - dt*a(k, l)*dphd(i, l)
                   rjac(m, cols(i, l)) = rjac(m, cols(i, l)) + dt*ahat(k, l)*dpphd(i, l)
                end do
            end do
        end do
        call compare('lobatto3 residual', fvec, ref)
        call compare('lobatto3 jacobian', reshape(jac, [n*n]), reshape(rjac, [n*n]))
    end subroutine check_lobatto

    subroutine check_rates()
        type(field_can_t) :: f
        real(dp) :: zdot(4), ref(4), hp

        call random_field(f)
        call gc_rates_canonical(f, zdot)
        hp = f%dH(1)/f%dpth(1)
        ref(1) = -(f%dH(2) - f%hth/f%hph*f%dH(3))/f%dpth(1)
        ref(2) = hp
        ref(3) = (f%vpar - hp*f%hth)/f%hph
        ref(4) = -(f%dH(3) - hp*f%dpth(3))
        call compare('guiding-centre rates', zdot, ref)
    end subroutine check_rates

end program test_orbit_systems_generated
