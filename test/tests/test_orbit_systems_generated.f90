program test_orbit_systems_generated
    ! The generated residuals and Jacobians of the implicit symplectic schemes
    ! against the hand-derived formulas they replaced (SIMPLE af9c3346), on
    ! random field values. Every field_can_t entry is drawn independently, so
    ! a missing, extra or misindexed term changes the result at O(1).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: field_can_t
    use orbit_symplectic_base, only: coeff_rk_gauss, coeff_rk_lobatto
    use orbit_systems_generated
    use field_can_mod, only: get_derivatives
    use simple_gpu, only: gpu_get_derivatives2_radial_mixed, gc_rates_boozer
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
        call check_radial_mixed()
        call check_rates_boozer()
        call check_geometry_transport()
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

    subroutine check_geometry_transport()
        type(field_can_t) :: f
        real(dp) :: state(3), delta(3), old(3, 5), expected(3, 5)
        real(dp) :: hessian(6, 5), actual(3, 5)
        integer :: q

        call random_vector(state)
        call random_vector(delta)
        do q = 1, 5
            old(:, q) = quadratic_gradient(state, q)
            expected(:, q) = quadratic_gradient(state + delta, q)
            hessian(:, q) = q*[2d0, 1d0, 2d0, 4d0, 3d0, 6d0]
        end do
        f%dAth = old(:, 1)
        f%dAph = old(:, 2)
        f%dhth = old(:, 3)
        f%dhph = old(:, 4)
        f%dBmod = old(:, 5)
        f%d2Ath = hessian(:, 1)
        f%d2Aph = hessian(:, 2)
        f%d2hth = hessian(:, 3)
        f%d2hph = hessian(:, 4)
        f%d2Bmod = hessian(:, 5)
        call gc_geometry_gradient_transport(delta, f, actual)
        call compare('transported geometry gradients', reshape(actual, [15]), &
                     reshape(expected, [15]))
    end subroutine check_geometry_transport

    !> Independent derivative of a quadratic field with every mixed term.
    function quadratic_gradient(state, q) result(gradient)
        real(dp), intent(in) :: state(3)
        integer, intent(in) :: q
        real(dp) :: gradient(3)

        gradient = [q*(2*state(1) + state(2) + 2*state(3)) + 1, &
                    q*(4*state(2) + state(1) + 3*state(3)) + q + 1, &
                    q*(6*state(3) + 2*state(1) + 3*state(2)) + q + 2]
    end function quadratic_gradient

    !> Base field members of field_can_t, drawn independently.
    subroutine random_base_field(f)
        type(field_can_t), intent(out) :: f
        integer :: k

        f%Ath = uniform(-1d0, 1d0)
        f%Aph = uniform(-1d0, 1d0)
        f%hth = uniform(-1d0, 1d0)
        f%hph = uniform(0.5d0, 1.5d0)
        f%Bmod = uniform(0.5d0, 1.5d0)
        f%Phie = uniform(-1d0, 1d0)
        do k = 1, 3
            f%dAth(k) = uniform(-1d0, 1d0)
            f%dAph(k) = uniform(-1d0, 1d0)
            f%dhth(k) = uniform(-1d0, 1d0)
            f%dhph(k) = uniform(-1d0, 1d0)
            f%dBmod(k) = uniform(-1d0, 1d0)
            f%dPhie(k) = uniform(-1d0, 1d0)
        end do
        do k = 1, 6
            f%d2Ath(k) = uniform(-1d0, 1d0)
            f%d2Aph(k) = uniform(-1d0, 1d0)
            f%d2hth(k) = uniform(-1d0, 1d0)
            f%d2hph(k) = uniform(-1d0, 1d0)
            f%d2Bmod(k) = uniform(-1d0, 1d0)
            f%d2Phie(k) = uniform(-1d0, 1d0)
        end do
        f%mu = uniform(0d0, 1d0)
        f%ro0 = uniform(0.5d0, 2d0)
        f%dvpar = 7d0
        f%dH = 7d0
        f%dpth = 7d0
        f%d2vpar = 7d0
        f%d2H = 7d0
        f%d2pth = 7d0
    end subroutine random_base_field

    subroutine compare_fields(name, new, old, order)
        character(*), intent(in) :: name
        type(field_can_t), intent(in) :: new, old
        integer, intent(in) :: order

        call compare(name//' values', [new%vpar, new%H, new%pth], &
                     [old%vpar, old%H, old%pth])
        if (order < 1) return
        call compare(name//' gradients', [new%dvpar, new%dH, new%dpth], &
                     [old%dvpar, old%dH, old%dpth])
        if (order < 2) return
        call compare(name//' Hessians', [new%d2vpar, new%d2H, new%d2pth], &
                     [old%d2vpar, old%d2H, old%d2pth])
    end subroutine compare_fields

    subroutine check_radial_mixed()
        type(field_can_t) :: f, new, old
        real(dp) :: pphi

        call random_base_field(f)
        pphi = uniform(-1d0, 1d0)
        new = f
        old = f
        call gpu_get_derivatives2_radial_mixed(new, pphi)
        call old_radial_mixed(old, pphi)
        call compare_fields('gpu_get_derivatives2_radial_mixed', new, old, 2)
    end subroutine check_radial_mixed

    !> The hand-written OpenACC subset of SIMPLE d5111a87, which omits
    !> d2Phie: SIMPLE's electric potential is linear in s.
    subroutine old_radial_mixed(f, pphi)
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: pphi

        call get_derivatives(f, pphi)
        f%d2vpar(1:3) = -f%d2Aph(1:3)/f%ro0 - f%d2hph(1:3)*f%vpar
        f%d2vpar(1) = f%d2vpar(1) - 2.0_dp*f%dhph(1)*f%dvpar(1)
        f%d2vpar(2) = f%d2vpar(2) - (f%dhph(1)*f%dvpar(2) + f%dhph(2)*f%dvpar(1))
        f%d2vpar(3) = f%d2vpar(3) - (f%dhph(1)*f%dvpar(3) + f%dhph(3)*f%dvpar(1))
        f%d2vpar(1:3) = f%d2vpar(1:3)/f%hph
        f%d2H(1:3) = f%vpar*f%d2vpar(1:3) + f%mu*f%d2Bmod(1:3)
        f%d2H(1) = f%d2H(1) + f%dvpar(1)**2
        f%d2H(2) = f%d2H(2) + f%dvpar(1)*f%dvpar(2)
        f%d2H(3) = f%d2H(3) + f%dvpar(1)*f%dvpar(3)
        f%d2pth(1:3) = f%d2vpar(1:3)*f%hth + f%vpar*f%d2hth(1:3) + f%d2Ath(1:3)/f%ro0
        f%d2pth(1) = f%d2pth(1) + 2.0_dp*f%dvpar(1)*f%dhth(1)
        f%d2pth(2) = f%d2pth(2) + f%dvpar(1)*f%dhth(2) + f%dvpar(2)*f%dhth(1)
        f%d2pth(3) = f%d2pth(3) + f%dvpar(1)*f%dhth(3) + f%dvpar(3)*f%dhth(1)
        f%d2vpar(7:9) = -f%dhph/f%hph**2
        f%d2H(7:9) = f%dvpar(1:3)/f%hph + f%vpar*f%d2vpar(7:9)
        f%d2pth(7:9) = f%dhth/f%hph + f%hth*f%d2vpar(7:9)
    end subroutine old_radial_mixed

    !> The generated Boozer rates against the hand-written OpenACC rates of
    !> SIMPLE d5111a87, including p_phi with v_par close to zero.
    subroutine check_rates_boozer()
        real(dp) :: pph, mu, ro0, torflux, aphi, daphi, btheta, dbtheta, bphi
        real(dp) :: dbphi, bmod, dbmod(3), rates(4), ref(4)
        integer :: near_zero

        do near_zero = 0, 1
            mu = uniform(0d0, 1d0)
            ro0 = uniform(0.5d0, 2d0)
            torflux = uniform(4d0, 6d0)
            aphi = uniform(-1d0, 1d0)
            daphi = uniform(-1d0, 1d0)
            btheta = uniform(-1d0, 1d0)
            dbtheta = uniform(-1d0, 1d0)
            bphi = uniform(0.5d0, 1.5d0)
            dbphi = uniform(-1d0, 1d0)
            bmod = uniform(0.5d0, 1.5d0)
            dbmod = [uniform(-1d0, 1d0), uniform(-1d0, 1d0), uniform(-1d0, 1d0)]
            pph = uniform(-1d0, 1d0)
            if (near_zero == 1) pph = aphi/ro0*(1d0 + 1d-9)
            call gc_rates_boozer(pph, mu, ro0, torflux, aphi, daphi, btheta, &
                                 dbtheta, bphi, dbphi, bmod, dbmod(1), dbmod(2), &
                                 dbmod(3), rates(1), rates(2), rates(3), rates(4))
            call old_rates_boozer(pph, mu, ro0, torflux, aphi, daphi, btheta, dbtheta, &
                                  bphi, dbphi, bmod, dbmod, ref)
            call compare('Boozer guiding-centre rates', rates, ref)
        end do
    end subroutine check_rates_boozer

    subroutine old_rates_boozer(pph, mu, ro0, torflux, aphi, daphi, btheta, dbtheta, &
                                bphi, dbphi, bmod, dbmod, dydt)
        real(dp), intent(in) :: pph, mu, ro0, torflux, aphi, daphi, btheta, dbtheta
        real(dp), intent(in) :: bphi, dbphi, bmod, dbmod(3)
        real(dp), intent(out) :: dydt(4)
        real(dp) :: hth, hph, vpar, hprime, inv_bmod, inv_bmod2, inv_hph, inv_dpth
        real(dp) :: dhth_r, dhth_phi, dhph_r, dhph_theta, dhph_phi
        real(dp) :: dvpar_r, dvpar_theta, dvpar_phi, dh_r, dh_theta, dh_phi
        real(dp) :: dpth_r, dpth_phi

        inv_bmod = 1.0_dp/bmod
        inv_bmod2 = inv_bmod*inv_bmod
        hth = btheta*inv_bmod
        hph = bphi*inv_bmod
        dhth_r = dbtheta*inv_bmod - btheta*dbmod(1)*inv_bmod2
        dhth_phi = -btheta*dbmod(3)*inv_bmod2
        dhph_r = dbphi*inv_bmod - bphi*dbmod(1)*inv_bmod2
        dhph_theta = -bphi*dbmod(2)*inv_bmod2
        dhph_phi = -bphi*dbmod(3)*inv_bmod2
        inv_hph = 1.0_dp/hph
        vpar = (pph - aphi/ro0)*inv_hph
        dvpar_r = -(daphi/ro0 + dhph_r*vpar)*inv_hph
        dvpar_theta = -dhph_theta*vpar*inv_hph
        dvpar_phi = -dhph_phi*vpar*inv_hph
        dh_r = vpar*dvpar_r + mu*dbmod(1)
        dh_theta = vpar*dvpar_theta + mu*dbmod(2)
        dh_phi = vpar*dvpar_phi + mu*dbmod(3)
        dpth_r = dvpar_r*hth + vpar*dhth_r + torflux/ro0
        dpth_phi = dvpar_phi*hth + vpar*dhth_phi
        inv_dpth = 1.0_dp/dpth_r
        hprime = dh_r*inv_dpth
        dydt(1) = -(dh_theta - hth*inv_hph*dh_phi)*inv_dpth
        dydt(2) = hprime
        dydt(3) = (vpar - hprime*hth)*inv_hph
        dydt(4) = -(dh_phi - hprime*dpth_phi)
    end subroutine old_rates_boozer

end program test_orbit_systems_generated
