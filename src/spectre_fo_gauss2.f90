module spectre_fo_gauss2
    !> Canonical FO in a smooth SPECTRE layer. q is stacked rho, theta, zeta;
    !> P = e_i.v + A_i/ro0 in cm, H = (P-A/ro0)^T g^{-1}(P-A/ro0)/2.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use magfie_sub, only: spectre_field, TESLA_TO_GAUSS, M_TO_CM
    use lapack_interfaces, only: dgesv
    use fortnum_linalg, only: inv3
    use spectre_basis, only: spectre_vecpot_t, eval_spectre_vector_potential, &
        get_cheby_d2, get_zernike_d2
    implicit none
    private
    real(dp), parameter :: A_SCALE = TESLA_TO_GAUSS*M_TO_CM**2
    real(dp), allocatable :: gauge(:, :, :, :)
    real(dp), parameter :: GL_A(2, 2) = reshape([0.25_dp, &
        0.25_dp + sqrt(3.0_dp)/6.0_dp, 0.25_dp - sqrt(3.0_dp)/6.0_dp, &
        0.25_dp], [2, 2])
    public :: gauss2_prepare, gauss2_geometry, gauss2_potential, gauss2_hamiltonian
    public :: gauss2_step, gauss2_encode, gauss2_decode, gauss2_gauge_zeta
contains
    subroutine gauss2_prepare
        !> Match the interface one-forms mode by mode. Their nonconstant part
        !> must be exact: m Delta A_zeta + n Delta A_theta = 0. Constants are
        !> gradients of linear functions on the unwrapped angular cover.
        real(dp), allocatable :: radial(:, :), cheby(:, :), zernike(:, :, :)
        real(dp) :: delta(2, 2), wave(2), c, defect, scale
        integer :: k, l, j, r, m, n

        if (.not. allocated(spectre_field)) &
            error stop "gauss2: initialize SPECTRE field before configuring"
        if (allocated(gauge)) deallocate(gauge)
        allocate(gauge(2, 2, spectre_field%data%mn, spectre_field%data%Mvol))
        gauge = 0.0_dp
        do k = 1, spectre_field%data%Mvol - 1
            l = spectre_field%data%Lrad(k)
            allocate(radial(0:l, spectre_field%data%mn))
            if (k == 1) then
                allocate(zernike(0:l, 0:spectre_field%data%Mpol, 0:2))
                call get_zernike_d2(1.0_dp, l, spectre_field%data%Mpol, zernike)
                do j = 1, spectre_field%data%mn
                    radial(:, j) = zernike(:, spectre_field%data%im(j), 0)
                end do
                deallocate(zernike)
            else
                allocate(cheby(0:l, 0:2))
                call get_cheby_d2(1.0_dp, l, cheby)
                do j = 1, spectre_field%data%mn
                    radial(:, j) = cheby(:, 0)
                end do
                deallocate(cheby)
            end if
            scale = max(maxval(abs(spectre_field%data%vol(k)%Ate)), &
                maxval(abs(spectre_field%data%vol(k)%Aze)), tiny(1.0_dp))
            do j = 1, spectre_field%data%mn
                delta(1, 1) = dot_product(radial(:, j), &
                    spectre_field%data%vol(k)%Ate(:, j))
                delta(1, 2) = dot_product(radial(:, j), &
                    spectre_field%data%vol(k)%Ato(:, j))
                delta(2, 1) = dot_product(radial(:, j), &
                    spectre_field%data%vol(k)%Aze(:, j))
                delta(2, 2) = dot_product(radial(:, j), &
                    spectre_field%data%vol(k)%Azo(:, j))
                delta(1, 1) = delta(1, 1) - spectre_field%data%vol(k + 1)%Ate(0, j)
                delta(1, 2) = delta(1, 2) - spectre_field%data%vol(k + 1)%Ato(0, j)
                delta(2, 1) = delta(2, 1) - spectre_field%data%vol(k + 1)%Aze(0, j)
                delta(2, 2) = delta(2, 2) - spectre_field%data%vol(k + 1)%Azo(0, j)
                m = spectre_field%data%im(j)
                n = spectre_field%data%in(j)
                wave = [real(m, dp), -real(n, dp)]
                defect = maxval(abs(m*delta(2, :) + n*delta(1, :)))
                if (defect > 1.0e-11_dp*scale) then
                    print *, 'Non-exact SPECTRE interface gauge: ', k, j, defect
                    error stop 'gauss2: interface normal fields disagree'
                end if
                if (m /= 0 .or. n /= 0) then
                    do r = 1, 2
                        c = dot_product(wave, delta(:, r))/dot_product(wave, wave)
                        delta(:, r) = c*wave
                    end do
                end if
                ! The inner Chebyshev boundary has only its degree-zero mode.
                gauge(:, :, j, k + 1) = gauge(:, :, j, k) + delta
            end do
            deallocate(radial)
        end do
    end subroutine gauss2_prepare

    pure subroutine switch_der(rho, k, width, s)
        real(dp), intent(in) :: rho, width
        integer, intent(in) :: k
        real(dp), intent(out) :: s(0:3)
        real(dp) :: t, phi, small, sp, ph1, ph2, ph3

        t = (rho - real(k, dp))/width + 0.5_dp
        s = 0.0_dp
        if (t <= 0.0_dp) return
        s(0) = 1.0_dp
        if (t >= 1.0_dp) return
        ! Smooth compact switch: logistic(-1/t + 1/(1-t)). Every derivative
        ! vanishes at the ends, avoiding repeated kicks from finite smoothness.
        phi = -1.0_dp/t + 1.0_dp/(1.0_dp - t)
        if (abs(phi) > 700.0_dp) then
            s(0) = merge(1.0_dp, 0.0_dp, phi > 0.0_dp)
            return
        end if
        small = exp(-abs(phi))
        if (phi < 0.0_dp) then
            s(0) = small/(1.0_dp + small)
        else
            s(0) = 1.0_dp/(1.0_dp + small)
        end if
        sp = small/(1.0_dp + small)**2
        ph1 = 1.0_dp/t**2 + 1.0_dp/(1.0_dp - t)**2
        ph2 = -2.0_dp/t**3 + 2.0_dp/(1.0_dp - t)**3
        ph3 = 6.0_dp/t**4 + 6.0_dp/(1.0_dp - t)**4
        s(1) = sp*ph1/width
        s(2) = sp*((1.0_dp - 2.0_dp*s(0))*ph1**2 + ph2)/width**2
        s(3) = sp*((1.0_dp - 6.0_dp*s(0) + 6.0_dp*s(0)**2)*ph1**3 &
            + 3.0_dp*(1.0_dp - 2.0_dp*s(0))*ph1*ph2 + ph3)/width**3
    end subroutine switch_der

    subroutine gauss2_potential(q, width, ro0, a, da, d2a)
        real(dp), intent(in) :: q(3), width, ro0
        real(dp), intent(out) :: a(3), da(3, 3), d2a(3, 3, 3)
        real(dp) :: al(3), dl(3, 3), hl(3, 3, 3)
        real(dp) :: ar(3), dr(3, 3), hr(3, 3, 3), s(0:3)
        integer :: k, owner, i, j

        owner = min(max(int(q(1)) + 1, 1), spectre_field%data%Mvol)
        k = nint(q(1))
        if (k > 0 .and. k < spectre_field%data%Mvol) then
            if (abs(q(1) - k) < width/2.0_dp) then
                call volume_potential(q, k, al, dl, hl)
                call volume_potential(q, k + 1, ar, dr, hr)
                call switch_der(q(1), k, width, s)
                a = al + s(0)*(ar - al)
                da = dl + s(0)*(dr - dl)
                da(:, 1) = da(:, 1) + s(1)*(ar - al)
                d2a = hl + s(0)*(hr - hl)
                do j = 1, 3
                    do i = 1, 3
                        if (i == 1) d2a(:, i, j) = d2a(:, i, j) &
                            + s(1)*(dr(:, j) - dl(:, j))
                        if (j == 1) d2a(:, i, j) = d2a(:, i, j) &
                            + s(1)*(dr(:, i) - dl(:, i))
                    end do
                end do
                d2a(:, 1, 1) = d2a(:, 1, 1) + s(2)*(ar - al)
                a = a*A_SCALE/ro0
                da = da*A_SCALE/ro0
                d2a = d2a*A_SCALE/ro0
                return
            end if
        end if
        call volume_potential(q, owner, a, da, d2a)
        a = a*A_SCALE/ro0
        da = da*A_SCALE/ro0
        d2a = d2a*A_SCALE/ro0
    end subroutine gauss2_potential

    subroutine volume_potential(q, owner, a, da, d2a)
        real(dp), intent(in) :: q(3)
        integer, intent(in) :: owner
        real(dp), intent(out) :: a(3), da(3, 3), d2a(3, 3, 3)
        type(spectre_vecpot_t) :: av
        real(dp) :: arg, c, s, wave(3), value(2), deriv(2), h(6, 2)
        integer :: i, j, ii
        integer, parameter :: slot(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

        call eval_spectre_vector_potential(spectre_field%data, owner, &
            2.0_dp*(q(1) - owner + 1) - 1.0_dp, q(2), q(3), av)
        a = [0.0_dp, av%Ath, av%Azt]
        da = 0.0_dp
        da(2, :) = av%dAth
        da(3, :) = av%dAzt
        da(:, 1) = 2.0_dp*da(:, 1)
        h(:, 1) = av%d2Ath
        h(:, 2) = av%d2Azt
        d2a = 0.0_dp
        do j = 1, 3
            do i = 1, 3
                d2a(2:3, i, j) = h(slot(i, j), :)* &
                    real(2**(merge(1, 0, i == 1) + merge(1, 0, j == 1)), dp)
            end do
        end do
        do ii = 1, spectre_field%data%mn
            wave = [0.0_dp, real(spectre_field%data%im(ii), dp), &
                -real(spectre_field%data%in(ii), dp)]
            arg = wave(2)*q(2) + wave(3)*q(3)
            c = cos(arg)
            s = sin(arg)
            value = gauge(:, 1, ii, owner)*c + gauge(:, 2, ii, owner)*s
            deriv = -gauge(:, 1, ii, owner)*s + gauge(:, 2, ii, owner)*c
            a(2:3) = a(2:3) + value
            do j = 2, 3
                da(2:3, j) = da(2:3, j) + wave(j)*deriv
                do i = 2, 3
                    d2a(2:3, i, j) = d2a(2:3, i, j) - wave(i)*wave(j)*value
                end do
            end do
        end do
    end subroutine volume_potential

    real(dp) function gauss2_gauge_zeta(q, owner, ro0) result(shift)
        real(dp), intent(in) :: q(3), ro0
        integer, intent(in) :: owner
        real(dp) :: phase
        integer :: j
        !> Diagnostic expressed in the original volume gauge.
        shift = 0.0_dp
        do j = 1, spectre_field%data%mn
            phase = spectre_field%data%im(j)*q(2) - spectre_field%data%in(j)*q(3)
            shift = shift + gauge(2, 1, j, owner)*cos(phase) &
                + gauge(2, 2, j, owner)*sin(phase)
        end do
        shift = shift*A_SCALE/ro0
    end function gauss2_gauge_zeta

    subroutine radial_geometry(rho, owner, m, inner, outer, coeff)
        real(dp), intent(in) :: rho, inner(2), outer(2)
        integer, intent(in) :: owner, m
        real(dp), intent(out) :: coeff(2, 0:3)
        real(dp) :: f(0:3)
        integer :: p, r, l

        f = 0.0_dp
        if (owner == 1) then
            p = merge(2, m, m == 0)
            do r = 0, min(3, p)
                f(r) = rho**(p - r)
                do l = 0, r - 1
                    f(r) = f(r)*real(p - l, dp)
                end do
            end do
        else
            f(0) = rho - real(owner - 1, dp)
            f(1) = 1.0_dp
        end if
        do r = 0, 3
            coeff(:, r) = (outer - inner)*f(r)
        end do
        coeff(:, 0) = coeff(:, 0) + inner
    end subroutine radial_geometry

    subroutine gauss2_geometry(q, width, x, e, metric, dg, d2g)
        !> Smooth the embedding, not just A: the native stacked chart is C0.
        !> SPEC interface Fourier coefficients and axis power law are unchanged
        !> outside the layer. Analytic geometry derivatives extend through order 3.
        real(dp), intent(in) :: q(3), width
        real(dp), intent(out) :: x(3), e(3, 3), metric(3, 3)
        real(dp), intent(out) :: dg(3, 3, 3), d2g(3, 3, 3, 3)
        real(dp) :: jets(0:3, 0:3, 0:3, 2), coeff(2, 0:3), cl(2, 0:3)
        real(dp) :: cr(2, 0:3), sw(0:3), arg, c(0:3), factor, wave(2)
        real(dp) :: d(2, 3), d2(2, 3, 3), d3(2, 3, 3, 3), rz(2), cz, sz
        integer :: owner, iface, ii, m, n, r, t, z, l, comp, v, i, j, k
        integer :: idx(3)
        integer, parameter :: binom(0:3, 0:3) = reshape([ &
            1, 0, 0, 0, 1, 1, 0, 0, 1, 2, 1, 0, 1, 3, 3, 1], [4, 4])
        logical :: blend

        owner = min(max(int(q(1)) + 1, 1), spectre_field%data%Mvol)
        iface = nint(q(1))
        blend = iface > 0 .and. iface < spectre_field%data%Mvol
        if (blend) blend = abs(q(1) - iface) < width/2.0_dp
        if (blend) call switch_der(q(1), iface, width, sw)
        jets = 0.0_dp
        do ii = 1, spectre_field%data%mn
            m = spectre_field%data%im(ii)
            n = spectre_field%data%in(ii)
            wave = [real(m, dp), -real(n, dp)]
            arg = wave(1)*q(2) + wave(2)*q(3)
            c = [cos(arg), -sin(arg), -cos(arg), sin(arg)]
            do comp = 1, 2
                v = owner
                if (blend) v = iface
                call geometry_coeff(v, comp, cl)
                coeff = cl
                if (blend) then
                    call geometry_coeff(iface + 1, comp, cr)
                    do r = 0, 3
                        do l = 0, r
                            coeff(:, r) = coeff(:, r) + real(binom(l, r), dp)* &
                                sw(l)*(cr(:, r - l) - cl(:, r - l))
                        end do
                    end do
                end if
                do z = 0, 3
                    do t = 0, 3 - z
                        do r = 0, 3 - z - t
                            factor = wave(1)**t*wave(2)**z
                            jets(r, t, z, comp) = jets(r, t, z, comp) + factor* &
                                (coeff(1, r)*c(mod(t + z, 4)) &
                                - coeff(2, r)*c(mod(t + z + 1, 4)))
                        end do
                    end do
                end do
            end do
        end do
        jets = jets*M_TO_CM
        rz = jets(0, 0, 0, :)
        do i = 1, 3
            idx = 0
            idx(i) = idx(i) + 1
            d(:, i) = jets(idx(1), idx(2), idx(3), :)
            do j = 1, 3
                idx(j) = idx(j) + 1
                d2(:, i, j) = jets(idx(1), idx(2), idx(3), :)
                do k = 1, 3
                    idx(k) = idx(k) + 1
                    d3(:, i, j, k) = jets(idx(1), idx(2), idx(3), :)
                    idx(k) = idx(k) - 1
                end do
                idx(j) = idx(j) - 1
            end do
        end do
        cz = cos(q(3))
        sz = sin(q(3))
        x = [rz(1)*cz, rz(1)*sz, rz(2)]
        do i = 1, 3
            e(:, i) = [d(1, i)*cz, d(1, i)*sz, d(2, i)]
        end do
        e(:, 3) = e(:, 3) + [-rz(1)*sz, rz(1)*cz, 0.0_dp]
        metric = matmul(transpose(d), d)
        metric(3, 3) = metric(3, 3) + rz(1)**2
        do k = 1, 3
            do j = 1, 3
                do i = 1, 3
                    dg(i, j, k) = dot_product(d2(:, i, k), d(:, j)) &
                        + dot_product(d(:, i), d2(:, j, k))
                    do l = 1, 3
                        d2g(i, j, k, l) = dot_product(d3(:, i, k, l), d(:, j)) &
                            + dot_product(d2(:, i, k), d2(:, j, l)) &
                            + dot_product(d2(:, i, l), d2(:, j, k)) &
                            + dot_product(d(:, i), d3(:, j, k, l))
                    end do
                end do
            end do
            dg(3, 3, k) = dg(3, 3, k) + 2.0_dp*rz(1)*d(1, k)
            do l = 1, 3
                d2g(3, 3, k, l) = d2g(3, 3, k, l) + &
                    2.0_dp*(d(1, k)*d(1, l) + rz(1)*d2(1, k, l))
            end do
        end do
    contains
        subroutine geometry_coeff(vol, component, cf)
            integer, intent(in) :: vol, component
            real(dp), intent(out) :: cf(2, 0:3)
            real(dp) :: inner(2), outer(2)
            if (component == 1) then
                inner = [spectre_field%data%Rbc(ii, vol - 1), &
                    spectre_field%data%Rbs(ii, vol - 1)]
                outer = [spectre_field%data%Rbc(ii, vol), &
                    spectre_field%data%Rbs(ii, vol)]
            else
                inner = [spectre_field%data%Zbc(ii, vol - 1), &
                    spectre_field%data%Zbs(ii, vol - 1)]
                outer = [spectre_field%data%Zbc(ii, vol), &
                    spectre_field%data%Zbs(ii, vol)]
            end if
            call radial_geometry(q(1), vol, m, inner, outer, cf)
        end subroutine geometry_coeff
    end subroutine gauss2_geometry

    pure function inverse3(a) result(b)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: b(3, 3)
        integer :: info
        call inv3(a, b, info)
        if (info /= 0) b = ieee_value(0.0_dp, ieee_quiet_nan)
    end function inverse3

    subroutine gauss2_encode(x, v, width, ro0, z, ok)
        real(dp), intent(in) :: x(3), v(3), width, ro0
        real(dp), intent(inout) :: z(6)
        logical, intent(out) :: ok
        real(dp) :: xc(3), e(3, 3), g(3, 3), dg(3, 3, 3), d2g(3, 3, 3, 3)
        real(dp) :: a(3), da(3, 3), d2a(3, 3, 3), dq(3)
        integer :: it

        ok = .false.
        do it = 1, 20
            call gauss2_geometry(z(1:3), width, xc, e, g, dg, d2g)
            dq = matmul(inverse3(e), x - xc)
            z(1:3) = z(1:3) + dq
            if (maxval(abs(dq)) < 1.0e-13_dp) then
                ok = .true.
                exit
            end if
        end do
        call gauss2_potential(z(1:3), width, ro0, a, da, d2a)
        z(4:6) = matmul(transpose(e), v) + a
    end subroutine gauss2_encode

    subroutine gauss2_decode(z, width, ro0, x, v)
        real(dp), intent(in) :: z(6), width, ro0
        real(dp), intent(out) :: x(3), v(3)
        real(dp) :: e(3, 3), g(3, 3), dg(3, 3, 3), d2g(3, 3, 3, 3)
        real(dp) :: a(3), da(3, 3), d2a(3, 3, 3)
        call gauss2_geometry(z(1:3), width, x, e, g, dg, d2g)
        call gauss2_potential(z(1:3), width, ro0, a, da, d2a)
        v = matmul(transpose(inverse3(e)), z(4:6) - a)
    end subroutine gauss2_decode

    subroutine gauss2_hamiltonian(z, width, ro0, energy, grad, hess)
        real(dp), intent(in) :: z(6), width, ro0
        real(dp), intent(out) :: energy, grad(6), hess(6, 6)
        real(dp) :: x(3), e(3, 3), metric(3, 3), dg(3, 3, 3), d2g(3, 3, 3, 3)
        real(dp) :: a(3), da(3, 3), d2a(3, 3, 3), w(3), v(3)
        real(dp) :: gi(3, 3), dgi(3, 3, 3), d2gi(3, 3), tmp(3, 3)
        real(dp) :: basis_norm(3), unit_metric(3, 3)
        integer :: i, j

        call gauss2_geometry(z(1:3), width, x, e, metric, dg, d2g)
        call gauss2_potential(z(1:3), width, ro0, a, da, d2a)
        ! The guard on an unscaled metric mistakes coordinate aspect ratio
        ! for singular geometry. Invert the Gram matrix of unit basis vectors.
        do i = 1, 3
            basis_norm(i) = sqrt(metric(i, i))
        end do
        do j = 1, 3
            do i = 1, 3
                unit_metric(i, j) = metric(i, j)/(basis_norm(i)*basis_norm(j))
            end do
        end do
        gi = inverse3(unit_metric)
        do j = 1, 3
            do i = 1, 3
                gi(i, j) = gi(i, j)/(basis_norm(i)*basis_norm(j))
            end do
        end do
        w = z(4:6) - a
        v = matmul(gi, w)
        energy = 0.5_dp*dot_product(w, v)
        grad(4:6) = v
        hess(4:6, 4:6) = gi
        do i = 1, 3
            dgi(:, :, i) = -matmul(gi, matmul(dg(:, :, i), gi))
            grad(i) = 0.5_dp*dot_product(w, matmul(dgi(:, :, i), w)) &
                - dot_product(da(:, i), v)
            hess(i, 4:6) = matmul(dgi(:, :, i), w) - matmul(gi, da(:, i))
            hess(4:6, i) = hess(i, 4:6)
        end do
        do j = 1, 3
            do i = 1, 3
                tmp = matmul(dg(:, :, j), matmul(gi, dg(:, :, i))) &
                    + matmul(dg(:, :, i), matmul(gi, dg(:, :, j))) - d2g(:, :, i, j)
                d2gi = matmul(gi, matmul(tmp, gi))
                hess(i, j) = 0.5_dp*dot_product(w, matmul(d2gi, w)) &
                    - dot_product(da(:, i), matmul(dgi(:, :, j), w)) &
                    - dot_product(da(:, j), matmul(dgi(:, :, i), w)) &
                    + dot_product(da(:, i), matmul(gi, da(:, j))) &
                    - dot_product(d2a(:, i, j), v)
            end do
        end do
    end subroutine gauss2_hamiltonian

    subroutine transformed(z, width, ro0, e0, length, gmin, f, df, g)
        real(dp), intent(in) :: z(6), width, ro0, e0, length, gmin
        real(dp), intent(out) :: f(6), df(6, 6), g
        real(dp) :: energy, grad(6), hess(6, 6), gq, gqq, d, expd
        integer :: k

        call gauss2_hamiltonian(z, width, ro0, energy, grad, hess)
        g = 1.0_dp
        gq = 0.0_dp
        gqq = 0.0_dp
        do k = 1, spectre_field%data%Mvol - 1
            d = z(1) - real(k, dp)
            expd = (1.0_dp - gmin)*exp(-(d/length)**2)
            ! Product of positive refinements also covers neighbouring layers.
            gqq = gqq*(1.0_dp - expd) + 4.0_dp*gq*d/length**2*expd &
                + g*(2.0_dp/length**2 - 4.0_dp*d*d/length**4)*expd
            gq = gq*(1.0_dp - expd) + g*2.0_dp*d/length**2*expd
            g = g*(1.0_dp - expd)
        end do
        hess = g*hess
        hess(1, :) = hess(1, :) + gq*grad
        hess(:, 1) = hess(:, 1) + gq*grad
        hess(1, 1) = hess(1, 1) + (energy - e0)*gqq
        grad = g*grad
        grad(1) = grad(1) + (energy - e0)*gq
        f(1:3) = grad(4:6)
        f(4:6) = -grad(1:3)
        df(1:3, :) = hess(4:6, :)
        df(4:6, :) = -hess(1:3, :)
    end subroutine transformed

    subroutine gauss2_step(z0, width, ro0, e0, length, gmin, h, z1, dt, ok, tangent)
        real(dp), intent(in) :: z0(6), width, ro0, e0, length, gmin, h
        real(dp), intent(out) :: z1(6), dt
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: tangent(6, 6)
        real(dp) :: y(6, 2), f(6, 2), df(6, 6, 2), g(2), jac(12, 12)
        real(dp) :: rhs(12, 6), scale(6), resid(12), update(12), identity(6, 6)
        real(dp) :: roundoff(12), residual_floor(12)
        integer :: it, i, j, info, piv(12)
        logical :: stages_converged

        z1 = z0
        dt = 0.0_dp
        identity = 0.0_dp
        do i = 1, 6
            identity(i, i) = 1.0_dp
        end do
        scale = [1.0_dp, 1.0_dp, 1.0_dp, max(abs(z0(4:6)), 100.0_dp)]
        call transformed(z0, width, ro0, e0, length, gmin, f(:, 1), df(:, :, 1), g(1))
        y(:, 1) = z0 + h*(0.5_dp - sqrt(3.0_dp)/6.0_dp)*f(:, 1)
        y(:, 2) = z0 + h*(0.5_dp + sqrt(3.0_dp)/6.0_dp)*f(:, 1)
        ok = .false.
        stages_converged = .false.
        do it = 1, 20
            do j = 1, 2
                call transformed(y(:, j), width, ro0, e0, length, gmin, &
                    f(:, j), df(:, :, j), g(j))
                do i = 1, 2
                    jac(6*i - 5:6*i, 6*j - 5:6*j) = -h*GL_A(i, j)*df(:, :, j)
                end do
            end do
            do i = 1, 12
                jac(i, i) = jac(i, i) + 1.0_dp
            end do
            do i = 1, 2
                resid(6*i - 5:6*i) = y(:, i) - z0 &
                    - h*(GL_A(i, 1)*f(:, 1) + GL_A(i, 2)*f(:, 2))
            end do
            if (.not. all(ieee_is_finite(resid))) return
            ! A rounding-sized perturbation of the stages is amplified by
            ! the force Jacobian, especially for large unwrapped angles. Use
            ! that backward-error floor, not a fixed absolute force tolerance.
            roundoff = 8.0_dp*epsilon(1.0_dp)* &
                max(abs(reshape(y, [12])), 1.0_dp)
            residual_floor = 2.0e-13_dp*[scale, scale] &
                + matmul(abs(jac), roundoff)
            if (stages_converged .or. all(abs(resid) <= residual_floor)) then
                ok = .true.
                exit
            end if
            update = -resid
            call dgesv(12, 1, jac, 12, piv, update, 12, info)
            if (info /= 0) return
            y(:, 1) = y(:, 1) + update(1:6)
            y(:, 2) = y(:, 2) + update(7:12)
            ! Inside a narrow layer the force residual has a rounding floor
            ! amplified by its analytic Jacobian. Convergence of the stages
            ! themselves is the meaningful Newton criterion there.
            stages_converged = maxval(abs(update(1:6))/scale) < 3.0e-12_dp &
                .and. maxval(abs(update(7:12))/scale) < 3.0e-12_dp
        end do
        z1 = z0 + 0.5_dp*h*(f(:, 1) + f(:, 2))
        dt = 0.5_dp*h*sum(g)
        if (.not. ok) return
        if (.not. present(tangent)) return
        ! Analytic derivative of the implicit stages, with E and h held fixed.
        rhs(1:6, :) = identity
        rhs(7:12, :) = identity
        call dgesv(12, 6, jac, 12, piv, rhs, 12, info)
        if (info /= 0) then
            ok = .false.
            return
        end if
        tangent = identity + 0.5_dp*h*(matmul(df(:, :, 1), rhs(1:6, :)) &
            + matmul(df(:, :, 2), rhs(7:12, :)))
    end subroutine gauss2_step
end module spectre_fo_gauss2
