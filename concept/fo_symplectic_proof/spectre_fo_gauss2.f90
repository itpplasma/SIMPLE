module spectre_fo_gauss2
    !> Canonical FO in a smooth SPECTRE layer. q is stacked rho, theta, zeta;
    !> P = e_i.v + A_i/ro0 in cm, H = (P-A/ro0)^T g^{-1}(P-A/ro0)/2.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use magfie_sub, only: spectre_field, TESLA_TO_GAUSS, M_TO_CM
    use lapack_interfaces, only: dgesv
    use fortnum_linalg, only: inv3
    use spectre_basis, only: spectre_vecpot_t, &
        get_cheby_d2, get_zernike_d2
    use fo_kernel, only: fo_monitor_factor, fo_jet_product, &
                         fo_radial, fo_metric_monitor, fo_radial_first, &
                         fo_metric_monitor_first
    use fo_joint_kernel, only: fo_joint_first, fo_joint_full, fo_joint_flow
    implicit none
    private
    real(dp), parameter :: A_SCALE = TESLA_TO_GAUSS*M_TO_CM**2
    real(dp), allocatable :: gauge(:, :, :, :)
    type :: potential_terms_t
        integer, allocatable :: offset(:), degree(:)
        real(dp), allocatable :: coeff(:, :, :)
        real(dp), allocatable :: polynomial(:, :, :, :)
    end type potential_terms_t
    type(potential_terms_t), allocatable :: potential_terms(:)
    logical, allocatable :: gauge_active(:)
    real(dp), allocatable :: angular_factor(:, :, :)
    integer :: harmonic_m = 0, harmonic_n = 0, harmonic_stride = 1
    integer, allocatable :: harmonic_order(:, :)
    real(dp), parameter :: GL_A(2, 2) = reshape([0.25_dp, &
        0.25_dp + sqrt(3.0_dp)/6.0_dp, 0.25_dp - sqrt(3.0_dp)/6.0_dp, &
        0.25_dp], [2, 2])
    public :: gauss2_prepare, gauss2_geometry, gauss2_potential, gauss2_hamiltonian
    public :: gauss2_entry_step, gauss2_eval, gauss2_step, gauss2_phases
    public :: gauss2_encode, gauss2_decode, gauss2_gauge_zeta
contains
    subroutine gauss2_prepare
        !> Match the interface one-forms mode by mode. Their nonconstant part
        !> must be exact: m Delta A_zeta + n Delta A_theta = 0. Constants are
        !> gradients of linear functions on the unwrapped angular cover.
        real(dp), allocatable :: radial(:, :), cheby(:, :), zernike(:, :, :)
        real(dp) :: delta(2, 2), wave(2), c, defect, scale
        real(dp) :: input(2, 2, 0:8), poly(2, 2, 0:8)
        logical :: axis
        integer :: k, l, j, r, m, n, t, z, term, divisor, remainder

        if (.not. allocated(spectre_field)) &
            error stop "gauss2: initialize SPECTRE field before configuring"
        if (allocated(gauge)) deallocate (gauge)
        allocate (gauge(2, 2, spectre_field%data%mn, spectre_field%data%Mvol))
        gauge = 0.0_dp
        if (allocated(angular_factor)) deallocate (angular_factor)
        allocate (angular_factor(0:3, 0:3, spectre_field%data%mn))
        harmonic_m = maxval(abs(spectre_field%data%im))
        harmonic_stride = 0
        do j = 1, spectre_field%data%mn
            divisor = abs(spectre_field%data%in(j))
            do while (divisor /= 0)
                remainder = mod(harmonic_stride, divisor)
                harmonic_stride = divisor
                divisor = remainder
            end do
        end do
        harmonic_stride = max(1, harmonic_stride)
        harmonic_n = maxval(abs(spectre_field%data%in))/harmonic_stride
        if (allocated(harmonic_order)) deallocate (harmonic_order)
        allocate (harmonic_order(2, spectre_field%data%mn))
        harmonic_order(1, :) = spectre_field%data%im
        harmonic_order(2, :) = spectre_field%data%in/harmonic_stride
        angular_factor = 0.0_dp
        do j = 1, spectre_field%data%mn
            do z = 0, 3
                do t = 0, 3 - z
                    angular_factor(t, z, j) = real(spectre_field%data%im(j), dp)**t* &
                                              real(-spectre_field%data%in(j), dp)**z
                end do
            end do
        end do
        do k = 1, spectre_field%data%Mvol - 1
            l = spectre_field%data%Lrad(k)
            allocate (radial(0:l, spectre_field%data%mn))
            if (k == 1) then
                allocate (zernike(0:l, 0:spectre_field%data%Mpol, 0:2))
                call get_zernike_d2(1.0_dp, l, spectre_field%data%Mpol, zernike)
                do j = 1, spectre_field%data%mn
                    radial(:, j) = zernike(:, spectre_field%data%im(j), 0)
                end do
                deallocate (zernike)
            else
                allocate (cheby(0:l, 0:2))
                call get_cheby_d2(1.0_dp, l, cheby)
                do j = 1, spectre_field%data%mn
                    radial(:, j) = cheby(:, 0)
                end do
                deallocate (cheby)
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
            deallocate (radial)
        end do
        if (allocated(gauge_active)) deallocate (gauge_active)
        allocate (gauge_active(spectre_field%data%Mvol))
        if (allocated(potential_terms)) deallocate (potential_terms)
        allocate (potential_terms(spectre_field%data%Mvol))
        do k = 1, spectre_field%data%Mvol
            gauge_active(k) = any(abs(gauge(:, :, :, k)) > 0.0_dp)
            l = spectre_field%data%Lrad(k)
            allocate (potential_terms(k)%offset(spectre_field%data%mn + 1))
            allocate (potential_terms(k)%degree((l + 1)*spectre_field%data%mn))
            allocate (potential_terms(k)%coeff(2, 2, (l + 1)*spectre_field%data%mn))
            term = 1
            do j = 1, spectre_field%data%mn
                potential_terms(k)%offset(j) = term
                do r = 0, l
                    if (max(abs(spectre_field%data%vol(k)%Ate(r, j)), &
                            abs(spectre_field%data%vol(k)%Ato(r, j)), &
                            abs(spectre_field%data%vol(k)%Aze(r, j)), &
                            abs(spectre_field%data%vol(k)%Azo(r, j))) <= 0.0_dp) cycle
                    potential_terms(k)%degree(term) = r
                    potential_terms(k)%coeff(:, 1, term) = &
                        [spectre_field%data%vol(k)%Ate(r, j), &
                         spectre_field%data%vol(k)%Aze(r, j)]
                    potential_terms(k)%coeff(:, 2, term) = &
                        [spectre_field%data%vol(k)%Ato(r, j), &
                         spectre_field%data%vol(k)%Azo(r, j)]
                    term = term + 1
                end do
            end do
            potential_terms(k)%offset(spectre_field%data%mn + 1) = term
            ! Re-express supported radial degrees exactly; larger orders keep
            ! the recurrence provider. Public direct-phase queries retain it too.
            if (l <= 8 .and. spectre_field%data%Mpol <= 8) then
                allocate (potential_terms(k)%polynomial(2, 2, 0:8, &
                                                        spectre_field%data%mn))
                axis = k == 1
                do j = 1, spectre_field%data%mn
                    m = spectre_field%data%im(j)
                    input = 0.0_dp
                    input(1, 1, 0:l) = spectre_field%data%vol(k)%Ate(:, j)
                    input(2, 1, 0:l) = spectre_field%data%vol(k)%Aze(:, j)
                    input(1, 2, 0:l) = spectre_field%data%vol(k)%Ato(:, j)
                    input(2, 2, 0:l) = spectre_field%data%vol(k)%Azo(:, j)
                    poly = 0.0_dp
                    include 'fo_poly_coeff.inc'
                    potential_terms(k)%polynomial(:, :, :, j) = poly
                end do
            end if
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

    subroutine gauss2_potential(q, width, ro0, a, da, d2a, first_order, phase)
        !> first_order skips the second derivatives: a and da are bitwise
        !> unchanged, d2a is incomplete and must not be used.
        real(dp), intent(in) :: q(3), width, ro0
        real(dp), intent(out) :: a(3), da(3, 3), d2a(3, 3, 3)
        real(dp), intent(in), optional :: phase(2, spectre_field%data%mn)
        logical, intent(in), optional :: first_order
        real(dp) :: al(3), dl(3, 3), hl(3, 3, 3)
        real(dp) :: ar(3), dr(3, 3), hr(3, 3, 3), s(0:3)
        integer :: k, owner, i, j
        logical :: second

        second = .true.
        if (present(first_order)) second = .not. first_order
        owner = min(max(int(q(1)) + 1, 1), spectre_field%data%Mvol)
        k = nint(q(1))
        if (k > 0 .and. k < spectre_field%data%Mvol) then
            if (abs(q(1) - k) < width/2.0_dp) then
                call volume_potential(q, k, al, dl, hl, second, phase)
                call volume_potential(q, k + 1, ar, dr, hr, second, phase)
                call switch_der(q(1), k, width, s)
                a = al + s(0)*(ar - al)
                da = dl + s(0)*(dr - dl)
                da(:, 1) = da(:, 1) + s(1)*(ar - al)
                a = a*A_SCALE/ro0
                da = da*A_SCALE/ro0
                if (.not. second) then
                    d2a = 0.0_dp
                    return
                end if
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
                d2a = d2a*A_SCALE/ro0
                return
            end if
        end if
        call volume_potential(q, owner, a, da, d2a, second, phase)
        a = a*A_SCALE/ro0
        da = da*A_SCALE/ro0
        d2a = d2a*A_SCALE/ro0
    end subroutine gauss2_potential

    subroutine volume_potential(q, owner, a, da, d2a, second, phase)
        real(dp), intent(in) :: q(3)
        integer, intent(in) :: owner
        real(dp), intent(out) :: a(3), da(3, 3), d2a(3, 3, 3)
        real(dp), intent(in), optional :: phase(2, spectre_field%data%mn)
        logical, intent(in) :: second
        ! Keep the gauge loop below unchanged: gfortran vectorizes it with
        ! libmvec sin/cos, and a different loop body changes the last bits.
        type(spectre_vecpot_t) :: av
        real(dp) :: arg, c, s, wave(3), value(2), deriv(2), h(6, 2)
        integer :: i, j, ii
        integer, parameter :: slot(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

        call prepared_potential(owner, &
                      2.0_dp*(q(1) - owner + 1) - 1.0_dp, q(2), q(3), av, second, phase)
        a = [0.0_dp, av%Ath, av%Azt]
        da = 0.0_dp
        da(2, :) = av%dAth
        da(3, :) = av%dAzt
        da(:, 1) = 2.0_dp*da(:, 1)
        h(:, 1) = av%d2Ath
        h(:, 2) = av%d2Azt
        d2a = 0.0_dp
        do j = 1, 3
            if (.not. second) exit
            do i = 1, 3
                d2a(2:3, i, j) = h(slot(i, j), :)* &
                    real(2**(merge(1, 0, i == 1) + merge(1, 0, j == 1)), dp)
            end do
        end do
        if (.not. gauge_active(owner)) return
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

    subroutine prepared_potential(owner, s, theta, zeta, av, second, phase)
        ! The radial basis and mode/degree summation order match libneo.
        ! Only exactly zero coefficient tuples are omitted; no spectral cutoff.
        integer, intent(in) :: owner
        real(dp), intent(in) :: s, theta, zeta
        type(spectre_vecpot_t), intent(out) :: av
        real(dp), intent(in), optional :: phase(2, spectre_field%data%mn)
        logical, intent(in) :: second
        real(dp) :: cheby(0:spectre_field%data%Lrad(owner), 0:2)
        real(dp) :: zernike(0:spectre_field%data%Lrad(owner), &
                            0:spectre_field%data%Mpol, 0:2)
        real(dp) :: tt(0:2), arg, carg, sarg, rmi, rni
        real(dp) :: values(2, 10), ce(2), co(2), pv(2), qv(2)
        real(dp) :: pe(2), po(2), increment(2)
        integer :: ii, term, ll, mi

        if (present(phase)) then
            if (allocated(potential_terms(owner)%polynomial)) then
                call polynomial_potential(owner, s, av, second, phase)
                return
            end if
        end if
        if (owner == 1) then
            call get_zernike_d2(0.5_dp*(s + 1.0_dp), &
               spectre_field%data%Lrad(owner), spectre_field%data%Mpol, zernike, second)
        else
            call get_cheby_d2(s, spectre_field%data%Lrad(owner), cheby, second)
        end if
        values = 0.0_dp
        do ii = 1, spectre_field%data%mn
            mi = spectre_field%data%im(ii)
            rmi = real(mi, dp)
            rni = real(spectre_field%data%in(ii), dp)
            if (present(phase)) then
                carg = phase(1, ii)
                sarg = phase(2, ii)
            else
                arg = rmi*theta - rni*zeta
                carg = cos(arg)
                sarg = sin(arg)
            end if
            do term = potential_terms(owner)%offset(ii), &
                potential_terms(owner)%offset(ii + 1) - 1
                ll = potential_terms(owner)%degree(term)
                if (owner == 1) then
                    tt = [zernike(ll, mi, 0), 0.5_dp*zernike(ll, mi, 1), &
                          0.25_dp*zernike(ll, mi, 2)]
                else
                    tt = cheby(ll, :)
                end if
                ce = potential_terms(owner)%coeff(:, 1, term)
                co = potential_terms(owner)%coeff(:, 2, term)
                if (second) then
                    include 'fo_potential_2.inc'
                else
                    include 'fo_potential_1.inc'
                end if
            end do
        end do
        av%Ath = values(1, 1)
        av%Azt = values(2, 1)
        av%dAth = values(1, 2:4)
        av%dAzt = values(2, 2:4)
        av%d2Ath = values(1, 5:10)
        av%d2Azt = values(2, 5:10)
    end subroutine prepared_potential

    subroutine polynomial_potential(owner, s, av, second, phase)
        integer, intent(in) :: owner
        real(dp), intent(in) :: s, phase(2, spectre_field%data%mn)
        type(spectre_vecpot_t), intent(out) :: av
        logical, intent(in) :: second
        real(dp) :: poly(2, 2, 0:8), radial(2, 2, 0:2), values(2, 10)
        real(dp) :: x, radial_scale, rmi, rni, cp, sp
        integer :: ii, m

        include 'fo_poly_point.inc'
        values = 0.0_dp
        do ii = 1, spectre_field%data%mn
            poly = potential_terms(owner)%polynomial(:, :, :, ii)
            cp = phase(1, ii)
            sp = phase(2, ii)
            rmi = real(spectre_field%data%im(ii), dp)
            rni = real(spectre_field%data%in(ii), dp)
            m = spectre_field%data%im(ii)
            if (owner == 1) then
                if (second) then
                    include 'fo_poly_axis_2.inc'
                else
                    include 'fo_poly_axis_1.inc'
                end if
            else
                if (second) then
                    include 'fo_poly_horner_2.inc'
                else
                    include 'fo_poly_horner_1.inc'
                end if
            end if
            if (second) then
                include 'fo_poly_angular_2.inc'
            else
                include 'fo_poly_angular_1.inc'
            end if
        end do
        av%Ath = values(1, 1)
        av%Azt = values(2, 1)
        av%dAth = values(1, 2:4)
        av%dAzt = values(2, 2:4)
        av%d2Ath = values(1, 5:10)
        av%d2Azt = values(2, 5:10)
    end subroutine polynomial_potential

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

    subroutine radial_profile(rho, owner, m, f)
        !> Radial interpolation profile of one Fourier mode and its derivatives.
        real(dp), intent(in) :: rho
        integer, intent(in) :: owner, m
        real(dp), intent(out) :: f(0:3)
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
    end subroutine radial_profile

    subroutine radial_geometry(f, inner, outer, coeff)
        real(dp), intent(in) :: f(0:3), inner(2), outer(2)
        real(dp), intent(out) :: coeff(2, 0:3)
        integer :: r

        do r = 0, 3
            coeff(:, r) = (outer - inner)*f(r)
        end do
        coeff(:, 0) = coeff(:, 0) + inner
    end subroutine radial_geometry

    subroutine gauss2_phases(q, phase)
        ! Harmonics are exact complex products; no truncation or interpolation.
        ! Addition/subtraction and emitted readback identities are proved by FortSym.
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: phase(2, spectre_field%data%mn)
        real(dp) :: poloidal(2, 0:harmonic_m), toroidal(2, 0:harmonic_n)
        real(dp) :: left(2), right(2), next(2), arg
        integer :: j, m, n

        if (harmonic_m + harmonic_n > spectre_field%data%mn) then
            ! Sparse high harmonics do not repay constructing the whole basis.
            do j = 1, spectre_field%data%mn
                arg = real(spectre_field%data%im(j), dp)*q(2) &
                      + real(-spectre_field%data%in(j), dp)*q(3)
                phase(:, j) = [cos(arg), sin(arg)]
            end do
            return
        end if
        include 'fo_phase_seed.inc'
        do j = 2, harmonic_m
            left = poloidal(:, j - 1)
            right = poloidal(:, 1)
            include 'fo_phase_sum.inc'
            poloidal(:, j) = next
        end do
        do j = 2, harmonic_n
            left = toroidal(:, j - 1)
            right = toroidal(:, 1)
            include 'fo_phase_sum.inc'
            toroidal(:, j) = next
        end do
        do j = 1, spectre_field%data%mn
            m = harmonic_order(1, j)
            n = harmonic_order(2, j)
            left = poloidal(:, abs(m))
            right = toroidal(:, abs(n))
            include 'fo_phase_orders.inc'
            include 'fo_phase_difference.inc'
            phase(:, j) = next
        end do
    end subroutine gauss2_phases

    subroutine gauss2_geometry(q, width, x, e, metric, dg, d2g, &
                               first_order, phase, basis_phase, embedding)
        !> Smooth the embedding, not just A: the native stacked chart is C0.
        !> SPEC interface Fourier coefficients and axis power law are unchanged
        !> outside the layer. Analytic geometry derivatives extend through order 3.
        !> first_order stops the embedding jets at order 2 and returns d2g = 0;
        !> x, e, metric and dg are bitwise unchanged.
        !> Internal embedding requests return only R/Z jets; the Cartesian
        !> and dense metric outputs are then undefined.
        real(dp), intent(in) :: q(3), width
        real(dp), intent(out) :: x(3), e(3, 3), metric(3, 3)
        real(dp), intent(out) :: dg(3, 3, 3), d2g(3, 3, 3, 3)
        real(dp), intent(out), optional :: phase(2, spectre_field%data%mn)
        real(dp), intent(in), optional :: basis_phase(2, spectre_field%data%mn)
        real(dp), intent(out), optional :: embedding(0:3, 0:3, 0:3, 2)
        logical, intent(in), optional :: first_order
        real(dp) :: jets(0:3, 0:3, 0:3, 2), coeff(2, 0:3), cl(2, 0:3)
        real(dp) :: cr(2, 0:3), sw(0:3), arg, c(0:3), jet_c, jet_s
        real(dp) :: fl(0:3), fr(0:3)
        real(dp) :: profiles(0:3, 0:spectre_field%data%Mpol, 2)
        real(dp) :: d(2, 3), d2(2, 3, 3), d3(2, 3, 3, 3), rz(2), cz, sz
        integer :: owner, iface, ii, m, r, l, comp, v, i, j, k
        integer :: idx(3), nord
        integer, parameter :: binom(0:3, 0:3) = reshape([ &
            1, 0, 0, 0, 1, 1, 0, 0, 1, 2, 1, 0, 1, 3, 3, 1], [4, 4])
        logical :: blend

        nord = 3
        if (present(first_order)) then
            if (first_order) nord = 2
        end if
        owner = min(max(int(q(1)) + 1, 1), spectre_field%data%Mvol)
        iface = nint(q(1))
        blend = iface > 0 .and. iface < spectre_field%data%Mvol
        if (blend) blend = abs(q(1) - iface) < width/2.0_dp
        if (blend) call switch_der(q(1), iface, width, sw)
        v = owner
        if (blend) v = iface
        do m = 0, spectre_field%data%Mpol
            call radial_profile(q(1), v, m, profiles(:, m, 1))
            if (blend) call radial_profile(q(1), iface + 1, m, profiles(:, m, 2))
        end do
        jets = 0.0_dp
        do ii = 1, spectre_field%data%mn
            m = spectre_field%data%im(ii)
            if (present(basis_phase)) then
                c = [basis_phase(1, ii), -basis_phase(2, ii), &
                     -basis_phase(1, ii), basis_phase(2, ii)]
            else
                arg = real(m, dp)*q(2) + real(-spectre_field%data%in(ii), dp)*q(3)
            c = [cos(arg), -sin(arg), -cos(arg), sin(arg)]
            end if
            if (present(phase)) phase(:, ii) = c(0:3:3)
            fl = profiles(:, m, 1)
            if (blend) fr = profiles(:, m, 2)
            do comp = 1, 2
                call geometry_coeff(v, comp, fl, cl)
                coeff = cl
                if (blend) then
                    call geometry_coeff(iface + 1, comp, fr, cr)
                    do r = 0, 3
                        do l = 0, r
                            coeff(:, r) = coeff(:, r) + real(binom(l, r), dp)* &
                                sw(l)*(cr(:, r - l) - cl(:, r - l))
                        end do
                    end do
                end if
                if (nord == 2) then
                    include 'fo_geometry_2.inc'
                else
                    include 'fo_geometry_3.inc'
                end if
            end do
        end do
        jets = jets*M_TO_CM
        if (present(embedding)) then
            embedding = jets
            return
        end if
        rz = jets(0, 0, 0, :)
        do i = 1, 3
            idx = 0
            idx(i) = idx(i) + 1
            d(:, i) = jets(idx(1), idx(2), idx(3), :)
            do j = 1, 3
                idx(j) = idx(j) + 1
                d2(:, i, j) = jets(idx(1), idx(2), idx(3), :)
                if (nord < 3) then
                    idx(j) = idx(j) - 1
                    cycle
                end if
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
        if (nord < 3) d2g = 0.0_dp
        do k = 1, 3
            do j = 1, 3
                do i = 1, 3
                    dg(i, j, k) = dot_product(d2(:, i, k), d(:, j)) &
                        + dot_product(d(:, i), d2(:, j, k))
                    if (nord < 3) cycle
                    do l = 1, 3
                        d2g(i, j, k, l) = dot_product(d3(:, i, k, l), d(:, j)) &
                            + dot_product(d2(:, i, k), d2(:, j, l)) &
                            + dot_product(d2(:, i, l), d2(:, j, k)) &
                            + dot_product(d(:, i), d3(:, j, k, l))
                    end do
                end do
            end do
            dg(3, 3, k) = dg(3, 3, k) + 2.0_dp*rz(1)*d(1, k)
            if (nord < 3) cycle
            do l = 1, 3
                d2g(3, 3, k, l) = d2g(3, 3, k, l) + &
                    2.0_dp*(d(1, k)*d(1, l) + rz(1)*d2(1, k, l))
            end do
        end do
    contains
        subroutine geometry_coeff(vol, component, f, cf)
            integer, intent(in) :: vol, component
            real(dp), intent(in) :: f(0:3)
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
            call radial_geometry(f, inner, outer, cf)
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

    subroutine gauss2_hamiltonian(z, width, ro0, energy, grad, hess, &
            radial_metric, radial_grad, radial_hess)
        real(dp), intent(in) :: z(6), width, ro0
        real(dp), intent(out) :: energy, grad(6), hess(6, 6)
        real(dp), intent(out), optional :: radial_metric, radial_grad(3)
        real(dp), intent(out), optional :: radial_hess(3, 3)
        real(dp) :: f(6), df(6, 6), rj(10)

        call gauss2_eval(z, width, ro0, energy, f, df, &
            grad=grad, hess=hess, radial_jet=rj)
        if (present(radial_metric)) radial_metric = rj(1)
        if (present(radial_grad)) radial_grad = rj(2:4)
        if (present(radial_hess)) then
            radial_hess = reshape([rj(5), rj(6), rj(8), rj(6), rj(7), &
                rj(9), rj(8), rj(9), rj(10)], [3, 3])
        end if
    end subroutine gauss2_hamiltonian

    subroutine gauss2_eval(z, width, ro0, energy, f, df, e0, length, gmin, &
                    radial_scale, kenergy, monitor, grad, hess, radial_jet, first_order)
        !> first_order skips the second-order field jets. energy, kenergy,
        !> monitor, f and grad depend only on value and first-order jets in the
        !> generated kernel. CSE can change their rounding; df, hess and the
        !> second-order radial_jet entries are then meaningless.
        real(dp), intent(in) :: z(6), width, ro0
        real(dp), intent(out) :: energy, f(6), df(6, 6)
        real(dp), intent(in), optional :: e0, length, gmin, radial_scale
        real(dp), intent(out), optional :: kenergy, monitor, grad(6), hess(6, 6)
        real(dp), intent(out), optional :: radial_jet(10)
        logical, intent(in), optional :: first_order
        real(dp) :: aj(3, 10), gj(6, 10), phase(2, spectre_field%data%mn)
        real(dp) :: jets(0:3, 0:3, 0:3, 2), raw_grad(6), raw_hess(6, 6), raw_energy
        real(dp) :: x(3), basis(3, 3), g(3, 3), dg(3, 3, 3), d2g(3, 3, 3, 3)
        logical :: first

        first = .false.
        if (present(first_order)) first = first_order
        call gauss2_phases(z(1:3), phase)
        call gauss2_geometry(z(1:3), width, x, basis, g, dg, d2g, &
                             first_order=first, basis_phase=phase, embedding=jets)
        call potential_jets(z(1:3), width, ro0, aj, first_order, phase)
        if (first) then
            call fo_joint_first(z(4:6), aj, jets, gj, raw_energy, raw_grad, raw_hess)
        else
            call fo_joint_full(z(4:6), aj, jets, gj, raw_energy, raw_grad, raw_hess)
        end if
        call jet_eval(z, gj, energy, f, df, e0, length, gmin, &
                      radial_scale, kenergy, monitor, grad, hess, radial_jet, &
                      first_order, raw_energy, raw_grad, raw_hess)
    end subroutine gauss2_eval

    subroutine potential_jets(q, width, ro0, aj, first_order, phase)
        !> Potential jets: value, rho, theta, zeta, rr, rt, tt, rz, tz, zz.
        real(dp), intent(in) :: q(3), width, ro0
        real(dp), intent(out) :: aj(3, 10)
        real(dp), intent(in), optional :: phase(2, spectre_field%data%mn)
        logical, intent(in), optional :: first_order
        real(dp) :: av(3), da(3, 3), d2a(3, 3, 3)
        integer :: i, j, k

        call gauss2_potential(q, width, ro0, av, da, d2a, first_order, phase)
        aj(:, 1) = av
        aj(:, 2:4) = da
        k = 4
        do j = 1, 3
            do i = 1, j
                k = k + 1
                aj(:, k) = d2a(:, i, j)
            end do
        end do
    end subroutine potential_jets

    subroutine jet_eval(z, gj, energy, f, df, e0, length, gmin, &
                        radial_scale, kenergy, monitor, grad, hess, radial_jet, &
                        first_order, raw_energy, raw_grad, raw_hess)
        real(dp), intent(in) :: raw_energy, raw_grad(6), raw_hess(6, 6)
        real(dp), intent(in) :: z(6), gj(6, 10)
        real(dp), intent(out) :: energy, f(6), df(6, 6)
        real(dp), intent(in), optional :: e0, length, gmin, radial_scale
        real(dp), intent(out), optional :: kenergy, monitor, grad(6), hess(6, 6)
        real(dp), intent(out), optional :: radial_jet(10)
        logical, intent(in), optional :: first_order
        logical :: first
        real(dp) :: ke, gh(6), hh(6, 6), tj(10), rhoj(3), rj(10)
        real(dp) :: factor(3), product(3), energy0
        integer :: k

        ! Normalized P=e_i.v+A/ro0 has unit mass and charge.
        energy0 = 0.0_dp
        if (present(e0)) energy0 = e0
        rhoj = [1.0_dp, 0.0_dp, 0.0_dp]
        if (present(gmin) .and. present(length)) then
            do k = 1, spectre_field%data%Mvol - 1
                call fo_monitor_factor(z(1) - real(k, dp), length, gmin, &
                    factor(1), factor(2), factor(3))
                call fo_jet_product(rhoj(1), rhoj(2), rhoj(3), &
                    factor(1), factor(2), factor(3), &
                    product(1), product(2), product(3))
                rhoj = product
            end do
        end if
        tj = 0.0_dp
        tj(1) = rhoj(1)
        tj(2) = rhoj(2)
        tj(5) = rhoj(3)
        first = .false.
        if (present(first_order)) first = first_order
        if (present(radial_jet)) then
            call fo_radial(gj, rj)
            radial_jet = rj
        end if
        if (present(radial_scale)) then
            if (radial_scale > 0.0_dp) then
                if (first) then
                    call fo_radial_first(gj, rj)
                    call fo_metric_monitor_first(rj, rhoj, radial_scale, tj)
                else
                call fo_radial(gj, rj)
                call fo_metric_monitor(rj, rhoj, radial_scale, tj)
            end if
        end if
        end if
        energy = raw_energy
        gh = raw_grad
        hh = raw_hess
        call fo_joint_flow(energy, energy0, tj, gh, hh, ke, f, df)
        if (present(kenergy)) kenergy = ke
        if (present(monitor)) monitor = tj(1)
        if (present(grad)) grad = gh
        if (present(hess)) hess = hh
    end subroutine jet_eval

    pure real(dp) function gauss2_entry_step(h_entry, width, length, gmin) result(h)
        real(dp), intent(in) :: h_entry, width, length, gmin

        ! At a sheet the rho factor is at most gmin. With the frozen metric
        ! reference, the on-shell radial speed is bounded by length/(16*h_entry)
        ! times this factor. Resolve the compact width as well as the monitor
        ! length. This entry choice stays fixed for the complete FO segment.
        h = 0.5_dp*h_entry*min(1.0_dp, width/(gmin*length))
    end function gauss2_entry_step

    subroutine gauss2_step(z0, width, ro0, e0, length, gmin, h, z1, dt, ok, &
                           tangent, radial_scale, tolerance, residual, start_energy)
        !> start_energy returns H(z0), also when the Newton solve fails.
        real(dp), intent(in) :: z0(6), width, ro0, e0, length, gmin, h
        real(dp), intent(out) :: z1(6), dt
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: tangent(6, 6), residual, start_energy
        real(dp), intent(in), optional :: radial_scale, tolerance
        real(dp) :: dy(6, 2), f(6, 2), df(6, 6, 2), monitor(2)
        real(dp) :: jac(12, 12), rhs(12, 6), r(6, 2), delta(12)
        real(dp) :: eye(6, 6), energy, radial_ref, tol, res
        integer :: i, it, piv(12), info
        logical :: full

        ok = .false.
        z1 = z0
        dt = 0.0_dp
        if (present(tangent)) tangent = 0.0_dp
        if (present(residual)) residual = huge(1.0_dp)
        radial_ref = 0.0_dp
        if (present(radial_scale)) radial_ref = radial_scale
        tol = 1.0e-12_dp
        if (present(tolerance)) tol = tolerance
        eye = 0.0_dp
        do i = 1, 6
            eye(i, i) = 1.0_dp
        end do
        ! Only f(z0) and H(z0) are used here, so the Hessian is not formed.
        ! Splitting this evaluation over threads costs more than it saves.
        call gauss2_eval(z0, width, ro0, energy, f(:, 1), df(:, :, 1), &
                         e0, length, gmin, radial_ref, first_order=.true.)
        if (present(start_energy)) start_energy = energy
        ! Solve for stage increments: subtracting two large absolute momenta
        ! needlessly adds a cancellation floor to the Newton residual.
        dy(:, 1) = h*(0.5_dp - sqrt(3.0_dp)/6.0_dp)*f(:, 1)
        dy(:, 2) = h*(0.5_dp + sqrt(3.0_dp)/6.0_dp)*f(:, 1)
        do it = 1, 20
            ! The Jacobian is needed only if this iterate is not accepted. A
            ! small-step predictor may already be, so evaluate to first order and
            ! add the Hessian at the same point only when the solve continues.
            ! The first/full oracle bounds CSE rounding in the retained values.
            full = present(tangent)
            call stage_evals(full)
            do i = 1, 2
                r(:, i) = dy(:, i) - h*(GL_A(i, 1)*f(:, 1) + GL_A(i, 2)*f(:, 2))
            end do
            res = maxval(abs(r))
            if (present(residual)) residual = res
            if (.not. all(ieee_is_finite(r))) return
            if (res <= tol) then
                ok = .true.
                exit
            end if
            if (.not. full) call stage_evals(.true.)
            delta = -reshape(r, [12])
            call dgesv(12, 1, jac, 12, piv, delta, 12, info)
            if (info /= 0) return
            dy = dy + reshape(delta, [6, 2])
        end do
        if (.not. ok) return
        z1 = z0 + 0.5_dp*h*(f(:, 1) + f(:, 2))
        dt = 0.5_dp*h*sum(monitor)
        if (.not. present(tangent)) return
        ! Derivative of the exact stage equations at the computed stages.
        ! The finite-difference oracle also checks the terminated algorithm.
        rhs(1:6, :) = eye
        rhs(7:12, :) = eye
        call dgesv(12, 6, jac, 12, piv, rhs, 12, info)
        if (info /= 0) then
            ok = .false.
            return
        end if
        tangent = eye + 0.5_dp*h*(matmul(df(:, :, 1), rhs(1:6, :)) &
            + matmul(df(:, :, 2), rhs(7:12, :)))
    contains
        subroutine stage_evals(hessian)
            !> These small evaluations run serially, including when the caller
            !> already distributes markers over an OpenMP team.
            logical, intent(in) :: hessian
            real(dp) :: y(6, 2), stage_energy
            integer :: j, k
            logical :: first
            first = .not. hessian
            y(:, 1) = z0 + dy(:, 1)
            y(:, 2) = z0 + dy(:, 2)
            do j = 1, 2
                call gauss2_eval(y(:, j), width, ro0, stage_energy, &
                                 f(:, j), df(:, :, j), e0, length, gmin, radial_ref, &
                                 monitor=monitor(j), first_order=first)
                if (hessian) then
                    do k = 1, 2
                        jac(6*k - 5:6*k, 6*j - 5:6*j) = -h*GL_A(k, j)*df(:, :, j)
                    end do
                end if
            end do
            if (.not. hessian) return
            do k = 1, 12
                jac(k, k) = jac(k, k) + 1.0_dp
            end do
        end subroutine stage_evals
    end subroutine gauss2_step

end module spectre_fo_gauss2
