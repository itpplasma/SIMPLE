module gc_second_order
    ! Consumer physics: periodic sigma4 from BT 2016 (114), theta=-phi.
    ! The symbolic integrand has trigonometric degree <=6: G1 has degree
    ! <=2, D1(Pn) <=3, D1^2(Pn) <=5, and F1a <=4. Sixteen equispaced
    ! samples recover all harmonics exactly in exact arithmetic. This is
    ! finite Fourier interpolation, not a fitted or orbit-averaged correction.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use gc_generated, only: map1, map2base, gauge4
    implicit none
    private
    public :: map2, map1plus
contains
    subroutine map1plus(p, result)
        ! First-order position truncation with second-order velocity/phase
        ! closure. Keeping the latter also makes inverse initialization
        ! accurate enough that B(X)/epsilon does not amplify an O(eps^2)
        ! slow-coordinate bias into an O(eps) phase error on O(1) times.
        real(dp), intent(in) :: p(130)
        real(dp), intent(out) :: result(42)
        real(dp) :: first(42), jac(6, 6), first_jac(6, 6)
        call map2(p, result)
        call map1(p, first)
        result(1:3) = first(1:3)
        jac = reshape(result(7:42), [6, 6])
        first_jac = reshape(first(7:42), [6, 6])
        jac(1:3, :) = first_jac(1:3, :)
        result(7:42) = reshape(jac, [36])
    end subroutine map1plus

    subroutine map2(p, result)
        real(dp), intent(in) :: p(130)
        real(dp), intent(out) :: result(42)
        integer, parameter :: ns = 16, nh = 6
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: local(130), base(126), q(13), samples(13, ns)
        real(dp) :: primitive(12), ca(12), sa(12), phase, weight, meanq
        real(dp) :: sig(2), dsig(2, 6), basis(6), jac(6, 6), dbasis(6, 6)
        integer :: j, k, a, offset
        call map2base(p, base)
        local = p
        !$omp parallel do default(none) shared(samples) firstprivate(local) private(j)
        do j = 1, ns
            local(6) = 2*pi*(j - 1)/ns
            call gauge4(local, samples(:, j))
        end do
        !$omp end parallel do
        ! Solvability is required; removing a nonzero mean would conceal a
        ! representation error. A small roundoff residual is permitted.
        meanq = sum(samples(1, :))/ns
        if (abs(meanq) > 2.0e-9_dp*(1 + maxval(abs(samples(1, :))))) then
            result = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        primitive = 0.0_dp
        do k = 1, nh
            ca = 0.0_dp
            sa = 0.0_dp
            do j = 1, ns
                phase = 2*pi*(j - 1)*k/ns
                weight = 2.0_dp/ns
                ca = ca + weight*cos(phase)*samples(2:13, j)
                sa = sa + weight*sin(phase)*samples(2:13, j)
            end do
            ! partial_phi sigma4=-qtheta; choose <sigma4>=0.
            primitive = primitive + (sa*cos(k*p(6)) - ca*sin(k*p(6)))/k
        end do
        sig = primitive(1:2)
        dsig(1, 1:5) = primitive(3:7)
        dsig(2, 1:5) = primitive(8:12)
        call gauge4(p, q)
        dsig(:, 6) = -q(2:3)
        result(1:6) = base(1:6)
        jac = reshape(base(7:42), [6, 6])
        do a = 1, 2
            offset = 42*a
            basis = base(offset + 1:offset + 6)
            dbasis = reshape(base(offset + 7:offset + 42), [6, 6])
            result(1:6) = result(1:6) + basis*sig(a)
            do j = 1, 6
                jac(:, j) = jac(:, j) + dbasis(:, j)*sig(a) + basis*dsig(a, j)
            end do
        end do
        result(7:42) = reshape(jac, [36])
    end subroutine map2
end module gc_second_order
