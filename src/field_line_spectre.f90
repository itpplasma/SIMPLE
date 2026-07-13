module field_line_spectre
    !> Symplectic field-line tracer for a single SPECTRE/SPEC volume.
    !>
    !> The magnetic field line is the flow of a one-degree-of-freedom canonical
    !> system with the toroidal angle zeta as time, the poloidal angle theta as
    !> coordinate, and the covariant vector-potential component psi = A_theta as
    !> the conjugate momentum. The radial coordinate s is recovered implicitly
    !> from psi at fixed (theta, zeta). The derived flow is
    !>     dtheta/dzeta = -dA_zeta/ds / dA_theta/ds
    !>     dpsi/dzeta   =  dA_zeta/dtheta - dA_theta/dtheta * dA_zeta/ds/dA_theta/ds
    !> (the explicit d/dzeta terms cancel exactly). One step advances the pair
    !> (theta, psi) by a semi-implicit symplectic Euler map. Both update signs
    !> set the field-line orientation, so a sign flip reverses iota; keep them
    !> as written.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spectre_reader, only: spectre_data_t
    use spectre_basis, only: spectre_vecpot_t, eval_spectre_vector_potential

    implicit none
    private

    public :: field_line_psi, field_line_step, spectre_rz

    integer, parameter :: max_newton = 20
    real(dp), parameter :: newton_tol = 1.0e-13_dp

contains

    !> Canonical momentum psi = A_theta at (lvol, s, theta, zeta).
    pure function field_line_psi(data, lvol, s, theta, zeta) result(psi)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, theta, zeta
        real(dp) :: psi

        type(spectre_vecpot_t) :: av

        call eval_spectre_vector_potential(data, lvol, s, theta, zeta, av)
        psi = av%Ath
    end function field_line_psi

    !> Advance (s, theta, zeta, psi) by one symplectic Euler step of size h.
    !> The radial coordinate s is solved by Newton from
    !>     F(s1) = A_theta(s1,theta0,zeta0) - psi0 - h*dpsi/dzeta(s1,theta0,zeta0)
    !> with the analytic F' built from the ss and st second derivatives; theta
    !> and psi then advance explicitly. s carries the previous radius as the
    !> Newton seed and returns the solved radius.
    subroutine field_line_step(data, lvol, h, s, theta, zeta, psi)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: h
        real(dp), intent(inout) :: s, theta, zeta, psi

        type(spectre_vecpot_t) :: av
        real(dp) :: theta0, zeta0, psi0, s1, ratio, flow, fun, dfun, dgds
        integer :: iter
        logical :: converged

        theta0 = theta
        zeta0 = zeta
        psi0 = psi
        s1 = s
        ratio = 0.0_dp
        fun = 0.0_dp
        converged = .false.

        do iter = 1, max_newton
            call eval_spectre_vector_potential(data, lvol, s1, theta0, zeta0, av)
            ratio = av%dAzt(1)/av%dAth(1)
            flow = av%dAzt(2) - av%dAth(2)*ratio
            fun = av%Ath - psi0 - h*flow
            if (abs(fun) < newton_tol) then
                converged = .true.
                exit
            end if
            dgds = av%d2Azt(2) &
                   - (av%d2Ath(2)*av%dAzt(1) + av%dAth(2)*av%d2Azt(1))/av%dAth(1) &
                   + av%dAth(2)*av%dAzt(1)*av%d2Ath(1)/av%dAth(1)**2
            dfun = av%dAth(1) - h*dgds
            s1 = s1 - fun/dfun
        end do

        if (.not. converged) then
            write (*, '(A,I0,A,ES13.5,A,ES13.5)') &
                'field_line_step: Newton failed to converge in volume ', lvol, &
                ', s = ', s1, ', |F| = ', abs(fun)
            error stop 'field_line_step: Newton did not converge'
        end if

        theta = theta0 - h*ratio
        psi = av%Ath
        zeta = zeta0 + h
        s = s1
    end subroutine field_line_step

    !> Cylindrical R, Z of a point in volume lvol from the SPEC interface Fourier
    !> blend, mirroring libneo_coordinates_spectre for output only: interior
    !> volumes blend interfaces lvol-1 and lvol linearly in s; the axis volume
    !> (lvol == 1) uses the sbar**m power law (sbar**2 for m = 0) with
    !> sbar = (s + 1)/2.
    pure subroutine spectre_rz(data, lvol, s, theta, zeta, R, Z)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: R, Z

        integer :: ii, m, n
        real(dp) :: arg, carg, sarg, sbar, fj, alss, blss
        real(dp) :: rc, rs, zc, zs

        R = 0.0_dp
        Z = 0.0_dp
        do ii = 1, data%mn
            m = data%im(ii)
            n = data%in(ii)
            arg = real(m, dp)*theta - real(n, dp)*zeta
            carg = cos(arg)
            sarg = sin(arg)
            if (lvol == 1) then
                sbar = 0.5_dp*(s + 1.0_dp)
                if (m == 0) then
                    fj = sbar*sbar
                else
                    fj = sbar**m
                end if
                rc = data%Rbc(ii, 0) + (data%Rbc(ii, 1) - data%Rbc(ii, 0))*fj
                rs = data%Rbs(ii, 0) + (data%Rbs(ii, 1) - data%Rbs(ii, 0))*fj
                zc = data%Zbc(ii, 0) + (data%Zbc(ii, 1) - data%Zbc(ii, 0))*fj
                zs = data%Zbs(ii, 0) + (data%Zbs(ii, 1) - data%Zbs(ii, 0))*fj
            else
                alss = 0.5_dp*(1.0_dp - s)
                blss = 0.5_dp*(1.0_dp + s)
                rc = alss*data%Rbc(ii, lvol - 1) + blss*data%Rbc(ii, lvol)
                rs = alss*data%Rbs(ii, lvol - 1) + blss*data%Rbs(ii, lvol)
                zc = alss*data%Zbc(ii, lvol - 1) + blss*data%Zbc(ii, lvol)
                zs = alss*data%Zbs(ii, lvol - 1) + blss*data%Zbs(ii, lvol)
            end if
            R = R + rc*carg + rs*sarg
            Z = Z + zc*carg + zs*sarg
        end do
    end subroutine spectre_rz

end module field_line_spectre
