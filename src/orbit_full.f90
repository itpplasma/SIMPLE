module orbit_full
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use util, only: pi, twopi, c, e_charge, p_mass

    implicit none

    integer, parameter :: MAX_NEWTON_ITER = 20
    real(dp), parameter :: NEWTON_RTOL = 1d-12
    real(dp), parameter :: NEWTON_ATOL = 1d-15
    real(dp), parameter :: FD_STEP = 1d-8

    type :: FullOrbitState
        real(dp) :: z(6)           ! (R, phi, Z, p_R, p_phi, p_Z) in cylindrical
        real(dp) :: z_gc(5)        ! (s, theta, phi, lambda, v) guiding-center
        real(dp) :: mu             ! Magnetic moment (constant for Pauli, 0 for particle)
        real(dp) :: m              ! Particle mass [g]
        real(dp) :: q              ! Particle charge [statcoulomb]
        real(dp) :: dt             ! Timestep [normalized]
        real(dp) :: v0             ! Reference velocity for normalization
        integer :: orbit_model     ! 1=PAULI_PARTICLE, 2=PARTICLE
    end type FullOrbitState

contains

    subroutine init_full_orbit_state(state, s, theta, phi, lambda, v, &
                                      orbit_model, mass_amu, charge_e, dt, v0)
        use simple_coordinates, only: transform_vmec_to_cyl
        use field_coils_cyl, only: evaluate_cyl

        type(FullOrbitState), intent(out) :: state
        real(dp), intent(in) :: s, theta, phi, lambda, v
        integer, intent(in) :: orbit_model
        real(dp), intent(in) :: mass_amu, charge_e, dt, v0

        real(dp) :: x_vmec(3), x_cyl(3), jac_vmec_cyl(3,3)
        real(dp) :: R, phi_cyl, Z
        real(dp) :: B_cyl(3), Bmod, A_cyl(3), gradB(3)
        real(dp) :: b_unit(3), e_perp(3), gradB_cross_b(3)
        real(dp) :: v_par, v_perp, v_cyl(3)
        real(dp) :: norm_gcb

        state%orbit_model = orbit_model
        state%m = mass_amu * p_mass
        state%q = charge_e * e_charge
        state%dt = dt
        state%v0 = v0

        state%z_gc = [s, theta, phi, lambda, v]

        x_vmec = [s**2, theta, phi]
        call transform_vmec_to_cyl(x_vmec, x_cyl, jac_vmec_cyl)

        R = x_cyl(1)
        phi_cyl = x_cyl(2)
        Z = x_cyl(3)

        call evaluate_cyl(R, phi_cyl, Z, A_cyl, B_cyl, Bmod, gradB)

        ! B_cyl from field_coils_cyl is (B_R, R*B_phi, B_Z) - covariant-like
        ! Convert to physical/contravariant for unit vector: (B_R, B_phi, B_Z)
        b_unit(1) = B_cyl(1) / Bmod
        b_unit(2) = B_cyl(2) / (R * Bmod)
        b_unit(3) = B_cyl(3) / Bmod

        v_par = lambda * v
        v_perp = v * sqrt(max(0d0, 1d0 - lambda**2))

        call compute_grad_B_cross_b(gradB, b_unit, gradB_cross_b)
        norm_gcb = sqrt(gradB_cross_b(1)**2 + gradB_cross_b(2)**2 + &
                        gradB_cross_b(3)**2)

        if (norm_gcb > 1d-10) then
            e_perp = gradB_cross_b / norm_gcb
        else
            call compute_fallback_perpendicular(b_unit, e_perp)
        endif

        ! v_cyl is physical velocity: (v_R, v_phi_physical, v_Z)
        ! where v_phi_physical = R * omega (linear velocity in phi direction)
        v_cyl = v_par * b_unit + v_perp * e_perp

        if (orbit_model == 1) then
            state%mu = state%m * v_perp**2 / (2d0 * Bmod)
        else
            state%mu = 0d0
        endif

        state%z(1) = R
        state%z(2) = phi_cyl
        state%z(3) = Z
        ! Canonical momenta: p_R = m*v_R + (q/c)*A_R
        !                    p_phi = m*R*v_phi_phys + (q/c)*R*A_phi = m*R*v_cyl(2) + (q/c)*A_cyl(2)
        !                    p_Z = m*v_Z + (q/c)*A_Z
        ! Note: A_cyl(2) = R*A_phi (covariant)
        state%z(4) = state%m * v_cyl(1) + (state%q / c) * A_cyl(1)
        state%z(5) = state%m * R * v_cyl(2) + (state%q / c) * A_cyl(2)
        state%z(6) = state%m * v_cyl(3) + (state%q / c) * A_cyl(3)

    end subroutine init_full_orbit_state


    subroutine compute_grad_B_cross_b(gradB, b_unit, result)
        real(dp), intent(in) :: gradB(3), b_unit(3)
        real(dp), intent(out) :: result(3)

        result(1) = gradB(2) * b_unit(3) - gradB(3) * b_unit(2)
        result(2) = gradB(3) * b_unit(1) - gradB(1) * b_unit(3)
        result(3) = gradB(1) * b_unit(2) - gradB(2) * b_unit(1)
    end subroutine compute_grad_B_cross_b


    subroutine compute_fallback_perpendicular(b_unit, e_perp)
        real(dp), intent(in) :: b_unit(3)
        real(dp), intent(out) :: e_perp(3)

        real(dp) :: temp(3), norm

        if (abs(b_unit(1)) < 0.9d0) then
            temp = [1d0, 0d0, 0d0]
        else
            temp = [0d0, 1d0, 0d0]
        endif

        e_perp(1) = temp(2) * b_unit(3) - temp(3) * b_unit(2)
        e_perp(2) = temp(3) * b_unit(1) - temp(1) * b_unit(3)
        e_perp(3) = temp(1) * b_unit(2) - temp(2) * b_unit(1)

        norm = sqrt(e_perp(1)**2 + e_perp(2)**2 + e_perp(3)**2)
        e_perp = e_perp / norm
    end subroutine compute_fallback_perpendicular


    subroutine timestep_full_orbit(state, ierr)
        type(FullOrbitState), intent(inout) :: state
        integer, intent(out) :: ierr

        real(dp) :: z_new(6), z_old(6), fvec(6), jac(6,6)
        real(dp) :: delta(6), res_norm, res_norm_old
        integer :: iter, ipiv(6), info
        integer :: i

        ierr = 0
        z_old = state%z

        z_new = z_old
        call explicit_euler_step(state, z_old, z_new)

        res_norm_old = 1d30
        do iter = 1, MAX_NEWTON_ITER
            call eval_midpoint_residual(state, z_old, z_new, fvec)

            res_norm = sqrt(sum(fvec**2))

            if (res_norm < NEWTON_ATOL) exit
            if (res_norm / (sqrt(sum(z_new**2)) + 1d-30) < NEWTON_RTOL) exit

            call eval_midpoint_jacobian(state, z_old, z_new, jac)

            delta = -fvec
            call dgesv(6, 1, jac, 6, ipiv, delta, 6, info)

            if (info /= 0) then
                ierr = 2
                return
            endif

            z_new = z_new + delta

            if (res_norm > res_norm_old * 2d0) then
                ierr = 3
                return
            endif
            res_norm_old = res_norm
        enddo

        if (iter > MAX_NEWTON_ITER) then
            ierr = 1
            return
        endif

        if (z_new(1) <= 0d0) then
            ierr = 4
            return
        endif

        do i = 1, 6
            if (ieee_is_nan(z_new(i))) then
                ierr = 5
                return
            endif
        enddo

        state%z = z_new

    end subroutine timestep_full_orbit


    subroutine explicit_euler_step(state, z_old, z_new)
        use field_coils_cyl, only: evaluate_cyl

        type(FullOrbitState), intent(in) :: state
        real(dp), intent(in) :: z_old(6)
        real(dp), intent(out) :: z_new(6)

        real(dp) :: R, phi, Z, p_R, p_phi, p_Z
        real(dp) :: A_cyl(3), B_cyl(3), Bmod, gradB(3), dA(3,3)
        real(dp) :: v_R, v_phi, v_Z, v_cyl(3)
        real(dp) :: dp_R, dp_phi, dp_Z
        integer :: i

        R = z_old(1)
        phi = z_old(2)
        Z = z_old(3)
        p_R = z_old(4)
        p_phi = z_old(5)
        p_Z = z_old(6)

        call evaluate_cyl(R, phi, Z, A_cyl, B_cyl, Bmod, gradB, dA)

        v_R = (p_R - (state%q / c) * A_cyl(1)) / state%m
        v_phi = (p_phi - (state%q / c) * A_cyl(2)) / (state%m * R**2)
        v_Z = (p_Z - (state%q / c) * A_cyl(3)) / state%m
        v_cyl = [v_R, v_phi, v_Z]

        dp_R = state%m * R * v_phi**2
        do i = 1, 3
            dp_R = dp_R + (state%q / c) * dA(i, 1) * v_cyl(i)
        enddo
        dp_R = dp_R - state%mu * gradB(1)

        dp_phi = 0d0
        do i = 1, 3
            dp_phi = dp_phi + (state%q / c) * dA(i, 2) * v_cyl(i)
        enddo
        dp_phi = dp_phi - state%mu * gradB(2)

        dp_Z = 0d0
        do i = 1, 3
            dp_Z = dp_Z + (state%q / c) * dA(i, 3) * v_cyl(i)
        enddo
        dp_Z = dp_Z - state%mu * gradB(3)

        z_new(1) = R + state%dt * v_R
        z_new(2) = phi + state%dt * v_phi
        z_new(3) = Z + state%dt * v_Z
        z_new(4) = p_R + state%dt * dp_R
        z_new(5) = p_phi + state%dt * dp_phi
        z_new(6) = p_Z + state%dt * dp_Z

    end subroutine explicit_euler_step


    subroutine eval_midpoint_residual(state, z_old, z_new, fvec)
        use field_coils_cyl, only: evaluate_cyl

        type(FullOrbitState), intent(in) :: state
        real(dp), intent(in) :: z_old(6), z_new(6)
        real(dp), intent(out) :: fvec(6)

        real(dp) :: state_mid(6)
        real(dp) :: R_mid, phi_mid, Z_cyl_mid
        real(dp) :: p_R_mid, p_phi_mid, p_Z_mid
        real(dp) :: A_cyl(3), B_cyl(3), Bmod, gradB(3), dA(3,3)
        real(dp) :: v_R, v_phi, v_Z, v_cyl(3)
        real(dp) :: dp_R, dp_phi, dp_Z
        integer :: i

        state_mid = 0.5d0 * (z_old + z_new)
        R_mid = state_mid(1)
        phi_mid = state_mid(2)
        Z_cyl_mid = state_mid(3)
        p_R_mid = state_mid(4)
        p_phi_mid = state_mid(5)
        p_Z_mid = state_mid(6)

        call evaluate_cyl(R_mid, phi_mid, Z_cyl_mid, A_cyl, B_cyl, Bmod, gradB, dA)

        v_R = (p_R_mid - (state%q / c) * A_cyl(1)) / state%m
        v_phi = (p_phi_mid - (state%q / c) * A_cyl(2)) / (state%m * R_mid**2)
        v_Z = (p_Z_mid - (state%q / c) * A_cyl(3)) / state%m
        v_cyl = [v_R, v_phi, v_Z]

        dp_R = state%m * R_mid * v_phi**2
        do i = 1, 3
            dp_R = dp_R + (state%q / c) * dA(i, 1) * v_cyl(i)
        enddo
        dp_R = dp_R - state%mu * gradB(1)

        dp_phi = 0d0
        do i = 1, 3
            dp_phi = dp_phi + (state%q / c) * dA(i, 2) * v_cyl(i)
        enddo
        dp_phi = dp_phi - state%mu * gradB(2)

        dp_Z = 0d0
        do i = 1, 3
            dp_Z = dp_Z + (state%q / c) * dA(i, 3) * v_cyl(i)
        enddo
        dp_Z = dp_Z - state%mu * gradB(3)

        fvec(1) = z_new(1) - z_old(1) - state%dt * v_R
        fvec(2) = z_new(2) - z_old(2) - state%dt * v_phi
        fvec(3) = z_new(3) - z_old(3) - state%dt * v_Z
        fvec(4) = z_new(4) - z_old(4) - state%dt * dp_R
        fvec(5) = z_new(5) - z_old(5) - state%dt * dp_phi
        fvec(6) = z_new(6) - z_old(6) - state%dt * dp_Z

    end subroutine eval_midpoint_residual


    subroutine eval_midpoint_jacobian(state, z_old, z_new, jac)
        type(FullOrbitState), intent(in) :: state
        real(dp), intent(in) :: z_old(6), z_new(6)
        real(dp), intent(out) :: jac(6,6)

        real(dp) :: z_pert(6), fvec0(6), fvec_pert(6)
        real(dp) :: h
        integer :: j

        call eval_midpoint_residual(state, z_old, z_new, fvec0)

        do j = 1, 6
            z_pert = z_new
            h = FD_STEP * max(abs(z_new(j)), 1d0)
            z_pert(j) = z_new(j) + h

            call eval_midpoint_residual(state, z_old, z_pert, fvec_pert)

            jac(:, j) = (fvec_pert - fvec0) / h
        enddo

    end subroutine eval_midpoint_jacobian


    subroutine convert_full_to_gc(state, s, theta, phi, lambda, v)
        use simple_coordinates, only: transform_cyl_to_vmec
        use field_coils_cyl, only: evaluate_cyl

        type(FullOrbitState), intent(in) :: state
        real(dp), intent(out) :: s, theta, phi, lambda, v

        real(dp) :: R, phi_cyl, Z
        real(dp) :: p_R, p_phi, p_Z
        real(dp) :: A_cyl(3), B_cyl(3), Bmod, gradB(3)
        real(dp) :: v_R, v_phi, v_Z
        real(dp) :: v_par, v_mag
        real(dp) :: x_cyl(3), x_vmec(3)
        integer :: ierr

        R = state%z(1)
        phi_cyl = state%z(2)
        Z = state%z(3)
        p_R = state%z(4)
        p_phi = state%z(5)
        p_Z = state%z(6)

        call evaluate_cyl(R, phi_cyl, Z, A_cyl, B_cyl, Bmod, gradB)

        v_R = (p_R - (state%q / c) * A_cyl(1)) / state%m
        v_phi = (p_phi - (state%q / c) * A_cyl(2)) / (state%m * R**2)
        v_Z = (p_Z - (state%q / c) * A_cyl(3)) / state%m

        v_mag = sqrt(v_R**2 + R**2 * v_phi**2 + v_Z**2)

        ! B_cyl is (B_R, R*B_phi, B_Z), v_phi is angular velocity
        ! v_par = v dot b = v_R*B_R + (R*v_phi)*(B_cyl(2)/R) + v_Z*B_Z
        v_par = (v_R * B_cyl(1) + v_phi * B_cyl(2) + v_Z * B_cyl(3)) / Bmod

        x_cyl = [R, phi_cyl, Z]
        call transform_cyl_to_vmec(x_cyl, x_vmec, ierr)

        if (ierr /= 0) then
            s = -1d0
            theta = 0d0
            phi = phi_cyl
        else
            s = x_vmec(1)
            theta = x_vmec(2)
            phi = x_vmec(3)
        endif

        if (v_mag > 1d-30) then
            lambda = v_par / v_mag
        else
            lambda = 0d0
        endif
        v = v_mag

    end subroutine convert_full_to_gc


    function compute_energy(state) result(H)
        use field_coils_cyl, only: evaluate_cyl

        type(FullOrbitState), intent(in) :: state
        real(dp) :: H

        real(dp) :: R, phi, Z, p_R, p_phi, p_Z
        real(dp) :: A_cyl(3), B_cyl(3), Bmod, gradB(3)
        real(dp) :: v_R, v_phi, v_Z, v_sq

        R = state%z(1)
        phi = state%z(2)
        Z = state%z(3)
        p_R = state%z(4)
        p_phi = state%z(5)
        p_Z = state%z(6)

        call evaluate_cyl(R, phi, Z, A_cyl, B_cyl, Bmod, gradB)

        v_R = (p_R - (state%q / c) * A_cyl(1)) / state%m
        v_phi = (p_phi - (state%q / c) * A_cyl(2)) / (state%m * R**2)
        v_Z = (p_Z - (state%q / c) * A_cyl(3)) / state%m

        v_sq = v_R**2 + R**2 * v_phi**2 + v_Z**2

        H = 0.5d0 * state%m * v_sq + state%mu * Bmod

    end function compute_energy

end module orbit_full
