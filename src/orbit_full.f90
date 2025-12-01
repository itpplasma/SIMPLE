module orbit_full
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use util, only: c, e_charge, p_mass
    use neo_biotsavart, only: coils_t, compute_vector_potential, &
                               compute_magnetic_field

    implicit none

    integer, parameter :: MAX_NEWTON_ITER = 20
    real(dp), parameter :: NEWTON_RTOL = 1d-12
    real(dp), parameter :: NEWTON_ATOL = 1d-15
    real(dp), parameter :: FD_STEP = 1d-8

    type :: FullOrbitState
        real(dp) :: z(6)        ! (x, y, z, p_x, p_y, p_z) in Cartesian
        real(dp) :: z_gc(5)     ! (s, theta, phi, lambda, v) guiding-center
        real(dp) :: mu          ! Magnetic moment (Pauli), 0 for full orbit
        real(dp) :: m           ! Particle mass [g]
        real(dp) :: q           ! Particle charge [statcoulomb]
        real(dp) :: dt          ! Timestep [normalized]
        real(dp) :: v0          ! Reference velocity for normalization
        integer :: orbit_model  ! 1=Pauli particle, 2=full orbit
        type(coils_t) :: coils  ! Coil geometry in CGS units
    end type FullOrbitState

contains

    subroutine init_full_orbit_state(state, s, theta, phi, lambda, v, &
                                      orbit_model, mass_amu, charge_e, dt, v0, &
                                      coils)
        use simple_coordinates, only: transform_vmec_to_cart

        type(FullOrbitState), intent(out) :: state
        real(dp), intent(in) :: s, theta, phi, lambda, v
        integer, intent(in) :: orbit_model
        real(dp), intent(in) :: mass_amu, charge_e, dt, v0
        type(coils_t), intent(in) :: coils

        real(dp) :: x_vmec(3), x_cart(3)
        real(dp) :: A_cart(3), B_cart(3), Bmod, gradB(3)
        real(dp) :: b_unit(3), e_perp(3), gradB_cross_b(3)
        real(dp) :: v_par, v_perp, v_vec(3)
        real(dp) :: norm_gcb

        state%orbit_model = orbit_model
        state%m = mass_amu * p_mass
        state%q = charge_e * e_charge
        state%dt = dt
        state%v0 = v0
        state%coils = coils

        state%z_gc = [s, theta, phi, lambda, v]

        x_vmec = [s**2, theta, phi]
        call transform_vmec_to_cart(x_vmec, x_cart)

        call evaluate_field_and_gradB(state, x_cart, A_cart, B_cart, Bmod, gradB)

        b_unit = B_cart / Bmod

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

        v_vec = v_par * b_unit + v_perp * e_perp

        if (orbit_model == 1) then
            state%mu = state%m * v_perp**2 / (2d0 * Bmod)
        else
            state%mu = 0d0
        endif

        state%z(1:3) = x_cart
        state%z(4:6) = state%m * v_vec + (state%q / c) * A_cart

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


    subroutine evaluate_field_and_gradB(state, x, A_cart, B_cart, Bmod, gradB)
        type(FullOrbitState), intent(in) :: state
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A_cart(3), B_cart(3), Bmod, gradB(3)

        real(dp) :: x_pert(3), Bp(3), Bm(3)
        real(dp) :: h, Bp_mod, Bm_mod
        integer :: i

        A_cart = compute_vector_potential(state%coils, x)
        B_cart = compute_magnetic_field(state%coils, x)
        Bmod = sqrt(B_cart(1)**2 + B_cart(2)**2 + B_cart(3)**2)

        do i = 1, 3
            x_pert = x
            h = max(1d-6 * max(abs(x(i)), 1d0), 1d-6)
            x_pert(i) = x(i) + h
            Bp = compute_magnetic_field(state%coils, x_pert)
            Bp_mod = sqrt(Bp(1)**2 + Bp(2)**2 + Bp(3)**2)

            x_pert(i) = x(i) - h
            Bm = compute_magnetic_field(state%coils, x_pert)
            Bm_mod = sqrt(Bm(1)**2 + Bm(2)**2 + Bm(3)**2)

            gradB(i) = (Bp_mod - Bm_mod) / (2d0 * h)
        enddo
    end subroutine evaluate_field_and_gradB


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

        do i = 1, 6
            if (ieee_is_nan(z_new(i))) then
                ierr = 5
                return
            endif
        enddo

        state%z = z_new

    end subroutine timestep_full_orbit


    subroutine explicit_euler_step(state, z_old, z_new)
        type(FullOrbitState), intent(in) :: state
        real(dp), intent(in) :: z_old(6)
        real(dp), intent(out) :: z_new(6)

        real(dp) :: x_old(3), p_old(3)
        real(dp) :: xdot(3), pdot(3)

        x_old = z_old(1:3)
        p_old = z_old(4:6)

        call compute_rhs(state, x_old, p_old, xdot, pdot)

        z_new(1:3) = x_old + state%dt * xdot
        z_new(4:6) = p_old + state%dt * pdot

    end subroutine explicit_euler_step


    subroutine eval_midpoint_residual(state, z_old, z_new, fvec)
        type(FullOrbitState), intent(in) :: state
        real(dp), intent(in) :: z_old(6), z_new(6)
        real(dp), intent(out) :: fvec(6)

        real(dp) :: x_old(3), p_old(3)
        real(dp) :: x_new(3), p_new(3)
        real(dp) :: x_mid(3), p_mid(3)
        real(dp) :: xdot(3), pdot(3)

        x_old = z_old(1:3)
        p_old = z_old(4:6)

        x_new = z_new(1:3)
        p_new = z_new(4:6)

        x_mid = 0.5d0 * (x_old + x_new)
        p_mid = 0.5d0 * (p_old + p_new)

        call compute_rhs(state, x_mid, p_mid, xdot, pdot)

        fvec(1:3) = x_new - x_old - state%dt * xdot
        fvec(4:6) = p_new - p_old - state%dt * pdot

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
        use simple_coordinates, only: transform_cart_to_cyl, &
                                      transform_cyl_to_vmec

        type(FullOrbitState), intent(in) :: state
        real(dp), intent(out) :: s, theta, phi, lambda, v

        real(dp) :: x_cart(3), p(3)
        real(dp) :: A_cart(3), B_cart(3), Bmod, gradB(3)
        real(dp) :: v_vec(3)
        real(dp) :: v_par, v_mag
        real(dp) :: x_cyl(3), x_vmec(3)
        integer :: ierr

        x_cart = state%z(1:3)
        p = state%z(4:6)

        call evaluate_field_and_gradB(state, x_cart, A_cart, B_cart, Bmod, gradB)

        v_vec = (p - (state%q / c) * A_cart) / state%m

        v_mag = sqrt(v_vec(1)**2 + v_vec(2)**2 + v_vec(3)**2)

        v_par = (v_vec(1) * B_cart(1) + v_vec(2) * B_cart(2) + &
                 v_vec(3) * B_cart(3)) / Bmod

        call transform_cart_to_cyl(x_cart, x_cyl)
        call transform_cyl_to_vmec(x_cyl, x_vmec, ierr)

        if (ierr /= 0) then
            s = -1d0
            theta = 0d0
            phi = x_cyl(2)
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
        type(FullOrbitState), intent(in) :: state
        real(dp) :: H

        real(dp) :: x_cart(3), p(3)
        real(dp) :: A_cart(3), B_cart(3), Bmod, gradB(3)
        real(dp) :: v_vec(3), v_sq

        x_cart = state%z(1:3)
        p = state%z(4:6)

        call evaluate_field_and_gradB(state, x_cart, A_cart, B_cart, Bmod, gradB)

        v_vec = (p - (state%q / c) * A_cart) / state%m

        v_sq = v_vec(1)**2 + v_vec(2)**2 + v_vec(3)**2

        H = 0.5d0 * state%m * v_sq + state%mu * Bmod

    end function compute_energy


    subroutine compute_rhs(state, x, p, xdot, pdot)
        type(FullOrbitState), intent(in) :: state
        real(dp), intent(in) :: x(3), p(3)
        real(dp), intent(out) :: xdot(3), pdot(3)

        real(dp) :: A_cart(3), B_cart(3), Bmod, gradB(3)
        real(dp) :: v(3), v_cross_B(3)
        real(dp) :: q_over_c

        call evaluate_field_and_gradB(state, x, A_cart, B_cart, Bmod, gradB)

        v = (p - (state%q / c) * A_cart) / state%m

        q_over_c = state%q / c

        call cross_product3(v, B_cart, v_cross_B)

        pdot = q_over_c * v_cross_B
        if (state%mu /= 0d0) then
            pdot = pdot - state%mu * gradB
        endif

        xdot = v
    end subroutine compute_rhs


    subroutine cross_product3(a, b, c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp), intent(out) :: c(3)

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end subroutine cross_product3

end module orbit_full
