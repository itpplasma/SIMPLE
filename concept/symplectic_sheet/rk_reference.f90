module rk_reference
    !> Reference integrator: adaptive Dormand-Prince 5(4) on the piecewise
    !> Hamiltonian flow of the sharp sheet, with the branch switched exactly at
    !> the crossing. Each crossing is located by false position on the
    !> length of a single Dormand-Prince step from the last accepted state.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, grad_hamiltonian, sheet_distance, &
        region_of
    implicit none
    private

    public :: rk_advance

    real(dp), parameter :: C2 = 1.0_dp/5, C3 = 3.0_dp/10, C4 = 4.0_dp/5, &
        C5 = 8.0_dp/9
    real(dp), parameter :: A21 = 1.0_dp/5
    real(dp), parameter :: A31 = 3.0_dp/40, A32 = 9.0_dp/40
    real(dp), parameter :: A41 = 44.0_dp/45, A42 = -56.0_dp/15, A43 = 32.0_dp/9
    real(dp), parameter :: A51 = 19372.0_dp/6561, A52 = -25360.0_dp/2187, &
        A53 = 64448.0_dp/6561, A54 = -212.0_dp/729
    real(dp), parameter :: A61 = 9017.0_dp/3168, A62 = -355.0_dp/33, &
        A63 = 46732.0_dp/5247, A64 = 49.0_dp/176, A65 = -5103.0_dp/18656
    real(dp), parameter :: B1 = 35.0_dp/384, B3 = 500.0_dp/1113, &
        B4 = 125.0_dp/192, B5 = -2187.0_dp/6784, B6 = 11.0_dp/84
    real(dp), parameter :: E1 = 71.0_dp/57600, E3 = -71.0_dp/16695, &
        E4 = 71.0_dp/1920, E5 = -17253.0_dp/339200, E6 = 22.0_dp/525, &
        E7 = -1.0_dp/40
    real(dp), parameter :: LANDING_TOL = 1.0e-14_dp, TOUCH_FRACTION = 1.0e-10_dp

contains

    subroutine rk_advance(par, region, z, t_span, tol, h_try, crossings)
        !> Advance z by physical time t_span. region is updated at crossings;
        !> h_try carries the step-size guess between calls.
        type(field_params_t), intent(in) :: par
        integer, intent(inout) :: region
        real(dp), intent(inout) :: z(4)
        real(dp), intent(in) :: t_span, tol
        real(dp), intent(inout) :: h_try
        integer, intent(inout) :: crossings

        real(dp) :: t, h, z_new(4), err, h_land
        logical :: landed

        t = 0.0_dp
        do while (t < t_span)
            h = min(h_try, t_span - t)
            call dp_step(par, region, z, h, z_new, err)
            err = err/tol
            if (err > 1.0_dp) then
                h_try = h*max(0.2_dp, 0.9_dp*err**(-0.2_dp))
                cycle
            end if
            if (region_of(par, z_new(1:2)) /= region) then
                call locate(par, region, z, h, h_land, z_new, landed)
                if (.not. landed) then
                    h_try = 0.5_dp*h
                    cycle
                end if
                if (h_land <= TOUCH_FRACTION*h) then
                    ! The orbit grazes the sheet at the step start: it belongs to
                    ! the other branch. Switch back and redo the step.
                    region = 3 - region
                    crossings = crossings - 1
                    cycle
                end if
                h = h_land
                region = 3 - region
                crossings = crossings + 1
            end if
            z = z_new
            t = t + h
            if (err > 0.0_dp) then
                h_try = h*min(5.0_dp, 0.9_dp*err**(-0.2_dp))
            else
                h_try = 5.0_dp*h
            end if
        end do
    end subroutine rk_advance

    subroutine locate(par, region, z, h, h_land, z_land, landed)
        !> Step length to the sheet by Illinois false position on single steps.
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z(4), h
        real(dp), intent(out) :: h_land, z_land(4)
        logical, intent(out) :: landed

        real(dp) :: a, b, fa, fb, c, fc, err
        integer :: it, side

        ! A start on the sheet right after a landing counts as the home side.
        a = 0.0_dp
        fa = merge(-1.0_dp, 1.0_dp, region == 1)* &
             max(abs(sheet_distance(par, z(1:2))), tiny(1.0_dp))
        b = h
        call dp_step(par, region, z, b, z_land, err)
        fb = sheet_distance(par, z_land(1:2))
        landed = .false.
        side = 0
        c = b
        do it = 1, 200
            c = (a*fb - b*fa)/(fb - fa)
            call dp_step(par, region, z, c, z_land, err)
            fc = sheet_distance(par, z_land(1:2))
            if (abs(fc) <= LANDING_TOL) then
                landed = .true.
                exit
            end if
            if (fc*fb > 0.0_dp) then
                b = c
                fb = fc
                if (side == -1) fa = 0.5_dp*fa
                side = -1
            else
                a = c
                fa = fc
                if (side == 1) fb = 0.5_dp*fb
                side = 1
            end if
        end do
        h_land = c
    end subroutine locate

    subroutine dp_step(par, region, z, h, z_new, err)
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z(4), h
        real(dp), intent(out) :: z_new(4), err

        real(dp) :: k1(4), k2(4), k3(4), k4(4), k5(4), k6(4), k7(4), e(4)

        k1 = rhs(z)
        k2 = rhs(z + h*A21*k1)
        k3 = rhs(z + h*(A31*k1 + A32*k2))
        k4 = rhs(z + h*(A41*k1 + A42*k2 + A43*k3))
        k5 = rhs(z + h*(A51*k1 + A52*k2 + A53*k3 + A54*k4))
        k6 = rhs(z + h*(A61*k1 + A62*k2 + A63*k3 + A64*k4 + A65*k5))
        z_new = z + h*(B1*k1 + B3*k3 + B4*k4 + B5*k5 + B6*k6)
        k7 = rhs(z_new)
        e = h*(E1*k1 + E3*k3 + E4*k4 + E5*k5 + E6*k6 + E7*k7)
        err = maxval(abs(e)/max(abs(z), abs(z_new), 1.0e-3_dp))

    contains

        function rhs(y) result(f)
            real(dp), intent(in) :: y(4)
            real(dp) :: f(4), g(4)

            call grad_hamiltonian(par, region, y, g)
            f = [g(3), g(4), -g(1), -g(2)]
        end function rhs

    end subroutine dp_step

end module rk_reference
