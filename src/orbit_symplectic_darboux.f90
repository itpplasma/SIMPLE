module orbit_symplectic_darboux
    !> Implicit midpoint rule in a global canonical chart that is regular at the
    !> magnetic axis (experimental, integmode = DARBOUX_MIDPOINT).
    !>
    !> SIMPLE's guiding-centre phase space carries the canonical pairs
    !> (theta, p_theta) and (phi, p_phi); the radial coordinate r only labels
    !> p_theta(r, theta, phi, p_phi) (Albert, Kasilov, Kernbichler, JCP 403,
    !> 109065 (2020)). Near the axis p_theta = A_theta/ro0 + h_theta v_par ~ s,
    !> so (theta, p_theta) are polar action-angle-like coordinates. With
    !>   J = sigma p_theta >= 0, alpha = sigma theta, sigma = sign(dp_theta/dr),
    !>   X = sqrt(2 J) cos(alpha), Y = sqrt(2 J) sin(alpha),
    !> one has dp_theta ^ dtheta = dJ ^ dalpha = dX ^ dY exactly, so
    !> (X, Y, phi, p_phi) is a Darboux chart with the axis at X = Y = 0. The
    !> Hamiltonian is smooth there whenever the field representation is smooth
    !> in the pseudo-Cartesian coordinates sqrt(s) (cos, sin) theta. The step is
    !> the implicit midpoint rule applied in this single global chart: there is
    !> no axis reflection, no radius floor and no chart handoff.
    !>
    !> Field data are still evaluated in (r, theta, phi): every canonical point
    !> is mapped back with theta = sigma atan2(Y, X) and a scalar Newton solve
    !> p_theta(r, theta, phi, p_phi) = sigma J (non-canonical quadrature). Only
    !> first field derivatives enter; the canonical gradient
    !>   dH/dX = X dH/dJ - Y (dH/dalpha)/(2 J),
    !>   dH/dY = Y dH/dJ + X (dH/dalpha)/(2 J)
    !> is a sum of products of regular factors, so no 1/rho cancellation occurs.
    !> The implicit stage is solved by fixed-point iteration, which converges
    !> because h |D(J grad H)| << 1 for the guiding-centre time step.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: pi, twopi
    use field_can_mod, only: field_can_t, get_derivatives, eval_field => evaluate
    use orbit_symplectic_base, only: symplectic_integrator_t, sympl_rmax, &
        SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_OUTSIDE_DOMAIN, &
        SYMPLECTIC_STEP_MAXITER

    implicit none

    private
    public :: orbit_timestep_sympl_darboux_midpoint, darboux_state_from_field, &
        darboux_eval, darboux_rhs, darboux_stats_reset, darboux_stats_t, darboux_stats

    integer, parameter :: maxit_stage = 100
    integer, parameter :: maxit_invert = 50

    type :: darboux_stats_t
        integer(8) :: steps = 0
        integer(8) :: stage_iterations = 0
        integer(8) :: invert_iterations = 0
        integer(8) :: failures = 0
        integer :: max_stage_iterations = 0
    end type darboux_stats_t

    type(darboux_stats_t) :: darboux_stats
    !$omp threadprivate(darboux_stats)

contains

    logical function darboux_debug()
        !> Prototype diagnostics: SIMPLE_DARBOUX_DEBUG=1 prints failed steps.
        integer :: length

        call get_environment_variable('SIMPLE_DARBOUX_DEBUG', length=length)
        darboux_debug = length > 0
    end function darboux_debug

    subroutine darboux_stats_reset()
        darboux_stats = darboux_stats_t()
    end subroutine darboux_stats_reset

    pure real(dp) function wrap_near(angle, reference)
        !> Representative of angle mod 2 pi closest to reference.
        real(dp), intent(in) :: angle, reference

        wrap_near = reference + modulo(angle - reference + pi, twopi) - pi
    end function wrap_near

    subroutine darboux_state_from_field(f, z, sigma, w)
        !> Canonical Darboux coordinates of the phase point whose field state f
        !> holds p_theta (f%pth) at non-canonical coordinates z.
        type(field_can_t), intent(in) :: f
        real(dp), intent(in) :: z(4)
        real(dp), intent(out) :: sigma, w(4)
        real(dp) :: rad

        sigma = sign(1.0_dp, f%dpth(1))
        rad = sqrt(2.0_dp*max(sigma*f%pth, 0.0_dp))
        w(1) = rad*cos(sigma*z(2))
        w(2) = rad*sin(sigma*z(2))
        w(3) = z(3)
        w(4) = z(4)
    end subroutine darboux_state_from_field

    subroutine darboux_rhs(f, w, sigma, rhs)
        !> Canonical vector field J grad H at Darboux point w from a field state
        !> f evaluated (with get_derivatives) at the matching (r, theta).
        type(field_can_t), intent(in) :: f
        real(dp), intent(in) :: w(4), sigma
        real(dp), intent(out) :: rhs(4)

        real(dp) :: x, y, jact, h_p, h_alpha, h_j, h_ph, h_pph

        call regular_xy(w, x, y, jact)
        ! Canonical derivatives at fixed p_theta, then to (J, alpha).
        h_p = f%dH(1)/f%dpth(1)
        h_alpha = sigma*(f%dH(2) - h_p*f%dpth(2))
        h_j = sigma*h_p
        h_ph = f%dH(3) - h_p*f%dpth(3)
        h_pph = f%dH(4) - h_p*f%dpth(4)

        ! dX/dt = -dH/dY, dY/dt = dH/dX, dphi/dt = dH/dp_phi, dp_phi/dt = -dH/dphi
        rhs(1) = -(y*h_j + x*h_alpha/(2.0_dp*jact))
        rhs(2) = x*h_j - y*h_alpha/(2.0_dp*jact)
        rhs(3) = h_pph
        rhs(4) = -h_ph
    end subroutine darboux_rhs

    pure subroutine regular_xy(w, x, y, jact)
        !> The chart is regular at the axis; nudge off the exact origin so that
        !> theta is defined. The resulting force error is O(1e-140).
        real(dp), intent(in) :: w(4)
        real(dp), intent(out) :: x, y, jact
        real(dp), parameter :: j_tiny = 1.0e-280_dp

        x = w(1)
        y = w(2)
        jact = 0.5_dp*(x*x + y*y)
        if (jact < j_tiny) then
            x = sqrt(2.0_dp*j_tiny)
            y = 0.0_dp
            jact = j_tiny
        end if
    end subroutine regular_xy

    subroutine darboux_eval(w, sigma, theta_ref, r, theta, f, rhs, status)
        !> Map Darboux point w to (r, theta) by Newton on p_theta(r) = sigma J
        !> (r on entry is the initial guess), evaluate the field there and
        !> return the canonical vector field.
        real(dp), intent(in) :: w(4), sigma, theta_ref
        real(dp), intent(inout) :: r
        real(dp), intent(out) :: theta
        type(field_can_t), intent(inout) :: f
        real(dp), intent(out) :: rhs(4)
        integer, intent(out) :: status

        real(dp) :: x, y, jact, dr
        integer :: it

        status = SYMPLECTIC_STEP_OK
        call regular_xy(w, x, y, jact)
        theta = wrap_near(sigma*atan2(y, x), theta_ref)
        if (.not. (r > 0.0_dp)) r = 1.0e-12_dp
        do it = 1, maxit_invert
            call eval_field(f, r, theta, w(3), 0)
            call get_derivatives(f, w(4))
            darboux_stats%invert_iterations = darboux_stats%invert_iterations + 1
            dr = (f%pth - sigma*jact)/f%dpth(1)
            ! Converged: the field at r is exact to the inversion roundoff.
            if (abs(dr) <= 4.0_dp*epsilon(1.0_dp)*r) exit
            if (r - dr <= 0.0_dp) then
                r = 0.5_dp*r
            else
                r = r - dr
            end if
            if (r > 2.0_dp*sympl_rmax) then
                status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
        end do
        if (it > maxit_invert) then
            if (darboux_debug()) print '(a,4es12.4)', 'invert fail', r, dr, &
                f%pth, sigma*jact
            status = SYMPLECTIC_STEP_MAXITER
            return
        end if
        call darboux_rhs(f, w, sigma, rhs)
    end subroutine darboux_eval

    subroutine darboux_midpoint_step(si, f, status)
        !> One implicit midpoint step w1 = w0 + h F((w0 + w1)/2) in the Darboux
        !> chart. The midpoint radius r_m is an extra unknown with the
        !> constraint p_theta(r_m, theta_m, phi_m, p_phi,m) = sigma J_m, solved
        !> jointly with the stage: one field evaluation per iteration.
        !> On failure si and f are left unchanged.
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: status

        type(field_can_t) :: fm
        real(dp) :: sigma, w0(4), wm(4), wm_new(4), w1(4), rhs(4), tol(4)
        real(dp) :: r_mid, th_mid, r_end, th_end, xy_scale, x, y, jact, dr
        real(dp) :: rad, wt(4), hist(11, 16)
        integer :: it

        status = SYMPLECTIC_STEP_OK
        ! The caller guarantees field data at si%z and the canonical p_theta in
        ! f%pth; derivatives are recomputed because initialisation, reseeding
        ! and collision kicks only refresh values (get_val).
        fm = f
        call get_derivatives(fm, si%z(4))
        fm%pth = f%pth
        call darboux_state_from_field(fm, si%z, sigma, w0)
        xy_scale = sqrt(2.0_dp*abs(fm%dpth(1))*max(si%z(1), 1.0e-8_dp))

        ! Explicit predictor from the field at the current point.
        call darboux_rhs(fm, w0, sigma, rhs)
        wm = w0 + 0.5_dp*si%dt*rhs
        r_mid = si%z(1)

        do it = 1, maxit_stage
            call regular_xy(wm, x, y, jact)
            th_mid = wrap_near(sigma*atan2(y, x), si%z(2))
            call eval_field(fm, r_mid, th_mid, wm(3), 0)
            call get_derivatives(fm, wm(4))
            darboux_stats%stage_iterations = darboux_stats%stage_iterations + 1
            dr = (fm%pth - sigma*jact)/fm%dpth(1)
            ! Force at the phase point actually evaluated, (r_m, theta_m):
            ! its Darboux radius sqrt(2 sigma p_theta(r_m)) replaces |w_m|.
            ! Mixing the field at r_m with an inconsistent |w_m| diverges
            ! near the axis, where the ratio of the two radii is unbounded.
            rad = sqrt(2.0_dp*max(sigma*fm%pth, 0.0_dp)/(2.0_dp*jact))
            wt = [x*rad, y*rad, wm(3), wm(4)]
            call darboux_rhs(fm, wt, sigma, rhs)
            wm_new = w0 + 0.5_dp*si%dt*rhs
            tol(1:2) = si%rtol*max(abs(wm_new(1:2)), xy_scale) + si%atol
            tol(3) = si%rtol*max(abs(wm_new(3)), 1.0_dp) + si%atol
            tol(4) = si%rtol*max(abs(wm_new(4)), abs(w0(4)), 1.0_dp) + si%atol
            ! The radius is converged when its residual is below what the
            ! (X, Y) tolerance implies for J = (X**2 + Y**2)/2.
            if (all(abs(wm_new - wm) <= tol) .and. abs(dr) <= si%rtol*r_mid &
                + (abs(x)*tol(1) + abs(y)*tol(2))/abs(fm%dpth(1))) then
                wm = wm_new
                exit
            end if
            if (it <= 16) hist(:, it) = [r_mid, dr, wm_new - wm, wm(1:2), jact, &
                fm%dpth(1), fm%pth]
            wm = wm_new
            if (r_mid - dr <= 0.0_dp) then
                r_mid = 0.5_dp*r_mid
            else
                r_mid = r_mid - dr
            end if
            if (r_mid > 2.0_dp*sympl_rmax) then
                status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
        end do
        darboux_stats%max_stage_iterations = &
            max(darboux_stats%max_stage_iterations, it)
        if (it > maxit_stage) then
            if (darboux_debug()) then
                print '(a,4es13.5)', 'stage fail w0=', w0
                do it = 1, 16
                    print '(a,i3,11es11.3)', 'it', it, hist(:, it)
                end do
            end if
            darboux_stats%failures = darboux_stats%failures + 1
            status = SYMPLECTIC_STEP_MAXITER
            return
        end if

        w1 = 2.0_dp*wm - w0
        r_end = 2.0_dp*r_mid - si%z(1)
        if (.not. (r_end > 0.0_dp)) r_end = r_mid
        call darboux_eval(w1, sigma, th_mid, r_end, th_end, fm, rhs, status)
        if (status /= SYMPLECTIC_STEP_OK) return
        if (r_end > sympl_rmax) then
            status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
            return
        end if

        si%z = [r_end, th_end, w1(3), w1(4)]
        f = fm
        ! Keep the canonical momentum exactly; the field state agrees with it
        ! to the inversion roundoff.
        f%pth = sigma*0.5_dp*(w1(1)**2 + w1(2)**2)
        si%pthold = f%pth
        darboux_stats%steps = darboux_stats%steps + 1
    end subroutine darboux_midpoint_step

    recursive subroutine orbit_timestep_sympl_darboux_midpoint(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field
        integer :: ktau

        ierr = SYMPLECTIC_STEP_OK
        do ktau = 1, si%ntau
            accepted_integrator = si
            accepted_field = f
            call darboux_midpoint_step(si, f, ierr)
            if (ierr /= SYMPLECTIC_STEP_OK .and. darboux_debug()) &
                print '(a,i3,a,4es13.5,a,2es12.4)', 'darboux fail status', ierr, &
                ' z=', accepted_integrator%z, ' pth,dpth_r=', &
                accepted_field%pth, accepted_field%dpth(1)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                return
            end if
        end do
    end subroutine orbit_timestep_sympl_darboux_midpoint

end module orbit_symplectic_darboux
