module orbit_symplectic_quasi

    use util, only: pi, twopi
    use field_can_mod, only: eval_field => evaluate, field_can_t, get_derivatives, &
                             get_derivatives2
use orbit_symplectic_base, only: symplectic_integrator_t, multistage_integrator_t, &
   orbit_timestep_quasi_i, coeff_rk_gauss, coeff_rk_lobatto, f_rk_lobatto, sympl_rmax, &
          SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_OUTSIDE_DOMAIN, SYMPLECTIC_STEP_MAXITER, &
                                     SYMPLECTIC_STEP_LINEAR_SOLVE
    use orbit_symplectic_base, only: S_MAX, symplectic_newton_warning_mode, &
                                     EXPL_IMPL_EULER, IMPL_EXPL_EULER
    use orbit_systems_generated, only: euler1_residual, euler2_residual, &
                 midpoint_residual, gauss1_residual, gauss2_residual, gauss3_residual, &
                                       gauss4_residual, gc_rates_canonical, &
                                       gc_rates_canonical_jacobian
use fortnum_multiroot, only: multiroot_hybrids, multiroot_hybrids_warm_t
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use vector_potentail_mod, only: torflux
    use diag_counters, only: count_event, EVT_R_NEGATIVE, EVT_WARNING_MAXIT_ACCEPT

implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)
save


type(field_can_t) :: f
type(symplectic_integrator_t) :: si
  !$omp threadprivate(f, si)



! Step size carried across macro-steps by the adaptive explicit integrators.
! Without this each macro-step restarts the step-size search from a default
! guess and spends most of its evaluations re-discovering a step it already
! knew, which inflates the evaluation count by more than an order of magnitude
! and would make any cost comparison meaningless.
!
! The carry belongs to one orbit, not to the thread: leaving it set when the
! next orbit starts makes that orbit's first macro-step depend on which orbit
! the thread happened to run before it, so the same particle gives different
! evaluation counts (and, through the controller's history, slightly different
! trajectories) depending on scheduling. orbit_sympl_init resets it.
    real(dp) :: explicit_h_carry = 0d0
    !$omp threadprivate(explicit_h_carry)


contains

    subroutine reset_explicit_step_carry
        explicit_h_carry = 0d0
    end subroutine reset_explicit_step_carry



    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Restore the chart invariants the rest of SIMPLE maintains inside a step, and
! report the status the caller's boundary bisection expects.
!
! The symplectic schemes reach both of these through their Newton solve, which
! inspects every iterate. The explicit drivers hand a whole substep to an
! adaptive integrator in one call and see nothing between its endpoints, so
! both have to be re-imposed here, at the same dtaumin granularity the
! symplectic schemes work at.
!
! Leaving the plasma. A marker that crossed the last closed flux surface used
! to be integrated onwards through extrapolated splines and reported as still
! confined: at a 1 s trace that put 44 of 1000 markers at s > 1, one as far out
! as s = 2.8, which is why the explicit modes found 76 losses where the
! symplectic schemes and RK4/5 both independently found about 207. Returning
! OUTSIDE_DOMAIN hands the step to locate_symplectic_boundary, which bisects on
! si%dt exactly as it does for the symplectic schemes.
!
! Crossing the axis. r is a radius, so r < 0 is not a position but a chart
! artefact of passing through the magnetic axis; every other integrator in
! SIMPLE continues on the opposite ray, (r, theta) -> (|r|, theta + pi), rather
! than the pre-#370 teleport to r = 0.01. Without it the next substep evaluates
! the field at a negative radius, outside the spline domain entirely.
!
! h_entry is the step-size carry as it stood before the substep. A step that
! does not complete must not leave its trial step size behind: the boundary
! bisection re-runs the step up to 64 times on ever smaller fractions of dt,
! and without this the carry would come back at bisection scale and throttle
! every later step of the orbit.
    subroutine explicit_step_finish(h_entry, ierr)
        real(dp), intent(in) :: h_entry
        integer, intent(out) :: ierr

        ! A non-finite state has to be caught before the comparisons below, because
        ! every comparison against NaN is false: NaN > sympl_rmax does not fire, so a
        ! NaN orbit would be waved through as a good step and reported as confined.
        ! The adaptive drivers usually fail first -- a NaN error norm never satisfies
        ! err <= 1, so the step is rejected until the size underflows -- but the
        ! fixed-step TDRK has no error estimate to fail, and nothing else on this path
        ! inspects the state.
        if (any(si%z(1:4) /= si%z(1:4))) then
            ierr = SYMPLECTIC_STEP_MAXITER
            explicit_h_carry = h_entry
            return
        end if

        if (si%z(1) < 0d0) then
            call count_event(EVT_R_NEGATIVE)
            si%z(1) = -si%z(1)
            si%z(2) = si%z(2) + pi
        end if

        ierr = SYMPLECTIC_STEP_OK
        if (si%z(1) > sympl_rmax) ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
        if (ierr /= SYMPLECTIC_STEP_OK) explicit_h_carry = h_entry
    end subroutine explicit_step_finish



  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine f_ode(tau, z, zdot)
  !

  real(dp), intent(in)  :: tau
  real(dp), intent(in)  :: z(:)
  real(dp), intent(out) :: zdot(:)

  call eval_field(f, z(1), z(2), z(3), 0)
  call get_derivatives(f, z(4))
        call gc_rates_canonical(f, zdot)

end subroutine f_ode

    !> Second time derivative F'(z) F(z), sharing the canonical RHS and its
    !> generated Jacobian so explicit integrators use the same equations.
    subroutine g_ode(tau, z, g)
        real(dp), intent(in)  :: tau
        real(dp), intent(in)  :: z(:)
        real(dp), intent(out) :: g(:)

        real(dp) :: zdot(4), dF(4, 4)

        call eval_field(f, z(1), z(2), z(3), 2)
        call get_derivatives2(f, z(4))
        call gc_rates_canonical(f, zdot)
        call gc_rates_canonical_jacobian(f, dF)

        g(1:4) = matmul(dF, zdot)
    end subroutine g_ode

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Adapter from g_ode to fortnum's ode_rhs2_t.
    subroutine g_ode_fortnum(t, y, g, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(in), optional :: ctx

        call g_ode(t, y, g)
    end subroutine g_ode_fortnum

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine orbit_timestep_rk45(ierr)
  !
  use odeint_allroutines_sub, only : odeint_allroutines
  integer, intent(out) :: ierr
  integer :: ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    call odeint_allroutines(si%z, 4, ktau*si%dt, (ktau+1)*si%dt, si%rtol, f_ode)
    ktau = ktau+1
  end do
end subroutine orbit_timestep_rk45

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Adapter from f_ode to fortnum's ode_rhs_t, which carries an optional context
! argument that this right-hand side does not need.
    subroutine f_ode_fortnum(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        call f_ode(t, y, dydt)
    end subroutine f_ode_fortnum

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Gauss-Radau (IAS15-style, 15th order) on the same canonical Hamiltonian the
! symplectic schemes use.
    subroutine orbit_timestep_radau15(ierr)
        !
        use fortnum_ode_gauss_radau, only: ode_integrate_radau
        use fortnum_ode, only: ode_problem_t, ode_solution_t
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK
        integer, intent(out) :: ierr
        integer :: ktau, nlast
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp) :: h_entry

        ierr = 0
        ktau = 0
        problem%rhs => f_ode_fortnum
        problem%rtol = si%rtol
        problem%atol = si%atol
        allocate (problem%y0(4))
        do while (ktau .lt. si%ntau)
            h_entry = explicit_h_carry
            problem%t0 = ktau*si%dt
            problem%t1 = (ktau + 1)*si%dt
            problem%y0 = si%z(1:4)
            problem%h0 = explicit_h_carry
            call ode_integrate_radau(problem, solution, status)
            if (status%code /= FORTNUM_OK) then
                ierr = SYMPLECTIC_STEP_MAXITER
                explicit_h_carry = h_entry
                return
            end if
            nlast = size(solution%t)
            si%z(1:4) = solution%y(:, nlast)
            if (nlast > 1) explicit_h_carry = solution%h(nlast)
            call explicit_step_finish(h_entry, ierr)
            if (ierr /= SYMPLECTIC_STEP_OK) return
            ktau = ktau + 1
        end do
        call sync_field_to_state
    end subroutine orbit_timestep_radau15

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Re-evaluate the field at the converged end-of-step state.
!
! Only z(1:4) are integrated; the caller reconstructs the parallel velocity as
! z(5) = f%vpar/(pabs*sqrt(2)), reading f%vpar off the module-global field_can_t
! that f_ode fills as a side effect. After a multi-stage step that side effect
! belongs to the last interior stage node, not to the end of the step, so z(5)
! carries an error the tolerance cannot remove: measured on the sweep,
! components 1-4 of the final state converged to 5e-9 while z(5) stalled at
! 3.3e-4 and got no better from rtol 1e-8 to 1e-10.
!
! One evaluation per macro-step, against the hundreds a step already costs.
    subroutine sync_field_to_state
        real(dp) :: zdot_discard(4)

        call f_ode(0d0, si%z(1:4), zdot_discard)
    end subroutine sync_field_to_state

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Gragg-Bulirsch-Stoer extrapolation on the same canonical Hamiltonian.
    subroutine orbit_timestep_gbs16(ierr)
        !
        use fortnum_ode_extrapolation, only: ode_integrate_gbs
        use fortnum_ode, only: ode_problem_t, ode_solution_t
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK
        integer, intent(out) :: ierr
        integer :: ktau, nlast
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp) :: h_entry

        ierr = 0
        ktau = 0
        problem%rhs => f_ode_fortnum
        problem%rtol = si%rtol
        problem%atol = si%atol
        allocate (problem%y0(4))
        do while (ktau .lt. si%ntau)
            h_entry = explicit_h_carry
            problem%t0 = ktau*si%dt
            problem%t1 = (ktau + 1)*si%dt
            problem%y0 = si%z(1:4)
            problem%h0 = explicit_h_carry
            call ode_integrate_gbs(problem, solution, status)
            if (status%code /= FORTNUM_OK) then
                ierr = SYMPLECTIC_STEP_MAXITER
                explicit_h_carry = h_entry
                return
            end if
            nlast = size(solution%t)
            si%z(1:4) = solution%y(:, nlast)
            if (nlast > 1) explicit_h_carry = solution%h(nlast)
            call explicit_step_finish(h_entry, ierr)
            if (ierr /= SYMPLECTIC_STEP_OK) return
            ktau = ktau + 1
        end do
        call sync_field_to_state
    end subroutine orbit_timestep_gbs16

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Two-derivative Runge-Kutta on the canonical Hamiltonian: the Runge-Kutta-
! Nystrom transplant.
!
! RKN integrates y'' = f(t,y) and does not apply to the guiding-centre system as
! written, which is first order. Differentiating once gives zdd = F'(z)F(z) =
! G(z), which is in special second-order form because zd is itself determined by
! z, and re-synchronising the velocity as F(z) at the start of every step turns
! the scheme into a special two-derivative Runge-Kutta method (Chan & Tsai,
! Numer. Algorithms 53 (2010) 171). Its order conditions coincide with the
! Nystrom ones, which is what lets an RKN tableau be used directly here.
!
! Because the velocity is recomputed rather than propagated, only the position
! order conditions are needed, so order 4 costs two stages instead of three.
!
! Fixed step, like the symplectic schemes: fortnum's TDRK has no embedded error
! estimate yet, so the accuracy parameter is the step count rather than a
! tolerance. npoiper2 therefore sets resolution here exactly as it does for the
! symplectic integrators, which also makes the two directly comparable.
    subroutine orbit_timestep_tdrk24(ierr)
        !
        use fortnum_ode_tdrk, only: tdrk_integrate_fixed
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK
        integer, intent(out) :: ierr
        integer :: ktau, nfev_f, nfev_g
        real(dp) :: yend(4)
        type(fortnum_status_t) :: status
        real(dp) :: h_entry

        ierr = 0
        ktau = 0
        do while (ktau .lt. si%ntau)
            h_entry = explicit_h_carry
            call tdrk_integrate_fixed(f_ode_fortnum, g_ode_fortnum, &
                         ktau*si%dt, (ktau + 1)*si%dt, si%z(1:4), 1, 4, .false., yend, &
                                      nfev_f, nfev_g, status)
            if (status%code /= FORTNUM_OK) then
                ierr = SYMPLECTIC_STEP_MAXITER
                explicit_h_carry = h_entry
                return
            end if
            si%z(1:4) = yend
            call explicit_step_finish(h_entry, ierr)
            if (ierr /= SYMPLECTIC_STEP_OK) return
            ktau = ktau + 1
        end do
        call sync_field_to_state
    end subroutine orbit_timestep_tdrk24

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Two-derivative Runge-Kutta with error control: the same method as
! orbit_timestep_tdrk24, driven by a tolerance instead of a step count.
!
! The embedded estimate needs the three-stage order-4 tableau. The two-stage
! one, which orbit_timestep_tdrk24 uses, cannot carry one: with c = (0, 1/2)
! the order-3 conditions already force the order-4 weights, so the two
! solutions coincide and their difference is identically zero. One extra G
! evaluation per step buys the estimate.
!
! Both modes are kept. Fixed step is what makes TDRK comparable to the
! symplectic schemes at matched resolution; tolerance control is what makes it
! comparable to Gauss-Radau, Bulirsch-Stoer and Cash-Karp at matched accuracy.
    subroutine orbit_timestep_tdrk24a(ierr)
        !
        use fortnum_ode_tdrk, only: tdrk_integrate_adaptive
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK
        integer, intent(out) :: ierr
        integer :: ktau, nfev_f, nfev_g, naccept, nreject
        real(dp) :: yend(4), hlast
        type(fortnum_status_t) :: status
        real(dp) :: h_entry

        ierr = 0
        ktau = 0
        do while (ktau .lt. si%ntau)
            h_entry = explicit_h_carry
            call tdrk_integrate_adaptive(f_ode_fortnum, g_ode_fortnum, &
                            ktau*si%dt, (ktau + 1)*si%dt, si%z(1:4), si%rtol, si%atol, &
                explicit_h_carry, yend, hlast, nfev_f, nfev_g, naccept, nreject, status)
            if (status%code /= FORTNUM_OK) then
                ierr = SYMPLECTIC_STEP_MAXITER
                explicit_h_carry = h_entry
                return
            end if
            si%z(1:4) = yend
            explicit_h_carry = hlast
            call explicit_step_finish(h_entry, ierr)
            if (ierr /= SYMPLECTIC_STEP_OK) return
            ktau = ktau + 1
        end do
        call sync_field_to_state
    end subroutine orbit_timestep_tdrk24a

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Cash-Karp RK5(4) on the 4D canonical chart.
!
! integmode = 0 already runs a Cash-Karp, but on the 5D drift-kinetic form with
! a different right-hand side and different coordinates, so its error and its
! field-evaluation count are not commensurable with the methods above. This
! mode puts the classical adaptive workhorse on exactly the same f_ode as the
! symplectic schemes, Gauss-Radau, Bulirsch-Stoer and TDRK, which is what makes
! the work-precision comparison a comparison of methods rather than of
! formulations.
    subroutine orbit_timestep_cashkarp45(ierr)
        !
        use fortnum_ode, only: ode_integrate, ode_problem_t, ode_workspace_t, &
                               ode_solution_t
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK
        integer, intent(out) :: ierr
        integer :: ktau, nlast
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp) :: h_entry

        ierr = 0
        ktau = 0
        problem%rhs => f_ode_fortnum
        problem%rtol = si%rtol
        problem%atol = si%atol
        allocate (problem%y0(4))
        do while (ktau .lt. si%ntau)
            h_entry = explicit_h_carry
            problem%t0 = ktau*si%dt
            problem%t1 = (ktau + 1)*si%dt
            problem%y0 = si%z(1:4)
            problem%h0 = explicit_h_carry
            call ode_integrate(problem, workspace, solution, status)
            if (status%code /= FORTNUM_OK) then
                ierr = SYMPLECTIC_STEP_MAXITER
                explicit_h_carry = h_entry
                return
            end if
            nlast = solution%nsteps + 1
            si%z(1:4) = solution%y(:, nlast)
            ! solution%h is trimmed to nsteps entries, one FEWER than solution%t and
            ! solution%y, because a step size belongs to the interval between two
            ! recorded points rather than to a point. ode_integrate_radau and
            ! ode_integrate_gbs allocate nsteps+1 for both, so the same expression is
            ! in bounds there; here nlast would read one past the end and carry
            ! whatever that memory held into the next macro-step's initial step size.
            if (solution%nsteps > 0) explicit_h_carry = solution%h(solution%nsteps)
            call explicit_step_finish(h_entry, ierr)
            if (ierr /= SYMPLECTIC_STEP_OK) return
            ktau = ktau + 1
        end do
        call sync_field_to_state
    end subroutine orbit_timestep_cashkarp45

end module orbit_symplectic_quasi
