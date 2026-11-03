module orbit_symplectic_base
use field_can_mod, only: eval_field => evaluate, field_can_t, get_val, get_derivatives, &
                             get_derivatives2

    implicit none

! Define real(dp) kind parameter
    integer, parameter :: dp = kind(1.0d0)

    logical, parameter :: extrap_field = .True.  ! do extrapolation after final iteration

    ! Upper radial bound for Newton trials: normalized-flux edge s = 1 by
    ! default; SPECTRE uses Mvol + 1, with physical loss localized at Mvol.
    real(dp) :: sympl_rmax = 1.0d0
    !$acc declare copyin(sympl_rmax)

    ! Integration methods
    integer, parameter :: RK45 = 0, EXPL_IMPL_EULER = 1, IMPL_EXPL_EULER = 2, &
             MIDPOINT = 3, GAUSS1 = 4, GAUSS2 = 5, GAUSS3 = 6, GAUSS4 = 7, &
             CASH_KARP = 8, DORMAND_PRINCE = 9, LOBATTO3 = 15

    ! High-order explicit methods on the 4D canonical chart, for benchmarking
    ! against the symplectic schemes above. They integrate the same canonical
    ! Hamiltonian through the same f_ode right-hand side, so a comparison per
    ! field evaluation is like for like -- unlike integmode = 0, which runs the
    ! 5D drift-kinetic form with a different right-hand side and coordinates.
    !
    ! RADAU15 is the first-order-system formulation of IAS15 (Rein & Spiegel
    ! 2015), the strongest published claim that a non-symplectic method beats
    ! symplectic ones on long integrations. GBS16 is Gragg-Bulirsch-Stoer
    ! extrapolation, the classical high-accuracy explicit method of celestial
    ! mechanics.
    ! TDRK24 is the Runge-Kutta-Nystrom transplant: an RKN tableau applied to
    ! the first-order guiding-centre system through zdd = F'(z)F(z), which is a
    ! special two-derivative Runge-Kutta method (Chan & Tsai 2010). Fixed step,
    ! so npoiper2 sets its resolution as it does for the symplectic schemes;
    ! TDRK24A is the same method under tolerance control, which is what makes it
    ! comparable to the adaptive methods rather than to the symplectic ones.
    !
    ! CASHKARP45 is the classical adaptive workhorse on this same chart.
    ! integmode = 0 also runs a Cash-Karp, but on the 5D drift-kinetic form with
    ! a different right-hand side and coordinates, so its accuracy and cost are
    ! not commensurable with anything here. CASH_KARP = 8 and DORMAND_PRINCE = 9
    ! above are a different thing again: they select the RK5(4) kernels of the
    ! GPU tracing paths (simple_gpu, simple_cuda_native). CASHKARP45 is the CPU
    ! method on the canonical chart and does not dispatch to them.
    integer, parameter :: RADAU15 = 20, GBS16 = 21, CASHKARP45 = 22, &
                          TDRK24 = 24, TDRK24A = 25
    integer, parameter :: SYMPLECTIC_STEP_OK = 0
    integer, parameter :: SYMPLECTIC_STEP_OUTSIDE_DOMAIN = 1
    integer, parameter :: SYMPLECTIC_STEP_MAXITER = 2
    integer, parameter :: SYMPLECTIC_STEP_LINEAR_SOLVE = 3
    integer, parameter :: SYMPLECTIC_STEP_BOUNDARY = 4
    integer, parameter :: SYMPLECTIC_STEP_EVENT_NOT_CONVERGED = 5
    integer, parameter :: SYMPLECTIC_STEP_BOUNDARY_LIMITED = 6
    ! Newton starting guess: z_n, or Lagrange extrapolation of the last two or
    ! three accepted states.
    integer, parameter :: PREDICT_NONE = 0, PREDICT_LINEAR = 1, &
                          PREDICT_QUADRATIC = 2
    integer :: sympl_predictor = PREDICT_QUADRATIC
    ! Newton stopping rule. STOP_CORRECTION accepts x_k once the correction
    ! x_k - x_(k-1) is below tolerance. STOP_RATE also accepts x_k when the
    ! contraction rate theta = |dx_k|/|dx_(k-1)| bounds the remaining error,
    ! theta/(1 - theta)*|dx_k| <= kappa*tolerance, which saves the iteration
    ! that would only confirm convergence.
    integer, parameter :: STOP_CORRECTION = 0, STOP_RATE = 1
    integer :: sympl_stop_rule = STOP_RATE
    real(dp) :: sympl_stop_kappa = 1d0
    !$acc declare copyin(sympl_predictor, sympl_stop_rule, sympl_stop_kappa)

    logical :: symplectic_newton_warning_mode = .true.
    real(dp), parameter :: symplectic_newton_warning_factor = 100.0_dp
    real(dp), parameter :: symplectic_spectre_newton_warning_factor = 10.0_dp
    real(dp) :: boundary_event_fraction_tolerance = -1d0
    real(dp) :: boundary_event_radial_tolerance = -1d0

    type :: symplectic_integrator_t
        real(dp) :: atol
        real(dp) :: rtol
        real(dp) :: warning_factor = symplectic_newton_warning_factor
        ! Method from orbit_sympl_init; composition stages dispatch on it
        ! because their Euler step and its adjoint alternate.
        integer :: mode = -1

        ! Current phase-space coordinates z and old pth
        real(dp), dimension(4) :: z  ! z = (r, th, ph, pphi)
        real(dp) :: pthold

        ! Timestep and variables from z0
        integer :: ntau
        real(dp) :: dt
        real(dp) :: pabs
        real(dp) :: last_step_fraction = 1d0
        real(dp) :: last_event_radial_residual = 0d0
        real(dp) :: last_event_fraction_width = 0d0

        ! Accepted-step history for the Newton starting guess: the last two
        ! increments z_n - z_(n-1), z_(n-1) - z_(n-2), and the state and step
        ! they belong to. The history is valid only while z and dt still match,
        ! so any external change of the state (chart switch, interface map,
        ! retry with a smaller step) discards it.
        integer :: nhist = 0
        real(dp) :: hist_dt = 0d0
        real(dp) :: hist_z(4) = 0d0
        real(dp) :: hist_dz(4, 2) = 0d0
    end type symplectic_integrator_t

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    ! Composition method with 2s internal stages according to Hairer, 2002 V.3.1
    !
    integer, parameter :: S_MAX = 32
    type :: multistage_integrator_t
        integer :: s
        real(dp) :: alpha(S_MAX), beta(S_MAX)
        type(symplectic_integrator_t) stages(2*S_MAX)
    end type multistage_integrator_t

    abstract interface
        subroutine orbit_timestep_sympl_i(si, f, ierr)
            import :: symplectic_integrator_t, field_can_t
            type(symplectic_integrator_t), intent(inout) :: si
            type(field_can_t), intent(inout) :: f
            integer, intent(out) :: ierr
        end subroutine orbit_timestep_sympl_i
    end interface

    abstract interface
        subroutine orbit_timestep_quasi_i(ierr)
            integer, intent(out) :: ierr
        end subroutine orbit_timestep_quasi_i
    end interface

contains

    pure logical function history_valid(si)
        !$acc routine seq
        ! The history continues the current state only if the last recorded step
        ! ended exactly here with the same step size.
        type(symplectic_integrator_t), intent(in) :: si

        history_valid = si%nhist > 0
        if (.not. history_valid) return
        history_valid = si%hist_dt == si%dt .and. all(si%hist_z == si%z)
    end function history_valid

    pure function predict_state(si, c) result(z)
        !$acc routine seq
        ! Lagrange extrapolation of the accepted states z_(n-2), z_(n-1), z_n to the
        ! fraction c of the coming step. Falls back to z_n without a valid history
        ! or when the guess leaves the radial domain.
        type(symplectic_integrator_t), intent(in) :: si
        real(dp), intent(in) :: c
        real(dp) :: z(4)

        z = si%z
        if (sympl_predictor == PREDICT_NONE) return
        if (.not. history_valid(si)) return
        z = si%z + c*si%hist_dz(:, 1)
        if (sympl_predictor /= PREDICT_LINEAR .and. si%nhist > 1) &
            z = z + 0.5d0*c*(c + 1d0)*(si%hist_dz(:, 1) - si%hist_dz(:, 2))
        if (.not. (z(1) > 0d0 .and. z(1) < sympl_rmax)) z = si%z
    end function predict_state

    pure subroutine record_history(si, z_start, keep)
        !$acc routine seq
        ! Append the accepted step z_start -> si%z to the history, or restart it
        ! when keep is false (chart switch).
        type(symplectic_integrator_t), intent(inout) :: si
        real(dp), intent(in) :: z_start(4)
        logical, intent(in) :: keep

        if (.not. keep) then
            si%nhist = 0
            return
        end if
        if (si%nhist > 0 .and. si%hist_dt == si%dt .and. all(si%hist_z == z_start)) then
            si%hist_dz(:, 2) = si%hist_dz(:, 1)
            si%nhist = min(si%nhist + 1, 2)
        else
            si%nhist = 1
        end if
        si%hist_dz(:, 1) = si%z - z_start
        si%hist_z = si%z
        si%hist_dt = si%dt
    end subroutine record_history

    pure logical function rate_converged(dnorm, dnorm_prev)
        !$acc routine seq
        ! Error bound for x_k from the contraction rate of the last two Newton
        ! corrections, both scaled by their tolerances (Hairer & Wanner II, IV.8).
        real(dp), intent(in) :: dnorm, dnorm_prev
        real(dp) :: theta

        rate_converged = .false.
        if (sympl_stop_rule == STOP_CORRECTION) return
        if (.not. (dnorm_prev > 0d0)) return
        theta = dnorm/dnorm_prev
        if (.not. (theta < 1d0)) return
        rate_converged = theta/(1d0 - theta)*dnorm <= sympl_stop_kappa
    end function rate_converged

    subroutine coeff_rk_gauss(n, a, b, c)
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(inout) :: a(n, n), b(n), c(n)

        if (n == 1) then
            a(1, 1) = 0.5d0
            b(1) = 1.0d0
            c(1) = 0.5d0
        elseif (n == 2) then
            a(1, 1) = 0.25d0
            a(1, 2) = -0.038675134594812d0
            a(2, 1) = 0.538675134594812d0
            a(2, 2) = 0.25d0

            b(1) = 0.5d0
            b(2) = 0.5d0

            c(1) = 0.211324865405187d0
            c(2) = 0.788675134594812d0
        elseif (n == 3) then
            a(1, 1) = 0.1388888888888889d0
            a(1, 2) = -0.03597666752493894d0
            a(1, 3) = 0.009789444015308318d0
            a(2, 1) = 0.3002631949808646d0
            a(2, 2) = 0.2222222222222222d0
            a(2, 3) = -0.022485417203086805d0
            a(3, 1) = 0.26798833376246944d0
            a(3, 2) = 0.48042111196938336d0
            a(3, 3) = 0.1388888888888889d0

            b(1) = 0.2777777777777778d0
            b(2) = 0.4444444444444444d0
            b(3) = 0.2777777777777778d0

            c(1) = 0.1127016653792583d0
            c(2) = 0.5d0
            c(3) = 0.8872983346207417d0
        elseif (n == 4) then  ! with help of coefficients from GeometricIntegrators.jl of Michael Kraus
            a(1, 1) = 0.086963711284363462428182d0
            a(1, 2) = -0.026604180084998794303397d0
            a(1, 3) = 0.012627462689404725035280d0
            a(1, 4) = -0.003555149685795683332096d0

            a(2, 1) = 0.188118117499868064967927d0
            a(2, 2) = 0.163036288715636523694030d0
            a(2, 3) = -0.027880428602470894855481d0
            a(2, 4) = 0.006735500594538155853808d0

            a(3, 1) = 0.167191921974188778543535d0
            a(3, 2) = 0.353953006033743966529670d0
            a(3, 3) = 0.163036288715636523694030d0
            a(3, 4) = -0.014190694931141143581010d0

            a(4, 1) = 0.177482572254522602550608d0
            a(4, 2) = 0.313445114741868369190314d0
            a(4, 3) = 0.352676757516271865977586d0
            a(4, 4) = 0.086963711284363462428182d0

            b(1) = 0.173927422568726924856364d0
            b(2) = 0.326072577431273047388061d0
            b(3) = 0.326072577431273047388061d0
            b(4) = 0.173927422568726924856364d0

            c(1) = 0.069431844202973713731097d0
            c(2) = 0.330009478207571871344328d0
            c(3) = 0.669990521792428128655672d0
            c(4) = 0.930568155797026341780054d0
        else
            ! not implemented
            a = 0d0
            b = 0d0
            c = 0d0
        end if
    end subroutine coeff_rk_gauss

    subroutine coeff_rk_lobatto(n, a, ahat, b, c)
        integer, intent(in) :: n
        real(dp), intent(inout) :: a(n, n), ahat(n, n), b(n), c(n)

        if (n == 3) then
            a(1, 1) = 0d0
            a(1, 2) = 0d0
            a(1, 3) = 0d0

            a(2, 1) = 0.20833333333333334d0
            a(2, 2) = 0.33333333333333333d0
            a(2, 3) = -0.041666666666666664d0

            a(3, 1) = 0.16666666666666667d0
            a(3, 2) = 0.66666666666666667d0
            a(3, 3) = 0.16666666666666667d0

            ahat(1, 1) = 0.16666666666666667d0
            ahat(1, 2) = -0.16666666666666667d0
            ahat(1, 3) = 0d0

            ahat(2, 1) = 0.16666666666666667d0
            ahat(2, 2) = 0.33333333333333333d0
            ahat(2, 3) = 0d0

            ahat(3, 1) = 0.16666666666666667d0
            ahat(3, 2) = 0.83333333333333333d0
            ahat(3, 3) = 0d0

            b(1) = 0.16666666666666667d0
            b(2) = 0.66666666666666667d0
            b(3) = 0.16666666666666667d0

            c(1) = 0d0
            c(2) = 0.5d0
            c(3) = 1.0d0

        else
            ! not implemented
            a = 0d0
            b = 0d0
            c = 0d0
        end if
    end subroutine coeff_rk_lobatto

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    ! Lobatto (IIIA)-(IIIB) Runge-Kutta method with s internal stages (n=4*s variables)
    !
    subroutine f_rk_lobatto(si, fs, s, x, fvec, jactype)
        !
        type(symplectic_integrator_t), intent(inout) :: si
        integer, intent(in) :: s
        type(field_can_t), intent(inout) :: fs(:)
        real(dp), intent(in) :: x(4*s)  ! = (rend, thend, phend, pphend)
        real(dp), intent(out) :: fvec(4*s)
        integer, intent(in) :: jactype  ! 0 = no second derivatives, 2 = second derivatives

        real(dp) :: a(s, s), ahat(s, s), b(s), c(s), Hprime(s)
        integer :: k, l  ! counters

        call coeff_rk_lobatto(s, a, ahat, b, c)

        call eval_field(fs(1), x(1), si%z(2), si%z(3), jactype)
        call stage_derivatives(fs(1), x(2))

        do k = 2, s
         call eval_field(fs(k), x(4*k - 3 - 2), x(4*k - 2 - 2), x(4*k - 1 - 2), jactype)
            call stage_derivatives(fs(k), x(4*k - 2))
        end do

        Hprime = fs%dH(1)/fs%dpth(1)

        fvec(1) = fs(1)%pth - si%pthold
        fvec(2) = x(2) - si%z(4)

        do l = 1, s
fvec(1) = fvec(1) + si%dt*ahat(1, l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
fvec(2) = fvec(2) + si%dt*ahat(1, l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
        end do

        do k = 2, s
            fvec(4*k - 3 - 2) = fs(k)%pth - si%pthold
            fvec(4*k - 2 - 2) = x(4*k - 2 - 2) - si%z(2)
            fvec(4*k - 1 - 2) = x(4*k - 1 - 2) - si%z(3)
            fvec(4*k - 2) = x(4*k - 2) - si%z(4)
        end do

        do l = 1, s
            do k = 2, s
      fvec(4*k-3-2) = fvec(4*k-3-2) + si%dt*ahat(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
      fvec(4*k-2-2) = fvec(4*k-2-2) - si%dt*a(k,l)*Hprime(l)                                         ! thdot
      fvec(4*k-1-2) = fvec(4*k-1-2) - si%dt*a(k,l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph     ! phdot
      fvec(4*k-2)   = fvec(4*k-2)   + si%dt*ahat(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
            end do
        end do

    contains

        subroutine stage_derivatives(fk, pphi)
            ! The Newton Jacobian needs the second derivatives of H and p_theta.
            type(field_can_t), intent(inout) :: fk
            real(dp), intent(in) :: pphi

            if (jactype == 2) then
                call get_derivatives2(fk, pphi)
            else
                call get_derivatives(fk, pphi)
            end if
        end subroutine stage_derivatives

    end subroutine f_rk_lobatto

end module orbit_symplectic_base
