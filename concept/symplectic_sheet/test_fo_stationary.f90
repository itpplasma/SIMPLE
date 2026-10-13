program test_fo_stationary
    !> Long-run stationarity of the microcanonical ensemble across the sheet.
    !>
    !> On the shell H = E at fixed p_phi the invariant measure
    !> dR dZ dp_R dp_Z delta(H - E) projects to the uniform density 2 pi dR dZ
    !> on the allowed region |p_phi - psi(r)| < R v of the poloidal plane. An
    !> ensemble drawn from it stays distributed so at every physical time under
    !> an exactly measure-preserving flow, so the fraction of markers inside the
    !> sheet must stay at the exact area ratio. A drift of that fraction is
    !> spurious accumulation across the sheet.
    !>
    !> Schemes: plain midpoint on the sharp sheet, and the Poincare-transformed
    !> midpoint on a smooth layer of width WIDTH.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use sheet_field, only: field_params_t, region_of, hamiltonian, minor_radius
    use fo_sheet, only: sheet_step, poincare_step, step_stats_t, &
        poincare_params_t, SCHEME_NAIVE
    use sheet_setup, only: shell_state
    implicit none

    integer, parameter :: SCHEME_POINCARE = 3
    integer, parameter :: NMARKER = 32768, NREC = 20, NGRID = 4000
    real(dp), parameter :: H_STEP = 0.1_dp, SPEED = 0.05_dp
    real(dp), parameter :: T_END = 1.0e4_dp, WIDTH = 2.5e-3_dp
    real(dp), parameter :: TWOPI = 6.283185307179586_dp, R_BOX = 1.0_dp

    type(field_params_t) :: par
    type(poincare_params_t) :: pp

    pp%g_min = 0.05_dp
    pp%length = 2.0_dp*WIDTH
    call run(SCHEME_NAIVE)
    par%width = WIDTH
    call run(SCHEME_POINCARE)

contains

    subroutine run(scheme)
        integer, intent(in) :: scheme

        real(dp) :: z(4), z_new(4), e0, max_de, frac(0:NREC), sigma, t, dt
        real(dp) :: t_rec, frac_exact
        integer :: im, region, inside(0:NREC), irec, failures
        integer(int64) :: nsteps
        type(step_stats_t) :: stats
        logical :: ok

        frac_exact = exact_inside_fraction()
        inside = 0
        max_de = 0.0_dp
        failures = 0
        nsteps = 0
        !$omp parallel do schedule(dynamic, 16) default(shared) &
        !$omp private(im, z, z_new, region, e0, stats, irec, t, dt, t_rec, ok) &
        !$omp reduction(+:inside, failures, nsteps) reduction(max:max_de)
        do im = 1, NMARKER
            call sample_marker(im, z)
            region = region_of(par, z(1:2))
            e0 = hamiltonian(par, region, z)
            stats = step_stats_t()
            if (region == 1) inside(0) = inside(0) + 1
            t = 0.0_dp
            irec = 1
            t_rec = T_END/real(NREC, dp)
            do while (irec <= NREC)
                if (scheme == SCHEME_POINCARE) then
                    call poincare_step(par, pp, z, e0, H_STEP, z_new, dt, ok)
                    if (.not. ok) failures = failures + 1
                    region = region_of(par, z_new(1:2))
                else
                    call sheet_step(par, scheme, region, z, H_STEP, z_new, stats)
                    dt = H_STEP
                end if
                z = z_new
                t = t + dt
                nsteps = nsteps + 1
                do while (irec <= NREC .and. t >= t_rec)
                    if (region == 1) inside(irec) = inside(irec) + 1
                    max_de = max(max_de, &
                        abs(hamiltonian(par, region, z)/e0 - 1.0_dp))
                    irec = irec + 1
                    t_rec = T_END*real(irec, dp)/real(NREC, dp)
                end do
            end do
            failures = failures + stats%failures
        end do
        !$omp end parallel do

        frac = real(inside, dp)/real(NMARKER, dp)
        sigma = sqrt(frac_exact*(1.0_dp - frac_exact)/real(NMARKER, dp))
        print '(/,A,I0,A,F8.5,A,I0,A,F6.2)', 'scheme ', scheme, &
            '  exact inside fraction ', frac_exact, '  solver failures ', &
            failures, '  steps per unit time ', &
            real(nsteps, dp)/(real(NMARKER, dp)*T_END)
        print '(A,ES10.2)', '  max |dH/H| ', max_de
        print '(A,F9.5,A,F9.5,A,F8.5)', '  inside fraction: start ', frac(0), &
            '  mean of last half ', sum(frac(NREC/2:NREC))/real(NREC/2 + 1, dp), &
            '  one-sample sigma ', sigma
        print '(A,*(F7.2))', '  (frac - start)/sigma per record:', &
            (frac(1:NREC) - frac(0))/sigma
    end subroutine run

    subroutine sample_marker(im, z)
        !> Rejection sampling of (R, Z) uniform on the allowed region and a
        !> uniform poloidal velocity direction, from a per-marker stream.
        integer, intent(in) :: im
        real(dp), intent(out) :: z(4)

        integer(int64) :: state
        integer :: k
        real(dp) :: x, y
        logical :: ok

        state = modulo(int(im, int64)*7919_int64 + 12345_int64, &
                       2147483646_int64) + 1
        do k = 1, 10
            x = uniform(state)
        end do
        do
            x = (2.0_dp*uniform(state) - 1.0_dp)*R_BOX
            y = (2.0_dp*uniform(state) - 1.0_dp)*R_BOX
            call shell_state(par, SPEED, hypot(x, y), atan2(y, x), &
                TWOPI*uniform(state), z, ok)
            if (ok) exit
        end do
    end subroutine sample_marker

    real(dp) function exact_inside_fraction() result(frac)
        !> Area of the allowed region inside the sheet over its total area,
        !> by midpoint quadrature on a fine grid over the bounding box.
        real(dp) :: x, y, z(4), d, area_in, area
        integer :: i, j
        logical :: ok

        d = 2.0_dp*R_BOX/real(NGRID, dp)
        area_in = 0.0_dp
        area = 0.0_dp
        do j = 1, NGRID
            y = -R_BOX + (real(j, dp) - 0.5_dp)*d
            do i = 1, NGRID
                x = -R_BOX + (real(i, dp) - 0.5_dp)*d
                call shell_state(par, SPEED, hypot(x, y), atan2(y, x), 0.0_dp, &
                    z, ok)
                if (.not. ok) cycle
                area = area + 1.0_dp
                if (minor_radius(par, z(1:2)) < par%r_sheet) &
                    area_in = area_in + 1.0_dp
            end do
        end do
        frac = area_in/area
    end function exact_inside_fraction

    real(dp) function uniform(state) result(u)
        !> Park-Miller minimal standard generator (48271, 2^31 - 1).
        integer(int64), intent(inout) :: state

        state = modulo(48271_int64*state, 2147483647_int64)
        u = real(state, dp)/2147483647.0_dp
    end function uniform

end program test_fo_stationary
