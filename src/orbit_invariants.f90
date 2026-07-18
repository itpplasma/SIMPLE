module orbit_invariants
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use field_can_mod, only: evaluate_field => evaluate, integ_to_ref
    use parmot_mod, only: rho0 => ro0
    use simple, only: tracer_t, encode_symplectic_state
    use util, only: twopi
    use velo_mod, only: isw_field_type
    use magfie_sub, only: TEST

    implicit none
    private

    integer, parameter, public :: INVARIANT_SUCCESS = 0
    integer, parameter, public :: INVARIANT_INVALID_INPUT = 1
    integer, parameter, public :: INVARIANT_NOT_AXISYMMETRIC = 2
    integer, parameter, public :: INVARIANT_NO_INTERSECTION = 3
    integer, parameter, public :: INVARIANT_ROOT_FAILURE = 4
    integer, parameter, public :: INVARIANT_CAPACITY = 5

    integer, parameter, public :: SECTION_EXTREMUM_UNKNOWN = 0
    integer, parameter, public :: SECTION_B_MINIMUM = 1
    integer, parameter, public :: SECTION_B_MAXIMUM = 2

    integer, parameter :: N_THETA_SCAN = 128
    integer, parameter :: N_RADIAL_SCAN = 257
    integer, parameter :: MAX_SECTION_BRANCHES = 16
    integer, parameter :: MAX_ROOT_ITERATIONS = 80
    real(dp), parameter :: R_MIN = 1.0e-6_dp
    real(dp), parameter :: UNIT_RADIAL_MAX = 1.0_dp - 1.0e-8_dp
    real(dp), parameter :: TEST_RADIAL_MAX = 0.5_dp - 1.0e-8_dp
    real(dp), parameter :: AXISYM_TOL = 1.0e-10_dp
    real(dp), parameter :: ROOT_TOL = 1.0e-11_dp

    type, public :: guiding_center_invariants_t
        real(dp) :: h0 = 0.0_dp
        real(dp) :: j_perp = 0.0_dp
        real(dp) :: p_phi = 0.0_dp
    end type guiding_center_invariants_t

    type, public :: invariant_start_t
        real(dp) :: state(5) = 0.0_dp
        integer :: sigma = 0
        integer :: section_branch = 0
        integer :: section_kind = SECTION_EXTREMUM_UNKNOWN
        real(dp) :: residual = 0.0_dp
    end type invariant_start_t

    type, public :: invariant_start_result_t
        integer :: status = INVARIANT_INVALID_INPUT
        type(invariant_start_t), allocatable :: starts(:)
    end type invariant_start_result_t

    public :: invariants_from_state, states_from_invariants
    public :: potato_to_simple_invariants, invariant_flux_convention

contains

    subroutine invariants_from_state(tracer, state, invariants, status)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(in) :: state(5)
        type(guiding_center_invariants_t), intent(out) :: invariants
        integer, intent(out) :: status

        type(tracer_t) :: work
        real(dp) :: axis_flux, outward_sign, psi_star_native, zsympl(4)

        invariants = guiding_center_invariants_t()
        status = INVARIANT_INVALID_INPUT
        if (.not. all(ieee_is_finite(state))) return
        if (rho0 <= 0.0_dp) return
        if (state(4) <= 0.0_dp) return
        if (abs(state(5)) > 1.0_dp) return

        work = tracer
        if (.not. field_is_axisymmetric(work)) then
            status = INVARIANT_NOT_AXISYMMETRIC
            return
        end if
        call flux_convention(work, axis_flux, outward_sign)
        call encode_symplectic_state(work%f, state, zsympl, invariants%h0, &
            invariants%j_perp, psi_star_native)
        invariants%p_phi = outward_sign*(psi_star_native - axis_flux)
        status = INVARIANT_SUCCESS
    end subroutine invariants_from_state

    pure subroutine potato_to_simple_invariants(potato, psi_axis, psi_edge, &
            v0_ratio, simple_invariants)
        !> Convert raw POTATO invariants to SIMPLE's portable convention.
        !>
        !> POTATO stores psi*=(c/q)P_phi in the equilibrium flux gauge. SIMPLE's
        !> public invariant uses an axis-zero, outward-positive flux gauge. v0_ratio
        !> is v0_POTATO/v0_SIMPLE; it is one when both runs use the same E_ref,
        !> particle mass, and charge state.
        type(guiding_center_invariants_t), intent(in) :: potato
        real(dp), intent(in) :: psi_axis, psi_edge, v0_ratio
        type(guiding_center_invariants_t), intent(out) :: simple_invariants

        real(dp) :: direction, velocity_scale2

        direction = sign(1.0_dp, psi_edge - psi_axis)
        velocity_scale2 = v0_ratio**2
        simple_invariants%h0 = velocity_scale2*potato%h0
        simple_invariants%j_perp = velocity_scale2*potato%j_perp
        simple_invariants%p_phi = direction*(potato%p_phi - psi_axis)
    end subroutine potato_to_simple_invariants

    subroutine invariant_flux_convention(tracer, axis_flux, outward_sign)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(out) :: axis_flux, outward_sign

        call flux_convention(tracer, axis_flux, outward_sign)
    end subroutine invariant_flux_convention

    subroutine states_from_invariants(tracer, invariants, result)
        type(tracer_t), intent(in) :: tracer
        type(guiding_center_invariants_t), intent(in) :: invariants
        type(invariant_start_result_t), intent(out) :: result

        type(invariant_start_t) :: candidate
        type(invariant_start_t), allocatable :: buffer(:)
        real(dp) :: residuals(N_RADIAL_SCAN, MAX_SECTION_BRANCHES)
        real(dp) :: thetas(N_RADIAL_SCAN, MAX_SECTION_BRANCHES)
        logical :: valid(N_RADIAL_SCAN, MAX_SECTION_BRANCHES)
        integer :: kinds(N_RADIAL_SCAN, MAX_SECTION_BRANCHES)
        real(dp) :: axis_flux, outward_sign, radius, radial_max, step_radius
        integer :: counts(N_RADIAL_SCAN)
        integer :: branch, branch_right, i, nfound, sigma, solve_status

        result%status = INVARIANT_INVALID_INPUT
        if (allocated(result%starts)) deallocate (result%starts)
        allocate (result%starts(0))
        if (.not. ieee_is_finite(invariants%h0)) return
        if (.not. ieee_is_finite(invariants%j_perp)) return
        if (.not. ieee_is_finite(invariants%p_phi)) return
        if (rho0 <= 0.0_dp) return
        if (invariants%h0 <= 0.0_dp) return
        if (invariants%j_perp < 0.0_dp) return
        if (.not. field_is_axisymmetric(tracer)) then
            result%status = INVARIANT_NOT_AXISYMMETRIC
            return
        end if
        call flux_convention(tracer, axis_flux, outward_sign)

        residuals = 0.0_dp
        thetas = 0.0_dp
        valid = .false.
        kinds = SECTION_EXTREMUM_UNKNOWN
        counts = 0
        radial_max = radial_upper_bound()
        step_radius = (radial_max - R_MIN)/real(N_RADIAL_SCAN - 1, dp)

        do i = 1, N_RADIAL_SCAN
            radius = R_MIN + real(i - 1, dp)*step_radius
            call section_geometry_at_radius(tracer, radius, thetas(i, :), &
                kinds(i, :), counts(i))
        end do

        allocate (buffer(2*MAX_SECTION_BRANCHES*N_RADIAL_SCAN))
        nfound = 0
        do sigma = -1, 1, 2
            do i = 1, N_RADIAL_SCAN
                radius = R_MIN + real(i - 1, dp)*step_radius
                call section_residuals_at_radius(tracer, invariants, axis_flux, &
                    outward_sign, radius, sigma, thetas(i, :), counts(i), &
                    residuals(i, :), valid(i, :))
            end do

            do i = 1, N_RADIAL_SCAN - 1
                do branch = 1, counts(i)
                    if (.not. valid(i, branch)) cycle
                    branch_right = matching_branch(thetas(i, branch), &
                        kinds(i, branch), thetas(i + 1, :), kinds(i + 1, :), &
                        valid(i + 1, :), counts(i + 1))
                    if (branch_right == 0) cycle
                    if (residuals(i, branch)*residuals(i + 1, branch_right) > &
                        0.0_dp) cycle

                    call refine_invariant_root(tracer, invariants, axis_flux, &
                        outward_sign, sigma, branch, &
                        R_MIN + real(i - 1, dp)*step_radius, &
                        R_MIN + real(i, dp)*step_radius, thetas(i, branch), &
                        thetas(i + 1, branch_right), candidate, solve_status)
                    if (solve_status /= INVARIANT_SUCCESS) cycle
                    if (is_duplicate(buffer, nfound, candidate)) cycle
                    nfound = nfound + 1
                    buffer(nfound) = candidate
                end do
            end do
        end do

        deallocate (result%starts)
        allocate (result%starts(nfound))
        if (nfound > 0) result%starts = buffer(1:nfound)
        if (nfound > 0) then
            result%status = INVARIANT_SUCCESS
        else
            result%status = INVARIANT_NO_INTERSECTION
        end if
    end subroutine states_from_invariants

    logical function field_is_axisymmetric(tracer)
        type(tracer_t), intent(in) :: tracer

        real(dp), parameter :: RADII(3) = [0.19_dp, 0.51_dp, 0.83_dp]
        real(dp), parameter :: THETAS(3) = [0.17_dp, 1.31_dp, 4.73_dp]
        real(dp), parameter :: PHI_TEST = 0.731_dp
        type(tracer_t) :: left, right
        real(dp) :: scale
        integer :: i

        field_is_axisymmetric = .false.
        do i = 1, size(RADII)
            left = tracer
            right = tracer
            call evaluate_field(left%f, RADII(i)*radial_upper_bound(), &
                THETAS(i), 0.0_dp, 0)
            call evaluate_field(right%f, RADII(i)*radial_upper_bound(), &
                THETAS(i), PHI_TEST, 0)
            scale = max(1.0_dp, abs(left%f%Bmod), abs(right%f%Bmod))
            if (abs(left%f%Bmod - right%f%Bmod) > AXISYM_TOL*scale) return
            if (abs(left%f%dBmod(3)) > AXISYM_TOL*scale) return
            if (abs(right%f%dBmod(3)) > AXISYM_TOL*scale) return
            scale = max(1.0_dp, abs(left%f%Aph), abs(right%f%Aph))
            if (abs(left%f%Aph - right%f%Aph) > AXISYM_TOL*scale) return
            if (abs(left%f%dAph(3)) > AXISYM_TOL*scale) return
            if (abs(right%f%dAph(3)) > AXISYM_TOL*scale) return
            scale = max(1.0_dp, abs(left%f%hph), abs(right%f%hph))
            if (abs(left%f%hph - right%f%hph) > AXISYM_TOL*scale) return
            if (abs(left%f%dhph(3)) > AXISYM_TOL*scale) return
            if (abs(right%f%dhph(3)) > AXISYM_TOL*scale) return
        end do
        field_is_axisymmetric = .true.
    end function field_is_axisymmetric

    subroutine section_geometry_at_radius(tracer, radius, thetas, kinds, count)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: thetas(MAX_SECTION_BRANCHES)
        integer, intent(out) :: kinds(MAX_SECTION_BRANCHES)
        integer, intent(out) :: count

        integer :: branch

        call find_section_extrema(tracer, radius, thetas, count)
        kinds = SECTION_EXTREMUM_UNKNOWN
        do branch = 1, count
            kinds(branch) = extremum_kind(tracer, radius, thetas(branch))
        end do
    end subroutine section_geometry_at_radius

    subroutine section_residuals_at_radius(tracer, invariants, axis_flux, &
            outward_sign, radius, sigma, &
            thetas, count, residuals, valid)
        type(tracer_t), intent(in) :: tracer
        type(guiding_center_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: axis_flux, outward_sign, radius
        integer, intent(in) :: sigma
        real(dp), intent(in) :: thetas(MAX_SECTION_BRANCHES)
        integer, intent(in) :: count
        real(dp), intent(out) :: residuals(MAX_SECTION_BRANCHES)
        logical, intent(out) :: valid(MAX_SECTION_BRANCHES)

        type(tracer_t) :: work
        real(dp) :: lambda2
        integer :: branch

        residuals = 0.0_dp
        valid = .false.
        do branch = 1, count
            work = tracer
            call evaluate_field(work%f, radius, thetas(branch), 0.0_dp, 0)
            if (work%f%Bmod <= 0.0_dp) cycle
            lambda2 = 1.0_dp - invariants%j_perp*work%f%Bmod/invariants%h0
            if (lambda2 < 0.0_dp) cycle
            residuals(branch) = outward_sign*(work%f%Aph - axis_flux + &
                rho0*sqrt(invariants%h0)*real(sigma, dp)* &
                sqrt(max(0.0_dp, lambda2))*work%f%hph) - invariants%p_phi
            valid(branch) = .true.
        end do
    end subroutine section_residuals_at_radius

    integer function matching_branch(theta, kind, candidates, candidate_kinds, &
            candidate_valid, count)
        real(dp), intent(in) :: theta
        integer, intent(in) :: kind, count
        real(dp), intent(in) :: candidates(MAX_SECTION_BRANCHES)
        integer, intent(in) :: candidate_kinds(MAX_SECTION_BRANCHES)
        logical, intent(in) :: candidate_valid(MAX_SECTION_BRANCHES)

        real(dp) :: distance, best_distance
        integer :: i

        matching_branch = 0
        best_distance = huge(1.0_dp)
        do i = 1, count
            if (.not. candidate_valid(i)) cycle
            if (candidate_kinds(i) /= kind) cycle
            distance = angular_distance(candidates(i), theta)
            if (distance >= best_distance) cycle
            matching_branch = i
            best_distance = distance
        end do
    end function matching_branch

    subroutine find_section_extrema(tracer, radius, roots, count)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: roots(MAX_SECTION_BRANCHES)
        integer, intent(out) :: count

        real(dp) :: dtheta, left, right, fleft, fright, root
        integer :: i

        roots = 0.0_dp
        count = 0
        dtheta = twopi/real(N_THETA_SCAN, dp)
        left = 0.5_dp*dtheta
        fleft = theta_derivative(tracer, radius, left)
        do i = 0, N_THETA_SCAN - 1
            left = (real(i, dp) + 0.5_dp)*dtheta
            right = left + dtheta
            fright = theta_derivative(tracer, radius, right)
            if (fleft*fright <= 0.0_dp) then
                call bisect_theta_extremum(tracer, radius, left, right, root)
                root = modulo(root, twopi)
                if (.not. root_is_duplicate(roots, count, root)) then
                    if (count == MAX_SECTION_BRANCHES) exit
                    count = count + 1
                    roots(count) = root
                end if
            end if
            fleft = fright
        end do
        call sort_roots(roots, count)
    end subroutine find_section_extrema

    subroutine bisect_theta_extremum(tracer, radius, left_in, right_in, root)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(in) :: radius, left_in, right_in
        real(dp), intent(out) :: root

        real(dp) :: left, right, middle, fleft, fmiddle
        integer :: iteration

        left = left_in
        right = right_in
        fleft = theta_derivative(tracer, radius, left)
        do iteration = 1, MAX_ROOT_ITERATIONS
            middle = 0.5_dp*(left + right)
            fmiddle = theta_derivative(tracer, radius, middle)
            if (abs(fmiddle) <= ROOT_TOL) exit
            if (right - left <= ROOT_TOL) exit
            if (fleft*fmiddle <= 0.0_dp) then
                right = middle
            else
                left = middle
                fleft = fmiddle
            end if
        end do
        root = 0.5_dp*(left + right)
    end subroutine bisect_theta_extremum

    real(dp) function theta_derivative(tracer, radius, theta)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(in) :: radius, theta

        type(tracer_t) :: work

        work = tracer
        call evaluate_field(work%f, radius, theta, 0.0_dp, 0)
        theta_derivative = work%f%dBmod(2)
    end function theta_derivative

    subroutine refine_invariant_root(tracer, invariants, axis_flux, outward_sign, &
            sigma, branch, &
            left_in, right_in, theta_left_in, &
            theta_right_in, solution, status)
        type(tracer_t), intent(in) :: tracer
        type(guiding_center_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: axis_flux, outward_sign
        integer, intent(in) :: sigma, branch
        real(dp), intent(in) :: left_in, right_in, theta_left_in, theta_right_in
        type(invariant_start_t), intent(out) :: solution
        integer, intent(out) :: status

        real(dp) :: left, right, middle, fleft, fmiddle
        real(dp) :: theta_left, theta_right, theta_middle, theta_guess
        real(dp) :: z_integ(5), z_ref(3)
        integer :: iteration, kind
        logical :: ok

        solution = invariant_start_t()
        status = INVARIANT_ROOT_FAILURE
        left = left_in
        right = right_in
        theta_left = theta_left_in
        theta_right = theta_right_in
        theta_guess = theta_left
        call residual_near_theta(tracer, invariants, axis_flux, outward_sign, &
            left, theta_guess, sigma, fleft, theta_left, kind, ok)
        if (.not. ok) return

        do iteration = 1, MAX_ROOT_ITERATIONS
            middle = 0.5_dp*(left + right)
            theta_middle = angular_midpoint(theta_left, theta_right)
            theta_guess = theta_middle
            call residual_near_theta(tracer, invariants, axis_flux, outward_sign, &
                middle, theta_guess, sigma, fmiddle, theta_middle, kind, ok)
            if (.not. ok) return
            if (abs(fmiddle) <= ROOT_TOL*max(1.0_dp, abs(invariants%p_phi))) exit
            if (right - left <= ROOT_TOL) exit
            if (fleft*fmiddle <= 0.0_dp) then
                right = middle
                theta_right = theta_middle
            else
                left = middle
                theta_left = theta_middle
                fleft = fmiddle
            end if
        end do

        z_integ(1:3) = [middle, theta_middle, 0.0_dp]
        z_integ(4) = sqrt(invariants%h0)
        call pitch_at_position(tracer, invariants, z_integ(1:3), sigma, &
            z_integ(5), ok)
        if (.not. ok) return
        call integ_to_ref(z_integ(1:3), z_ref)
        solution%state(1:3) = z_ref
        solution%state(4:5) = z_integ(4:5)
        solution%sigma = sigma
        solution%section_branch = branch
        solution%section_kind = kind
        solution%residual = fmiddle
        status = INVARIANT_SUCCESS
    end subroutine refine_invariant_root

    subroutine residual_near_theta(tracer, invariants, axis_flux, outward_sign, &
            radius, theta_guess, sigma, &
            residual, theta, kind, ok)
        type(tracer_t), intent(in) :: tracer
        type(guiding_center_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: axis_flux, outward_sign, radius, theta_guess
        integer, intent(in) :: sigma
        real(dp), intent(out) :: residual, theta
        integer, intent(out) :: kind
        logical, intent(out) :: ok

        type(tracer_t) :: work
        real(dp) :: roots(MAX_SECTION_BRANCHES), distance, best_distance
        real(dp) :: lambda2
        integer :: count, i, best

        residual = 0.0_dp
        theta = 0.0_dp
        kind = SECTION_EXTREMUM_UNKNOWN
        ok = .false.
        call find_section_extrema(tracer, radius, roots, count)
        if (count == 0) return
        best = 1
        best_distance = angular_distance(roots(1), theta_guess)
        do i = 2, count
            distance = angular_distance(roots(i), theta_guess)
            if (distance < best_distance) then
                best = i
                best_distance = distance
            end if
        end do
        theta = roots(best)
        work = tracer
        call evaluate_field(work%f, radius, theta, 0.0_dp, 0)
        lambda2 = 1.0_dp - invariants%j_perp*work%f%Bmod/invariants%h0
        if (lambda2 < 0.0_dp) return
        residual = outward_sign*(work%f%Aph - axis_flux + &
            rho0*sqrt(invariants%h0)*real(sigma, dp)*sqrt(lambda2)*work%f%hph) - &
            invariants%p_phi
        kind = extremum_kind(tracer, radius, theta)
        ok = .true.
    end subroutine residual_near_theta

    integer function extremum_kind(tracer, radius, theta)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(in) :: radius, theta

        real(dp), parameter :: ANGLE_STEP = 1.0e-5_dp
        real(dp) :: slope

        slope = theta_derivative(tracer, radius, theta + ANGLE_STEP) - &
            theta_derivative(tracer, radius, theta - ANGLE_STEP)
        extremum_kind = SECTION_EXTREMUM_UNKNOWN
        if (slope > 0.0_dp) extremum_kind = SECTION_B_MINIMUM
        if (slope < 0.0_dp) extremum_kind = SECTION_B_MAXIMUM
    end function extremum_kind

    subroutine pitch_at_position(tracer, invariants, position, sigma, pitch, ok)
        type(tracer_t), intent(in) :: tracer
        type(guiding_center_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: position(3)
        integer, intent(in) :: sigma
        real(dp), intent(out) :: pitch
        logical, intent(out) :: ok

        type(tracer_t) :: work
        real(dp) :: lambda2

        work = tracer
        call evaluate_field(work%f, position(1), position(2), position(3), 0)
        lambda2 = 1.0_dp - invariants%j_perp*work%f%Bmod/invariants%h0
        ok = lambda2 >= 0.0_dp
        if (ok) then
            pitch = real(sigma, dp)*sqrt(max(0.0_dp, lambda2))
        else
            pitch = 0.0_dp
        end if
    end subroutine pitch_at_position

    subroutine flux_convention(tracer, axis_flux, outward_sign)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(out) :: axis_flux, outward_sign

        type(tracer_t) :: axis_field, edge_field
        real(dp) :: delta_flux

        axis_field = tracer
        edge_field = tracer
        call evaluate_field(axis_field%f, R_MIN, 0.0_dp, 0.0_dp, 0)
        call evaluate_field(edge_field%f, radial_upper_bound(), 0.0_dp, 0.0_dp, 0)
        axis_flux = axis_field%f%Aph
        delta_flux = edge_field%f%Aph - axis_flux
        outward_sign = sign(1.0_dp, delta_flux)
    end subroutine flux_convention

    pure real(dp) function radial_upper_bound()
        if (isw_field_type == TEST) then
            radial_upper_bound = TEST_RADIAL_MAX
        else
            radial_upper_bound = UNIT_RADIAL_MAX
        end if
    end function radial_upper_bound

    logical function root_is_duplicate(roots, count, root)
        real(dp), intent(in) :: roots(MAX_SECTION_BRANCHES), root
        integer, intent(in) :: count

        integer :: i

        root_is_duplicate = .false.
        do i = 1, count
            if (angular_distance(roots(i), root) < 1.0e-8_dp) then
                root_is_duplicate = .true.
                return
            end if
        end do
    end function root_is_duplicate

    logical function is_duplicate(buffer, count, candidate)
        type(invariant_start_t), intent(in) :: buffer(:), candidate
        integer, intent(in) :: count

        integer :: i

        is_duplicate = .false.
        do i = 1, count
            if (buffer(i)%sigma /= candidate%sigma) cycle
            if (abs(buffer(i)%state(1) - candidate%state(1)) > 1.0e-7_dp) cycle
            if (angular_distance(buffer(i)%state(2), candidate%state(2)) > &
                1.0e-7_dp) cycle
            is_duplicate = .true.
            return
        end do
    end function is_duplicate

    subroutine sort_roots(roots, count)
        real(dp), intent(inout) :: roots(MAX_SECTION_BRANCHES)
        integer, intent(in) :: count

        real(dp) :: value
        integer :: i, j

        do i = 2, count
            value = roots(i)
            j = i - 1
            do while (j >= 1)
                if (roots(j) <= value) exit
                roots(j + 1) = roots(j)
                j = j - 1
            end do
            roots(j + 1) = value
        end do
    end subroutine sort_roots

    pure real(dp) function angular_distance(left, right)
        real(dp), intent(in) :: left, right

        angular_distance = abs(modulo(left - right + 0.5_dp*twopi, twopi) - &
            0.5_dp*twopi)
    end function angular_distance

    pure real(dp) function angular_midpoint(left, right)
        real(dp), intent(in) :: left, right

        angular_midpoint = modulo(left + 0.5_dp* &
            (modulo(right - left + 0.5_dp*twopi, twopi) - 0.5_dp*twopi), twopi)
    end function angular_midpoint

end module orbit_invariants
