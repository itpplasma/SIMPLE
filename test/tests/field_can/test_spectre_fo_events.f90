program test_spectre_fo_events
    ! Near-tangent cross-and-return: both endpoints belong to volume 1,
    ! but intermediate points belong to volume 2. Axisymmetry supplies an
    ! independent canonical-momentum oracle, without a GC exit projection.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_sub, only: spectre_field, set_magfie_spectre_field, &
        M_TO_CM, TESLA_TO_GAUSS
    use libneo_coordinates, only: spectre_coordinate_system_t, &
        make_spectre_coordinate_system
    use field_spectre, only: spectre_field_t, create_spectre_field
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_configure, &
        spectre_fo_advance_until_exit, spectre_fo_canonical_pzeta, &
        spectre_fo_enter, SPECTRE_FO_OK, SPECTRE_FO_LOSS, &
        SPECTRE_FO_GYROPHASE_FAIL
    implicit none

    type(spectre_field_t) :: sf
    type(spectre_fo_state_t) :: state
    character(len=1024) :: h5file
    real(dp), parameter :: RO = 1.0e4_dp
    real(dp) :: u(3), x(3), b(3), bmod, n(3), tangent(3), dt, used, y(5)
    real(dp) :: p0, p1, defect, acov(3), hcov(3), scale, xref(3), vref(3)
    integer :: status, owner, j, visits, cross_entries(2)
    logical :: exited

    call get_command_argument(1, h5file)
    call create_spectre_field(sf, trim(h5file), status)
    if (status /= 0) error stop 'fixture initialization failed'
    call set_magfie_spectre_field(sf)
    call spectre_fo_configure( 1.0e-6_dp)

    u = [1.0_dp - 1.0e-12_dp, 0.0_dp, 0.0_dp]
    call sample(u, x, b, bmod)
    n = [1.0_dp, 0.0_dp, 0.0_dp]
    tangent = cross(n, b)
    tangent = tangent/norm2(tangent)
    state = spectre_fo_state_t()
    state%x = x - 0.0008_dp*(RO/bmod)*n
    state%v = 0.05_dp*n + sqrt(1.0_dp - 0.05_dp**2)*tangent
    state%u = u
    state%owner = 1
    state%iface = 1
    state%active = .true.
    state%exit_rho = 10.0_dp
    dt = 0.1_dp*RO/bmod
    state%h_rk = dt
    call spectre_fo_canonical_pzeta(state, RO, p0, status)
    if (status /= SPECTRE_FO_OK) error stop 'initial invariant sample failed'
    state%entry_pzeta = p0
    state%entry_owner = 1
    scale = norm2(state%v)*norm2(state%x)
    xref = state%x
    vref = state%v
    call reference_trace(xref, vref, dt, 10000)
    call spectre_fo_advance_until_exit(state, dt, RO, used, y, owner, &
        exited, status)
    if (status /= SPECTRE_FO_OK .or. exited) error stop 'incomplete tangent trace'
    if (abs(used/dt - 1.0_dp) > 1.0e-12_dp) error stop 'missing tangent time'
    if (state%owner /= 1) error stop 'tangent orbit did not return'
    print '(A,ES12.4)', 'cross-and-return velocity error against RK4: ', &
        norm2(state%v - vref)
    if (norm2(state%v - vref) > 2.0e-6_dp) &
        error stop 'missed near-tangent trajectory impulse'
    call spectre_fo_canonical_pzeta(state, RO, p1, status)
    if (status /= SPECTRE_FO_OK) error stop 'final invariant sample failed'
    defect = abs(p1 - p0)/scale
    print '(A,ES12.4)', 'cross-and-return momentum defect: ', defect
    if (defect > 1.0e-8_dp) error stop 'missed near-tangent sheet impulse'

    ! A separately sampled trace witnesses the excursion: it really crosses
    ! both ways, rather than satisfying the invariant by staying in volume 1.
    call sample(u, x, b, bmod)
    state%x = x - 0.0008_dp*(RO/bmod)*n
    state%v = 0.05_dp*n + sqrt(1.0_dp - 0.05_dp**2)*tangent
    state%u = u
    state%owner = 1
    state%h_rk = dt
    visits = 0
    do j = 1, 20
        call spectre_fo_advance_until_exit(state, dt/20, RO, used, y, owner, &
            exited, status)
        if (status /= SPECTRE_FO_OK .or. exited) error stop 'sampled trace failed'
        if (state%owner == 2) visits = visits + 1
    end do
    if (visits == 0 .or. state%owner /= 1) error stop 'missing excursion witness'

    ! Entry at the sheet from either gauge, including Larmor samples on the
    ! other side: the recorded residual must equal actual P in the entry gauge.
    cross_entries = 0
    do owner = 1, 2
        do j = 1, 16
            y = [1.0_dp, 8.0_dp*atan(1.0_dp)*real(mod(j - 1, 8), dp)/8, &
                0.0_dp, 1.0_dp, merge(0.5_dp, -0.5_dp, j <= 8)]
            call spectre_fo_enter(state, y, owner, 1, RO, status)
            if (owner == 1 .and. j == 1) then
                ! This no-root fixture exceeds the finite defect bound.
                if (status /= SPECTRE_FO_GYROPHASE_FAIL .or. state%active) &
                    error stop 'unbounded no-root entry was accepted'
                cycle
            end if
            if (status /= SPECTRE_FO_OK) error stop 'sheet entry failed'
            call spectre_fo_canonical_pzeta(state, RO, p1, status)
            if (status /= SPECTRE_FO_OK) error stop 'entry invariant sample failed'
            if (state%owner /= owner) then
                cross_entries(owner) = cross_entries(owner) + 1
                u = y(1:3)
                u(1) = 1.0_dp - 1.0e-12_dp
                call spectre_field%evaluate(u, acov, hcov, bmod)
                p0 = acov(3)
                u(1) = 1.0_dp
                call spectre_field%evaluate(u, acov, hcov, bmod)
                p0 = (acov(3) - p0)*TESLA_TO_GAUSS*M_TO_CM**2/RO
                p1 = p1 - real(state%owner - owner, dp)*p0
            end if
            defect = p1 - state%entry_pzeta
            if (abs(defect - state%entry_pzeta_defect) > 1.0e-10_dp) &
                error stop 'entry residual uses inconsistent gauges'
            if (abs(defect) > 1.0e-3_dp*max(abs(state%entry_pzeta), 1.0_dp)) &
                error stop 'unbounded entry momentum defect'
        end do
        print '(A,I0,A,I0)', 'entry owner ', owner, ' cross-gauge entries: ', &
            cross_entries(owner)
    end do
    if (any(cross_entries == 0)) error stop 'missing cross-gauge entry coverage'

    ! An outer-domain trial stage must retry deterministically. Only an
    ! actually located terminal point can report physical loss.
    u = [1.99_dp, 0.0_dp, 0.0_dp]
    call sample(u, x, b, bmod)
    state = spectre_fo_state_t()
    state%x = x
    state%v = n
    state%u = u
    state%owner = 2
    state%iface = 1
    state%active = .true.
    state%exit_rho = 10.0_dp
    dt = 2.0_dp*RO/bmod
    state%h_rk = dt
    call spectre_fo_advance_until_exit(state, dt, RO, used, y, owner, &
        exited, status)
    if (status /= SPECTRE_FO_LOSS) error stop 'outer-domain retry did not find loss'
    if (used < 0.0_dp .or. used > dt) error stop 'invalid loss time budget'
    if (state%u(1) < 2.0_dp - 2.0e-12_dp) error stop 'unlocalized outer loss'
    call zero_time_field_jump
    call inward_thin_volume
    call hidden_opposite_crossing
    print '(A)', 'PASS'

contains

    subroutine zero_time_field_jump
        ! A zero-time switch into a stronger field must recompute the gyro cap.
        ! Loose RK tolerance isolates that bound; independent RK4 supplies the
        ! trajectory oracle, without sharing the trial or event algorithm.
        type(spectre_field_t) :: jump
        type(spectre_fo_state_t) :: orbit
        real(dp) :: u(3), x(3), b(3), bmod, dt, used, y(5), pz, xr(3), vr(3)
        integer :: status, owner
        logical :: exited

        jump = sf
        jump%data%vol(2)%Ate = 20.0_dp*jump%data%vol(2)%Ate
        jump%data%vol(2)%Aze = 20.0_dp*jump%data%vol(2)%Aze
        jump%data%vol(2)%Ato = 20.0_dp*jump%data%vol(2)%Ato
        jump%data%vol(2)%Azo = 20.0_dp*jump%data%vol(2)%Azo
        call set_magfie_spectre_field(jump)
        call spectre_fo_configure( 0.1_dp)
        u = [1.0_dp - 1.0e-12_dp, 0.0_dp, 0.0_dp]
        call sample(u, x, b, bmod)
        u(1) = 1.0_dp
        orbit = spectre_fo_state_t()
        call spectre_field%coords%evaluate_cart(u, orbit%x)
        orbit%x = orbit%x*M_TO_CM
        orbit%v = [1.0_dp, 0.0_dp, 0.0_dp]
        orbit%u = u
        orbit%owner = 1
        orbit%iface = 1
        orbit%active = .true.
        orbit%exit_rho = 10.0_dp
        call spectre_fo_canonical_pzeta(orbit, RO, pz, status)
        if (status /= SPECTRE_FO_OK) error stop 'field-jump initial sample failed'
        orbit%entry_owner = 1
        orbit%entry_pzeta = pz
        dt = 0.1_dp*RO/bmod
        orbit%h_rk = dt
        ! A nanometre outward offset gives the oracle the outgoing side
        ! at t=0, rather than the inverse's arbitrary exact-sheet convention.
        xr = orbit%x + [1.0e-7_dp, 0.0_dp, 0.0_dp]
        vr = orbit%v
        call reference_trace(xr, vr, dt, 10000)
        call spectre_fo_advance_until_exit(orbit, dt, RO, used, y, owner, &
            exited, status)
        if (status /= SPECTRE_FO_OK .or. exited) error stop 'field-jump trace failed'
        if (abs(used/dt - 1.0_dp) > 1.0e-12_dp) error stop 'missing field-jump time'
        print '(A,ES12.4)', 'zero-time field-jump velocity error against RK4: ', &
            norm2(orbit%v - vr)
        if (norm2(orbit%v - vr) > 1.0e-7_dp) &
            error stop 'gyro cap lost after zero-time owner change'
        call set_magfie_spectre_field(sf)
        call spectre_fo_configure( 1.0e-6_dp)
    end subroutine zero_time_field_jump

    subroutine inward_thin_volume
        ! Split the fixture's outer volume with a very thin intermediate
        ! chart. Its coefficients scale with its width to keep |B| finite.
        ! A long inward trial reaches volume 1 through a reversed-field thin
        ! volume 2. Independent RK4 must see its impulse at the first sheet.
        type(spectre_fo_state_t) :: orbit
        real(dp) :: u(3), x(3), b(3), bmod, pz, dt, used, y(5), xr(3), vr(3)
        integer :: status, owner
        logical :: exited

        call set_thin_field(1.0e-4_dp)
        call spectre_fo_configure( 1.0e-6_dp)
        u = [2.00001_dp, 0.0_dp, 0.0_dp]
        call sample(u, x, b, bmod)
        orbit = spectre_fo_state_t()
        orbit%x = x
        orbit%v = [-1.0_dp, 0.0_dp, 0.0_dp]
        orbit%u = u
        orbit%owner = 3
        orbit%iface = 2
        orbit%active = .true.
        orbit%exit_rho = 10.0_dp
        call spectre_fo_canonical_pzeta(orbit, RO, pz, status)
        if (status /= SPECTRE_FO_OK) error stop 'thin-volume initial sample failed'
        orbit%entry_owner = 3
        orbit%entry_pzeta = pz
        dt = 0.1_dp*RO/bmod
        orbit%h_rk = dt
        xr = orbit%x
        vr = orbit%v
        call reference_trace(xr, vr, dt, 100000)
        call spectre_fo_advance_until_exit(orbit, dt, RO, used, y, owner, &
            exited, status)
        if (status /= SPECTRE_FO_OK .or. exited) &
            error stop 'thin-volume traversal failed'
        if (abs(used/dt - 1.0_dp) > 1.0e-12_dp) error stop 'missing inward time'
        if (orbit%owner /= 1) error stop 'inward orbit did not reach volume 1'
        print '(A,ES12.4)', 'thin-volume velocity error against RK4: ', &
            norm2(orbit%v - vr)
        if (norm2(orbit%v - vr) > 2.0e-6_dp) &
            error stop 'inward trial skipped thin-volume impulse'
    end subroutine inward_thin_volume

    subroutine hidden_opposite_crossing
        ! An early inward excursion must precede the home-field continuation's
        ! later outward landing in a thin volume. Check the first sheet's
        ! impulse with independent geometric-field RK4.
        type(spectre_fo_state_t) :: orbit
        real(dp) :: u(3), x(3), b(3), bmod, dt, used, y(5), pz, xr(3), vr(3)
        real(dp) :: tangent(3), normal(3)
        integer :: status, owner
        logical :: exited

        call set_thin_field(1.0e-5_dp)
        call spectre_fo_configure( 1.0e-6_dp)
        u = [1.025_dp, 0.0_dp, 0.0_dp]
        call sample(u, x, b, bmod)
        normal = [1.0_dp, 0.0_dp, 0.0_dp]
        tangent = cross(b, normal)
        tangent = tangent/norm2(tangent)
        orbit = spectre_fo_state_t()
        orbit%x = x
        orbit%v = -0.01_dp*normal + sqrt(1.0_dp - 0.01_dp**2)*tangent
        orbit%u = u
        orbit%owner = 2
        orbit%iface = 1
        orbit%active = .true.
        orbit%exit_rho = 10.0_dp
        call spectre_fo_canonical_pzeta(orbit, RO, pz, status)
        if (status /= SPECTRE_FO_OK) error stop 'hidden-pulse initial sample failed'
        orbit%entry_owner = 2
        orbit%entry_pzeta = pz
        dt = 0.1_dp*RO/bmod
        orbit%h_rk = dt
        xr = orbit%x
        vr = orbit%v
        call reference_trace(xr, vr, dt, 100000)
        call spectre_fo_advance_until_exit(orbit, dt, RO, used, y, owner, &
            exited, status)
        if (status /= SPECTRE_FO_OK .or. exited) error stop 'hidden-pulse trace failed'
        if (abs(used/dt - 1.0_dp) > 1.0e-12_dp) error stop 'missing hidden-pulse time'
        print '(A,ES12.4)', 'landing-trial hidden-pulse velocity error against RK4: ', &
            norm2(orbit%v - vr)
        if (orbit%owner /= 1 .or. norm2(orbit%v - vr) > 2.0e-6_dp) &
            error stop 'landing trial skipped the first inward pulse'
    end subroutine hidden_opposite_crossing

    subroutine set_thin_field(width)
        real(dp), intent(in) :: width
        type(spectre_field_t) :: thin

        thin = sf
        thin%data%Mvol = 3
        thin%data%Nvol = 3
        thin%data%Lrad = [sf%data%Lrad, sf%data%Lrad(2)]
        thin%data%vol = [sf%data%vol, sf%data%vol(2)]
        call split_geometry(sf%data%Rbc, thin%data%Rbc, width)
        call split_geometry(sf%data%Rbs, thin%data%Rbs, width)
        call split_geometry(sf%data%Zbc, thin%data%Zbc, width)
        call split_geometry(sf%data%Zbs, thin%data%Zbs, width)
        call linear_sheet_coefficients(sf%data%vol(2)%Ate, &
            thin%data%vol(2)%Ate, width)
        call linear_sheet_coefficients(sf%data%vol(2)%Aze, &
            thin%data%vol(2)%Aze, width)
        call linear_sheet_coefficients(sf%data%vol(2)%Ato, &
            thin%data%vol(2)%Ato, width)
        call linear_sheet_coefficients(sf%data%vol(2)%Azo, &
            thin%data%vol(2)%Azo, width)
        thin%data%vol(2)%Ate = -2.0_dp*thin%data%vol(2)%Ate
        thin%data%vol(2)%Aze = -2.0_dp*thin%data%vol(2)%Aze
        thin%data%vol(2)%Ato = -2.0_dp*thin%data%vol(2)%Ato
        thin%data%vol(2)%Azo = -2.0_dp*thin%data%vol(2)%Azo
        call make_spectre_coordinate_system(thin%coords, thin%data)
        call set_magfie_spectre_field(thin)
    end subroutine set_thin_field

    subroutine linear_sheet_coefficients(original, linear, width)
        ! Keep the outer volume's inner-edge dA/ds constant through the thin
        ! chart before reversing it, giving a finite sheet-force contrast.
        ! Chebyshev T_l'(-1)=(-1)^(l-1) l^2.
        real(dp), intent(in) :: original(0:, :)
        real(dp), intent(out) :: linear(0:, :)
        real(dp), intent(in) :: width
        integer :: l

        linear = 0.0_dp
        do l = 1, size(original, 1) - 1
            linear(1, :) = linear(1, :) &
                + WIDTH*real((-1)**(l - 1)*l*l, dp)*original(l, :)
        end do
        linear(0, :) = linear(1, :)
    end subroutine linear_sheet_coefficients

    subroutine split_geometry(original, split, width)
        real(dp), intent(in) :: original(:, 0:)
        real(dp), allocatable, intent(out) :: split(:, :)
        real(dp), intent(in) :: width

        allocate (split(size(original, 1), 0:3))
        split(:, 0:1) = original(:, 0:1)
        split(:, 2) = (1.0_dp - WIDTH)*original(:, 1) + WIDTH*original(:, 2)
        split(:, 3) = original(:, 2)
    end subroutine split_geometry

    subroutine reference_trace(x, v, dt, nsteps)
        ! Independent fixed-step RK4, always sampling the geometric field.
        ! Fine fixed steps resolve the sheet-force quadrature independently;
        ! this oracle has no owner continuation, landing, or GC projection.
        real(dp), intent(inout) :: x(3), v(3)
        real(dp), intent(in) :: dt
        integer, intent(in) :: nsteps
        real(dp) :: h, kx(3, 4), kv(3, 4)
        integer :: i

        h = dt/real(nsteps, dp)
        do i = 1, nsteps
            kx(:, 1) = v
            call acceleration(x, v, kv(:, 1))
            kx(:, 2) = v + 0.5_dp*h*kv(:, 1)
            call acceleration(x + 0.5_dp*h*kx(:, 1), kx(:, 2), kv(:, 2))
            kx(:, 3) = v + 0.5_dp*h*kv(:, 2)
            call acceleration(x + 0.5_dp*h*kx(:, 2), kx(:, 3), kv(:, 3))
            kx(:, 4) = v + h*kv(:, 3)
            call acceleration(x + h*kx(:, 3), kx(:, 4), kv(:, 4))
            x = x + h*(kx(:, 1) + 2*kx(:, 2) + 2*kx(:, 3) + kx(:, 4))/6
            v = v + h*(kv(:, 1) + 2*kv(:, 2) + 2*kv(:, 3) + kv(:, 4))/6
        end do
    end subroutine reference_trace

    subroutine acceleration(x, v, dv)
        real(dp), intent(in) :: x(3), v(3)
        real(dp), intent(out) :: dv(3)
        real(dp) :: cyl(3), u(3), xu(3), b(3), bmod
        integer :: ierr

        cyl = [sqrt(x(1)**2 + x(2)**2)/M_TO_CM, atan2(x(2), x(1)), &
            x(3)/M_TO_CM]
        select type (cs => spectre_field%coords)
        type is (spectre_coordinate_system_t)
            call cs%from_cyl(cyl, u, ierr)
        class default
            error stop 'unexpected reference coordinates'
        end select
        if (ierr /= 0) error stop 'reference location failed'
        call sample(u, xu, b, bmod)
        dv = (bmod/RO)*cross(v, b)
    end subroutine acceleration

    subroutine sample(u, x, b, bmod)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3), b(3), bmod
        real(dp) :: a(3), hc(3), sqgb(3), g(3, 3), gi(3, 3), sg, e(3, 3)

        call spectre_field%coords%evaluate_cart(u, x)
        call spectre_field%coords%covariant_basis(u, e)
        call spectre_field%coords%metric_tensor(u, g, gi, sg)
        call spectre_field%evaluate(u, a, hc, bmod, sqgb)
        b = matmul(e, sqgb/sg)/bmod
        x = x*M_TO_CM
        bmod = bmod*TESLA_TO_GAUSS
    end subroutine sample

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross
end program test_spectre_fo_events
