program test_magfie_refcoords_splines
    !> Verify magfie_refcoords assembly against finite-difference oracles.
    !>
    !> Oracle sources:
    !> - Bmod, hcov: splined_field_t%evaluate
    !> - bder, dhcov: 5-point central finite differences of evaluate()
    !> - sqrtg, ginv: finite differences of coordinate_system_t%evaluate_cart
    !>
    !> Target:
    !> - magfie_refcoords (via init_magfie(REFCOORDS)) for a splined field, and
    !>   hcurl of its finite-difference path for the unsplined VMEC field
    !>   against the curl definition

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec, tracer_t, init_params
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use field_coils, only: coils_field_t, create_coils_field
    use field_vmec, only: vmec_field_t, create_vmec_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use magfie_sub, only: init_magfie, magfie, REFCOORDS, BOOZER, &
                          set_magfie_refcoords_field, axis_radial_coordinate
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system, UNKNOWN, RHO_TOR
    use alpha_lifetime_sub, only: velo_can, velo_axis_regularized, orbit_timestep_axis
    use finite_check, only: is_finite

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: f
    type(vmec_field_t) :: vmec_field
    real(dp) :: dummy

    integer, parameter :: n_points = 240
    real(dp), parameter :: h_rho = 1.0e-4_dp
    real(dp), parameter :: h_ang = 1.0e-5_dp

    real(dp), parameter :: tol_bmod_rel = 1.0e-12_dp
    real(dp), parameter :: tol_bmod_abs = 1.0e-12_dp
    real(dp), parameter :: tol_sqrtg_rel = 5.0e-7_dp
    real(dp), parameter :: tol_sqrtg_abs = 1.0e-10_dp
    real(dp), parameter :: tol_bder_rel = 5.0e-8_dp
    real(dp), parameter :: tol_bder_abs = 1.0e-10_dp
    real(dp), parameter :: tol_hcov_rel = 1.0e-12_dp
    real(dp), parameter :: tol_hcov_abs = 1.0e-12_dp
    real(dp), parameter :: tol_hctr_rel = 1.0e-8_dp
    real(dp), parameter :: tol_hctr_abs = 5.0e-11_dp
    real(dp), parameter :: tol_hcurl_rel = 1.0e-6_dp
    real(dp), parameter :: tol_hcurl_abs = 1.0e-12_dp

    real(dp) :: phi_period
    integer :: n_failed, i, axis_error
    logical :: radial_is_rho

    n_failed = 0

    ! Derivative consistency needs only a coarse VMEC coordinate fixture.
    ! Keep quintic splines; multiplier 1 bounds the coordinate arrays to 0.43 GiB.
    call init_vmec('wout_ncsx.nc', 5, 5, 1, dummy)
    call init_reference_coordinates('wout_ncsx.nc')

    call create_coils_field('coils.simple', raw_coils)
    ! Resolve curl(h) near the axis: 16 radial samples introduce a spurious
    ! B_parallel_star zero just outside the default reconstruction disc.
    call create_splined_field(raw_coils, ref_coords, f, n_r=64, n_th=17, n_phi=16)

    call set_magfie_refcoords_field(f)
    call init_magfie(REFCOORDS)
    call axis_radial_coordinate(radial_is_rho, axis_error)
    if (axis_error /= 0 .or. .not. radial_is_rho) &
        error stop 'splined VMEC RK chart must use rho'
    ! Retained reference storage must not select the chart of another backend.
    call init_magfie(BOOZER)
    call axis_radial_coordinate(radial_is_rho, axis_error)
    if (axis_error /= 0 .or. radial_is_rho) &
        error stop 'inactive reference field changed the RK axis chart'
    call init_magfie(REFCOORDS)
    call check_active_rho_rhs(n_failed)
    call check_default_axis_disc(n_failed)

    phi_period = twopi/real(nper, dp)

    do i = 1, n_points
        call check_point(i, phi_period, f, n_failed)
    end do

    call create_vmec_field(vmec_field)
    call set_magfie_refcoords_field(vmec_field)
    call axis_radial_coordinate(radial_is_rho, axis_error)
    if (axis_error /= 0 .or. radial_is_rho) &
        error stop 'unsplined VMEC RK chart must use flux'
    do i = 1, 20
        call check_unsplined_hcurl(i, phi_period, vmec_field, n_failed)
    end do
    call check_geometric_chart_rk(n_failed)

    if (n_failed > 0) then
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' mismatch(es) in magfie_refcoords test'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'PASSED: magfie_refcoords matches FD oracles'
    print *, '================================'

contains

    subroutine check_geometric_chart_rk(n_failed)
        integer, intent(inout) :: n_failed
        class(coordinate_system_t), allocatable :: chart
        type(splined_field_t) :: chart_field
        real(dp) :: initial(5), flux_result(5), geometric_result(5)
        integer :: status

        call make_chartmap_coordinate_system(chart, 'wout_ncsx.chartmap.nc')
        call create_splined_field(raw_coils, chart, chart_field, &
                                  n_r=16, n_th=17, n_phi=16)
        initial = [0.45_dp, 0.7_dp, 0.1_dp, 1.0_dp, -0.3_dp]
        call init_magfie(REFCOORDS)
        select type (coords => chart_field%coords)
        class is (chartmap_coordinate_system_t)
            coords%rho_convention = RHO_TOR
            call set_magfie_refcoords_field(chart_field)
            flux_result = initial
            call orbit_timestep_axis(flux_result, 1.0e-3_dp, 1.0e-3_dp, &
                                     1.0e-10_dp, status)
            if (status /= 0) error stop 'Reference chart RK step failed'
            coords%rho_convention = UNKNOWN
            call set_magfie_refcoords_field(chart_field)
            geometric_result = initial
            call orbit_timestep_axis(geometric_result, 1.0e-3_dp, 1.0e-3_dp, &
                                     1.0e-10_dp, status)
            if (status /= 0 .or. .not. all(is_finite(geometric_result))) then
                print *, 'FAIL: geometric disk chart rejected by RK'
                n_failed = n_failed + 1
            else if (maxval(abs(geometric_result - flux_result)) > 1.0e-12_dp &
                     .or. maxval(abs(geometric_result - initial)) <= 1.0e-12_dp) then
                print *, 'FAIL: geometric radius changed the physical RK step'
                n_failed = n_failed + 1
            end if
        class default
            error stop 'Expected chartmap reference coordinates'
        end select
    end subroutine check_geometric_chart_rk

    subroutine check_active_rho_rhs(n_failed)
        integer, intent(inout) :: n_failed
        type(tracer_t) :: orbit
        real(dp), parameter :: native_floor = 1.0e-6_dp, theta = 0.7_dp
        real(dp) :: zpolar(5), zaxis(5), vpolar(5), expected(5), actual(5)
        real(dp) :: axis_value(5), probe(5), after(5), after_cart(2)
        real(dp) :: rho, angle, scale, spread, step
        integer :: ray, status

        ! Reuse the active splined coil field and existing VMEC normalization.
        ! No new field initialization is needed for this physical RHS oracle.
        call init_params(orbit, 2, 4, 3.5e6_dp, 256, 1, 1.0e-12_dp)
        rho = 0.005_dp
        zpolar = [rho, theta, 0.1_dp, 1.0_dp, -0.3_dp]
        zaxis = [rho*cos(theta), rho*sin(theta), zpolar(3:5)]
        call velo_can(0.0_dp, zpolar, vpolar)
        ! Independent native-rho chain rule, not the production map helper.
        expected(1) = vpolar(1)*cos(theta) - rho*sin(theta)*vpolar(2)
        expected(2) = vpolar(1)*sin(theta) + rho*cos(theta)*vpolar(2)
        expected(3:5) = vpolar(3:5)
        call velo_axis_regularized(0.0_dp, zaxis, actual, native_floor)
        if (.not. all(is_finite([expected, actual]))) &
            error stop 'Active REF rho RHS produced nonfinite phase velocity'
        scale = max(maxval(abs(expected)), tiny(1.0_dp))
        if (maxval(abs(actual - expected)) > 1.0e-11_dp*scale) then
            print *, 'FAIL: active REF rho regularization changed outer polar RHS'
            n_failed = n_failed + 1
        end if

        ! A short near-axis push must return native rho and reproduce the
        ! independently evaluated Cartesian tangent to first order. Choose a
        ! small physical displacement rather than a device-dependent time.
        scale = max(maxval(abs(expected(1:2))), tiny(1.0_dp))
        step = min(1.0e-5_dp*rho/scale, &
                   1.0e-5_dp/max(maxval(abs(vpolar(2:5))), tiny(1.0_dp)))
        after = zpolar
        call orbit_timestep_axis(after, step, step, 1.0e-12_dp, status, &
                                 axis_floor=native_floor, near_axis_limit=0.01_dp)
        if (status /= 0 .or. .not. all(is_finite(after))) &
            error stop 'Active REF rho near-axis push failed'
        after_cart = [after(1)*cos(after(2)), after(1)*sin(after(2))]
        if (maxval(abs((after_cart - zaxis(1:2))/step - expected(1:2))) > &
            5.0e-3_dp*scale) then
            print *, 'FAIL: active REF near-axis push lost native rho tangent'
            n_failed = n_failed + 1
        end if

        ! Immediately inside the native-rho handoff, recover the polar field
        ! rather than a second chart's radial units.
        rho = (1.0_dp - 1.0e-6_dp)*native_floor
        zpolar(1) = rho
        zaxis(1:2) = [rho*cos(theta), rho*sin(theta)]
        call velo_can(0.0_dp, zpolar, vpolar)
        expected(1) = vpolar(1)*cos(theta) - rho*sin(theta)*vpolar(2)
        expected(2) = vpolar(1)*sin(theta) + rho*cos(theta)*vpolar(2)
        expected(3:5) = vpolar(3:5)
        call velo_axis_regularized(0.0_dp, zaxis, actual, native_floor)
        if (.not. all(is_finite([expected, actual]))) &
            error stop 'Active REF rho handoff produced nonfinite velocity'
        scale = max(maxval(abs(expected)), tiny(1.0_dp))
        if (maxval(abs(actual - expected)) > 1.0e-6_dp*scale) then
            print *, 'FAIL: active REF rho reconstruction has a handoff seam'
            n_failed = n_failed + 1
        end if

        ! Check the actual reconstruction's axis limit using the same active
        ! rho field. Toroidal position and particle state remain fixed.
        probe = [0.0_dp, 0.0_dp, zpolar(3:5)]
        call velo_axis_regularized(0.0_dp, probe, axis_value, native_floor)
        if (.not. all(is_finite(axis_value))) &
            error stop 'Active REF rho reconstruction has no finite axis value'
        spread = 0.0_dp
        do ray = 0, 7
            angle = twopi*real(ray, dp)/8.0_dp
            probe(1:2) = 1.0e-14_dp*[cos(angle), sin(angle)]
            call velo_axis_regularized(0.0_dp, probe, actual, native_floor)
            if (.not. all(is_finite(actual))) &
                error stop 'Active REF rho axis limit depends on a singular ray'
            spread = max(spread, maxval(abs(actual - axis_value)))
        end do
        scale = max(maxval(abs(axis_value)), 1.0e-10_dp)
        if (spread > 1.0e-6_dp*scale) then
            print *, 'FAIL: active REF rho axis limit depends on approach ray'
            n_failed = n_failed + 1
        end if
    end subroutine check_active_rho_rhs

    subroutine check_default_axis_disc(n_failed)
        !> The default regularization disc and near-axis zone are defined in flux
        !> units (1e-8, 0.01) and must keep their Cartesian size (1e-4, 0.1) on a
        !> rho-native chart. A step that starts within 1e-6 of the axis crosses it
        !> at Cartesian radius ~1e-6 and used to exhaust the RK step budget with
        !> the unconverted defaults (ierr = 2). The reference integrates the same
        !> orbit in many short steps with the equivalent explicit native limits.
        integer, intent(inout) :: n_failed
        type(tracer_t) :: orbit
        real(dp), parameter :: rho0(2) = [1.0e-8_dp, 1.0e-6_dp]
        real(dp), parameter :: pitch(3) = [0.3_dp, -0.6_dp, 0.05_dp]
        real(dp), parameter :: tol = 1.0e-9_dp
        real(dp) :: z(5), zref(5), dist, dtau
        integer :: irho, ipitch, istep, status, status_ref
        integer, parameter :: n_sub = 40

        call init_params(orbit, 2, 4, 3.5e6_dp, 256, 1, 1.0e-10_dp)
        dtau = orbit%dtau
        do irho = 1, size(rho0)
            do ipitch = 1, size(pitch)
                z = [rho0(irho), 0.4_dp, 0.1_dp, 1.0_dp, pitch(ipitch)]
                zref = z
                call orbit_timestep_axis(z, dtau, dtau, 1.0e-10_dp, status)
                status_ref = 0
                do istep = 1, n_sub
                    if (status_ref == 0) &
                        call orbit_timestep_axis(zref, dtau/n_sub, dtau/n_sub, &
                                                 1.0e-10_dp, status_ref, &
                                                 axis_floor=1.0e-4_dp, &
                                                 near_axis_limit=0.1_dp)
                end do
                if (status /= 0 .or. status_ref /= 0 .or. &
                    .not. all(is_finite([z, zref]))) then
                    print *, 'FAIL: default axis disc orbit failed', irho, ipitch, &
                        status, status_ref
                    n_failed = n_failed + 1
                    cycle
                end if
                dist = hypot(z(1)*cos(z(2)) - zref(1)*cos(zref(2)), &
                             z(1)*sin(z(2)) - zref(1)*sin(zref(2)))
                if (dist > tol .or. maxval(abs(z(3:5) - zref(3:5))) > tol) then
                    print *, 'FAIL: default axis disc differs from explicit disc', &
                        irho, ipitch, dist
                    n_failed = n_failed + 1
                end if
            end do
        end do
    end subroutine check_default_axis_disc

    subroutine check_point(i, phi_period, f, n_failed)
        integer, intent(in) :: i
        real(dp), intent(in) :: phi_period
        type(splined_field_t), intent(in) :: f
        integer, intent(inout) :: n_failed

        real(dp) :: x(3)
        real(dp) :: bmod_t, sqrtg_t, bder_t(3), hcov_t(3), hctr_t(3), hcurl_t(3)
        real(dp) :: B0, h0(3)
        real(dp) :: bder_fd(3), hcov_fd(3), hctr_fd(3), hcurl_fd(3)
        real(dp) :: dhcov_fd(3, 3)
        real(dp) :: ginv_fd(3, 3), sqrtg_fd

        x(1) = 0.25_dp + 0.55_dp*real(mod(13*i, 997), dp)/997.0_dp
        x(2) = 0.1_dp + (twopi - 0.2_dp)*real(mod(29*i, 991), dp)/991.0_dp
        x(3) = 0.1_dp + (phi_period - 0.2_dp)*real(mod(31*i, 983), dp)/983.0_dp

        call oracle_field_fd(f, x, phi_period, B0, h0, bder_fd, dhcov_fd)
        hcov_fd = h0

        call metric_from_cart_fd(ref_coords, x, phi_period, ginv_fd, sqrtg_fd)
        hctr_fd = matmul(ginv_fd, hcov_fd)
        call compute_hcurl(sqrtg_fd, dhcov_fd, hcurl_fd)

        call magfie(x, bmod_t, sqrtg_t, bder_t, hcov_t, hctr_t, hcurl_t)

        call check_scalar('bmod', x, bmod_t, B0, tol_bmod_abs, tol_bmod_rel, n_failed)
        call check_scalar('sqrtg', x, sqrtg_t, sqrtg_fd, tol_sqrtg_abs, &
                          tol_sqrtg_rel, n_failed)
        call check_vec('bder', x, bder_t, bder_fd, tol_bder_abs, tol_bder_rel, n_failed)
        call check_vec('hcov', x, hcov_t, hcov_fd, tol_hcov_abs, tol_hcov_rel, n_failed)
        call check_vec('hctr', x, hctr_t, hctr_fd, tol_hctr_abs, tol_hctr_rel, n_failed)
        call check_vec('hcurl', x, hcurl_t, hcurl_fd, tol_hcurl_abs, tol_hcurl_rel, &
                       n_failed)
    end subroutine check_point

    subroutine oracle_field_fd(f, x, phi_period, B0, h0, bder, dh)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3), phi_period
        real(dp), intent(out) :: B0, h0(3)
        real(dp), intent(out) :: bder(3), dh(3, 3)

        real(dp) :: Acov0(3), Bp1, Bm1, Bp2, Bm2
        real(dp) :: hp1(3), hm1(3), hp2(3), hm2(3)
        integer :: k
        real(dp) :: step, dBdx

        call f%evaluate(x, Acov0, h0, B0)

        do k = 1, 3
            call eval_shift(f, x, k, +1, phi_period, Bp1, hp1)
            call eval_shift(f, x, k, -1, phi_period, Bm1, hm1)
            call eval_shift(f, x, k, +2, phi_period, Bp2, hp2)
            call eval_shift(f, x, k, -2, phi_period, Bm2, hm2)

            if (k == 1) then
                step = h_rho
            else
                step = h_ang
            end if

            dBdx = (-Bp2 + 8.0_dp*Bp1 - 8.0_dp*Bm1 + Bm2)/(12.0_dp*step)
            bder(k) = dBdx/max(B0, 1.0e-30_dp)
            dh(k, :) = (-hp2 + 8.0_dp*hp1 - 8.0_dp*hm1 + hm2)/(12.0_dp*step)
        end do
    end subroutine oracle_field_fd

    subroutine check_unsplined_hcurl(i, phi_period, fld, n_failed)
        !> hcurl^i = eps^ijk d_j h_k / sqrt(g) spelled out with central
        !> differences of the field's own evaluate, so the oracle shares no
        !> index convention with magfie.
        integer, intent(in) :: i
        real(dp), intent(in) :: phi_period
        type(vmec_field_t), intent(in) :: fld
        integer, intent(inout) :: n_failed

        real(dp), parameter :: h_fd = 1.0e-5_dp
        real(dp) :: x(3), xp(3), xm(3), Acov(3), hp(3), hm(3), bp, bm
        real(dp) :: bmod, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)
        real(dp) :: dh(3, 3), hcurl_fd(3)
        integer :: j

        x(1) = 0.25_dp + 0.5_dp*real(mod(13*i, 97), dp)/97.0_dp
        x(2) = 0.1_dp + (twopi - 0.2_dp)*real(mod(29*i, 89), dp)/89.0_dp
        x(3) = 0.1_dp + (phi_period - 0.2_dp)*real(mod(31*i, 83), dp)/83.0_dp

        call magfie(x, bmod, sqrtg, bder, hcov, hctr, hcurl)
        do j = 1, 3
            xp = x
            xm = x
            xp(j) = x(j) + h_fd
            xm(j) = x(j) - h_fd
            call fld%evaluate(xp, Acov, hp, bp)
            call fld%evaluate(xm, Acov, hm, bm)
            dh(:, j) = (hp - hm)/(2.0_dp*h_fd)  ! dh(k, j) = d h_k / d x_j
        end do
        hcurl_fd(1) = (dh(3, 2) - dh(2, 3))/sqrtg
        hcurl_fd(2) = (dh(1, 3) - dh(3, 1))/sqrtg
        hcurl_fd(3) = (dh(2, 1) - dh(1, 2))/sqrtg

        call check_vec('hcurl unsplined', x, hcurl, hcurl_fd, tol_hcurl_abs, &
                       1.0e-4_dp, n_failed)
    end subroutine check_unsplined_hcurl

    subroutine eval_shift(f, x, k, m, phi_period, Bmod, hcov)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k, m
        real(dp), intent(in) :: phi_period
        real(dp), intent(out) :: Bmod, hcov(3)

        real(dp) :: x_shift(3)
        real(dp) :: Acov(3)

        x_shift = x
        select case (k)
        case (1)
            x_shift(1) = x(1) + real(m, dp)*h_rho
        case (2)
            x_shift(2) = wrap_angle(x(2) + real(m, dp)*h_ang, twopi)
        case (3)
            x_shift(3) = wrap_angle(x(3) + real(m, dp)*h_ang, phi_period)
        case default
            error stop 'eval_shift: invalid coordinate index'
        end select

        call f%evaluate(x_shift, Acov, hcov, Bmod)
    end subroutine eval_shift

    subroutine metric_from_cart_fd(cs, x, phi_period, ginv, sqrtg)
        class(*), intent(in) :: cs
        real(dp), intent(in) :: x(3), phi_period
        real(dp), intent(out) :: ginv(3, 3), sqrtg

        real(dp) :: u0(3), u_p1(3), u_m1(3), u_p2(3), u_m2(3)
        real(dp) :: xc_p1(3), xc_m1(3), xc_p2(3), xc_m2(3)
        real(dp) :: e1(3), e2(3), e3(3)
        real(dp) :: g(3, 3), det

        u0 = [x(1)**2, x(2), x(3)]

        u_p1 = [(x(1) + h_rho)**2, x(2), x(3)]
        u_m1 = [(x(1) - h_rho)**2, x(2), x(3)]
        u_p2 = [(x(1) + 2.0_dp*h_rho)**2, x(2), x(3)]
        u_m2 = [(x(1) - 2.0_dp*h_rho)**2, x(2), x(3)]
        call eval_cart(cs, u_p1, xc_p1)
        call eval_cart(cs, u_m1, xc_m1)
        call eval_cart(cs, u_p2, xc_p2)
        call eval_cart(cs, u_m2, xc_m2)
        e1 = (-xc_p2 + 8.0_dp*xc_p1 - 8.0_dp*xc_m1 + xc_m2)/(12.0_dp*h_rho)

        u_p1 = [x(1)**2, wrap_angle(x(2) + h_ang, twopi), x(3)]
        u_m1 = [x(1)**2, wrap_angle(x(2) - h_ang, twopi), x(3)]
        u_p2 = [x(1)**2, wrap_angle(x(2) + 2.0_dp*h_ang, twopi), x(3)]
        u_m2 = [x(1)**2, wrap_angle(x(2) - 2.0_dp*h_ang, twopi), x(3)]
        call eval_cart(cs, u_p1, xc_p1)
        call eval_cart(cs, u_m1, xc_m1)
        call eval_cart(cs, u_p2, xc_p2)
        call eval_cart(cs, u_m2, xc_m2)
        e2 = (-xc_p2 + 8.0_dp*xc_p1 - 8.0_dp*xc_m1 + xc_m2)/(12.0_dp*h_ang)

        u_p1 = [x(1)**2, x(2), wrap_angle(x(3) + h_ang, phi_period)]
        u_m1 = [x(1)**2, x(2), wrap_angle(x(3) - h_ang, phi_period)]
        u_p2 = [x(1)**2, x(2), wrap_angle(x(3) + 2.0_dp*h_ang, phi_period)]
        u_m2 = [x(1)**2, x(2), wrap_angle(x(3) - 2.0_dp*h_ang, phi_period)]
        call eval_cart(cs, u_p1, xc_p1)
        call eval_cart(cs, u_m1, xc_m1)
        call eval_cart(cs, u_p2, xc_p2)
        call eval_cart(cs, u_m2, xc_m2)
        e3 = (-xc_p2 + 8.0_dp*xc_p1 - 8.0_dp*xc_m1 + xc_m2)/(12.0_dp*h_ang)

        sqrtg = dot_product(e1, cross_product(e2, e3))

        g(1, 1) = dot_product(e1, e1)
        g(1, 2) = dot_product(e1, e2)
        g(1, 3) = dot_product(e1, e3)
        g(2, 1) = g(1, 2)
        g(2, 2) = dot_product(e2, e2)
        g(2, 3) = dot_product(e2, e3)
        g(3, 1) = g(1, 3)
        g(3, 2) = g(2, 3)
        g(3, 3) = dot_product(e3, e3)

        det = g(1, 1)*(g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2)) &
              - g(1, 2)*(g(2, 1)*g(3, 3) - g(2, 3)*g(3, 1)) &
              + g(1, 3)*(g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))

        ginv(1, 1) = (g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2))/det
        ginv(1, 2) = (g(1, 3)*g(3, 2) - g(1, 2)*g(3, 3))/det
        ginv(1, 3) = (g(1, 2)*g(2, 3) - g(1, 3)*g(2, 2))/det
        ginv(2, 1) = (g(2, 3)*g(3, 1) - g(2, 1)*g(3, 3))/det
        ginv(2, 2) = (g(1, 1)*g(3, 3) - g(1, 3)*g(3, 1))/det
        ginv(2, 3) = (g(1, 3)*g(2, 1) - g(1, 1)*g(2, 3))/det
        ginv(3, 1) = (g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))/det
        ginv(3, 2) = (g(1, 2)*g(3, 1) - g(1, 1)*g(3, 2))/det
        ginv(3, 3) = (g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1))/det
    end subroutine metric_from_cart_fd

    subroutine eval_cart(cs, u, xcart)
        class(*), intent(in) :: cs
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: xcart(3)

        select type (cs)
        class is (coordinate_system_t)
            call cs%evaluate_cart(u, xcart)
        class default
            error stop 'eval_cart: invalid coordinate_system_t'
        end select
    end subroutine eval_cart

    subroutine compute_hcurl(sqrtg, dh, hcurl)
        real(dp), intent(in) :: sqrtg
        real(dp), intent(in) :: dh(3, 3)
        real(dp), intent(out) :: hcurl(3)

        if (abs(sqrtg) <= 0.0_dp) error stop 'compute_hcurl: sqrtg must be nonzero'

        hcurl(1) = (dh(2, 3) - dh(3, 2))/sqrtg
        hcurl(2) = (dh(3, 1) - dh(1, 3))/sqrtg
        hcurl(3) = (dh(1, 2) - dh(2, 1))/sqrtg
    end subroutine compute_hcurl

    subroutine check_scalar(name, x, a, b, tol_abs, tol_rel, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3), a, b, tol_abs, tol_rel
        integer, intent(inout) :: n_failed
        real(dp) :: abs_err

        abs_err = abs(a - b)
        if (abs_err > tol_abs + tol_rel*abs(b)) then
            print *, 'FAIL: ', trim(name), ' at x=(rho,th,ph)=', x
            print *, '  got=', a, ' expected=', b, ' abs=', abs_err
            n_failed = n_failed + 1
        end if
    end subroutine check_scalar

    subroutine check_vec(name, x, a, b, tol_abs, tol_rel, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3), a(3), b(3), tol_abs, tol_rel
        integer, intent(inout) :: n_failed
        integer :: k
        real(dp) :: abs_err

        do k = 1, 3
            abs_err = abs(a(k) - b(k))
            if (abs_err > tol_abs + tol_rel*abs(b(k))) then
                print *, 'FAIL: ', trim(name), '(', k, ') at x=(rho,th,ph)=', x
                print *, '  got=', a(k), ' expected=', b(k), ' abs=', abs_err
                n_failed = n_failed + 1
            end if
        end do
    end subroutine check_vec

    pure function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

    pure real(dp) function wrap_angle(angle, period) result(wrapped)
        real(dp), intent(in) :: angle, period

        wrapped = modulo(angle, period)
        if (wrapped < 0.0_dp) wrapped = wrapped + period
    end function wrap_angle

end program test_magfie_refcoords_splines
