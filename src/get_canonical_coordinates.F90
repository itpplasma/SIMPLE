module exchange_get_cancoord_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    logical, public :: onlytheta
    real(dp), public :: vartheta_c, varphi_c, sqg, aiota, Bcovar_vartheta, &
        Bcovar_varphi, A_theta, A_phi, theta, Bctrvr_vartheta, Bctrvr_varphi

    !$omp threadprivate(onlytheta, vartheta_c, varphi_c, sqg, aiota)
    !$omp threadprivate(Bcovar_vartheta, Bcovar_varphi, A_theta, A_phi)
    !$omp threadprivate(theta, Bctrvr_vartheta, Bctrvr_varphi)

end module exchange_get_cancoord_mod


module get_can_sub
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spl_three_to_five_sub
    use stencil_utils
    use field, only: magnetic_field_t, vmec_field_t, create_vmec_field, field_clone
    use field_newton, only: newton_theta_from_canonical
    use interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                           construct_batch_splines_1d, construct_batch_splines_3d, &
                           evaluate_batch_splines_1d_der2, &
                           evaluate_batch_splines_1d_der3, &
                           evaluate_batch_splines_3d_der, &
                           evaluate_batch_splines_3d_der2, &
                           destroy_batch_splines_1d, destroy_batch_splines_3d

    implicit none
    private

    public :: get_canonical_coordinates, get_canonical_coordinates_with_field
    public :: splint_can_coord
    public :: can_to_vmec, vmec_to_can, vmec_to_cyl
    public :: deallocate_can_coord
    public :: reset_canflux_batch_splines

    ! Constants
    real(dp), parameter :: TWOPI = 2.0_dp*3.14159265358979_dp

    ! Module variable to store the field for use in subroutines
    class(magnetic_field_t), allocatable :: current_field
    !$omp threadprivate(current_field)

    ! Batch spline for A_phi (vector potential)
    type(BatchSplineData1D), save :: aphi_batch_spline
    logical, save :: aphi_batch_spline_ready = .false.

    ! Batch spline for G_c (generating function)
    type(BatchSplineData3D), save :: G_batch_spline
    logical, save :: G_batch_spline_ready = .false.

    ! Batch spline for sqg_c, B_vartheta_c, B_varphi_c (3 quantities)
    type(BatchSplineData3D), save :: sqg_Bt_Bp_batch_spline
    logical, save :: sqg_Bt_Bp_batch_spline_ready = .false.

contains


subroutine get_canonical_coordinates_with_field(field)
    implicit none

    class(magnetic_field_t), intent(in) :: field

    ! Store field in module variable for use in nested subroutines
    call field_clone(field, current_field)

    call reset_canflux_batch_splines

    ! Call the actual implementation
    call get_canonical_coordinates_impl

end subroutine get_canonical_coordinates_with_field


subroutine get_canonical_coordinates
    ! Backward compatibility wrapper - uses VMEC field by default
    type(vmec_field_t) :: vmec_field

    call create_vmec_field(vmec_field)
    call get_canonical_coordinates_with_field(vmec_field)
end subroutine get_canonical_coordinates


subroutine reset_canflux_batch_splines
    if (aphi_batch_spline_ready) then
        call destroy_batch_splines_1d(aphi_batch_spline)
        aphi_batch_spline_ready = .false.
    end if
    if (G_batch_spline_ready) then
        call destroy_batch_splines_3d(G_batch_spline)
        G_batch_spline_ready = .false.
    end if
    if (sqg_Bt_Bp_batch_spline_ready) then
        call destroy_batch_splines_3d(sqg_Bt_Bp_batch_spline)
        sqg_Bt_Bp_batch_spline_ready = .false.
    end if
end subroutine reset_canflux_batch_splines


subroutine get_canonical_coordinates_impl
    use canonical_coordinates_mod, only: ns_c, n_theta_c, n_phi_c, &
                                         hs_c, h_theta_c, h_phi_c, &
                                         ns_s_c, ns_tp_c, &
                                         nh_stencil, G_c, sqg_c, &
                                         B_vartheta_c, B_varphi_c
    use vector_potentail_mod, only: ns, hs
    use exchange_get_cancoord_mod, only: vartheta_c, varphi_c, sqg, aiota, &
                                         Bcovar_vartheta, Bcovar_varphi, &
                                         onlytheta
    use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, ns_s, ns_tp
    use odeint_allroutines_sub, only: odeint_allroutines

    implicit none

    real(dp), parameter :: relerr = 1.0e-10_dp
    integer :: i_theta, i_phi, i_sten, ndim, is_beg
    integer, dimension(:), allocatable :: ipoi_t, ipoi_p
    real(dp), dimension(:), allocatable :: y, dy
    real(dp) :: dstencil_theta(-nh_stencil:nh_stencil), &
        dstencil_phi(-nh_stencil:nh_stencil)

    real(dp) :: r, r1, r2, G_beg, dG_c_dt, dG_c_dp
    integer :: is
    integer :: i_ctr

    ns_c = ns
    n_theta_c = n_theta
    n_phi_c = n_phi
    h_theta_c = h_theta
    h_phi_c = h_phi
    hs_c = hs

    ! Initialize derivative stencils using stencil_utils module
    call init_derivative_stencil(nh_stencil, h_theta_c, dstencil_theta)
    call init_derivative_stencil(nh_stencil, h_phi_c, dstencil_phi)

    allocate(ipoi_t(1 - nh_stencil:n_theta_c + nh_stencil))
    allocate(ipoi_p(1 - nh_stencil:n_phi_c + nh_stencil))

    do i_theta = 1, n_theta_c
        ipoi_t(i_theta) = i_theta
    end do

    do i_phi = 1, n_phi_c
        ipoi_p(i_phi) = i_phi
    end do

    do i_sten = 1, nh_stencil
        ipoi_t(1 - i_sten) = ipoi_t(n_theta - i_sten)
        ipoi_t(n_theta_c + i_sten) = ipoi_t(1 + i_sten)
        ipoi_p(1 - i_sten) = ipoi_p(n_phi_c - i_sten)
        ipoi_p(n_phi_c + i_sten) = ipoi_p(1 + i_sten)
    end do

    allocate(G_c(ns_c, n_theta_c, n_phi_c))
    allocate(sqg_c(ns_c, n_theta_c, n_phi_c))
    allocate(B_vartheta_c(ns_c, n_theta_c, n_phi_c))
    allocate(B_varphi_c(ns_c, n_theta_c, n_phi_c))

    onlytheta = .false.
    ndim = 1
    is_beg = 1
    G_beg = 1.0e-8_dp

    i_ctr = 0
!$omp parallel private(y, dy, i_theta, i_phi, is, r1, r2, r, dG_c_dt, dG_c_dp)
!$omp critical
    allocate(y(ndim), dy(ndim))
!$omp end critical

!$omp do
    do i_theta = 1, n_theta_c
!$omp critical
        i_ctr = i_ctr + 1
        call print_progress('integrate ODE: ', i_ctr, n_theta_c)
!$omp end critical
        vartheta_c = h_theta_c*real(i_theta - 1, dp)
        do i_phi = 1, n_phi_c
            varphi_c = h_phi_c*real(i_phi - 1, dp)

            G_c(is_beg, i_theta, i_phi) = G_beg
            y(1) = G_beg

            do is = is_beg - 1, 2, -1
                r1 = hs_c*real(is, dp)
                r2 = hs_c*real(is - 1, dp)

                call odeint_allroutines(y, ndim, r1, r2, relerr, rhs_cancoord)

                G_c(is, i_theta, i_phi) = y(1)
            end do

            y(1) = G_beg

            do is = is_beg + 1, ns_c
                r1 = hs_c*real(is - 2, dp)
                r2 = hs_c*real(is - 1, dp)
                if (is == 2) r1 = 1.0e-8_dp

                call odeint_allroutines(y, ndim, r1, r2, relerr, rhs_cancoord)

                G_c(is, i_theta, i_phi) = y(1)
            end do
        end do
    end do
!$omp end do

    i_ctr = 0
!$omp barrier
!$omp do
    do i_theta = 1, n_theta_c
!$omp critical
        i_ctr = i_ctr + 1
        call print_progress('compute components: ', i_ctr, n_theta_c)
!$omp end critical
        vartheta_c = h_theta_c*real(i_theta - 1, dp)
        do i_phi = 1, n_phi_c
            varphi_c = h_phi_c*real(i_phi - 1, dp)

            do is = 2, ns_c
                r = hs_c*real(is - 1, dp)
                y(1) = G_c(is, i_theta, i_phi)

                call rhs_cancoord(r, y, dy)

                dG_c_dt = sum(dstencil_theta*G_c(is, &
                    ipoi_t(i_theta - nh_stencil:i_theta + nh_stencil), i_phi))
                dG_c_dp = sum(dstencil_phi*G_c(is, i_theta, &
                    ipoi_p(i_phi - nh_stencil:i_phi + nh_stencil)))
                sqg_c(is, i_theta, i_phi) = sqg*(1.0_dp + aiota*dG_c_dt + dG_c_dp)
                B_vartheta_c(is, i_theta, i_phi) = Bcovar_vartheta + &
                    (aiota*Bcovar_vartheta + Bcovar_varphi)*dG_c_dt
                B_varphi_c(is, i_theta, i_phi) = Bcovar_varphi + &
                    (aiota*Bcovar_vartheta + Bcovar_varphi)*dG_c_dp
            end do
            ! Extrapolate on-axis point (is=1) with parabola
            sqg_c(1, i_theta, i_phi) = 3.0_dp*(sqg_c(2, i_theta, i_phi) &
                - sqg_c(3, i_theta, i_phi)) + sqg_c(4, i_theta, i_phi)
            B_vartheta_c(1, i_theta, i_phi) = 0.0_dp
            B_varphi_c(1, i_theta, i_phi) = 3.0_dp*(B_varphi_c(2, i_theta, i_phi) &
                - B_varphi_c(3, i_theta, i_phi)) + B_varphi_c(4, i_theta, i_phi)
        end do
    end do
!$omp end do

!$omp critical
    deallocate(y, dy)
!$omp end critical
!$omp end parallel

    ns_s_c = ns_s
    ns_tp_c = ns_tp

    onlytheta = .true.

    ! Build batch splines from computed grids
    call build_canflux_aphi_batch_spline
    call build_canflux_G_batch_spline
    call build_canflux_sqg_Bt_Bp_batch_spline

    deallocate(ipoi_t, ipoi_p, sqg_c, B_vartheta_c, B_varphi_c, G_c)

end subroutine get_canonical_coordinates_impl


subroutine build_canflux_aphi_batch_spline
    use vector_potentail_mod, only: ns, hs, sA_phi
    use new_vmec_stuff_mod, only: ns_A

    integer :: order

    if (aphi_batch_spline_ready) then
        call destroy_batch_splines_1d(aphi_batch_spline)
        aphi_batch_spline_ready = .false.
    end if

    order = ns_A
    if (order < 3 .or. order > 5) then
        error stop "build_canflux_aphi_batch_spline: spline order must be 3..5"
    end if

    aphi_batch_spline%order = order
    aphi_batch_spline%num_points = ns
    aphi_batch_spline%periodic = .false.
    aphi_batch_spline%x_min = 0.0_dp
    aphi_batch_spline%h_step = hs
    aphi_batch_spline%num_quantities = 1

    allocate(aphi_batch_spline%coeff(1, 0:order, ns))
    aphi_batch_spline%coeff(1, 0:order, :) = sA_phi(1:order + 1, :)

    aphi_batch_spline_ready = .true.
end subroutine build_canflux_aphi_batch_spline


subroutine build_canflux_G_batch_spline
    use canonical_coordinates_mod, only: ns_c, n_theta_c, n_phi_c, &
                                         hs_c, h_theta_c, h_phi_c, &
                                         ns_s_c, ns_tp_c, G_c

    integer :: order(3)
    real(dp) :: x_min(3), x_max(3)
    logical :: periodic(3)
    real(dp), allocatable :: y_batch(:, :, :, :)

    if (G_batch_spline_ready) then
        call destroy_batch_splines_3d(G_batch_spline)
        G_batch_spline_ready = .false.
    end if

    order = [ns_s_c, ns_tp_c, ns_tp_c]
    if (any(order < 3) .or. any(order > 5)) then
        error stop "build_canflux_G_batch_spline: spline order must be 3..5"
    end if

    x_min = [0.0_dp, 0.0_dp, 0.0_dp]
    x_max(1) = hs_c*real(ns_c - 1, dp)
    x_max(2) = h_theta_c*real(n_theta_c - 1, dp)
    x_max(3) = h_phi_c*real(n_phi_c - 1, dp)

    periodic = [.false., .true., .true.]

    allocate(y_batch(ns_c, n_theta_c, n_phi_c, 1))
    y_batch(:, :, :, 1) = G_c

    call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, &
                                    G_batch_spline)
    G_batch_spline_ready = .true.
    deallocate(y_batch)
end subroutine build_canflux_G_batch_spline


subroutine build_canflux_sqg_Bt_Bp_batch_spline
    use canonical_coordinates_mod, only: ns_c, n_theta_c, n_phi_c, &
                                         hs_c, h_theta_c, h_phi_c, &
                                         ns_s_c, ns_tp_c, &
                                         sqg_c, B_vartheta_c, B_varphi_c

    integer :: order(3)
    real(dp) :: x_min(3), x_max(3)
    logical :: periodic(3)
    real(dp), allocatable :: y_batch(:, :, :, :)

    if (sqg_Bt_Bp_batch_spline_ready) then
        call destroy_batch_splines_3d(sqg_Bt_Bp_batch_spline)
        sqg_Bt_Bp_batch_spline_ready = .false.
    end if

    order = [ns_s_c, ns_tp_c, ns_tp_c]
    if (any(order < 3) .or. any(order > 5)) then
        error stop "build_canflux_sqg_Bt_Bp_batch_spline: spline order must be 3..5"
    end if

    x_min = [0.0_dp, 0.0_dp, 0.0_dp]
    x_max(1) = hs_c*real(ns_c - 1, dp)
    x_max(2) = h_theta_c*real(n_theta_c - 1, dp)
    x_max(3) = h_phi_c*real(n_phi_c - 1, dp)

    periodic = [.false., .true., .true.]

    allocate(y_batch(ns_c, n_theta_c, n_phi_c, 3))
    y_batch(:, :, :, 1) = sqg_c
    y_batch(:, :, :, 2) = B_vartheta_c
    y_batch(:, :, :, 3) = B_varphi_c

    call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, &
                                    sqg_Bt_Bp_batch_spline)
    sqg_Bt_Bp_batch_spline_ready = .true.
    deallocate(y_batch)
end subroutine build_canflux_sqg_Bt_Bp_batch_spline


subroutine rhs_cancoord(r, y, dy)
    use exchange_get_cancoord_mod, only: vartheta_c, varphi_c, sqg, aiota, &
                                         Bcovar_vartheta, Bcovar_varphi, &
                                         theta, onlytheta
    use spline_vmec_sub
#ifdef GVEC_AVAILABLE
    use vmec_field_adapter
#else
    use vmec_field_eval
#endif

    implicit none

    real(dp), intent(in) :: r
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: dy(:)

    real(dp), parameter :: epserr = 1.0e-14_dp
    integer :: iter
    real(dp) :: s, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, &
        alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, Bcovar_r
    logical :: converged

    real(dp) :: vartheta, daiota_ds, deltheta

    s = r**2

    if (allocated(current_field)) then
        call vmec_iota_interpolate_with_field(current_field, s, aiota, daiota_ds)
    else
        call vmec_iota_interpolate(s, aiota, daiota_ds)
    end if

    vartheta = vartheta_c + aiota*y(1)
    varphi = varphi_c + y(1)

    ! Newton iteration to find field-specific theta from canonical theta
    if (allocated(current_field)) then
        theta = vartheta
        call newton_theta_from_canonical(current_field, s, vartheta, varphi, &
                                         theta, converged)
        if (.not. converged) then
            print *, 'WARNING: Newton iteration failed in rhs_cancoord'
        end if
    else
        theta = vartheta
        do iter = 1, 100
            call vmec_lambda_interpolate(s, theta, varphi, alam, dl_dt)
            deltheta = (vartheta - theta - alam)/(1.0_dp + dl_dt)
            theta = theta + deltheta
            if (abs(deltheta) < epserr) exit
        end do
    end if

    if (onlytheta) return

    if (allocated(current_field)) then
        call vmec_field_evaluate_with_field(current_field, s, theta, varphi, &
            A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
            sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
            Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    else
        call vmec_field_evaluate(s, theta, varphi, A_theta, A_phi, &
            dA_theta_ds, dA_phi_ds, aiota, &
            sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
            Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    end if

    dy(1) = -(Bcovar_r + daiota_ds*Bcovar_vartheta*y(1)) / &
            (aiota*Bcovar_vartheta + Bcovar_varphi)
    dy(1) = 2.0_dp*r*dy(1)

end subroutine rhs_cancoord


subroutine print_progress(message, progress, total)
    character(*), intent(in) :: message
    integer, intent(in) :: progress, total

    write(*, '(A, I4, A, I4)', advance='no') message, progress, ' of ', total

    if (progress < total) then
        write(*, '(A)', advance="no") char(13)
    else
        write(*, *)
    end if
end subroutine print_progress


subroutine splint_can_coord(fullset, mode_secders, r, vartheta_c, varphi_c, &
                            A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                            d2A_phi_dr2, d3A_phi_dr3, &
                            sqg_c, dsqg_c_dr, dsqg_c_dt, dsqg_c_dp, &
                            B_vartheta_c, dB_vartheta_c_dr, &
                            dB_vartheta_c_dt, dB_vartheta_c_dp, &
                            B_varphi_c, dB_varphi_c_dr, &
                            dB_varphi_c_dt, dB_varphi_c_dp, &
                            d2sqg_rr, d2sqg_rt, d2sqg_rp, &
                            d2sqg_tt, d2sqg_tp, d2sqg_pp, &
                            d2bth_rr, d2bth_rt, d2bth_rp, &
                            d2bth_tt, d2bth_tp, d2bth_pp, &
                            d2bph_rr, d2bph_rt, d2bph_rp, &
                            d2bph_tt, d2bph_tp, d2bph_pp, G_c)

    use vector_potentail_mod, only: torflux
    use new_vmec_stuff_mod, only: nper
    use chamb_mod, only: rnegflag
    use diag_mod, only: icounter

    implicit none

    logical, intent(in) :: fullset
    integer, intent(in) :: mode_secders
    real(dp), intent(in) :: r
    real(dp), intent(in) :: vartheta_c, varphi_c

    real(dp), intent(out) :: A_phi, A_theta, dA_phi_dr, dA_theta_dr
    real(dp), intent(out) :: d2A_phi_dr2, d3A_phi_dr3
    real(dp), intent(out) :: sqg_c, dsqg_c_dr, dsqg_c_dt, dsqg_c_dp
    real(dp), intent(out) :: B_vartheta_c, dB_vartheta_c_dr
    real(dp), intent(out) :: dB_vartheta_c_dt, dB_vartheta_c_dp
    real(dp), intent(out) :: B_varphi_c, dB_varphi_c_dr
    real(dp), intent(out) :: dB_varphi_c_dt, dB_varphi_c_dp
    real(dp), intent(out) :: d2sqg_rr, d2sqg_rt, d2sqg_rp
    real(dp), intent(out) :: d2sqg_tt, d2sqg_tp, d2sqg_pp
    real(dp), intent(out) :: d2bth_rr, d2bth_rt, d2bth_rp
    real(dp), intent(out) :: d2bth_tt, d2bth_tp, d2bth_pp
    real(dp), intent(out) :: d2bph_rr, d2bph_rt, d2bph_rp
    real(dp), intent(out) :: d2bph_tt, d2bph_tp, d2bph_pp
    real(dp), intent(out) :: G_c

    real(dp) :: r_eval
    real(dp) :: rho_tor, drhods, drhods2, d2rhods2m
    real(dp) :: x_eval(3)
    real(dp) :: y_eval(3), dy_eval(3, 3), d2y_eval(6, 3)
    real(dp) :: y_G(1), dy_G(3, 1)
    real(dp) :: y1d(1), dy1d(1), d2y1d(1)
    real(dp) :: theta_wrapped, phi_wrapped
    real(dp) :: qua, dqua_dr, dqua_dt, dqua_dp
    real(dp) :: d2qua_dr2, d2qua_drdt, d2qua_drdp
    real(dp) :: d2qua_dt2, d2qua_dtdp, d2qua_dp2

!$omp atomic
    icounter = icounter + 1
    r_eval = r
    if (r_eval <= 0.0_dp) then
        rnegflag = .true.
        r_eval = abs(r_eval)
    end if

    A_theta = torflux*r_eval
    dA_theta_dr = torflux

    ! Interpolate A_phi using batch spline (1D)
    if (.not. aphi_batch_spline_ready) then
        error stop "splint_can_coord: Aphi batch spline not initialized"
    end if

    if (mode_secders > 0) then
        ! Need third derivative - use der3 which computes all derivatives in one pass
        block
            real(dp) :: d3y1d(1)
            call evaluate_batch_splines_1d_der3(aphi_batch_spline, r_eval, &
                                                y1d, dy1d, d2y1d, d3y1d)
            d3A_phi_dr3 = d3y1d(1)
        end block
    else
        call evaluate_batch_splines_1d_der2(aphi_batch_spline, r_eval, &
                                            y1d, dy1d, d2y1d)
        d3A_phi_dr3 = 0.0_dp
    end if
    A_phi = y1d(1)
    dA_phi_dr = dy1d(1)
    d2A_phi_dr2 = d2y1d(1)

    ! Prepare coordinates for 3D interpolation
    rho_tor = sqrt(r_eval)
    theta_wrapped = modulo(vartheta_c, TWOPI)
    phi_wrapped = modulo(varphi_c, TWOPI/real(nper, dp))

    x_eval(1) = rho_tor
    x_eval(2) = theta_wrapped
    x_eval(3) = phi_wrapped

    ! Chain rule coefficients for rho -> s conversion
    ! rho = sqrt(s), drho/ds = 0.5/rho, d2rho/ds2 = -0.25/rho^3
    drhods = 0.5_dp/rho_tor
    drhods2 = drhods**2
    d2rhods2m = drhods2/rho_tor  ! -d2rho/ds2 (negated for chain rule)

    ! Interpolate G if needed
    if (fullset) then
        if (.not. G_batch_spline_ready) then
            error stop "splint_can_coord: G batch spline not initialized"
        end if
        call evaluate_batch_splines_3d_der(G_batch_spline, x_eval, y_G, dy_G)
        G_c = y_G(1)
    else
        G_c = 0.0_dp
    end if

    ! Interpolate sqg, B_vartheta, B_varphi (3 quantities)
    if (.not. sqg_Bt_Bp_batch_spline_ready) then
        error stop "splint_can_coord: sqg_Bt_Bp batch spline not initialized"
    end if

    if (mode_secders == 2) then
        call evaluate_batch_splines_3d_der2(sqg_Bt_Bp_batch_spline, x_eval, &
                                            y_eval, dy_eval, d2y_eval)

        ! Extract sqg_c (quantity 1)
        qua = y_eval(1)
        dqua_dr = dy_eval(1, 1)
        dqua_dt = dy_eval(2, 1)
        dqua_dp = dy_eval(3, 1)
        d2qua_dr2 = d2y_eval(1, 1)
        d2qua_drdt = d2y_eval(2, 1)
        d2qua_drdp = d2y_eval(3, 1)
        d2qua_dt2 = d2y_eval(4, 1)
        d2qua_dtdp = d2y_eval(5, 1)
        d2qua_dp2 = d2y_eval(6, 1)

        d2qua_dr2 = d2qua_dr2*drhods2 - dqua_dr*d2rhods2m
        dqua_dr = dqua_dr*drhods
        d2qua_drdt = d2qua_drdt*drhods
        d2qua_drdp = d2qua_drdp*drhods

        sqg_c = qua
        dsqg_c_dr = dqua_dr
        dsqg_c_dt = dqua_dt
        dsqg_c_dp = dqua_dp
        d2sqg_rr = d2qua_dr2
        d2sqg_rt = d2qua_drdt
        d2sqg_rp = d2qua_drdp
        d2sqg_tt = d2qua_dt2
        d2sqg_tp = d2qua_dtdp
        d2sqg_pp = d2qua_dp2

        ! Extract B_vartheta_c (quantity 2)
        qua = y_eval(2)
        dqua_dr = dy_eval(1, 2)
        dqua_dt = dy_eval(2, 2)
        dqua_dp = dy_eval(3, 2)
        d2qua_dr2 = d2y_eval(1, 2)
        d2qua_drdt = d2y_eval(2, 2)
        d2qua_drdp = d2y_eval(3, 2)
        d2qua_dt2 = d2y_eval(4, 2)
        d2qua_dtdp = d2y_eval(5, 2)
        d2qua_dp2 = d2y_eval(6, 2)

        d2qua_dr2 = d2qua_dr2*drhods2 - dqua_dr*d2rhods2m
        dqua_dr = dqua_dr*drhods
        d2qua_drdt = d2qua_drdt*drhods
        d2qua_drdp = d2qua_drdp*drhods

        B_vartheta_c = qua
        dB_vartheta_c_dr = dqua_dr
        dB_vartheta_c_dt = dqua_dt
        dB_vartheta_c_dp = dqua_dp
        d2bth_rr = d2qua_dr2
        d2bth_rt = d2qua_drdt
        d2bth_rp = d2qua_drdp
        d2bth_tt = d2qua_dt2
        d2bth_tp = d2qua_dtdp
        d2bth_pp = d2qua_dp2

        ! Extract B_varphi_c (quantity 3)
        qua = y_eval(3)
        dqua_dr = dy_eval(1, 3)
        dqua_dt = dy_eval(2, 3)
        dqua_dp = dy_eval(3, 3)
        d2qua_dr2 = d2y_eval(1, 3)
        d2qua_drdt = d2y_eval(2, 3)
        d2qua_drdp = d2y_eval(3, 3)
        d2qua_dt2 = d2y_eval(4, 3)
        d2qua_dtdp = d2y_eval(5, 3)
        d2qua_dp2 = d2y_eval(6, 3)

        d2qua_dr2 = d2qua_dr2*drhods2 - dqua_dr*d2rhods2m
        dqua_dr = dqua_dr*drhods
        d2qua_drdt = d2qua_drdt*drhods
        d2qua_drdp = d2qua_drdp*drhods

        B_varphi_c = qua
        dB_varphi_c_dr = dqua_dr
        dB_varphi_c_dt = dqua_dt
        dB_varphi_c_dp = dqua_dp
        d2bph_rr = d2qua_dr2
        d2bph_rt = d2qua_drdt
        d2bph_rp = d2qua_drdp
        d2bph_tt = d2qua_dt2
        d2bph_tp = d2qua_dtdp
        d2bph_pp = d2qua_dp2

    else
        ! mode_secders == 0 or 1: only need first derivatives
        call evaluate_batch_splines_3d_der(sqg_Bt_Bp_batch_spline, x_eval, &
                                           y_eval, dy_eval)

        ! Extract sqg_c (quantity 1)
        sqg_c = y_eval(1)
        dsqg_c_dr = dy_eval(1, 1)*drhods
        dsqg_c_dt = dy_eval(2, 1)
        dsqg_c_dp = dy_eval(3, 1)

        ! Extract B_vartheta_c (quantity 2)
        B_vartheta_c = y_eval(2)
        dB_vartheta_c_dr = dy_eval(1, 2)*drhods
        dB_vartheta_c_dt = dy_eval(2, 2)
        dB_vartheta_c_dp = dy_eval(3, 2)

        ! Extract B_varphi_c (quantity 3)
        B_varphi_c = y_eval(3)
        dB_varphi_c_dr = dy_eval(1, 3)*drhods
        dB_varphi_c_dt = dy_eval(2, 3)
        dB_varphi_c_dp = dy_eval(3, 3)

        ! Zero out second derivatives if not needed
        d2sqg_rr = 0.0_dp
        d2sqg_rt = 0.0_dp
        d2sqg_rp = 0.0_dp
        d2sqg_tt = 0.0_dp
        d2sqg_tp = 0.0_dp
        d2sqg_pp = 0.0_dp
        d2bth_rr = 0.0_dp
        d2bth_rt = 0.0_dp
        d2bth_rp = 0.0_dp
        d2bth_tt = 0.0_dp
        d2bth_tp = 0.0_dp
        d2bth_pp = 0.0_dp
        d2bph_rr = 0.0_dp
        d2bph_rt = 0.0_dp
        d2bph_rp = 0.0_dp
        d2bph_tt = 0.0_dp
        d2bph_tp = 0.0_dp
        d2bph_pp = 0.0_dp

        if (mode_secders == 1) then
            ! Compute only d2/dr2 derivatives
            call evaluate_batch_splines_3d_der2(sqg_Bt_Bp_batch_spline, x_eval, &
                                                y_eval, dy_eval, d2y_eval)
            d2sqg_rr = d2y_eval(1, 1)*drhods2 - dy_eval(1, 1)*d2rhods2m
            d2bth_rr = d2y_eval(1, 2)*drhods2 - dy_eval(1, 2)*d2rhods2m
            d2bph_rr = d2y_eval(1, 3)*drhods2 - dy_eval(1, 3)*d2rhods2m
        end if
    end if

end subroutine splint_can_coord


subroutine can_to_vmec(r, vartheta_c_in, varphi_c_in, theta_vmec, varphi_vmec)
    use exchange_get_cancoord_mod, only: vartheta_c, varphi_c, theta

    implicit none

    real(dp), intent(in) :: r, vartheta_c_in, varphi_c_in
    real(dp), intent(out) :: theta_vmec, varphi_vmec

    logical :: fullset
    integer :: mode_secders
    real(dp) :: r_local
    real(dp) :: A_phi, A_theta, dA_phi_dr, dA_theta_dr, d2A_phi_dr2, d3A_phi_dr3
    real(dp) :: sqg_c, dsqg_c_dr, dsqg_c_dt, dsqg_c_dp
    real(dp) :: B_vartheta_c, dB_vartheta_c_dr, dB_vartheta_c_dt, dB_vartheta_c_dp
    real(dp) :: B_varphi_c, dB_varphi_c_dr, dB_varphi_c_dt, dB_varphi_c_dp
    real(dp) :: G_c
    real(dp) :: d2sqg_rr, d2sqg_rt, d2sqg_rp, d2sqg_tt, d2sqg_tp, d2sqg_pp
    real(dp) :: d2bth_rr, d2bth_rt, d2bth_rp, d2bth_tt, d2bth_tp, d2bth_pp
    real(dp) :: d2bph_rr, d2bph_rt, d2bph_rp, d2bph_tt, d2bph_tp, d2bph_pp
    real(dp), dimension(1) :: y, dy

    fullset = .true.
    mode_secders = 0
    r_local = r

    call splint_can_coord(fullset, mode_secders, r_local, vartheta_c_in, &
        varphi_c_in, A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
        d2A_phi_dr2, d3A_phi_dr3, &
        sqg_c, dsqg_c_dr, dsqg_c_dt, dsqg_c_dp, &
        B_vartheta_c, dB_vartheta_c_dr, dB_vartheta_c_dt, dB_vartheta_c_dp, &
        B_varphi_c, dB_varphi_c_dr, dB_varphi_c_dt, dB_varphi_c_dp, &
        d2sqg_rr, d2sqg_rt, d2sqg_rp, d2sqg_tt, d2sqg_tp, d2sqg_pp, &
        d2bth_rr, d2bth_rt, d2bth_rp, d2bth_tt, d2bth_tp, d2bth_pp, &
        d2bph_rr, d2bph_rt, d2bph_rp, d2bph_tt, d2bph_tp, d2bph_pp, G_c)

    vartheta_c = vartheta_c_in
    varphi_c = varphi_c_in
    y(1) = G_c

    ! Transform from r (toroidal flux) to rho_tor for ODE integration
    call rhs_cancoord(sqrt(r_local), y, dy)

    theta_vmec = theta
    varphi_vmec = varphi_c_in + G_c

end subroutine can_to_vmec


subroutine deallocate_can_coord
    call reset_canflux_batch_splines
end subroutine deallocate_can_coord


subroutine vmec_to_can(r, theta, varphi, vartheta_c, varphi_c)
    ! Input : r,theta,varphi      - VMEC coordinates
    ! Output: vartheta_c,varphi_c - canonical coordinates

    use spline_vmec_sub
    use new_vmec_stuff_mod, only: nper
    use vector_potentail_mod, only: torflux
    use chamb_mod, only: rnegflag
#ifdef GVEC_AVAILABLE
    use vmec_field_adapter
#else
    use vmec_field_eval
#endif

    implicit none

    real(dp), parameter :: epserr = 1.0e-14_dp
    integer, parameter :: niter = 100
    integer :: iter
    real(dp), intent(in) :: r, theta, varphi
    real(dp), intent(out) :: vartheta_c, varphi_c
    real(dp) :: delthe, delphi, alam, dl_dt, vartheta
    real(dp) :: rho_tor, x_eval(3), y_G(1), dy_G(3, 1)
    real(dp) :: G_c, dG_c_dt, dG_c_dp, aiota
    real(dp) :: ts, ps, dts_dtc, dts_dpc, dps_dtc, dps_dpc, det
    real(dp) :: y1d(1), dy1d(1), d2y1d(1), dA_phi_dr, dA_theta_dr
    real(dp) :: r_local

    r_local = r
    if (r_local <= 0.0_dp) then
        rnegflag = .true.
        r_local = abs(r_local)
    end if

    if (allocated(current_field)) then
        call vmec_lambda_interpolate_with_field(current_field, r_local, theta, &
                                                varphi, alam, dl_dt)
    else
        call vmec_lambda_interpolate(r_local, theta, varphi, alam, dl_dt)
    end if

    vartheta = theta + alam

    vartheta_c = vartheta
    varphi_c = varphi

    ! Get iota from A_phi interpolation
    dA_theta_dr = torflux
    call evaluate_batch_splines_1d_der2(aphi_batch_spline, r_local, y1d, &
                                        dy1d, d2y1d)
    dA_phi_dr = dy1d(1)
    aiota = -dA_phi_dr/dA_theta_dr

    do iter = 1, niter
        rho_tor = sqrt(r_local)
        x_eval(1) = rho_tor
        x_eval(2) = modulo(vartheta_c, TWOPI)
        x_eval(3) = modulo(varphi_c, TWOPI/real(nper, dp))

        call evaluate_batch_splines_3d_der(G_batch_spline, x_eval, y_G, dy_G)
        G_c = y_G(1)
        dG_c_dt = dy_G(2, 1)
        dG_c_dp = dy_G(3, 1)

        ts = vartheta_c + aiota*G_c - vartheta
        ps = varphi_c + G_c - varphi
        dts_dtc = 1.0_dp + aiota*dG_c_dt
        dts_dpc = aiota*dG_c_dp
        dps_dtc = dG_c_dt
        dps_dpc = 1.0_dp + dG_c_dp
        det = 1.0_dp + aiota*dG_c_dt + dG_c_dp

        delthe = (ps*dts_dpc - ts*dps_dpc)/det
        delphi = (ts*dps_dtc - ps*dts_dtc)/det

        vartheta_c = vartheta_c + delthe
        varphi_c = varphi_c + delphi
        if (abs(delthe) + abs(delphi) < epserr) exit
    end do

end subroutine vmec_to_can


subroutine vmec_to_cyl(s, theta, varphi, Rcyl, Zcyl)
    use spline_vmec_sub
#ifdef GVEC_AVAILABLE
    use vmec_field_adapter
#else
    use vmec_field_eval
#endif

    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: Rcyl, Zcyl

    real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
        R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp

    if (allocated(current_field)) then
        call vmec_data_interpolate_with_field(current_field, s, theta, varphi, &
            A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
            R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
            dl_ds, dl_dt, dl_dp)
    else
        call vmec_data_interpolate(s, theta, varphi, A_phi, A_theta, &
            dA_phi_ds, dA_theta_ds, aiota, &
            R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
            dl_ds, dl_dt, dl_dp)
    end if

    Rcyl = R
    Zcyl = Z
end subroutine vmec_to_cyl

end module get_can_sub
