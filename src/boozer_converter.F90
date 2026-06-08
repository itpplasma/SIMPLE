module boozer_sub
    use spl_three_to_five_sub
    use interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                           construct_batch_splines_1d, construct_batch_splines_3d, &
                           evaluate_batch_splines_1d_der2, &
                           evaluate_batch_splines_1d_der3, &
                           evaluate_batch_splines_3d_der, &
                           evaluate_batch_splines_3d_der2, &
                           destroy_batch_splines_1d, destroy_batch_splines_3d
    use field, only: magnetic_field_t, field_clone
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    ! Public API
    public :: get_boozer_coordinates, get_boozer_coordinates_with_field
    public :: splint_boozer_coord
    public :: vmec_to_boozer, boozer_to_vmec
    public :: delthe_delphi_BV
    public :: reset_boozer_batch_splines
    public :: load_boozer_from_chartmap
    public :: export_boozer_chartmap

    ! Constants
    real(dp), parameter :: TWOPI = 2.0_dp*3.14159265358979_dp

    ! Field storage for nested subroutine calls
    class(magnetic_field_t), allocatable :: current_field
!$omp threadprivate(current_field)

    ! Batch spline data for Bmod and B_r interpolation
    type(BatchSplineData3D), allocatable :: bmod_br_batch_spline
    logical, save :: bmod_br_batch_spline_ready = .false.
    real(dp), allocatable, save :: bmod_grid(:, :, :)
    real(dp), allocatable, save :: br_grid(:, :, :)

    ! Batch spline for A_phi (vector potential)
    type(BatchSplineData1D), allocatable :: aphi_batch_spline
    logical, save :: aphi_batch_spline_ready = .false.
    ! Canonical A_phi sampled natively on the Boozer rho grid in
    ! compute_boozer_data (a flux function, like s_Bcovar_tp_B). export_boozer_chartmap
    ! writes this directly so the chartmap A_phi shares the rho abscissa of
    ! B_theta/B_phi -- no resampling of the s-grid A_phi spline.
    real(dp), allocatable, save :: aphi_rho(:)

    ! Batch spline for B_theta, B_phi covariant components
    type(BatchSplineData1D), allocatable :: bcovar_tp_batch_spline
    logical, save :: bcovar_tp_batch_spline_ready = .false.

    ! Batch splines for coordinate transformations (VMEC <-> Boozer)
    type(BatchSplineData3D), save :: delt_delp_V_batch_spline
    logical, save :: delt_delp_V_batch_spline_ready = .false.
    real(dp), allocatable, save :: delt_delp_V_grid(:, :, :, :)

    type(BatchSplineData3D), save :: delt_delp_B_batch_spline
    logical, save :: delt_delp_B_batch_spline_ready = .false.
    real(dp), allocatable, save :: delt_delp_B_grid(:, :, :, :)

    type :: boozer_state_t
        real(dp) :: torflux = 0.0_dp
        integer :: nper = 1
        logical :: use_B_r = .false.
        integer :: bmod_br_num_quantities = 0
        logical :: aphi_over_rho = .false.
    end type boozer_state_t

    ! Device-accessible Boozer runtime state shared by host and OpenACC code.
    ! sync_boozer_state copies the libneo globals into this object before GPU use.
    type(boozer_state_t), save :: boozer_state
    public :: sync_boozer_state, boozer_state

    ! Module statics referenced by splint_boozer_coord inside OpenACC device
    ! kernels need device symbols. declare create + update device (in
    ! sync_boozer_state / spline construction) keeps them coherent.
    !$acc declare create(boozer_state)
    !$acc declare create(bmod_br_batch_spline, aphi_batch_spline, bcovar_tp_batch_spline)

contains

    !> Copy scalar Boozer field parameters from the libneo modules into the
    !> shared runtime state. Host-only; call after field init and before any
    !> OpenACC tracing kernel.
    subroutine sync_boozer_state()
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper
        use boozer_coordinates_mod, only: use_B_r
        boozer_state%torflux = torflux
        boozer_state%nper = nper
        boozer_state%use_B_r = use_B_r
        !$acc update device(boozer_state)
        !$acc update device(bmod_br_batch_spline, aphi_batch_spline, bcovar_tp_batch_spline)
    end subroutine sync_boozer_state

    !> Initialize Boozer coordinates using given magnetic field
    subroutine get_boozer_coordinates_with_field(field)

        class(magnetic_field_t), intent(in) :: field

        ! Store field in module variable for use in nested subroutines
        call field_clone(field, current_field)
        call reset_boozer_batch_splines

        ! Call the actual implementation
        call get_boozer_coordinates_impl

    end subroutine get_boozer_coordinates_with_field

    !> Initialize Boozer coordinates using VMEC field (backward compatibility)
    subroutine get_boozer_coordinates
        use field, only: vmec_field_t, create_vmec_field

        type(vmec_field_t) :: vmec_field
        call create_vmec_field(vmec_field)
        call get_boozer_coordinates_with_field(vmec_field)

    end subroutine get_boozer_coordinates

    subroutine get_boozer_coordinates_impl

        use vector_potentail_mod, only: ns, hs
        use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, ns_s, ns_tp
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B

        implicit none

        ns_s_B = ns_s
        ns_tp_B = ns_tp
        ns_B = ns
        n_theta_B = n_theta
        n_phi_B = n_phi

        hs_B = hs*real(ns - 1, dp)/real(ns_B - 1, dp)
        h_theta_B = h_theta*real(n_theta - 1, dp)/real(n_theta_B - 1, dp)
        h_phi_B = h_phi*real(n_phi - 1, dp)/real(n_phi_B - 1, dp)

        call compute_boozer_data

        call build_boozer_aphi_batch_spline
        call build_boozer_bcovar_tp_batch_spline
        call build_boozer_bmod_br_batch_spline
        call build_boozer_delt_delp_batch_splines

    end subroutine get_boozer_coordinates_impl

    subroutine splint_boozer_coord(r, vartheta_B, varphi_B, mode_secders, &
                                   A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                   d2A_phi_dr2, d3A_phi_dr3, &
                                   B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                                   B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                                   Bmod_B, dBmod_B, d2Bmod_B, &
                                   B_r, dB_r, d2B_r)

        implicit none
        !$acc routine seq

        integer, intent(in) :: mode_secders

        real(dp), intent(in) :: r, vartheta_B, varphi_B
        real(dp), intent(out) :: A_phi, A_theta, dA_phi_dr, dA_theta_dr
        real(dp), intent(out) :: d2A_phi_dr2, d3A_phi_dr3
        real(dp), intent(out) :: B_vartheta_B, dB_vartheta_B, d2B_vartheta_B
        real(dp), intent(out) :: B_varphi_B, dB_varphi_B, d2B_varphi_B
        real(dp), intent(out) :: Bmod_B, B_r
        real(dp), intent(out) :: dBmod_B(3), dB_r(3)
        real(dp), intent(out) :: d2Bmod_B(6), d2B_r(6)

        real(dp) :: r_eval, rho_tor, drhods, drhods2, d2rhods2m
        real(dp) :: qua, dqua_dr, dqua_dt, dqua_dp
        real(dp) :: d2qua_dr2, d2qua_drdt, d2qua_drdp, d2qua_dt2, &
                    d2qua_dtdp, d2qua_dp2
        real(dp) :: x_eval(3), y_eval(2), dy_eval(3, 2), d2y_eval(6, 2)
        real(dp) :: theta_wrapped, phi_wrapped
        real(dp) :: y1d(2), dy1d(2), d2y1d(2)
        real(dp) :: d3y1d(1)

        ! Negative r only occurs transiently inside the Newton solve; the
        ! symplectic integrator clamps the radial coordinate itself, so the
        ! abs() below is the only handling needed. The dodiag evaluation
        ! counter and the chamb_mod rnegflag are host-only diagnostics that are
        ! omitted here (threadprivate variables are not allowed in acc routines).
        r_eval = abs(r)

        A_theta = boozer_state%torflux*r_eval
        dA_theta_dr = boozer_state%torflux

        ! Interpolate A_phi. VMEC-derived Boozer data is tabulated on a uniform-s
        ! grid (evaluate at s directly). A chartmap tabulates A_phi on the file's
        ! uniform rho grid, so evaluate at rho and chain-rule to s, exactly like
        ! Bmod/B_theta/B_phi below. The chartmap branch uses local chain-rule
        ! coefficients so the VMEC-s branch keeps its original instruction order.
        ! Guards (error stop) are omitted so this routine is callable from OpenACC
        ! device kernels (error stop is not supported there).
        if (boozer_state%aphi_over_rho) then
            ! A_phi(rho) -> s-derivatives by the chain rule for rho = sqrt(s):
            !   rho'   = 1/(2 rho)              [drds]
            !   rho''  = -1/(4 rho^3)           [-d2rds2m]
            !   rho''' = 3/(8 rho^5)            [d3rds3]
            ! f'  = g' rho'
            ! f'' = g'' rho'^2 + g' rho''
            ! f'''= g''' rho'^3 + 3 g'' rho' rho'' + g' rho'''
            block
                real(dp) :: rho_a, drds, drds2, d2rds2m, drds3, d3rds3, d3y1d(1)
                rho_a = sqrt(r_eval)
                drds = 0.5_dp/rho_a
                drds2 = drds**2
                d2rds2m = drds2/rho_a  ! -d2rho/ds2
                if (mode_secders > 0) then
                    call evaluate_batch_splines_1d_der3(aphi_batch_spline, rho_a, &
                                                        y1d(1:1), dy1d(1:1), &
                                                        d2y1d(1:1), d3y1d)
                    drds3 = drds*drds2  ! (drho/ds)^3 = 1/(8 rho^3)
                    d3rds3 = 3.0_dp/(8.0_dp*rho_a**5)  ! d3rho/ds3
                    d3A_phi_dr3 = d3y1d(1)*drds3 &
                                  - 3.0_dp*d2y1d(1)*drds*d2rds2m &
                                  + dy1d(1)*d3rds3
                else
                    call evaluate_batch_splines_1d_der2(aphi_batch_spline, rho_a, &
                                                        y1d(1:1), dy1d(1:1), d2y1d(1:1))
                    d3A_phi_dr3 = 0.0_dp
                end if
                A_phi = y1d(1)
                dA_phi_dr = dy1d(1)*drds
                d2A_phi_dr2 = d2y1d(1)*drds2 - dy1d(1)*d2rds2m
            end block
        else
            if (mode_secders > 0) then
                ! Need third derivative - use der3 which computes all in one pass
                block
                    real(dp) :: d3y1d(1)
                    call evaluate_batch_splines_1d_der3(aphi_batch_spline, r_eval, &
                                                        y1d(1:1), dy1d(1:1), &
                                                        d2y1d(1:1), d3y1d)
                    d3A_phi_dr3 = d3y1d(1)
                end block
            else
                call evaluate_batch_splines_1d_der2(aphi_batch_spline, r_eval, &
                                                    y1d(1:1), dy1d(1:1), d2y1d(1:1))
                d3A_phi_dr3 = 0.0_dp
            end if
            A_phi = y1d(1)
            dA_phi_dr = dy1d(1)
            d2A_phi_dr2 = d2y1d(1)
        end if

        ! Interpolation of mod-B (and B_r if use_B_r)
        rho_tor = sqrt(r_eval)
        theta_wrapped = modulo(vartheta_B, TWOPI)
        phi_wrapped = modulo(varphi_B, TWOPI/real(boozer_state%nper, dp))

        x_eval(1) = rho_tor
        x_eval(2) = theta_wrapped
        x_eval(3) = phi_wrapped

        ! Chain rule coefficients for rho -> s conversion
        drhods = 0.5_dp/rho_tor
        drhods2 = drhods**2
        d2rhods2m = drhods2/rho_tor  ! -d2rho/ds2 (negative of second derivative)

        if (mode_secders == 2) then
            call evaluate_batch_splines_3d_der2(bmod_br_batch_spline, x_eval, &
                                                y_eval(1:boozer_state%bmod_br_num_quantities), &
                                                dy_eval(:, 1:boozer_state%bmod_br_num_quantities), &
                                                d2y_eval(:, 1:boozer_state%bmod_br_num_quantities))

            ! Extract Bmod (quantity 1)
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

            Bmod_B = qua

            dBmod_B(1) = dqua_dr
            dBmod_B(2) = dqua_dt
            dBmod_B(3) = dqua_dp

            d2Bmod_B(1) = d2qua_dr2
            d2Bmod_B(2) = d2qua_drdt
            d2Bmod_B(3) = d2qua_drdp
            d2Bmod_B(4) = d2qua_dt2
            d2Bmod_B(5) = d2qua_dtdp
            d2Bmod_B(6) = d2qua_dp2

            ! Extract B_r (quantity 2, if present)
            if (boozer_state%use_B_r) then
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

                B_r = qua*drhods

                dB_r(1) = dqua_dr*drhods - qua*d2rhods2m
                dB_r(2) = dqua_dt*drhods
                dB_r(3) = dqua_dp*drhods

                d2B_r(1) = d2qua_dr2*drhods - 2.0_dp*dqua_dr*d2rhods2m + &
                           qua*drhods*(3.0_dp/4.0_dp)/r_eval**2
                d2B_r(2) = d2qua_drdt*drhods - dqua_dt*d2rhods2m
                d2B_r(3) = d2qua_drdp*drhods - dqua_dp*d2rhods2m
                d2B_r(4) = d2qua_dt2*drhods
                d2B_r(5) = d2qua_dtdp*drhods
                d2B_r(6) = d2qua_dp2*drhods
            else
                B_r = 0.0_dp
                dB_r = 0.0_dp
                d2B_r = 0.0_dp
            end if
        else
            call evaluate_batch_splines_3d_der(bmod_br_batch_spline, x_eval, &
                                               y_eval(1:boozer_state%bmod_br_num_quantities), &
                                               dy_eval(:, 1:boozer_state%bmod_br_num_quantities))

            Bmod_B = y_eval(1)
            dBmod_B(1) = dy_eval(1, 1)*drhods
            dBmod_B(2) = dy_eval(2, 1)
            dBmod_B(3) = dy_eval(3, 1)

            d2Bmod_B = 0.0_dp

            if (mode_secders == 1) then
                call evaluate_batch_splines_3d_der2(bmod_br_batch_spline, x_eval, &
                                                    y_eval(1:boozer_state%bmod_br_num_quantities), &
                                                    dy_eval(:, &
                                                            1:boozer_state%bmod_br_num_quantities), &
                                                    d2y_eval(:, &
                                                             1:boozer_state%bmod_br_num_quantities))
                d2Bmod_B(1) = d2y_eval(1, 1)*drhods2 - dy_eval(1, 1)*d2rhods2m
            end if

            if (boozer_state%use_B_r) then
                qua = y_eval(2)
                dqua_dr = dy_eval(1, 2)
                dqua_dt = dy_eval(2, 2)
                dqua_dp = dy_eval(3, 2)

                dqua_dr = dqua_dr*drhods
                B_r = qua*drhods

                dB_r(1) = dqua_dr*drhods - qua*d2rhods2m
                dB_r(2) = dqua_dt*drhods
                dB_r(3) = dqua_dp*drhods

                d2B_r = 0.0_dp
                if (mode_secders == 1) then
                    d2qua_dr2 = d2y_eval(1, 2)*drhods2 - dy_eval(1, 2)*d2rhods2m
                    d2B_r(1) = d2qua_dr2*drhods - 2.0_dp*dqua_dr*d2rhods2m + &
                               qua*drhods*(3.0_dp/4.0_dp)/r_eval**2
                end if
            else
                B_r = 0.0_dp
                dB_r = 0.0_dp
                d2B_r = 0.0_dp
            end if
        end if

        ! Interpolation of B_\vartheta and B_\varphi (flux functions)
        call evaluate_batch_splines_1d_der2(bcovar_tp_batch_spline, rho_tor, y1d, &
                                            dy1d, d2y1d)
        B_vartheta_B = y1d(1)
        dB_vartheta_B = dy1d(1)
        B_varphi_B = y1d(2)
        dB_varphi_B = dy1d(2)
        dB_vartheta_B = dB_vartheta_B*drhods
        dB_varphi_B = dB_varphi_B*drhods
        if (mode_secders > 0) then
            d2B_vartheta_B = d2y1d(1)*drhods2 - dy1d(1)*d2rhods2m
            d2B_varphi_B = d2y1d(2)*drhods2 - dy1d(2)*d2rhods2m
        else
            d2B_vartheta_B = 0.0_dp
            d2B_varphi_B = 0.0_dp
        end if

    end subroutine splint_boozer_coord

!> Computes delta_vartheta = vartheta_B - theta_V and delta_varphi = varphi_B - varphi_V
    !> and their first derivatives over angles.
    !> isw=0: given as functions of VMEC coordinates (r, vartheta, varphi)
    !> isw=1: given as functions of Boozer coordinates (r, vartheta, varphi)
    subroutine delthe_delphi_BV(isw, r, vartheta, varphi, deltheta_BV, delphi_BV, &
                                ddeltheta_BV, ddelphi_BV)
        use boozer_coordinates_mod, only: use_del_tp_B
        use chamb_mod, only: rnegflag

        integer, intent(in) :: isw
        real(dp), intent(in) :: r, vartheta, varphi
        real(dp), intent(out) :: deltheta_BV, delphi_BV
        real(dp), dimension(2), intent(out) :: ddeltheta_BV, ddelphi_BV

        real(dp) :: rho_tor, x_eval(3), y_eval(2), dy_eval(3, 2)
        real(dp) :: r_local

        r_local = r
        if (r_local <= 0.0_dp) then
            rnegflag = .true.
            r_local = abs(r_local)
        end if

        rho_tor = sqrt(r_local)
        x_eval(1) = rho_tor
        x_eval(2) = vartheta
        x_eval(3) = varphi

        if (isw .eq. 0) then
            if (.not. delt_delp_V_batch_spline_ready) then
                error stop "delthe_delphi_BV: V batch spline not initialized"
            end if
            call evaluate_batch_splines_3d_der(delt_delp_V_batch_spline, x_eval, &
                                               y_eval, dy_eval)
        elseif (isw .eq. 1) then
            if (.not. use_del_tp_B) then
                print *, 'delthe_delphi_BV : Boozer data is not loaded'
                return
            end if
            if (.not. delt_delp_B_batch_spline_ready) then
                error stop "delthe_delphi_BV: B batch spline not initialized"
            end if
            call evaluate_batch_splines_3d_der(delt_delp_B_batch_spline, x_eval, &
                                               y_eval, dy_eval)
        else
            print *, 'delthe_delphi_BV : unknown value of switch isw'
            return
        end if

        deltheta_BV = y_eval(1)
        delphi_BV = y_eval(2)

        ddeltheta_BV(1) = dy_eval(2, 1)
        ddelphi_BV(1) = dy_eval(2, 2)

        ddeltheta_BV(2) = dy_eval(3, 1)
        ddelphi_BV(2) = dy_eval(3, 2)

    end subroutine delthe_delphi_BV

!> Convert VMEC coordinates (r, theta, varphi) to Boozer coordinates (vartheta_B,
!> varphi_B)
    subroutine vmec_to_boozer(r, theta, varphi, vartheta_B, varphi_B)
        use new_vmec_stuff_mod, only: nper

        real(dp), intent(in) :: r, theta, varphi
        real(dp), intent(out) :: vartheta_B, varphi_B
        real(dp) :: deltheta_BV, delphi_BV
        real(dp), dimension(2) :: ddeltheta_BV, ddelphi_BV

        call delthe_delphi_BV(0, r, theta, varphi, deltheta_BV, delphi_BV, &
                              ddeltheta_BV, ddelphi_BV)

        vartheta_B = modulo(theta + deltheta_BV, TWOPI)
        varphi_B = modulo(varphi + delphi_BV, TWOPI/real(nper, dp))

    end subroutine vmec_to_boozer

!> Convert Boozer coordinates (r, vartheta_B, varphi_B) to VMEC coordinates (theta,
!> varphi)
    subroutine boozer_to_vmec(r, vartheta_B, varphi_B, theta, varphi)
        use boozer_coordinates_mod, only: use_del_tp_B

        real(dp), intent(in) :: r, vartheta_B, varphi_B
        real(dp), intent(out) :: theta, varphi

        real(dp), parameter :: epserr = 1.0e-14_dp
        integer, parameter :: niter = 100

        integer :: iter
        real(dp) :: deltheta_BV, delphi_BV
        real(dp) :: f1, f2, f11, f12, f21, f22, delthe, delphi, det
        real(dp), dimension(2) :: ddeltheta_BV, ddelphi_BV

        if (use_del_tp_B) then

            call delthe_delphi_BV(1, r, vartheta_B, varphi_B, deltheta_BV, delphi_BV, &
                                  ddeltheta_BV, ddelphi_BV)

            theta = vartheta_B - deltheta_BV
            varphi = varphi_B - delphi_BV
        else
            theta = vartheta_B
            varphi = varphi_B
        end if

! Newton method:

        do iter = 1, niter

            call delthe_delphi_BV(0, r, theta, varphi, deltheta_BV, delphi_BV, &
                                  ddeltheta_BV, ddelphi_BV)

            f1 = theta + deltheta_BV - vartheta_B
            f2 = varphi + delphi_BV - varphi_B
            f11 = 1.0_dp + ddeltheta_BV(1)
            f12 = ddeltheta_BV(2)
            f21 = ddelphi_BV(1)
            f22 = 1.0_dp + ddelphi_BV(2)

            det = f11*f22 - f12*f21
            delthe = (f2*f12 - f1*f22)/det
            delphi = (f1*f21 - f2*f11)/det

            theta = theta + delthe
            varphi = varphi + delphi
            if (abs(delthe) + abs(delphi) .lt. epserr) exit
        end do

!  theta=modulo(theta,TWOPI)
!  varphi=modulo(varphi,TWOPI/real(nper, dp))

    end subroutine boozer_to_vmec

    subroutine compute_boozer_data
        ! Computes Boozer coordinate transformations and magnetic field data
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_Bcovar_tp_B, &
                                          use_B_r, use_del_tp_B
        use binsrc_sub, only: binsrc
        use plag_coeff_sub, only: plag_coeff
        use spline_vmec_sub
        use vmec_field_eval

        implicit none

        real(dp), parameter :: s_min = 1.0e-6_dp, rho_min = sqrt(s_min)

        integer :: i, i_rho, i_theta, i_phi, npoilag, nder, nshift
        integer :: ibeg, iend, nqua
        real(dp) :: s, theta, varphi, A_theta, A_phi
        real(dp) :: dA_theta_ds, dA_phi_ds, aiota
        real(dp) :: sqg, alam, dl_ds, dl_dt, dl_dp
        real(dp) :: Bctrvr_vartheta, Bctrvr_varphi
        real(dp) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi
        real(dp) :: Bcovar_vartheta_B, Bcovar_varphi_B
        real(dp) :: denomjac, G00, Gbeg, aper
        real(dp) :: per_theta, per_phi, gridcellnum
        real(dp), allocatable :: wint_t(:), wint_p(:), theta_V(:), theta_B(:)
        real(dp), allocatable :: phi_V(:), phi_B(:), aiota_arr(:), rho_tor(:)
        real(dp), allocatable :: Bcovar_theta_V(:, :), Bcovar_varphi_V(:, :)
        real(dp), allocatable :: bmod_Vg(:, :), alam_2D(:, :)
        real(dp), allocatable :: deltheta_BV_Vg(:, :), delphi_BV_Vg(:, :)
        real(dp), allocatable :: splcoe_t(:, :)
        real(dp), allocatable :: splcoe_p(:, :), coef(:, :)
        real(dp), allocatable :: perqua_t(:, :), perqua_p(:, :)
        real(dp), allocatable :: perqua_2D(:, :, :), Gfunc(:, :, :)
        real(dp), allocatable :: Bcovar_symfl(:, :, :, :)

        nqua = 6
        gridcellnum = real((n_theta_B - 1)*(n_phi_B - 1), dp)

        npoilag = ns_tp_B + 1
        nder = 0
        nshift = npoilag/2

        print *, 'Transforming to Boozer coordinates'

        if (use_B_r) then
            print *, 'B_r is computed'
        else
            print *, 'B_r is not computed'
        end if

        G00 = 0.0_dp

        allocate (rho_tor(ns_B))
        if (allocated(aphi_rho)) deallocate (aphi_rho)
        allocate (aphi_rho(ns_B))
        allocate (aiota_arr(1))
        allocate (Gfunc(1, 1, 1))
        allocate (Bcovar_symfl(1, 1, 1, 1))
        if (use_B_r) then
            deallocate (aiota_arr, Gfunc, Bcovar_symfl)
            allocate (aiota_arr(ns_B))
            allocate (Gfunc(ns_B, n_theta_B, n_phi_B))
            allocate (Bcovar_symfl(3, ns_B, n_theta_B, n_phi_B))
        end if

        allocate (Bcovar_theta_V(n_theta_B, n_phi_B))
        allocate (Bcovar_varphi_V(n_theta_B, n_phi_B))
        allocate (bmod_Vg(n_theta_B, n_phi_B))
        allocate (alam_2D(n_theta_B, n_phi_B))
        allocate (deltheta_BV_Vg(n_theta_B, n_phi_B))
        allocate (delphi_BV_Vg(n_theta_B, n_phi_B))
        allocate (wint_t(0:ns_tp_B), wint_p(0:ns_tp_B))
        allocate (coef(0:nder, npoilag))
        allocate (theta_V(2 - n_theta_B:2*n_theta_B - 1))
        allocate (theta_B(2 - n_theta_B:2*n_theta_B - 1))
        allocate (phi_V(2 - n_phi_B:2*n_phi_B - 1))
        allocate (phi_B(2 - n_phi_B:2*n_phi_B - 1))
        allocate (perqua_t(nqua, 2 - n_theta_B:2*n_theta_B - 1))
        allocate (perqua_p(nqua, 2 - n_phi_B:2*n_phi_B - 1))
        allocate (perqua_2D(nqua, n_theta_B, n_phi_B))

        allocate (splcoe_t(0:ns_tp_B, n_theta_B))
        allocate (splcoe_p(0:ns_tp_B, n_phi_B))

! allocate data arrays for Boozer data:
        if (.not. allocated(s_Bcovar_tp_B)) &
            allocate (s_Bcovar_tp_B(2, ns_s_B + 1, ns_B))

        ! Allocate module-level grids
        call ensure_grid_3d(bmod_grid, ns_B, n_theta_B, n_phi_B)
        if (use_B_r) call ensure_grid_3d(br_grid, ns_B, n_theta_B, n_phi_B)
        call ensure_grid_4d(delt_delp_V_grid, ns_B, n_theta_B, n_phi_B, 2)
        if (use_del_tp_B) call ensure_grid_4d(delt_delp_B_grid, ns_B, n_theta_B, &
                                              n_phi_B, 2)

        do i = 0, ns_tp_B
            wint_t(i) = h_theta_B**(i + 1)/real(i + 1, dp)
            wint_p(i) = h_phi_B**(i + 1)/real(i + 1, dp)
        end do

        ! Set theta_V and phi_V linear, with value 0 at index 1 and stepsize h.
        ! Then expand this in both directions beyond 1:n_theta_B.
        do i_theta = 1, n_theta_B
            theta_V(i_theta) = real(i_theta - 1, dp)*h_theta_B
        end do
        per_theta = real(n_theta_B - 1, dp)*h_theta_B
        theta_V(2 - n_theta_B:0) = theta_V(1:n_theta_B - 1) - per_theta
        theta_V(n_theta_B + 1:2*n_theta_B - 1) = theta_V(2:n_theta_B) + per_theta

        do i_phi = 1, n_phi_B
            phi_V(i_phi) = real(i_phi - 1, dp)*h_phi_B
        end do
        per_phi = real(n_phi_B - 1, dp)*h_phi_B
        phi_V(2 - n_phi_B:0) = phi_V(1:n_phi_B - 1) - per_phi
        phi_V(n_phi_B + 1:2*n_phi_B - 1) = phi_V(2:n_phi_B) + per_phi

        do i_rho = 1, ns_B
            rho_tor(i_rho) = max(real(i_rho - 1, dp)*hs_B, rho_min)
            s = rho_tor(i_rho)**2

            do i_theta = 1, n_theta_B
                theta = real(i_theta - 1, dp)*h_theta_B
                do i_phi = 1, n_phi_B
                    varphi = real(i_phi - 1, dp)*h_phi_B

                    if (allocated(current_field)) then
                        call vmec_field_evaluate_with_field(current_field, &
                                                            s, theta, varphi, &
                                                            A_theta, &
                                                            A_phi, &
                                                            dA_theta_ds, &
                                                            dA_phi_ds, &
                                                            aiota, &
                                                            sqg, alam, dl_ds, &
                                                            dl_dt, dl_dp, &
                                                            Bctrvr_vartheta, &
                                                            Bctrvr_varphi, &
                                                            Bcovar_r, Bcovar_vartheta, &
                                                            Bcovar_varphi)
                    else
                        call vmec_field_evaluate(s, theta, varphi, &
                                                 A_theta, A_phi, dA_theta_ds, &
                                                 dA_phi_ds, aiota, &
                                                 sqg, alam, dl_ds, &
                                                 dl_dt, dl_dp, &
                                                 Bctrvr_vartheta, &
                                                 Bctrvr_varphi, &
                                                 Bcovar_r, Bcovar_vartheta, &
                                                 Bcovar_varphi)
                    end if

                    alam_2D(i_theta, i_phi) = alam
                    bmod_Vg(i_theta, i_phi) = &
                        sqrt(Bctrvr_vartheta*Bcovar_vartheta &
                             + Bctrvr_varphi*Bcovar_varphi)
                    Bcovar_theta_V(i_theta, i_phi) = Bcovar_vartheta*(1.0_dp + dl_dt)
                    Bcovar_varphi_V(i_theta, i_phi) = &
                        Bcovar_varphi + Bcovar_vartheta*dl_dp
                    perqua_2D(4, i_theta, i_phi) = Bcovar_r
                    perqua_2D(5, i_theta, i_phi) = Bcovar_vartheta
                    perqua_2D(6, i_theta, i_phi) = Bcovar_varphi
                end do
            end do

! covariant components $B_\vartheta$ and $B_\varphi$ of Boozer coordinates:
            Bcovar_vartheta_B = sum(Bcovar_theta_V(2:n_theta_B, 2:n_phi_B))/gridcellnum
            Bcovar_varphi_B = sum(Bcovar_varphi_V(2:n_theta_B, 2:n_phi_B))/gridcellnum
            s_Bcovar_tp_B(1, 1, i_rho) = Bcovar_vartheta_B
            s_Bcovar_tp_B(2, 1, i_rho) = Bcovar_varphi_B
            ! A_phi is a flux function; store it natively on the rho grid for
            ! export_boozer_chartmap (same abscissa as B_theta/B_phi above).
            aphi_rho(i_rho) = A_phi

            denomjac = 1.0_dp/(aiota*Bcovar_vartheta_B + Bcovar_varphi_B)
            Gbeg = G00 + Bcovar_vartheta_B*denomjac*alam_2D(1, 1)

            splcoe_t(0, :) = Bcovar_theta_V(:, 1)

            call spl_per(ns_tp_B, n_theta_B, h_theta_B, splcoe_t)

            delphi_BV_Vg(1, 1) = 0.0_dp
            do i_theta = 1, n_theta_B - 1
                delphi_BV_Vg(i_theta + 1, 1) = &
                    delphi_BV_Vg(i_theta, 1) &
                    + sum(wint_t*splcoe_t(:, i_theta))
            end do
            ! Remove linear increasing component from delphi_BV_Vg
            aper = (delphi_BV_Vg(n_theta_B, 1) &
                    - delphi_BV_Vg(1, 1))/real(n_theta_B - 1, dp)
            do i_theta = 2, n_theta_B
                delphi_BV_Vg(i_theta, 1) = &
                    delphi_BV_Vg(i_theta, 1) - aper*real(i_theta - 1, dp)
            end do

            do i_theta = 1, n_theta_B
                splcoe_p(0, :) = Bcovar_varphi_V(i_theta, :)

                call spl_per(ns_tp_B, n_phi_B, h_phi_B, splcoe_p)

                do i_phi = 1, n_phi_B - 1
                    delphi_BV_Vg(i_theta, i_phi + 1) = &
                        delphi_BV_Vg(i_theta, i_phi) &
                        + sum(wint_p*splcoe_p(:, i_phi))
                end do
                aper = (delphi_BV_Vg(i_theta, n_phi_B) &
                        - delphi_BV_Vg(i_theta, 1))/real(n_phi_B - 1, dp)
                do i_phi = 2, n_phi_B
                    delphi_BV_Vg(i_theta, i_phi) = &
                        delphi_BV_Vg(i_theta, i_phi) &
                        - aper*real(i_phi - 1, dp)
                end do
            end do

! difference between Boozer and VMEC toroidal angle,
! $\Delta \varphi_{BV}=\varphi_B-\varphi=G$:
            delphi_BV_Vg = denomjac*delphi_BV_Vg + Gbeg
! difference between Boozer and VMEC poloidal angle,
! $\Delta \vartheta_{BV}=\vartheta_B-\theta$:
            deltheta_BV_Vg = aiota*delphi_BV_Vg + alam_2D

            delt_delp_V_grid(i_rho, :, :, 1) = deltheta_BV_Vg
            delt_delp_V_grid(i_rho, :, :, 2) = delphi_BV_Vg

! At this point, all quantities are specified on
! equidistant grid in VMEC angles $(\theta,\varphi)$

! Re-interpolate to equidistant grid in $(\vartheta_B,\varphi)$:

            do i_phi = 1, n_phi_B
                perqua_t(1, 1:n_theta_B) = deltheta_BV_Vg(:, i_phi)
                perqua_t(2, 1:n_theta_B) = delphi_BV_Vg(:, i_phi)
                perqua_t(3, 1:n_theta_B) = bmod_Vg(:, i_phi)
                perqua_t(4:6, 1:n_theta_B) = perqua_2D(4:6, :, i_phi)
                ! Extend range of theta values
                perqua_t(:, 2 - n_theta_B:0) = perqua_t(:, 1:n_theta_B - 1)
                perqua_t(:, n_theta_B + 1:2*n_theta_B - 1) = perqua_t(:, 2:n_theta_B)
                theta_B = theta_V + perqua_t(1, :)
                do i_theta = 1, n_theta_B

                    call binsrc(theta_B, 2 - n_theta_B, 2*n_theta_B - 1, &
                                theta_V(i_theta), i)

                    ibeg = i - nshift
                    iend = ibeg + ns_tp_B

                    call plag_coeff(npoilag, nder, theta_V(i_theta), &
                                    theta_B(ibeg:iend), coef)

                    perqua_2D(:, i_theta, i_phi) = matmul(perqua_t(:, ibeg:iend), &
                                                          coef(0, :))
                end do
            end do

! End re-interpolate to equidistant grid in $(\vartheta_B,\varphi)$

! Re-interpolate to equidistant grid in $(\vartheta_B,\varphi_B)$:

            do i_theta = 1, n_theta_B
                perqua_p(:, 1:n_phi_B) = perqua_2D(:, i_theta, :)
                perqua_p(:, 2 - n_phi_B:0) = perqua_p(:, 1:n_phi_B - 1)
                ! Extend range of phi values
                perqua_p(:, n_phi_B + 1:2*n_phi_B - 1) = perqua_p(:, 2:n_phi_B)
                phi_B = phi_V + perqua_p(2, :)
                do i_phi = 1, n_phi_B

                    call binsrc(phi_B, 2 - n_phi_B, 2*n_phi_B - 1, phi_V(i_phi), i)

                    ibeg = i - nshift
                    iend = ibeg + ns_tp_B

                    call plag_coeff(npoilag, nder, phi_V(i_phi), phi_B(ibeg:iend), coef)

                    perqua_2D(:, i_theta, i_phi) = matmul(perqua_p(:, ibeg:iend), &
                                                          coef(0, :))
                end do
            end do

            if (use_del_tp_B) then
                delt_delp_B_grid(i_rho, :, :, 1) = perqua_2D(1, :, :)
                delt_delp_B_grid(i_rho, :, :, 2) = perqua_2D(2, :, :)
            end if
            bmod_grid(i_rho, :, :) = perqua_2D(3, :, :)

! End re-interpolate to equidistant grid in $(\vartheta_B,\varphi_B)$

            if (use_B_r) then
                aiota_arr(i_rho) = aiota
                Gfunc(i_rho, :, :) = perqua_2D(2, :, :)
! covariant components $B_k$ in symmetry flux coordinates on equidistant grid of
! Boozer coordinates:
                Bcovar_symfl(:, i_rho, :, :) = perqua_2D(4:6, :, :)
            end if

        end do

        if (use_B_r) then
            call compute_br_from_symflux(rho_tor, aiota_arr, Gfunc, Bcovar_symfl)
            deallocate (aiota_arr, Gfunc, Bcovar_symfl)
        end if

        deallocate (Bcovar_theta_V, Bcovar_varphi_V, bmod_Vg, alam_2D, &
                    deltheta_BV_Vg, delphi_BV_Vg, &
                    wint_t, wint_p, coef, theta_V, theta_B, phi_V, phi_B, &
                    perqua_t, perqua_p, perqua_2D)

        print *, 'done'

    end subroutine compute_boozer_data

    !> Compute radial covariant magnetic field B_rho from symmetry flux coordinates
    subroutine compute_br_from_symflux(rho_tor, aiota_arr, Gfunc, Bcovar_symfl)
        use boozer_coordinates_mod, only: ns_B, n_phi_B
        use plag_coeff_sub, only: plag_coeff

        real(dp), intent(in) :: rho_tor(:)
        real(dp), intent(in) :: aiota_arr(:)
        real(dp), intent(in) :: Gfunc(:, :, :)
        real(dp), intent(in) :: Bcovar_symfl(:, :, :, :)

        integer, parameter :: NPOILAG = 5
        integer, parameter :: NDER = 1

        integer :: i_rho, i_phi, ibeg, iend, nshift
        real(dp) :: coef(0:NDER, NPOILAG)

        nshift = (NPOILAG - 1)/2

        do i_rho = 1, ns_B
            ibeg = i_rho - nshift
            iend = ibeg + NPOILAG - 1
            if (ibeg < 1) then
                ibeg = 1
                iend = ibeg + NPOILAG - 1
            else if (iend > ns_B) then
                iend = ns_B
                ibeg = iend - NPOILAG + 1
            end if

            call plag_coeff(NPOILAG, NDER, rho_tor(i_rho), rho_tor(ibeg:iend), coef)

            ! Compute B_rho (we spline covariant component B_rho instead of B_s)
            do i_phi = 1, n_phi_B
                br_grid(i_rho, :, i_phi) = &
                    2.0_dp*rho_tor(i_rho)*Bcovar_symfl(1, i_rho, :, i_phi) &
                    - matmul(coef(1, :)*aiota_arr(ibeg:iend), Gfunc(ibeg:iend, &
                                                                    :, i_phi)) &
                    *Bcovar_symfl(2, i_rho, :, i_phi) &
                    - matmul(coef(1, :), Gfunc(ibeg:iend, :, i_phi)) &
                    *Bcovar_symfl(3, i_rho, :, i_phi)
            end do
        end do

    end subroutine compute_br_from_symflux

    !> Ensure 3D grid is allocated with correct dimensions
    subroutine ensure_grid_3d(grid, n1, n2, n3)
        real(dp), allocatable, intent(inout) :: grid(:, :, :)
        integer, intent(in) :: n1, n2, n3

        if (.not. allocated(grid)) then
            allocate (grid(n1, n2, n3))
        else if (any(shape(grid) /= [n1, n2, n3])) then
            deallocate (grid)
            allocate (grid(n1, n2, n3))
        end if
    end subroutine ensure_grid_3d

    !> Ensure 4D grid is allocated with correct dimensions
    subroutine ensure_grid_4d(grid, n1, n2, n3, n4)
        real(dp), allocatable, intent(inout) :: grid(:, :, :, :)
        integer, intent(in) :: n1, n2, n3, n4

        if (.not. allocated(grid)) then
            allocate (grid(n1, n2, n3, n4))
        else if (any(shape(grid) /= [n1, n2, n3, n4])) then
            deallocate (grid)
            allocate (grid(n1, n2, n3, n4))
        end if
    end subroutine ensure_grid_4d

    subroutine reset_boozer_batch_splines
        boozer_state = boozer_state_t()
        !$acc update device(boozer_state)
        if (aphi_batch_spline_ready) then
            call destroy_batch_splines_1d(aphi_batch_spline)
            aphi_batch_spline_ready = .false.
        end if
        if (bcovar_tp_batch_spline_ready) then
            call destroy_batch_splines_1d(bcovar_tp_batch_spline)
            bcovar_tp_batch_spline_ready = .false.
        end if
        if (bmod_br_batch_spline_ready) then
            call destroy_batch_splines_3d(bmod_br_batch_spline)
            bmod_br_batch_spline_ready = .false.
            boozer_state%bmod_br_num_quantities = 0
        end if
        if (allocated(bmod_grid)) deallocate (bmod_grid)
        if (allocated(br_grid)) deallocate (br_grid)
        if (delt_delp_V_batch_spline_ready) then
            call destroy_batch_splines_3d(delt_delp_V_batch_spline)
            delt_delp_V_batch_spline_ready = .false.
        end if
        if (allocated(delt_delp_V_grid)) deallocate (delt_delp_V_grid)
        if (delt_delp_B_batch_spline_ready) then
            call destroy_batch_splines_3d(delt_delp_B_batch_spline)
            delt_delp_B_batch_spline_ready = .false.
        end if
        if (allocated(delt_delp_B_grid)) deallocate (delt_delp_B_grid)

        ! Descriptors are allocatable so that, under -gpu=mem:managed/separate,
        ! the module common is not emitted as a value-initialized device global
        ! (which trips an nvfortran NVVM codegen bug). Keep them allocated for the
        ! program lifetime; destroy_batch_splines_* above frees only %coeff.
        if (.not. allocated(aphi_batch_spline)) allocate (aphi_batch_spline)
        if (.not. allocated(bcovar_tp_batch_spline)) allocate (bcovar_tp_batch_spline)
        if (.not. allocated(bmod_br_batch_spline)) allocate (bmod_br_batch_spline)
    end subroutine reset_boozer_batch_splines

    subroutine load_boozer_from_chartmap(filename)
        !> Populate module-level Boozer batch splines from an extended chartmap
        !> NetCDF file, bypassing the VMEC-based compute_boozer_data path.
        use vector_potentail_mod, only: torflux, ns, hs
        use new_vmec_stuff_mod, only: nper, rmajor, ns_A, ns_s, ns_tp, &
                                      vmec_B_scale, vmec_RZ_scale
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, &
                                          n_phi_B, hs_B, h_theta_B, h_phi_B, &
                                          use_B_r, use_del_tp_B
        use boozer_chartmap_io, only: boozer_chartmap_data_t, read_boozer_chartmap

        character(len=*), intent(in) :: filename

        type(boozer_chartmap_data_t) :: d
        real(dp), allocatable :: y_aphi(:, :), y_bcovar(:, :), y_bmod(:, :, :, :)
        real(dp) :: s_min, s_max
        real(dp) :: b_scale, rz_scale, covar_scale, flux_scale
        integer :: spline_order
        integer :: order_3d(3)
        logical :: periodic_3d(3)
        real(dp) :: x_min_3d(3), x_max_3d(3)

        call reset_boozer_batch_splines

        ! Single shared parse: base-unit arrays on the endpoint-included field grid.
        call read_boozer_chartmap(filename, d)

        ! Apply the VMEC scaling knobs so a chartmap behaves like a VMEC run
        ! (matches boozer_chartmap_field_t and test_chartmap_scaling). Base files
        ! are exported at scale 1, so this is a no-op by default.
        b_scale = vmec_B_scale
        rz_scale = vmec_RZ_scale
        covar_scale = b_scale*rz_scale
        flux_scale = covar_scale*rz_scale
        d%A_phi = flux_scale*d%A_phi
        d%B_theta = covar_scale*d%B_theta
        d%B_phi = covar_scale*d%B_phi
        d%Bmod = b_scale*d%Bmod

        ! Set global parameters used by splint_boozer_coord
        torflux = flux_scale*d%torflux
        nper = d%nfp
        rmajor = d%rmajor*rz_scale

        ! Set boozer_coordinates_mod parameters
        ns_s_B = 5
        ns_tp_B = 5
        ns_B = d%n_rho
        n_theta_B = d%n_theta
        n_phi_B = d%n_phi
        hs_B = d%rho(2) - d%rho(1)
        h_theta_B = d%h_theta
        h_phi_B = d%h_phi
        use_B_r = .false.
        use_del_tp_B = .false.

        ! Set vector_potentail_mod parameters for A_phi spline
        s_min = d%rho_min**2
        s_max = d%rho_max**2
        ns = d%n_rho
        hs = (s_max - s_min) / real(d%n_rho - 1, dp)
        ns_A = 5

        ! Build A_phi batch spline over the file's uniform rho grid. The chartmap
        ! tabulates A_phi against rho (like B_theta/B_phi/Bmod), so splint_boozer_coord
        ! evaluates it at rho and chain-rules to s.
        spline_order = ns_A
        allocate (y_aphi(d%n_rho, 1))
        y_aphi(:, 1) = d%A_phi
        call construct_batch_splines_1d(d%rho(1), d%rho(d%n_rho), y_aphi, spline_order, &
                                        .false., aphi_batch_spline)
        aphi_batch_spline_ready = .true.
        boozer_state%aphi_over_rho = .true.
        deallocate (y_aphi)

        ! Build B_theta, B_phi batch spline over rho_tor
        spline_order = ns_s_B
        allocate (y_bcovar(d%n_rho, 2))
        y_bcovar(:, 1) = d%B_theta
        y_bcovar(:, 2) = d%B_phi
        call construct_batch_splines_1d(d%rho(1), d%rho(d%n_rho), y_bcovar, spline_order, &
                                        .false., bcovar_tp_batch_spline)
        bcovar_tp_batch_spline_ready = .true.
        deallocate (y_bcovar)

        ! Build Bmod 3D batch spline over (rho_tor, theta_B, phi_B)
        ! Uses endpoint-included field grid matching original VMEC Boozer splines
        order_3d = [ns_s_B, ns_tp_B, ns_tp_B]
        periodic_3d = [.false., .true., .true.]
        x_min_3d = [d%rho(1), 0.0_dp, 0.0_dp]
        x_max_3d = [d%rho(d%n_rho), h_theta_B * real(d%n_theta - 1, dp), &
                     h_phi_B * real(d%n_phi - 1, dp)]

        allocate (y_bmod(d%n_rho, d%n_theta, d%n_phi, 1))
        y_bmod(:, :, :, 1) = d%Bmod
        call construct_batch_splines_3d(x_min_3d, x_max_3d, y_bmod, order_3d, &
                                        periodic_3d, bmod_br_batch_spline)
        bmod_br_batch_spline_ready = .true.
        boozer_state%bmod_br_num_quantities = 1
        deallocate (y_bmod)

        print *, 'Loaded Boozer splines from chartmap: ', trim(filename)
        print *, '  nfp=', d%nfp, ' ns=', d%n_rho, ' ntheta_field=', &
                 d%n_theta, ' nphi_field=', d%n_phi
        print *, '  torflux=', torflux

        ! The chartmap loader builds the batch splines inline (not via
        ! build_boozer_bmod_br_batch_spline), so refresh the device-accessible
        ! field-parameter mirrors here too. splint_boozer_coord reads them on
        ! both host and device.
        call sync_boozer_state

    end subroutine load_boozer_from_chartmap

    subroutine export_boozer_chartmap(filename)
        !> Export Boozer coordinate data computed by get_boozer_coordinates()
        !> to an extended chartmap NetCDF file. Must be called after
        !> get_boozer_coordinates() and while VMEC splines are still active
        !> (needed for X, Y, Z geometry evaluation).
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper, rmajor
        use boozer_coordinates_mod, only: ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_Bcovar_tp_B
        use spline_vmec_sub, only: splint_vmec_data
        use netcdf

        character(len=*), intent(in) :: filename

        integer :: ncid, status
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: dim_theta_field, dim_zeta_field
        integer :: var_rho, var_theta, var_zeta
        integer :: var_x, var_y, var_z
        integer :: var_aphi, var_btheta, var_bphi, var_bmod, var_nfp
        integer :: i_rho, i_theta, i_phi
        integer :: n_theta_out, n_phi_out
        real(dp) :: rho_tor, s, theta_B, phi_B, theta_V, phi_V
        real(dp) :: R, Zval, alam
        real(dp) :: A_phi_dum, A_theta_dum, dA_phi_ds, dA_theta_ds, aiota
        real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
        real(dp) :: dl_ds, dl_dt, dl_dp
        real(dp), parameter :: rho_min = sqrt(1.0e-6_dp)
        real(dp), allocatable :: rho_arr(:), theta_arr(:), zeta_arr(:)
        real(dp), allocatable :: A_phi_arr(:), B_theta_arr(:), B_phi_arr(:)
        real(dp), allocatable :: x_arr(:, :, :), y_arr(:, :, :), z_arr(:, :, :)

        ! Chartmap geometry requires endpoint-excluded angular grids (libneo validator).
        ! Boozer field data uses endpoint-included grids (matching original splines).
        n_theta_out = n_theta_B - 1
        n_phi_out = n_phi_B - 1

        allocate (rho_arr(ns_B))
        allocate (theta_arr(n_theta_out), zeta_arr(n_phi_out))
        allocate (A_phi_arr(ns_B), B_theta_arr(ns_B), B_phi_arr(ns_B))
        allocate (x_arr(ns_B, n_theta_out, n_phi_out))
        allocate (y_arr(ns_B, n_theta_out, n_phi_out))
        allocate (z_arr(ns_B, n_theta_out, n_phi_out))

        ! Radial grid
        do i_rho = 1, ns_B
            rho_arr(i_rho) = max(real(i_rho - 1, dp) * hs_B, rho_min)
        end do
        ! Angular grids (endpoint excluded, for chartmap geometry)
        do i_theta = 1, n_theta_out
            theta_arr(i_theta) = real(i_theta - 1, dp) * h_theta_B
        end do
        do i_phi = 1, n_phi_out
            zeta_arr(i_phi) = real(i_phi - 1, dp) * h_phi_B
        end do

        ! Extract 1D profiles as genuine functions of the rho grid. A_phi, like
        ! B_theta/B_phi, was sampled natively on the rho grid in compute_boozer_data,
        ! so all three share the rho abscissa the loader reads them back on.
        if (.not. allocated(aphi_rho)) then
            error stop "export_boozer_chartmap: aphi_rho not populated"
        end if
        do i_rho = 1, ns_B
            A_phi_arr(i_rho) = aphi_rho(i_rho)
            B_theta_arr(i_rho) = s_Bcovar_tp_B(1, 1, i_rho)
            B_phi_arr(i_rho) = s_Bcovar_tp_B(2, 1, i_rho)
        end do

        ! Compute X, Y, Z geometry on the Boozer grid (endpoint-excluded)
        do i_phi = 1, n_phi_out
            do i_theta = 1, n_theta_out
                do i_rho = 1, ns_B
                    rho_tor = rho_arr(i_rho)
                    s = rho_tor**2
                    theta_B = theta_arr(i_theta)
                    phi_B = zeta_arr(i_phi)

                    ! Convert Boozer angles to VMEC angles
                    call boozer_to_vmec(s, theta_B, phi_B, theta_V, phi_V)

                    ! Evaluate VMEC geometry at (s, theta_V, phi_V)
                    call splint_vmec_data(s, theta_V, phi_V, &
                        A_phi_dum, A_theta_dum, dA_phi_ds, dA_theta_ds, aiota, &
                        R, Zval, alam, dR_ds, dR_dt, dR_dp, &
                        dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

                    x_arr(i_rho, i_theta, i_phi) = R * cos(phi_V)
                    y_arr(i_rho, i_theta, i_phi) = R * sin(phi_V)
                    z_arr(i_rho, i_theta, i_phi) = Zval
                end do
            end do
        end do

        ! Write NetCDF file
        status = nf90_create(trim(filename), nf90_clobber, ncid)
        call nc_assert(status, "create " // trim(filename))

        ! Dimensions: geometry uses endpoint-excluded grids (chartmap validator),
        ! field data uses endpoint-included grids (exact spline reproduction)
        call nc_assert(nf90_def_dim(ncid, "rho", ns_B, dim_rho), "def_dim rho")
        call nc_assert(nf90_def_dim(ncid, "theta", n_theta_out, dim_theta), &
                        "def_dim theta")
        call nc_assert(nf90_def_dim(ncid, "zeta", n_phi_out, dim_zeta), &
                        "def_dim zeta")
        call nc_assert(nf90_def_dim(ncid, "theta_field", n_theta_B, dim_theta_field), &
                        "def_dim theta_field")
        call nc_assert(nf90_def_dim(ncid, "zeta_field", n_phi_B, dim_zeta_field), &
                        "def_dim zeta_field")

        ! Coordinate variables
        call nc_assert(nf90_def_var(ncid, "rho", nf90_double, [dim_rho], var_rho), &
                        "def_var rho")
        call nc_assert(nf90_def_var(ncid, "theta", nf90_double, [dim_theta], &
                        var_theta), "def_var theta")
        call nc_assert(nf90_def_var(ncid, "zeta", nf90_double, [dim_zeta], &
                        var_zeta), "def_var zeta")

        ! Geometry (NF90 reverses dims: Fortran (rho,theta,zeta) -> NetCDF (zeta,theta,rho))
        call nc_assert(nf90_def_var(ncid, "x", nf90_double, &
                        [dim_rho, dim_theta, dim_zeta], var_x), "def_var x")
        call nc_assert(nf90_put_att(ncid, var_x, "units", "cm"), "att x units")
        call nc_assert(nf90_def_var(ncid, "y", nf90_double, &
                        [dim_rho, dim_theta, dim_zeta], var_y), "def_var y")
        call nc_assert(nf90_put_att(ncid, var_y, "units", "cm"), "att y units")
        call nc_assert(nf90_def_var(ncid, "z", nf90_double, &
                        [dim_rho, dim_theta, dim_zeta], var_z), "def_var z")
        call nc_assert(nf90_put_att(ncid, var_z, "units", "cm"), "att z units")

        ! Boozer field data
        call nc_assert(nf90_def_var(ncid, "A_phi", nf90_double, [dim_rho], var_aphi), &
                        "def_var A_phi")
        call nc_assert(nf90_def_var(ncid, "B_theta", nf90_double, [dim_rho], &
                        var_btheta), "def_var B_theta")
        call nc_assert(nf90_def_var(ncid, "B_phi", nf90_double, [dim_rho], &
                        var_bphi), "def_var B_phi")
        call nc_assert(nf90_def_var(ncid, "Bmod", nf90_double, &
                        [dim_rho, dim_theta_field, dim_zeta_field], var_bmod), &
                        "def_var Bmod")
        call nc_assert(nf90_def_var(ncid, "num_field_periods", nf90_int, var_nfp), &
                        "def_var nfp")

        ! Global attributes
        call nc_assert(nf90_put_att(ncid, nf90_global, "rho_convention", "rho_tor"), &
                        "att rho_convention")
        call nc_assert(nf90_put_att(ncid, nf90_global, "zeta_convention", "boozer"), &
                        "att zeta_convention")
        call nc_assert(nf90_put_att(ncid, nf90_global, "rho_lcfs", rho_arr(ns_B)), &
                        "att rho_lcfs")
        call nc_assert(nf90_put_att(ncid, nf90_global, "boozer_field", 1), &
                        "att boozer_field")
        call nc_assert(nf90_put_att(ncid, nf90_global, "torflux", torflux), &
                        "att torflux")
        ! Major radius (metres, base units) so the chartmap reader can restore
        ! new_vmec_stuff_mod::rmajor; stevvo/params_init need it for dtaumin.
        call nc_assert(nf90_put_att(ncid, nf90_global, "rmajor", rmajor), &
                        "att rmajor")

        call nc_assert(nf90_enddef(ncid), "enddef")

        ! Write data
        call nc_assert(nf90_put_var(ncid, var_rho, rho_arr), "put rho")
        call nc_assert(nf90_put_var(ncid, var_theta, theta_arr), "put theta")
        call nc_assert(nf90_put_var(ncid, var_zeta, zeta_arr), "put zeta")
        call nc_assert(nf90_put_var(ncid, var_x, x_arr), "put x")
        call nc_assert(nf90_put_var(ncid, var_y, y_arr), "put y")
        call nc_assert(nf90_put_var(ncid, var_z, z_arr), "put z")
        call nc_assert(nf90_put_var(ncid, var_aphi, A_phi_arr), "put A_phi")
        call nc_assert(nf90_put_var(ncid, var_btheta, B_theta_arr), "put B_theta")
        call nc_assert(nf90_put_var(ncid, var_bphi, B_phi_arr), "put B_phi")
        call nc_assert(nf90_put_var(ncid, var_bmod, bmod_grid), "put Bmod")
        call nc_assert(nf90_put_var(ncid, var_nfp, nper), "put nfp")

        call nc_assert(nf90_close(ncid), "close")

        print *, 'Exported Boozer chartmap to ', trim(filename)
        print *, '  nfp=', nper, ' ns=', ns_B, ' ntheta=', n_theta_out, &
                 ' nphi=', n_phi_out
        print *, '  torflux=', torflux

    contains

        subroutine nc_assert(stat, loc)
            integer, intent(in) :: stat
            character(len=*), intent(in) :: loc
            if (stat /= nf90_noerr) then
                print *, "export_boozer_chartmap: NetCDF error at ", trim(loc), &
                         ": ", trim(nf90_strerror(stat))
                error stop
            end if
        end subroutine nc_assert

    end subroutine export_boozer_chartmap

    subroutine build_boozer_aphi_batch_spline
        use vector_potentail_mod, only: ns, hs, sA_phi
        use new_vmec_stuff_mod, only: ns_A

        integer :: order

        if (aphi_batch_spline_ready) then
            call destroy_batch_splines_1d(aphi_batch_spline)
            aphi_batch_spline_ready = .false.
        end if

        order = ns_A
        if (order < 3 .or. order > 5) then
            error stop "build_boozer_aphi_batch_spline: spline order must be 3..5"
        end if

        aphi_batch_spline%order = order
        aphi_batch_spline%num_points = ns
        aphi_batch_spline%periodic = .false.
        aphi_batch_spline%x_min = 0.0_dp
        aphi_batch_spline%h_step = hs
        aphi_batch_spline%num_quantities = 1

        allocate (aphi_batch_spline%coeff(1, 0:order, ns))
        aphi_batch_spline%coeff(1, 0:order, :) = sA_phi(1:order + 1, :)

        aphi_batch_spline_ready = .true.
        boozer_state%aphi_over_rho = .false.  ! VMEC A_phi lives on the uniform-s grid
    end subroutine build_boozer_aphi_batch_spline

    subroutine build_boozer_bcovar_tp_batch_spline
        use boozer_coordinates_mod, only: ns_s_B, ns_B, hs_B, s_Bcovar_tp_B

        integer :: order
        real(dp) :: x_min, x_max
        real(dp), allocatable :: y_batch(:, :)

        if (bcovar_tp_batch_spline_ready) then
            call destroy_batch_splines_1d(bcovar_tp_batch_spline)
            bcovar_tp_batch_spline_ready = .false.
        end if

        order = ns_s_B
        if (order < 3 .or. order > 5) then
            error stop "build_boozer_bcovar_tp_batch_spline: spline order must be 3..5"
        end if

        x_min = 0.0_dp
        x_max = hs_B*real(ns_B - 1, dp)

        allocate (y_batch(ns_B, 2))
        y_batch(:, 1) = s_Bcovar_tp_B(1, 1, :)
        y_batch(:, 2) = s_Bcovar_tp_B(2, 1, :)

        call construct_batch_splines_1d(x_min, x_max, y_batch, order, .false., &
                                        bcovar_tp_batch_spline)
        bcovar_tp_batch_spline_ready = .true.
        deallocate (y_batch)
    end subroutine build_boozer_bcovar_tp_batch_spline

    subroutine build_boozer_bmod_br_batch_spline
        ! Combined Bmod + Br batch spline (1 or 2 quantities depending on use_B_r)
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, use_B_r

        real(dp) :: x_min(3), x_max(3)
        real(dp), allocatable :: y_batch(:, :, :, :)
        integer :: order(3), nq
        logical :: periodic(3)

        if (.not. allocated(bmod_grid)) then
            error stop "build_boozer_bmod_br_batch_spline: bmod_grid not allocated"
        end if
        if (use_B_r .and. .not. allocated(br_grid)) then
            error stop "build_boozer_bmod_br_batch_spline: br_grid not allocated"
        end if

        if (bmod_br_batch_spline_ready) then
            call destroy_batch_splines_3d(bmod_br_batch_spline)
            bmod_br_batch_spline_ready = .false.
            boozer_state%bmod_br_num_quantities = 0
        end if

        order = [ns_s_B, ns_tp_B, ns_tp_B]
        if (any(order < 3) .or. any(order > 5)) then
            error stop "build_boozer_bmod_br_batch_spline: spline order must be 3..5"
        end if

        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max(1) = hs_B*real(ns_B - 1, dp)
        x_max(2) = h_theta_B*real(n_theta_B - 1, dp)
        x_max(3) = h_phi_B*real(n_phi_B - 1, dp)

        periodic = [.false., .true., .true.]

        ! Determine number of quantities: 1 (Bmod only) or 2 (Bmod + Br)
        if (use_B_r) then
            nq = 2
        else
            nq = 1
        end if

        allocate (y_batch(ns_B, n_theta_B, n_phi_B, nq))
        y_batch(:, :, :, 1) = bmod_grid(:, :, :)
        if (use_B_r) then
            y_batch(:, :, :, 2) = br_grid(:, :, :)
        end if

        call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, &
                                        bmod_br_batch_spline)
        bmod_br_batch_spline_ready = .true.
        boozer_state%bmod_br_num_quantities = nq
        deallocate (y_batch)

        ! All scalar field parameters are final once the Bmod/B_r spline is
        ! built; refresh the device-accessible mirrors used by splint_boozer_coord
        ! (also keeps the host CPU path correct, since splint now reads them).
        call sync_boozer_state
    end subroutine build_boozer_bmod_br_batch_spline

    subroutine build_boozer_delt_delp_batch_splines
        ! Use the simple 4D grids populated in compute_boozer_data
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, use_del_tp_B

        integer :: order(3)
        real(dp) :: x_min(3), x_max(3)
        logical :: periodic(3)
        real(dp), allocatable :: y_batch(:, :, :, :)

        if (.not. allocated(delt_delp_V_grid)) then
            error stop &
                "build_boozer_delt_delp_batch_splines: delt_delp_V_grid not allocated"
        end if

        if (delt_delp_V_batch_spline_ready) then
            call destroy_batch_splines_3d(delt_delp_V_batch_spline)
            delt_delp_V_batch_spline_ready = .false.
        end if

        order = [ns_s_B, ns_tp_B, ns_tp_B]
        if (any(order < 3) .or. any(order > 5)) then
            error stop "build_boozer_delt_delp_batch_splines: order must be 3..5"
        end if

        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max(1) = hs_B*real(ns_B - 1, dp)
        x_max(2) = h_theta_B*real(n_theta_B - 1, dp)
        x_max(3) = h_phi_B*real(n_phi_B - 1, dp)

        periodic = [.false., .true., .true.]

        ! Build spline directly from 4D grid (already populated in compute_boozer_data)
        allocate (y_batch(ns_B, n_theta_B, n_phi_B, 2))
        y_batch(:, :, :, 1) = delt_delp_V_grid(:, :, :, 1)
        y_batch(:, :, :, 2) = delt_delp_V_grid(:, :, :, 2)

        call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, &
                                        delt_delp_V_batch_spline)
        delt_delp_V_batch_spline_ready = .true.

        if (use_del_tp_B) then
            if (.not. allocated(delt_delp_B_grid)) then
       error stop "build_boozer_delt_delp_batch_splines: delt_delp_B_grid not allocated"
            end if

            if (delt_delp_B_batch_spline_ready) then
                call destroy_batch_splines_3d(delt_delp_B_batch_spline)
                delt_delp_B_batch_spline_ready = .false.
            end if

            y_batch(:, :, :, 1) = delt_delp_B_grid(:, :, :, 1)
            y_batch(:, :, :, 2) = delt_delp_B_grid(:, :, :, 2)

            call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, &
                                            delt_delp_B_batch_spline)
            delt_delp_B_batch_spline_ready = .true.
        end if

        deallocate (y_batch)
    end subroutine build_boozer_delt_delp_batch_splines

end module boozer_sub
