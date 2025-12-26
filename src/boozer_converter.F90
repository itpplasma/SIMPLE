module boozer_sub
    use boozer_data_mod, only: boozer_set_current_field, compute_boozer_data, &
                               reset_boozer_batch_splines, &
                               build_boozer_aphi_batch_spline, &
                               build_boozer_bcovar_tp_batch_spline, &
                               build_boozer_bmod_br_batch_spline, &
                               build_boozer_delt_delp_batch_splines, &
                               aphi_batch_spline, aphi_batch_spline_ready, &
                               bcovar_tp_batch_spline, bcovar_tp_batch_spline_ready, &
                               bmod_br_batch_spline, bmod_br_batch_spline_ready, &
                               bmod_br_num_quantities, &
                               delt_delp_V_batch_spline, &
                               delt_delp_V_batch_spline_ready, &
                               delt_delp_B_batch_spline, delt_delp_B_batch_spline_ready
    use interpolate, only: evaluate_batch_splines_1d_der2, &
                           evaluate_batch_splines_1d_der3, &
                           evaluate_batch_splines_1d_many_der2, &
                           evaluate_batch_splines_1d_many_der3, &
                           evaluate_batch_splines_3d_der, &
                           evaluate_batch_splines_3d_der2, &
                           evaluate_batch_splines_3d_many_der2
    use field, only: magnetic_field_t
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    ! Public API
    public :: get_boozer_coordinates, get_boozer_coordinates_with_field
    public :: splint_boozer_coord, splint_boozer_coord_many
    public :: vmec_to_boozer, boozer_to_vmec
    public :: delthe_delphi_BV
    public :: reset_boozer_batch_splines

    ! Constants
    real(dp), parameter :: TWOPI = 2.0_dp*3.14159265358979_dp

contains

    !> Initialize Boozer coordinates using given magnetic field
	    subroutine get_boozer_coordinates_with_field(field)

	        class(magnetic_field_t), intent(in) :: field

	        ! Store field in module variable for use in nested subroutines
	        call boozer_set_current_field(field)
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
                                          hs_B, h_theta_B, h_phi_B, use_B_r

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

    subroutine splint_boozer_coord(r, vartheta_B, varphi_B, &
                                   A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                   d2A_phi_dr2, d3A_phi_dr3, &
                                   B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                                   B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                                   Bmod_B, dBmod_B, d2Bmod_B, &
                                   B_r, dB_r, d2B_r)

        use boozer_coordinates_mod, only: use_B_r
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper
        use chamb_mod, only: rnegflag
        use diag_mod, only: icounter

        implicit none

        integer, parameter :: mode_secders = 1

        real(dp) :: r, vartheta_B, varphi_B, &
                    A_phi, A_theta, dA_phi_dr, dA_theta_dr, d2A_phi_dr2, d3A_phi_dr3, &
                    B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                    B_varphi_B, dB_varphi_B, d2B_varphi_B, Bmod_B, B_r
        real(dp), dimension(3) :: dBmod_B, dB_r
        real(dp), dimension(6) :: d2Bmod_B, d2B_r

        real(dp) :: rho_tor, drhods, drhods2, d2rhods2m
        real(dp) :: qua, dqua_dr, dqua_dt, dqua_dp
        real(dp) :: d2qua_dr2, d2qua_drdt, d2qua_drdp, d2qua_dt2, &
                    d2qua_dtdp, d2qua_dp2
        real(dp) :: x_eval(3), y_eval(2), dy_eval(3, 2), d2y_eval(6, 2)
        real(dp) :: theta_wrapped, phi_wrapped
        real(dp) :: y1d(2), dy1d(2), d2y1d(2)

!$omp atomic
        icounter = icounter + 1
        if (r .le. 0.0_dp) then
            rnegflag = .true.
            r = abs(r)
        end if

        A_theta = torflux*r
        dA_theta_dr = torflux

        ! Interpolate A_phi over s (batch spline 1D)
        if (.not. aphi_batch_spline_ready) then
            error stop "splint_boozer_coord: Aphi batch spline not initialized"
        end if

        if (mode_secders > 0) then
            ! Need third derivative - use der3 which computes all in one pass
            block
                real(dp) :: d3y1d(1)
                call evaluate_batch_splines_1d_der3(aphi_batch_spline, r, &
                                                    y1d(1:1), dy1d(1:1), &
                                                    d2y1d(1:1), d3y1d)
                d3A_phi_dr3 = d3y1d(1)
            end block
        else
            call evaluate_batch_splines_1d_der2(aphi_batch_spline, r, y1d(1:1), &
                                                dy1d(1:1), d2y1d(1:1))
            d3A_phi_dr3 = 0.0_dp
        end if
        A_phi = y1d(1)
        dA_phi_dr = dy1d(1)
        d2A_phi_dr2 = d2y1d(1)

        ! Interpolation of mod-B (and B_r if use_B_r)
        rho_tor = sqrt(r)
        theta_wrapped = modulo(vartheta_B, TWOPI)
        phi_wrapped = modulo(varphi_B, TWOPI/real(nper, dp))

        if (.not. bmod_br_batch_spline_ready) then
            error stop "splint_boozer_coord: Bmod/Br batch spline not initialized"
        end if

        x_eval(1) = rho_tor
        x_eval(2) = theta_wrapped
        x_eval(3) = phi_wrapped
        call evaluate_batch_splines_3d_der2(bmod_br_batch_spline, x_eval, &
                                            y_eval(1:bmod_br_num_quantities), &
                                            dy_eval(:, 1:bmod_br_num_quantities), &
                                            d2y_eval(:, 1:bmod_br_num_quantities))

        ! Chain rule coefficients for rho -> s conversion
        drhods = 0.5_dp/rho_tor
        drhods2 = drhods**2
        d2rhods2m = drhods2/rho_tor  ! -d2rho/ds2 (negative of second derivative)

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
        if (use_B_r) then
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
                       qua*drhods*(3.0_dp/4.0_dp)/r**2
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

        ! Interpolation of B_\vartheta and B_\varphi (flux functions)
        if (.not. bcovar_tp_batch_spline_ready) then
            error stop "splint_boozer_coord: Bcovar_tp batch spline not initialized"
        end if

        call evaluate_batch_splines_1d_der2(bcovar_tp_batch_spline, rho_tor, y1d, &
                                            dy1d, d2y1d)
        B_vartheta_B = y1d(1)
        dB_vartheta_B = dy1d(1)
        d2B_vartheta_B = d2y1d(1)
        B_varphi_B = y1d(2)
        dB_varphi_B = dy1d(2)
        d2B_varphi_B = d2y1d(2)

        d2B_vartheta_B = d2B_vartheta_B*drhods2 - dB_vartheta_B*d2rhods2m
        d2B_varphi_B = d2B_varphi_B*drhods2 - dB_varphi_B*d2rhods2m
        dB_vartheta_B = dB_vartheta_B*drhods
        dB_varphi_B = dB_varphi_B*drhods

    end subroutine splint_boozer_coord

    subroutine splint_boozer_coord_many(npts, r, vartheta_B, varphi_B, &
                                        A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                        d2A_phi_dr2, d3A_phi_dr3, &
                                        B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                                        B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                                        Bmod_B, dBmod_B, d2Bmod_B)
        use boozer_coordinates_mod, only: use_B_r
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper
        use chamb_mod, only: rnegflag
        use diag_mod, only: icounter

        integer, intent(in) :: npts
        real(dp), intent(inout) :: r(npts)
        real(dp), intent(in) :: vartheta_B(npts), varphi_B(npts)
        real(dp), intent(out) :: A_phi(npts), A_theta(npts)
        real(dp), intent(out) :: dA_phi_dr(npts), dA_theta_dr(npts)
        real(dp), intent(out) :: d2A_phi_dr2(npts), d3A_phi_dr3(npts)
        real(dp), intent(out) :: B_vartheta_B(npts), dB_vartheta_B(npts)
        real(dp), intent(out) :: d2B_vartheta_B(npts)
        real(dp), intent(out) :: B_varphi_B(npts), dB_varphi_B(npts), d2B_varphi_B(npts)
        real(dp), intent(out) :: Bmod_B(npts)
        real(dp), intent(out) :: dBmod_B(3, npts), d2Bmod_B(6, npts)

        integer :: ipt, n_quantities
        real(dp) :: rho_tor(npts), drhods(npts), drhods2(npts), d2rhods2m(npts)
        real(dp) :: theta_wrapped(npts), phi_wrapped(npts)
        real(dp) :: x_eval(3, npts)
        real(dp) :: y_aphi(1, npts), dy_aphi(1, npts), d2y_aphi(1, npts)
        real(dp) :: d3y_aphi(1, npts)
        real(dp) :: y_bmod(2, npts), dy_bmod(3, 2, npts), d2y_bmod(6, 2, npts)
        real(dp) :: y_btp(2, npts), dy_btp(2, npts), d2y_btp(2, npts)

        if (.not. aphi_batch_spline_ready) then
            error stop "splint_boozer_coord_many: Aphi batch spline not initialized"
        end if
        if (.not. bmod_br_batch_spline_ready) then
            error stop "splint_boozer_coord_many: Bmod/Br batch spline not initialized"
        end if
        if (.not. bcovar_tp_batch_spline_ready) then
            error stop "splint_boozer_coord_many: Bcovar_tp spline not initialized"
        end if

        do ipt = 1, npts
            if (r(ipt) <= 0.0_dp) then
                rnegflag = .true.
                r(ipt) = abs(r(ipt))
            end if
        end do

        A_theta = torflux*r
        dA_theta_dr = torflux

        call evaluate_batch_splines_1d_many_der3(aphi_batch_spline, r, &
                                                 y_aphi, dy_aphi, d2y_aphi, d3y_aphi)
        A_phi = y_aphi(1, :)
        dA_phi_dr = dy_aphi(1, :)
        d2A_phi_dr2 = d2y_aphi(1, :)
        d3A_phi_dr3 = d3y_aphi(1, :)

        rho_tor = sqrt(r)
        do ipt = 1, npts
            theta_wrapped(ipt) = modulo(vartheta_B(ipt), TWOPI)
            phi_wrapped(ipt) = modulo(varphi_B(ipt), TWOPI/real(nper, dp))
        end do

        x_eval(1, :) = rho_tor
        x_eval(2, :) = theta_wrapped
        x_eval(3, :) = phi_wrapped

        n_quantities = bmod_br_num_quantities
        call evaluate_batch_splines_3d_many_der2(bmod_br_batch_spline, x_eval, &
                                                 y_bmod(1:n_quantities, :), &
                                                 dy_bmod(:, 1:n_quantities, :), &
                                                 d2y_bmod(:, 1:n_quantities, :))

        drhods = 0.5_dp/rho_tor
        drhods2 = drhods**2
        d2rhods2m = drhods2/rho_tor

        Bmod_B = y_bmod(1, :)
        dBmod_B(1, :) = dy_bmod(1, 1, :)*drhods
        dBmod_B(2, :) = dy_bmod(2, 1, :)
        dBmod_B(3, :) = dy_bmod(3, 1, :)

        d2Bmod_B(1, :) = d2y_bmod(1, 1, :)*drhods2 - dy_bmod(1, 1, :)*d2rhods2m
        d2Bmod_B(2, :) = d2y_bmod(2, 1, :)*drhods
        d2Bmod_B(3, :) = d2y_bmod(3, 1, :)*drhods
        d2Bmod_B(4, :) = d2y_bmod(4, 1, :)
        d2Bmod_B(5, :) = d2y_bmod(5, 1, :)
        d2Bmod_B(6, :) = d2y_bmod(6, 1, :)

        call evaluate_batch_splines_1d_many_der2(bcovar_tp_batch_spline, rho_tor, &
                                                 y_btp, dy_btp, d2y_btp)
        B_vartheta_B = y_btp(1, :)
        dB_vartheta_B = dy_btp(1, :)*drhods
        d2B_vartheta_B = d2y_btp(1, :)*drhods2 - dy_btp(1, :)*d2rhods2m
        B_varphi_B = y_btp(2, :)
        dB_varphi_B = dy_btp(2, :)*drhods
        d2B_varphi_B = d2y_btp(2, :)*drhods2 - dy_btp(2, :)*d2rhods2m

!$omp atomic
        icounter = icounter + npts

    end subroutine splint_boozer_coord_many

    !> Computes delta_vartheta = vartheta_B - theta_V
    !> and delta_varphi = varphi_B - varphi_V
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

    !> Convert VMEC coordinates (r, theta, varphi)
    !> to Boozer coordinates (vartheta_B, varphi_B)
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

    !> Convert Boozer coordinates (r, vartheta_B, varphi_B)
    !> to VMEC coordinates (theta, varphi)
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



end module boozer_sub
