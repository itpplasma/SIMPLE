module boozer_sub
    use spl_three_to_five_sub
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                           evaluate_batch_splines_3d_der2, destroy_batch_splines_3d
    use field, only : magnetic_field_t
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    real(dp), parameter :: twopi = 2.0_dp*3.14159265358979_dp

    ! Module variable to store the field for use in subroutines
    class(magnetic_field_t), allocatable :: current_field
    !$omp threadprivate(current_field)

    type(BatchSplineData3D), save :: bmod_batch_spline
    logical, save :: bmod_batch_spline_ready = .false.
    real(dp), allocatable, save :: bmod_grid(:, :, :)

    type(BatchSplineData3D), save :: br_batch_spline
    logical, save :: br_batch_spline_ready = .false.
    real(dp), allocatable, save :: br_grid(:, :, :)

contains

    subroutine get_boozer_coordinates_with_field(field)

! Field-agnostic version that accepts a magnetic_field_t object

        implicit none

        class(magnetic_field_t), intent(in) :: field

        ! Store field in module variable for use in nested subroutines
        if (allocated(current_field)) deallocate (current_field)
        allocate (current_field, source=field)
        call reset_boozer_batch_splines

        ! Call the actual implementation
        call get_boozer_coordinates_impl

    end subroutine get_boozer_coordinates_with_field

    subroutine get_boozer_coordinates

! Backward compatibility wrapper - uses VMEC field by default

        use field, only: vmec_field_t

        call get_boozer_coordinates_with_field(vmec_field_t())

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

        hs_B = hs*dble(ns - 1)/dble(ns_B - 1)
        h_theta_B = h_theta*dble(n_theta - 1)/dble(n_theta_B - 1)
        h_phi_B = h_phi*dble(n_phi - 1)/dble(n_phi_B - 1)

        call compute_boozer_data

        call spline_boozer_data
        call build_boozer_bmod_batch_spline
        if (use_B_r) call build_boozer_br_batch_spline

    end subroutine get_boozer_coordinates_impl

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine splint_boozer_coord(r, vartheta_B, varphi_B, &
                                   A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                   d2A_phi_dr2, d3A_phi_dr3, &
                                   B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                                   B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                                   Bmod_B, dBmod_B, d2Bmod_B, &
                                   B_r, dB_r, d2B_r)

        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_Bcovar_tp_B, &
                                          derf1, derf2, derf3, &
                                          use_B_r
        use vector_potentail_mod, only: ns, hs, torflux, sA_phi
        use new_vmec_stuff_mod, only: nper, ns_A
        use chamb_mod, only: rnegflag
        use diag_mod, only: icounter

        implicit none

        integer, parameter :: mode_secders = 1

        integer :: nstp, ns_A_p1, ns_s_p1
        integer :: k, is, i_theta, i_phi

        real(dp) :: r, vartheta_B, varphi_B, &
            A_phi, A_theta, dA_phi_dr, dA_theta_dr, d2A_phi_dr2, d3A_phi_dr3, &
            B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
            B_varphi_B, dB_varphi_B, d2B_varphi_B, Bmod_B, B_r
        real(dp), dimension(3) :: dBmod_B, dB_r
        real(dp), dimension(6) :: d2Bmod_B, d2B_r

        real(dp) :: ds, dtheta, dphi, rho_tor, drhods, drhods2, d2rhods2m
        real(dp) :: qua, dqua_dr, dqua_dt, dqua_dp
        real(dp) :: d2qua_dr2, d2qua_drdt, d2qua_drdp, d2qua_dt2, &
            d2qua_dtdp, d2qua_dp2
        real(dp) :: x_eval(3), y_eval(1), dy_eval(3, 1), d2y_eval(6, 1)
        real(dp) :: theta_wrapped, phi_wrapped

!$omp atomic
        icounter = icounter + 1
        if (r .le. 0.d0) then
            rnegflag = .true.
            r = abs(r)
        end if

        A_theta = torflux*r
        dA_theta_dr = torflux

        call normalize_angular_coordinates(vartheta_B, varphi_B, n_theta_B, n_phi_B, &
                                           h_theta_B, h_phi_B, &
                                           i_theta, i_phi, dtheta, dphi)

! Begin interpolation of vector potentials over $s$

        ds = r/hs
        is = max(0, min(ns - 1, int(ds)))
        ds = (ds - dble(is))*hs
        is = is + 1

        ns_A_p1 = ns_A + 1
        A_phi = sA_phi(ns_A_p1, is)
        dA_phi_dr = A_phi*derf1(ns_A_p1)
        d2A_phi_dr2 = A_phi*derf2(ns_A_p1)

        do k = ns_A, 3, -1
            A_phi = sA_phi(k, is) + ds*A_phi
            dA_phi_dr = sA_phi(k, is)*derf1(k) + ds*dA_phi_dr
            d2A_phi_dr2 = sA_phi(k, is)*derf2(k) + ds*d2A_phi_dr2
        end do

        A_phi = sA_phi(1, is) + ds*(sA_phi(2, is) + ds*A_phi)
        dA_phi_dr = sA_phi(2, is) + ds*dA_phi_dr

        if (mode_secders .gt. 0) then
            d3A_phi_dr3 = sA_phi(ns_A_p1, is)*derf3(ns_A_p1)

            do k = ns_A, 4, -1
                d3A_phi_dr3 = sA_phi(k, is)*derf3(k) + ds*d3A_phi_dr3
            end do

        end if

! End interpolation of vector potentials over $s$

!--------------------------------

        rho_tor = sqrt(r)
        ds = rho_tor/hs_B
        is = max(0, min(ns_B - 1, int(ds)))
        ds = (ds - dble(is))*hs_B
        is = is + 1

        nstp = ns_tp_B + 1
        ns_s_p1 = ns_s_B + 1

!--------------------------------
! Interpolation of mod-B:
!--------------------------------

        theta_wrapped = modulo(vartheta_B, twopi)
        phi_wrapped = modulo(varphi_B, twopi/dble(nper))

        if (.not. bmod_batch_spline_ready) then
            error stop "splint_boozer_coord: Bmod batch spline not initialized"
        end if

        x_eval(1) = rho_tor
        x_eval(2) = theta_wrapped
        x_eval(3) = phi_wrapped
        call evaluate_batch_splines_3d_der2(bmod_batch_spline, x_eval, y_eval, &
                                            dy_eval, d2y_eval)

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

! Coversion coefficients for derivatives over s
        drhods = 0.5d0/rho_tor
        drhods2 = drhods**2
        ! $-\dr^2 \rho / \rd s^2$ (second derivative with minus sign)
        d2rhods2m = drhods2/rho_tor

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

!--------------------------------
! End Interpolation of mod-B
!--------------------------------
! Interpolation of B_\vartheta and B_\varphi:
!--------------------------------

        B_vartheta_B = s_Bcovar_tp_B(1, ns_s_p1, is)
        dB_vartheta_B = B_vartheta_B*derf1(ns_s_p1)
        d2B_vartheta_B = B_vartheta_B*derf2(ns_s_p1)
        B_varphi_B = s_Bcovar_tp_B(2, ns_s_p1, is)
        dB_varphi_B = B_varphi_B*derf1(ns_s_p1)
        d2B_varphi_B = B_varphi_B*derf2(ns_s_p1)

        do k = ns_s_B, 3, -1
            B_vartheta_B = s_Bcovar_tp_B(1, k, is) + ds*B_vartheta_B
            dB_vartheta_B = s_Bcovar_tp_B(1, k, is)*derf1(k) + ds*dB_vartheta_B
            d2B_vartheta_B = s_Bcovar_tp_B(1, k, is)*derf2(k) + ds*d2B_vartheta_B
            B_varphi_B = s_Bcovar_tp_B(2, k, is) + ds*B_varphi_B
            dB_varphi_B = s_Bcovar_tp_B(2, k, is)*derf1(k) + ds*dB_varphi_B
            d2B_varphi_B = s_Bcovar_tp_B(2, k, is)*derf2(k) + ds*d2B_varphi_B
        end do

        B_vartheta_B = s_Bcovar_tp_B(1, 1, is) &
            + ds*(s_Bcovar_tp_B(1, 2, is) + ds*B_vartheta_B)
        dB_vartheta_B = s_Bcovar_tp_B(1, 2, is) + ds*dB_vartheta_B
        B_varphi_B = s_Bcovar_tp_B(2, 1, is) &
            + ds*(s_Bcovar_tp_B(2, 2, is) + ds*B_varphi_B)
        dB_varphi_B = s_Bcovar_tp_B(2, 2, is) + ds*dB_varphi_B

        d2B_vartheta_B = d2B_vartheta_B*drhods2 - dB_vartheta_B*d2rhods2m
        d2B_varphi_B = d2B_varphi_B*drhods2 - dB_varphi_B*d2rhods2m
        dB_vartheta_B = dB_vartheta_B*drhods
        dB_varphi_B = dB_varphi_B*drhods

!--------------------------------
! End interpolation of B_\vartheta and B_\varphi
!--------------------------------
! Interpolation of B_r:
!--------------------------------
! Note that splined quantity is $B_\rho$, not $B_s$

        if (use_B_r) then

            if (.not. br_batch_spline_ready) then
                error stop "splint_boozer_coord: Br batch spline not initialized"
            end if

            x_eval(1) = rho_tor
            x_eval(2) = theta_wrapped
            x_eval(3) = phi_wrapped
            call evaluate_batch_splines_3d_der2(br_batch_spline, x_eval, y_eval, &
                                                dy_eval, d2y_eval)

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

! Convert coefficients for derivatives over s
            drhods = 0.5d0/rho_tor
            drhods2 = drhods**2
            ! $-\dr^2 \rho / \rd s^2$ (second derivative with minus sign)
            d2rhods2m = drhods2/rho_tor

            d2qua_dr2 = d2qua_dr2*drhods2 - dqua_dr*d2rhods2m
            dqua_dr = dqua_dr*drhods
            d2qua_drdt = d2qua_drdt*drhods
            d2qua_drdp = d2qua_drdp*drhods

            B_r = qua*drhods

            dB_r(1) = dqua_dr*drhods - qua*d2rhods2m
            dB_r(2) = dqua_dt*drhods
            dB_r(3) = dqua_dp*drhods

            d2B_r(1) = d2qua_dr2*drhods - 2.d0*dqua_dr*d2rhods2m + &
                       qua*drhods*(3.d0/4.d0)/r**2
            d2B_r(2) = d2qua_drdt*drhods - dqua_dt*d2rhods2m
            d2B_r(3) = d2qua_drdp*drhods - dqua_dp*d2rhods2m
            d2B_r(4) = d2qua_dt2*drhods
            d2B_r(5) = d2qua_dtdp*drhods
            d2B_r(6) = d2qua_dp2*drhods

        else
            B_r = 0.d0
            dB_r = 0.d0
            d2B_r = 0.d0
        end if
!--------------------------------
! End interpolation of B_r
!--------------------------------
    end subroutine splint_boozer_coord

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine delthe_delphi_BV(isw, r, vartheta, varphi, deltheta_BV, delphi_BV, &
                                ddeltheta_BV, ddelphi_BV)

! Computes $\Delta \vartheta = \vartheta_B - \theta_V$ and
! $\Delta \varphi = \varphi_B - \varphi_V$
! and their first derivatives over angles for two cases:
! isw=0 - if they are given as functions of VMEC coordinates (r,vartheta,varphi)
! isw=1 - if they are given as functions of Boozer coordinates (r,vartheta,varphi)

        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_delt_delp_V, s_delt_delp_B, &
                                          ns_max, derf1, &
                                          use_del_tp_B
        use chamb_mod, only: rnegflag

        implicit none

        integer :: isw
        real(dp) :: r, vartheta, varphi, deltheta_BV, delphi_BV
        real(dp), dimension(2) :: ddeltheta_BV, ddelphi_BV

        integer, parameter :: n_qua = 2
        integer :: nstp, ns_s_p1
        integer :: k, is, i_theta, i_phi

        real(dp) :: ds, dtheta, dphi, rho_tor

        real(dp), dimension(n_qua) :: qua, dqua_dt, dqua_dp
        real(dp), dimension(n_qua, ns_max) :: sp_all, dsp_all_dt
        real(dp), dimension(n_qua, ns_max, ns_max) :: stp_all

        sp_all = 0.0_dp
        dsp_all_dt = 0.0_dp
        stp_all = 0.0_dp

        if (r .le. 0.d0) then
            rnegflag = .true.
            r = abs(r)
        end if

        call normalize_angular_coordinates(vartheta, varphi, n_theta_B, n_phi_B, &
                                           h_theta_B, h_phi_B, &
                                           i_theta, i_phi, dtheta, dphi)

!--------------------------------

        rho_tor = sqrt(r)
        ds = rho_tor/hs_B
        is = max(0, min(ns_B - 1, int(ds)))
        ds = (ds - dble(is))*hs_B
        is = is + 1

        nstp = ns_tp_B + 1
        ns_s_p1 = ns_s_B + 1

!--------------------------------

! Begin interpolation of all over $rho$

        if (isw .eq. 0) then
            stp_all(:, 1:nstp, 1:nstp) = &
                s_delt_delp_V(:, ns_s_p1, :, :, is, i_theta, i_phi)

            do k = ns_s_B, 1, -1
                stp_all(:, 1:nstp, 1:nstp) = &
                    s_delt_delp_V(:, k, :, :, is, i_theta, i_phi) &
                    + ds*stp_all(:, 1:nstp, 1:nstp)
            end do
        elseif (isw .eq. 1) then
            if (.not. use_del_tp_B) then
                print *, 'delthe_delphi_BV : Boozer data is not loaded'
                return
            end if
            stp_all(:, 1:nstp, 1:nstp) = &
                s_delt_delp_B(:, ns_s_p1, :, :, is, i_theta, i_phi)

            do k = ns_s_B, 1, -1
                stp_all(:, 1:nstp, 1:nstp) = &
                    s_delt_delp_B(:, k, :, :, is, i_theta, i_phi) &
                    + ds*stp_all(:, 1:nstp, 1:nstp)
            end do
        else
            print *, 'delthe_delphi_BV : unknown value of switch isw'
            return
        end if

! End interpolation of all over $rho$
!-------------------------------
! Begin interpolation of all over $\theta$

        sp_all(:, 1:nstp) = stp_all(:, nstp, 1:nstp)
        dsp_all_dt(:, 1:nstp) = sp_all(:, 1:nstp)*derf1(nstp)

        do k = ns_tp_B, 2, -1
            sp_all(:, 1:nstp) = stp_all(:, k, 1:nstp) + dtheta*sp_all(:, 1:nstp)
            dsp_all_dt(:, 1:nstp) = stp_all(:, k, 1:nstp)*derf1(k) &
                + dtheta*dsp_all_dt(:, 1:nstp)
        end do

        sp_all(:, 1:nstp) = stp_all(:, 1, 1:nstp) + dtheta*sp_all(:, 1:nstp)

! End interpolation of all over $\theta$
!--------------------------------
! Begin interpolation of all over $\varphi$

        qua = sp_all(:, nstp)
        dqua_dt = dsp_all_dt(:, nstp)
        dqua_dp = qua*derf1(nstp)

        do k = ns_tp_B, 2, -1
            qua = sp_all(:, k) + dphi*qua
            dqua_dt = dsp_all_dt(:, k) + dphi*dqua_dt
            dqua_dp = sp_all(:, k)*derf1(k) + dphi*dqua_dp
        end do

        qua = sp_all(:, 1) + dphi*qua
        dqua_dt = dsp_all_dt(:, 1) + dphi*dqua_dt

! End interpolation of all over $\varphi$

        deltheta_BV = qua(1)
        delphi_BV = qua(2)

        ddeltheta_BV(1) = dqua_dt(1)
        ddelphi_BV(1) = dqua_dt(2)

        ddeltheta_BV(2) = dqua_dp(1)
        ddelphi_BV(2) = dqua_dp(2)

    end subroutine delthe_delphi_BV

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vmec_to_boozer(r, theta, varphi, vartheta_B, varphi_B)

! Input : r,theta,varphi      - VMEC coordinates
! Output: vartheta_B,varphi_B - Boozer coordinates

        use new_vmec_stuff_mod, only: nper

        implicit none


        real(dp) :: r, theta, varphi, vartheta_B, varphi_B
        real(dp) :: deltheta_BV, delphi_BV
        real(dp), dimension(2) :: ddeltheta_BV, ddelphi_BV

        call delthe_delphi_BV(0, r, theta, varphi, deltheta_BV, delphi_BV, &
                              ddeltheta_BV, ddelphi_BV)

        vartheta_B = modulo(theta + deltheta_BV, twopi)
        varphi_B = modulo(varphi + delphi_BV, twopi/dble(nper))

    end subroutine vmec_to_boozer

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine boozer_to_vmec(r, vartheta_B, varphi_B, theta, varphi)

! Input : r,vartheta_B,varphi_B - Boozer coordinates
! Output: theta,varphi          - VMEC coordinates

        use boozer_coordinates_mod, only: use_del_tp_B

        implicit none

        real(dp), parameter :: epserr = 1.d-14
        integer, parameter :: niter = 100

        integer :: iter
        real(dp) :: r, theta, varphi, vartheta_B, varphi_B
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
            f11 = 1.d0 + ddeltheta_BV(1)
            f12 = ddeltheta_BV(2)
            f21 = ddelphi_BV(1)
            f22 = 1.d0 + ddelphi_BV(2)

            det = f11*f22 - f12*f21
            delthe = (f2*f12 - f1*f22)/det
            delphi = (f1*f21 - f2*f11)/det

            theta = theta + delthe
            varphi = varphi + delphi
            if (abs(delthe) + abs(delphi) .lt. epserr) exit
        end do

!  theta=modulo(theta,twopi)
!  varphi=modulo(varphi,twopi/dble(nper))

    end subroutine boozer_to_vmec

    subroutine compute_boozer_data
        ! Computes Boozer coordinate transformations and magnetic field data
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_Bcovar_tp_B, &
                                          s_Bmod_B, s_Bcovar_r_B, &
                                          s_delt_delp_V, s_delt_delp_B, &
                                          use_B_r, use_del_tp_B
        use binsrc_sub, only: binsrc
        use plag_coeff_sub, only: plag_coeff
        use spline_vmec_sub
#ifdef GVEC_AVAILABLE
        use vmec_field_adapter
#else
        use vmec_field_eval
#endif

        implicit none

        real(dp), parameter :: s_min = 1.d-6, rho_min = sqrt(s_min)

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
        gridcellnum = dble((n_theta_B - 1)*(n_phi_B - 1))

        npoilag = ns_tp_B + 1
        nder = 0
        nshift = npoilag/2

        print *, 'Transforming to Boozer coordinates'

        if (use_B_r) then
            print *, 'B_r is computed'
        else
            print *, 'B_r is not computed'
        end if

        G00 = 0.d0

        allocate (rho_tor(ns_B))
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

! allocate spline coefficients for Boozer data:
        if (.not. allocated(s_Bcovar_tp_B)) &
            allocate (s_Bcovar_tp_B(2, ns_s_B + 1, ns_B))
        if (.not. allocated(s_Bmod_B)) &
            allocate (s_Bmod_B(ns_s_B + 1, ns_tp_B + 1, ns_tp_B + 1, &
                ns_B, n_theta_B, n_phi_B))
        if (.not. allocated(bmod_grid)) then
            allocate (bmod_grid(ns_B, n_theta_B, n_phi_B))
        else if (any(shape(bmod_grid) /= [ns_B, n_theta_B, n_phi_B])) then
            deallocate (bmod_grid)
            allocate (bmod_grid(ns_B, n_theta_B, n_phi_B))
        end if
        if (use_B_r .and. .not. allocated(s_Bcovar_r_B)) &
            allocate (s_Bcovar_r_B(ns_s_B + 1, ns_tp_B + 1, &
                ns_tp_B + 1, ns_B, n_theta_B, n_phi_B))
        if (.not. allocated(s_delt_delp_V)) &
            allocate (s_delt_delp_V(2, ns_s_B + 1, ns_tp_B + 1, &
                ns_tp_B + 1, ns_B, n_theta_B, n_phi_B))
        if (use_del_tp_B .and. .not. allocated(s_delt_delp_B)) &
            allocate (s_delt_delp_B(2, ns_s_B + 1, ns_tp_B + 1, &
                ns_tp_B + 1, ns_B, n_theta_B, n_phi_B))

        do i = 0, ns_tp_B
            wint_t(i) = h_theta_B**(i + 1)/dble(i + 1)
            wint_p(i) = h_phi_B**(i + 1)/dble(i + 1)
        end do

        ! Set theta_V and phi_V linear, with value 0 at index 1 and stepsize h.
        ! Then expand this in both directions beyond 1:n_theta_B.
        do i_theta = 1, n_theta_B
            theta_V(i_theta) = dble(i_theta - 1)*h_theta_B
        end do
        per_theta = dble(n_theta_B - 1)*h_theta_B
        theta_V(2 - n_theta_B:0) = theta_V(1:n_theta_B - 1) - per_theta
        theta_V(n_theta_B + 1:2*n_theta_B - 1) = theta_V(2:n_theta_B) + per_theta

        do i_phi = 1, n_phi_B
            phi_V(i_phi) = dble(i_phi - 1)*h_phi_B
        end do
        per_phi = dble(n_phi_B - 1)*h_phi_B
        phi_V(2 - n_phi_B:0) = phi_V(1:n_phi_B - 1) - per_phi
        phi_V(n_phi_B + 1:2*n_phi_B - 1) = phi_V(2:n_phi_B) + per_phi

        do i_rho = 1, ns_B
            rho_tor(i_rho) = max(dble(i_rho - 1)*hs_B, rho_min)
            s = rho_tor(i_rho)**2

            do i_theta = 1, n_theta_B
                theta = dble(i_theta - 1)*h_theta_B
                do i_phi = 1, n_phi_B
                    varphi = dble(i_phi - 1)*h_phi_B

                    if (allocated(current_field)) then
               call vmec_field_evaluate_with_field(current_field, &
                    s, theta, varphi, A_theta, A_phi, &
                    dA_theta_ds, dA_phi_ds, aiota, &
                    sqg, alam, dl_ds, dl_dt, dl_dp, &
                    Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
                    else
                        call vmec_field_evaluate(s, theta, varphi, &
                            A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                            sqg, alam, dl_ds, dl_dt, dl_dp, &
                            Bctrvr_vartheta, Bctrvr_varphi, &
                            Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
                    end if

                    alam_2D(i_theta, i_phi) = alam
                    bmod_Vg(i_theta, i_phi) = &
                        sqrt(Bctrvr_vartheta*Bcovar_vartheta &
                        + Bctrvr_varphi*Bcovar_varphi)
                    Bcovar_theta_V(i_theta, i_phi) = Bcovar_vartheta*(1.d0 + dl_dt)
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

            denomjac = 1.d0/(aiota*Bcovar_vartheta_B + Bcovar_varphi_B)
            Gbeg = G00 + Bcovar_vartheta_B*denomjac*alam_2D(1, 1)

            splcoe_t(0, :) = Bcovar_theta_V(:, 1)

            call spl_per(ns_tp_B, n_theta_B, h_theta_B, splcoe_t)

            delphi_BV_Vg(1, 1) = 0.d0
            do i_theta = 1, n_theta_B - 1
                delphi_BV_Vg(i_theta + 1, 1) = &
                    delphi_BV_Vg(i_theta, 1) &
                    + sum(wint_t*splcoe_t(:, i_theta))
            end do
            ! Remove linear increasing component from delphi_BV_Vg
            aper = (delphi_BV_Vg(n_theta_B, 1) &
                - delphi_BV_Vg(1, 1))/dble(n_theta_B - 1)
            do i_theta = 2, n_theta_B
                delphi_BV_Vg(i_theta, 1) = &
                    delphi_BV_Vg(i_theta, 1) - aper*dble(i_theta - 1)
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
                    - delphi_BV_Vg(i_theta, 1))/dble(n_phi_B - 1)
                do i_phi = 2, n_phi_B
                    delphi_BV_Vg(i_theta, i_phi) = &
                        delphi_BV_Vg(i_theta, i_phi) &
                        - aper*dble(i_phi - 1)
                end do
            end do

! difference between Boozer and VMEC toroidal angle,
! $\Delta \varphi_{BV}=\varphi_B-\varphi=G$:
            delphi_BV_Vg = denomjac*delphi_BV_Vg + Gbeg
! difference between Boozer and VMEC poloidal angle,
! $\Delta \vartheta_{BV}=\vartheta_B-\theta$:
            deltheta_BV_Vg = aiota*delphi_BV_Vg + alam_2D

            s_delt_delp_V(1, 1, 1, 1, i_rho, :, :) = deltheta_BV_Vg
            s_delt_delp_V(2, 1, 1, 1, i_rho, :, :) = delphi_BV_Vg

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

                    call binsrc(theta_B, 2 - n_theta_B, 2*n_theta_B - 1, theta_V(i_theta), i)

                    ibeg = i - nshift
                    iend = ibeg + ns_tp_B

                    call plag_coeff(npoilag, nder, theta_V(i_theta), theta_B(ibeg:iend), coef)

                    perqua_2D(:, i_theta, i_phi) = matmul(perqua_t(:, ibeg:iend), coef(0, :))
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

                    perqua_2D(:, i_theta, i_phi) = matmul(perqua_p(:, ibeg:iend), coef(0, :))
                end do
            end do

            if (use_del_tp_B) s_delt_delp_B(:, 1, 1, 1, i_rho, :, :) = perqua_2D(1:2, :, :)
            s_Bmod_B(1, 1, 1, i_rho, :, :) = perqua_2D(3, :, :)
            bmod_grid(i_rho, :, :) = perqua_2D(3, :, :)

! End re-interpolate to equidistant grid in $(\vartheta_B,\varphi_B)$

            if (use_B_r) then
                aiota_arr(i_rho) = aiota
                Gfunc(i_rho, :, :) = perqua_2D(2, :, :)
! covariant components $B_k$ in symmetry flux coordinates on equidistant grid of Boozer coordinates:
                Bcovar_symfl(:, i_rho, :, :) = perqua_2D(4:6, :, :)
            end if

        end do

! Compute radial covariant magnetic field component in Boozer coordinates

        if (use_B_r) then
            if (allocated(coef)) deallocate (coef)
            nder = 1
            npoilag = 5
            nshift = npoilag/2
            allocate (coef(0:nder, npoilag))

            do i_rho = 1, ns_B
                ibeg = i_rho - nshift
                iend = ibeg + npoilag - 1
                if (ibeg .lt. 1) then
                    ibeg = 1
                    iend = ibeg + npoilag - 1
                elseif (iend .gt. ns_B) then
                    iend = ns_B
                    ibeg = iend - npoilag + 1
                end if

                call plag_coeff(npoilag, nder, rho_tor(i_rho), rho_tor(ibeg:iend), coef)

! We spline covariant component $B_\rho$ instead of $B_s$:
                do i_phi = 1, n_phi_B
                    s_Bcovar_r_B(1, 1, 1, i_rho, :, i_phi) = &
                        2.d0*rho_tor(i_rho) &
                        *Bcovar_symfl(1, i_rho, :, i_phi) &
                        - matmul(coef(1, :)*aiota_arr(ibeg:iend), &
                            Gfunc(ibeg:iend, :, i_phi)) &
                        *Bcovar_symfl(2, i_rho, :, i_phi) &
                        - matmul(coef(1, :), Gfunc(ibeg:iend, :, i_phi)) &
                        *Bcovar_symfl(3, i_rho, :, i_phi)
                end do

            end do
            deallocate (aiota_arr, Gfunc, Bcovar_symfl)

            if (.not. allocated(br_grid)) then
                allocate (br_grid(ns_B, n_theta_B, n_phi_B))
            else if (any(shape(br_grid) /= [ns_B, n_theta_B, n_phi_B])) then
                deallocate (br_grid)
                allocate (br_grid(ns_B, n_theta_B, n_phi_B))
            end if

            br_grid(:, :, :) = s_Bcovar_r_B(1, 1, 1, :, :, :)
        end if

! End compute radial covariant magnetic field component in Boozer coordinates

        deallocate (Bcovar_theta_V, Bcovar_varphi_V, bmod_Vg, alam_2D, &
                    deltheta_BV_Vg, delphi_BV_Vg, &
                    wint_t, wint_p, coef, theta_V, theta_B, phi_V, phi_B, &
                    perqua_t, perqua_p, perqua_2D)

        print *, 'done'

    end subroutine compute_boozer_data

    subroutine spline_boozer_data
        ! Splines Boozer coordinate data over angles and radial direction
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_Bcovar_tp_B, &
                                          s_Bmod_B, s_Bcovar_r_B, &
                                          s_delt_delp_V, s_delt_delp_B, &
                                          ns_max, derf1, derf2, derf3, &
                                          use_B_r, use_del_tp_B
        use spline_vmec_sub

        implicit none

        integer :: i_rho, i_theta, i_phi, i_qua, k, ist, isp
        real(dp), allocatable :: splcoe_r(:, :), splcoe_t(:, :), splcoe_p(:, :)

        print *, 'Splining Boozer data'

        if (use_del_tp_B) then
            print *, 'Delta theta and Delta phi in Boozer coordinates are splined'
        else
            print *, 'Delta theta and Delta phi in Boozer coordinates are not splined'
        end if

        allocate (splcoe_t(0:ns_tp_B, n_theta_B))
        allocate (splcoe_p(0:ns_tp_B, n_phi_B))

! splining over $\varphi$:

        do i_rho = 1, ns_B
            do i_theta = 1, n_theta_B

                do i_qua = 1, 2
                    splcoe_p(0, :) = s_delt_delp_V(i_qua, 1, 1, 1, i_rho, i_theta, :)

                    call spl_per(ns_tp_B, n_phi_B, h_phi_B, splcoe_p)

                    do k = 1, ns_tp_B
                        s_delt_delp_V(i_qua, 1, 1, k + 1, i_rho, i_theta, :) = splcoe_p(k, :)
                    end do

                    if (use_del_tp_B) then
                        splcoe_p(0, :) = s_delt_delp_B(i_qua, 1, 1, 1, i_rho, i_theta, :)

                        call spl_per(ns_tp_B, n_phi_B, h_phi_B, splcoe_p)

                        do k = 1, ns_tp_B
                            s_delt_delp_B(i_qua, 1, 1, k + 1, i_rho, i_theta, :) = splcoe_p(k, :)
                        end do
                    end if
                end do

                splcoe_p(0, :) = s_Bmod_B(1, 1, 1, i_rho, i_theta, :)

                call spl_per(ns_tp_B, n_phi_B, h_phi_B, splcoe_p)

                do k = 1, ns_tp_B
                    s_Bmod_B(1, 1, k + 1, i_rho, i_theta, :) = splcoe_p(k, :)
                end do

                if (use_B_r) then
                    splcoe_p(0, :) = s_Bcovar_r_B(1, 1, 1, i_rho, i_theta, :)

                    call spl_per(ns_tp_B, n_phi_B, h_phi_B, splcoe_p)

                    do k = 1, ns_tp_B
                        s_Bcovar_r_B(1, 1, k + 1, i_rho, i_theta, :) = splcoe_p(k, :)
                    end do
                end if

            end do
        end do

! splining over $\vartheta$:

        do i_rho = 1, ns_B
            do i_phi = 1, n_phi_B
                do isp = 1, ns_tp_B + 1

                    do i_qua = 1, 2
                        splcoe_t(0, :) = s_delt_delp_V(i_qua, 1, 1, isp, i_rho, :, i_phi)

                        call spl_per(ns_tp_B, n_theta_B, h_theta_B, splcoe_t)

                        do k = 1, ns_tp_B
                            s_delt_delp_V(i_qua, 1, k + 1, isp, i_rho, :, i_phi) = splcoe_t(k, :)
                        end do

                        if (use_del_tp_B) then
                            splcoe_t(0, :) = s_delt_delp_B(i_qua, 1, 1, isp, i_rho, :, i_phi)

                            call spl_per(ns_tp_B, n_theta_B, h_theta_B, splcoe_t)

                            do k = 1, ns_tp_B
                                s_delt_delp_B(i_qua, 1, k + 1, isp, i_rho, :, i_phi) = splcoe_t(k, :)
                            end do
                        end if
                    end do

                    splcoe_t(0, :) = s_Bmod_B(1, 1, isp, i_rho, :, i_phi)

                    call spl_per(ns_tp_B, n_theta_B, h_theta_B, splcoe_t)

                    do k = 1, ns_tp_B
                        s_Bmod_B(1, k + 1, isp, i_rho, :, i_phi) = splcoe_t(k, :)
                    end do

                    if (use_B_r) then
                        splcoe_t(0, :) = s_Bcovar_r_B(1, 1, isp, i_rho, :, i_phi)

                        call spl_per(ns_tp_B, n_theta_B, h_theta_B, splcoe_t)

                        do k = 1, ns_tp_B
                            s_Bcovar_r_B(1, k + 1, isp, i_rho, :, i_phi) = splcoe_t(k, :)
                        end do
                    end if

                end do
            end do
        end do

! splining over $\rho$:

        allocate (splcoe_r(0:ns_s_B, ns_B))

        do i_qua = 1, 2
            splcoe_r(0, :) = s_Bcovar_tp_B(i_qua, 1, :)

            call spl_reg(ns_s_B, ns_B, hs_B, splcoe_r)

            do k = 1, ns_s_B
                s_Bcovar_tp_B(i_qua, k + 1, :) = splcoe_r(k, :)
            end do
        end do

        do i_theta = 1, n_theta_B
            do i_phi = 1, n_phi_B
                do ist = 1, ns_tp_B + 1
                    do isp = 1, ns_tp_B + 1

                        do i_qua = 1, 2
                            splcoe_r(0, :) = s_delt_delp_V(i_qua, 1, ist, isp, :, i_theta, i_phi)

                            call spl_reg(ns_s_B, ns_B, hs_B, splcoe_r)

                            do k = 1, ns_s_B
                                s_delt_delp_V(i_qua, k + 1, ist, isp, :, i_theta, i_phi) = splcoe_r(k, :)
                            end do

                            if (use_del_tp_B) then
                                splcoe_r(0, :) = s_delt_delp_B(i_qua, 1, ist, isp, :, i_theta, i_phi)

                                call spl_reg(ns_s_B, ns_B, hs_B, splcoe_r)

                                do k = 1, ns_s_B
                                    s_delt_delp_B(i_qua, k + 1, ist, isp, :, i_theta, i_phi) = splcoe_r(k, :)
                                end do
                            end if
                        end do

                        splcoe_r(0, :) = s_Bmod_B(1, ist, isp, :, i_theta, i_phi)

                        call spl_reg(ns_s_B, ns_B, hs_B, splcoe_r)

                        do k = 1, ns_s_B
                            s_Bmod_B(k + 1, ist, isp, :, i_theta, i_phi) = splcoe_r(k, :)
                        end do

                        if (use_B_r) then
                            splcoe_r(0, :) = s_Bcovar_r_B(1, ist, isp, :, i_theta, i_phi)

                            call spl_reg(ns_s_B, ns_B, hs_B, splcoe_r)

                            do k = 1, ns_s_B
                                s_Bcovar_r_B(k + 1, ist, isp, :, i_theta, i_phi) = splcoe_r(k, :)
                            end do
                        end if

                    end do
                end do
            end do
        end do

        deallocate (splcoe_r, splcoe_t, splcoe_p)

        do k = 1, ns_max
            derf1(k) = dble(k - 1)
            derf2(k) = dble((k - 1)*(k - 2))
            derf3(k) = dble((k - 1)*(k - 2)*(k - 3))
        end do

        print *, 'done'

    end subroutine spline_boozer_data

    subroutine reset_boozer_batch_splines
        if (bmod_batch_spline_ready) then
            call destroy_batch_splines_3d(bmod_batch_spline)
            bmod_batch_spline_ready = .false.
        end if
        if (allocated(bmod_grid)) then
            deallocate (bmod_grid)
        end if
        if (br_batch_spline_ready) then
            call destroy_batch_splines_3d(br_batch_spline)
            br_batch_spline_ready = .false.
        end if
        if (allocated(br_grid)) then
            deallocate (br_grid)
        end if
    end subroutine reset_boozer_batch_splines

    subroutine build_boozer_bmod_batch_spline
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B

        real(dp) :: x_min(3), x_max(3)
        real(dp), allocatable :: y_batch(:, :, :, :)
        integer :: order(3)
        logical :: periodic(3)

        if (.not. allocated(bmod_grid)) then
            error stop "build_boozer_bmod_batch_spline: bmod_grid not allocated"
        end if

        if (bmod_batch_spline_ready) then
            call destroy_batch_splines_3d(bmod_batch_spline)
            bmod_batch_spline_ready = .false.
        end if

        order = [ns_s_B, ns_tp_B, ns_tp_B]
        if (any(order < 3) .or. any(order > 5)) then
            error stop "build_boozer_bmod_batch_spline: spline order must be 3..5"
        end if

        x_min = [0.0d0, 0.0d0, 0.0d0]
        x_max(1) = hs_B*dble(ns_B - 1)
        x_max(2) = h_theta_B*dble(n_theta_B - 1)
        x_max(3) = h_phi_B*dble(n_phi_B - 1)

        periodic = [.false., .true., .true.]

        allocate (y_batch(ns_B, n_theta_B, n_phi_B, 1))
        y_batch(:, :, :, 1) = bmod_grid(:, :, :)

        call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, &
                                        bmod_batch_spline)
        bmod_batch_spline_ready = .true.
    end subroutine build_boozer_bmod_batch_spline

    subroutine build_boozer_br_batch_spline
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, use_B_r

        real(dp) :: x_min(3), x_max(3)
        real(dp), allocatable :: y_batch(:, :, :, :)
        integer :: order(3)
        logical :: periodic(3)

        if (.not. use_B_r) then
            return
        end if

        if (.not. allocated(br_grid)) then
            error stop "build_boozer_br_batch_spline: br_grid not allocated"
        end if

        if (br_batch_spline_ready) then
            call destroy_batch_splines_3d(br_batch_spline)
            br_batch_spline_ready = .false.
        end if

        order = [ns_s_B, ns_tp_B, ns_tp_B]
        if (any(order < 3) .or. any(order > 5)) then
            error stop "build_boozer_br_batch_spline: spline order must be 3..5"
        end if

        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max(1) = hs_B*dble(ns_B - 1)
        x_max(2) = h_theta_B*dble(n_theta_B - 1)
        x_max(3) = h_phi_B*dble(n_phi_B - 1)

        periodic = [.false., .true., .true.]

        allocate (y_batch(ns_B, n_theta_B, n_phi_B, 1))
        y_batch(:, :, :, 1) = br_grid(:, :, :)

        call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, &
                                        br_batch_spline)
        br_batch_spline_ready = .true.
    end subroutine build_boozer_br_batch_spline

    subroutine normalize_angular_coordinates(vartheta, varphi, n_theta_B, n_phi_B, &
                                             h_theta_B, h_phi_B, &
                                             i_theta, i_phi, dtheta, dphi)
        use new_vmec_stuff_mod, only: nper
        
        implicit none
        
        real(dp), intent(in) :: vartheta, varphi
        integer, intent(in) :: n_theta_B, n_phi_B
        real(dp), intent(in) :: h_theta_B, h_phi_B
        integer, intent(out) :: i_theta, i_phi
        real(dp), intent(out) :: dtheta, dphi
        
        
        dtheta = modulo(vartheta, twopi)/h_theta_B
        i_theta = max(0, min(n_theta_B - 1, int(dtheta)))
        dtheta = (dtheta - dble(i_theta))*h_theta_B
        i_theta = i_theta + 1
        
        dphi = modulo(varphi, twopi/dble(nper))/h_phi_B
        i_phi = max(0, min(n_phi_B - 1, int(dphi)))
        dphi = (dphi - dble(i_phi))*h_phi_B
        i_phi = i_phi + 1
        
    end subroutine normalize_angular_coordinates


end module boozer_sub
