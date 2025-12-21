module boozer_sub
    use spl_three_to_five_sub
    use interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                           construct_batch_splines_1d, construct_batch_splines_3d, &
                           evaluate_batch_splines_1d_der2, &
                           evaluate_batch_splines_3d_der, evaluate_batch_splines_3d_der2, &
                           destroy_batch_splines_1d, destroy_batch_splines_3d
    use field, only : magnetic_field_t
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    real(dp), parameter :: twopi = 2.0_dp*3.14159265358979_dp

    ! Module variable to store the field for use in subroutines
    class(magnetic_field_t), allocatable :: current_field
    !$omp threadprivate(current_field)

    ! Combined Bmod + Br batch spline (1 or 2 quantities depending on use_B_r)
    type(BatchSplineData3D), save :: bmod_br_batch_spline
    logical, save :: bmod_br_batch_spline_ready = .false.
    integer, save :: bmod_br_num_quantities = 0
    real(dp), allocatable, save :: bmod_grid(:, :, :)
    real(dp), allocatable, save :: br_grid(:, :, :)

    type(BatchSplineData1D), save :: aphi_batch_spline
    logical, save :: aphi_batch_spline_ready = .false.

    type(BatchSplineData1D), save :: bcovar_tp_batch_spline
    logical, save :: bcovar_tp_batch_spline_ready = .false.

    ! Combined delta_theta + delta_phi grids (2 quantities each)
    type(BatchSplineData3D), save :: delt_delp_V_batch_spline
    logical, save :: delt_delp_V_batch_spline_ready = .false.
    real(dp), allocatable, save :: delt_delp_V_grid(:, :, :, :)

    type(BatchSplineData3D), save :: delt_delp_B_batch_spline
    logical, save :: delt_delp_B_batch_spline_ready = .false.
    real(dp), allocatable, save :: delt_delp_B_grid(:, :, :, :)

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

        call build_boozer_aphi_batch_spline
        call build_boozer_bcovar_tp_batch_spline
        call build_boozer_bmod_br_batch_spline
        call build_boozer_delt_delp_batch_splines

    end subroutine get_boozer_coordinates_impl

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine splint_boozer_coord(r, vartheta_B, varphi_B, &
                                   A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                   d2A_phi_dr2, d3A_phi_dr3, &
                                   B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                                   B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                                   Bmod_B, dBmod_B, d2Bmod_B, &
                                   B_r, dB_r, d2B_r)

        use boozer_coordinates_mod, only: ns_s_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          use_B_r
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper
        use chamb_mod, only: rnegflag
        use diag_mod, only: icounter

        implicit none

        integer, parameter :: mode_secders = 1

        integer :: ns_s_p1
        integer :: is, i_theta, i_phi

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
        real(dp) :: x_eval(3), y_eval(2), dy_eval(3, 2), d2y_eval(6, 2)
        real(dp) :: theta_wrapped, phi_wrapped
        real(dp) :: y1d(2), dy1d(2), d2y1d(2)

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

! Interpolation of A_phi over s (batch spline 1D)
        if (.not. aphi_batch_spline_ready) then
            error stop "splint_boozer_coord: Aphi batch spline not initialized"
        end if

        call evaluate_batch_splines_1d_der2(aphi_batch_spline, r, y1d(1:1), &
                                            dy1d(1:1), d2y1d(1:1))
        A_phi = y1d(1)
        dA_phi_dr = dy1d(1)
        d2A_phi_dr2 = d2y1d(1)
        if (mode_secders .gt. 0) then
            call evaluate_batch_splines_1d_der3_single(aphi_batch_spline, r, d3A_phi_dr3)
        else
            d3A_phi_dr3 = 0.0_dp
        end if

!--------------------------------

        rho_tor = sqrt(r)
        ds = rho_tor/hs_B
        is = max(0, min(ns_B - 1, int(ds)))
        ds = (ds - dble(is))*hs_B
        is = is + 1

        ns_s_p1 = ns_s_B + 1

!--------------------------------
! Interpolation of mod-B (and B_r if use_B_r):
!--------------------------------

        theta_wrapped = modulo(vartheta_B, twopi)
        phi_wrapped = modulo(varphi_B, twopi/dble(nper))

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

! Coversion coefficients for derivatives over s
        drhods = 0.5d0/rho_tor
        drhods2 = drhods**2
        ! $-\dr^2 \rho / \rd s^2$ (second derivative with minus sign)
        d2rhods2m = drhods2/rho_tor

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
! End Interpolation of mod-B and B_r
!--------------------------------
! Interpolation of B_\vartheta and B_\varphi (flux functions, batch spline 1D):
!--------------------------------

        if (.not. bcovar_tp_batch_spline_ready) then
            error stop "splint_boozer_coord: Bcovar_tp batch spline not initialized"
        end if

        call evaluate_batch_splines_1d_der2(bcovar_tp_batch_spline, rho_tor, y1d, dy1d, d2y1d)
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

!--------------------------------
! End interpolation of B_\vartheta and B_\varphi
!--------------------------------
    end subroutine splint_boozer_coord

    subroutine evaluate_batch_splines_1d_der3_single(spl, x, d3y)
        type(BatchSplineData1D), intent(in) :: spl
        real(dp), intent(in) :: x
        real(dp), intent(out) :: d3y

        real(dp) :: x_norm, x_local, xj
        integer :: interval_index, k_power
        integer :: N

        N = spl%order
        if (N < 3) then
            d3y = 0.0_dp
            return
        end if

        if (spl%periodic) then
            xj = modulo(x - spl%x_min, spl%h_step*(spl%num_points - 1)) + spl%x_min
        else
            xj = x
        end if

        x_norm = (xj - spl%x_min)/spl%h_step
        interval_index = max(0, min(spl%num_points - 2, int(x_norm)))
        x_local = (x_norm - dble(interval_index))*spl%h_step

        d3y = 0.0_dp
        if (N >= 3) then
            d3y = N*(N - 1)*(N - 2)*spl%coeff(1, N, interval_index + 1)
            do k_power = N - 1, 3, -1
                d3y = k_power*(k_power - 1)*(k_power - 2)* &
                      spl%coeff(1, k_power, interval_index + 1) + x_local*d3y
            end do
        end if
    end subroutine evaluate_batch_splines_1d_der3_single

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine delthe_delphi_BV(isw, r, vartheta, varphi, deltheta_BV, delphi_BV, &
                                ddeltheta_BV, ddelphi_BV)

! Computes $\Delta \vartheta = \vartheta_B - \theta_V$ and
! $\Delta \varphi = \varphi_B - \varphi_V$
! and their first derivatives over angles for two cases:
! isw=0 - if they are given as functions of VMEC coordinates (r,vartheta,varphi)
! isw=1 - if they are given as functions of Boozer coordinates (r,vartheta,varphi)

        use boozer_coordinates_mod, only: use_del_tp_B
        use chamb_mod, only: rnegflag

        implicit none

        integer, intent(in) :: isw
        real(dp), intent(inout) :: r
        real(dp), intent(in) :: vartheta, varphi
        real(dp), intent(out) :: deltheta_BV, delphi_BV
        real(dp), dimension(2), intent(out) :: ddeltheta_BV, ddelphi_BV

        real(dp) :: rho_tor, x_eval(3), y_eval(2), dy_eval(3, 2)

        if (r .le. 0.0_dp) then
            rnegflag = .true.
            r = abs(r)
        end if

        rho_tor = sqrt(r)
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

! allocate data arrays for Boozer data:
        if (.not. allocated(s_Bcovar_tp_B)) &
            allocate (s_Bcovar_tp_B(2, ns_s_B + 1, ns_B))

        ! Allocate simple 3D grid for Bmod
        if (.not. allocated(bmod_grid)) then
            allocate (bmod_grid(ns_B, n_theta_B, n_phi_B))
        else if (any(shape(bmod_grid) /= [ns_B, n_theta_B, n_phi_B])) then
            deallocate (bmod_grid)
            allocate (bmod_grid(ns_B, n_theta_B, n_phi_B))
        end if

        ! Allocate simple 3D grid for B_r (if needed)
        if (use_B_r) then
            if (.not. allocated(br_grid)) then
                allocate (br_grid(ns_B, n_theta_B, n_phi_B))
            else if (any(shape(br_grid) /= [ns_B, n_theta_B, n_phi_B])) then
                deallocate (br_grid)
                allocate (br_grid(ns_B, n_theta_B, n_phi_B))
            end if
        end if

        ! Allocate simple 4D grids for delta_theta/delta_phi
        if (.not. allocated(delt_delp_V_grid)) then
            allocate (delt_delp_V_grid(ns_B, n_theta_B, n_phi_B, 2))
        else if (any(shape(delt_delp_V_grid) /= [ns_B, n_theta_B, n_phi_B, 2])) then
            deallocate (delt_delp_V_grid)
            allocate (delt_delp_V_grid(ns_B, n_theta_B, n_phi_B, 2))
        end if

        if (use_del_tp_B) then
            if (.not. allocated(delt_delp_B_grid)) then
                allocate (delt_delp_B_grid(ns_B, n_theta_B, n_phi_B, 2))
            else if (any(shape(delt_delp_B_grid) /= [ns_B, n_theta_B, n_phi_B, 2])) then
                deallocate (delt_delp_B_grid)
                allocate (delt_delp_B_grid(ns_B, n_theta_B, n_phi_B, 2))
            end if
        end if

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

            if (use_del_tp_B) then
                delt_delp_B_grid(i_rho, :, :, 1) = perqua_2D(1, :, :)
                delt_delp_B_grid(i_rho, :, :, 2) = perqua_2D(2, :, :)
            end if
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
! Write directly to br_grid

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
                    br_grid(i_rho, :, i_phi) = &
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
        end if

! End compute radial covariant magnetic field component in Boozer coordinates

        deallocate (Bcovar_theta_V, Bcovar_varphi_V, bmod_Vg, alam_2D, &
                    deltheta_BV_Vg, delphi_BV_Vg, &
                    wint_t, wint_p, coef, theta_V, theta_B, phi_V, phi_B, &
                    perqua_t, perqua_p, perqua_2D)

        print *, 'done'

    end subroutine compute_boozer_data

    subroutine reset_boozer_batch_splines
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
            bmod_br_num_quantities = 0
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
    end subroutine reset_boozer_batch_splines

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
        x_max = hs_B*dble(ns_B - 1)

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
            bmod_br_num_quantities = 0
        end if

        order = [ns_s_B, ns_tp_B, ns_tp_B]
        if (any(order < 3) .or. any(order > 5)) then
            error stop "build_boozer_bmod_br_batch_spline: spline order must be 3..5"
        end if

        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max(1) = hs_B*dble(ns_B - 1)
        x_max(2) = h_theta_B*dble(n_theta_B - 1)
        x_max(3) = h_phi_B*dble(n_phi_B - 1)

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
        bmod_br_num_quantities = nq
        deallocate (y_batch)
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
            error stop "build_boozer_delt_delp_batch_splines: delt_delp_V_grid not allocated"
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
        x_max(1) = hs_B*dble(ns_B - 1)
        x_max(2) = h_theta_B*dble(n_theta_B - 1)
        x_max(3) = h_phi_B*dble(n_phi_B - 1)

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
