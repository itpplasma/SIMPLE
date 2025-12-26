module boozer_data_mod
    use spl_three_to_five_sub
    use interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                           construct_batch_splines_1d, construct_batch_splines_3d, &
                           destroy_batch_splines_1d, destroy_batch_splines_3d
    use field, only: magnetic_field_t, field_clone
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: boozer_set_current_field
    public :: compute_boozer_data
    public :: reset_boozer_batch_splines
    public :: build_boozer_aphi_batch_spline
    public :: build_boozer_bcovar_tp_batch_spline
    public :: build_boozer_bmod_br_batch_spline
    public :: build_boozer_delt_delp_batch_splines

    public :: current_field
    public :: aphi_batch_spline, aphi_batch_spline_ready
    public :: bcovar_tp_batch_spline, bcovar_tp_batch_spline_ready
    public :: bmod_br_batch_spline, bmod_br_batch_spline_ready
    public :: bmod_br_num_quantities
    public :: delt_delp_V_batch_spline, delt_delp_V_batch_spline_ready
    public :: delt_delp_B_batch_spline, delt_delp_B_batch_spline_ready

    real(dp), parameter :: TWOPI = 2.0_dp*3.14159265358979_dp

    class(magnetic_field_t), allocatable :: current_field
!$omp threadprivate(current_field)

    type(BatchSplineData3D), save :: bmod_br_batch_spline
    logical, save :: bmod_br_batch_spline_ready = .false.
    integer, save :: bmod_br_num_quantities = 0
    real(dp), allocatable, save :: bmod_grid(:, :, :)
    real(dp), allocatable, save :: br_grid(:, :, :)

    type(BatchSplineData1D), save :: aphi_batch_spline
    logical, save :: aphi_batch_spline_ready = .false.

    type(BatchSplineData1D), save :: bcovar_tp_batch_spline
    logical, save :: bcovar_tp_batch_spline_ready = .false.

    type(BatchSplineData3D), save :: delt_delp_V_batch_spline
    logical, save :: delt_delp_V_batch_spline_ready = .false.
    real(dp), allocatable, save :: delt_delp_V_grid(:, :, :, :)

    type(BatchSplineData3D), save :: delt_delp_B_batch_spline
    logical, save :: delt_delp_B_batch_spline_ready = .false.
    real(dp), allocatable, save :: delt_delp_B_grid(:, :, :, :)

contains

    subroutine boozer_set_current_field(field)
        class(magnetic_field_t), intent(in) :: field

        call field_clone(field, current_field)
    end subroutine boozer_set_current_field

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
        if (use_del_tp_B) then
            call ensure_grid_4d(delt_delp_B_grid, ns_B, n_theta_B, n_phi_B, 2)
        end if

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
                                                            A_theta, A_phi, &
                                                            dA_theta_ds, &
                                                            dA_phi_ds, aiota, &
                                                            sqg, alam, dl_ds, &
                                                            dl_dt, dl_dp, &
                                                            Bctrvr_vartheta, &
                                                            Bctrvr_varphi, &
                                                            Bcovar_r, &
                                                            Bcovar_vartheta, &
                                                            Bcovar_varphi)
                    else
                        call vmec_field_evaluate(s, theta, varphi, &
                                                 A_theta, A_phi, dA_theta_ds, &
                                                 dA_phi_ds, aiota, &
                                                 sqg, alam, dl_ds, dl_dt, dl_dp, &
                                                 Bctrvr_vartheta, Bctrvr_varphi, &
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
! covariant components $B_k$ in symmetry flux coordinates
! on equidistant grid of Boozer coordinates:
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
    subroutine compute_br_from_symflux(rho_tor, aiota_arr, Gfunc, Bcovar_symfl)
        use boozer_coordinates_mod, only: ns_B, n_theta_B, n_phi_B
        use plag_coeff_sub, only: plag_coeff

        real(dp), intent(in) :: rho_tor(:)
        real(dp), intent(in) :: aiota_arr(:)
        real(dp), intent(in) :: Gfunc(:, :, :)
        real(dp), intent(in) :: Bcovar_symfl(:, :, :, :)

        integer, parameter :: NPOILAG = 5
        integer, parameter :: NDER = 1

        integer :: i_rho, i_phi, ibeg, iend, nshift
        real(dp) :: coef(0:NDER, NPOILAG)

        nshift = NPOILAG/2

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
                    - matmul(coef(1, :)*aiota_arr(ibeg:iend), &
                             Gfunc(ibeg:iend, :, i_phi)) * &
                      Bcovar_symfl(2, i_rho, :, i_phi) &
                    - matmul(coef(1, :), Gfunc(ibeg:iend, :, i_phi)) &
                      *Bcovar_symfl(3, i_rho, :, i_phi)
            end do
        end do

    end subroutine compute_br_from_symflux
    subroutine ensure_grid_3d(grid, n1, n2, n3)
        real(dp), allocatable, intent(inout) :: grid(:, :, :)
        integer, intent(in) :: n1, n2, n3

        if (.not. allocated(grid)) then
            allocate(grid(n1, n2, n3))
        else if (any(shape(grid) /= [n1, n2, n3])) then
            deallocate(grid)
            allocate(grid(n1, n2, n3))
        end if
    end subroutine ensure_grid_3d
    subroutine ensure_grid_4d(grid, n1, n2, n3, n4)
        real(dp), allocatable, intent(inout) :: grid(:, :, :, :)
        integer, intent(in) :: n1, n2, n3, n4

        if (.not. allocated(grid)) then
            allocate(grid(n1, n2, n3, n4))
        else if (any(shape(grid) /= [n1, n2, n3, n4])) then
            deallocate(grid)
            allocate(grid(n1, n2, n3, n4))
        end if
    end subroutine ensure_grid_4d
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
            bmod_br_num_quantities = 0
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

end module boozer_data_mod
