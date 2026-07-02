program test_chartmap_aphi_abscissa
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use boozer_chartmap_io, only: boozer_chartmap_data_t, read_boozer_chartmap
    use boozer_chartmap, only: load_boozer_from_chartmap
    use boozer_sub, only: splint_boozer_coord, reset_boozer_batch_splines
    implicit none

    real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: aphi_c0 = 0.2_dp, aphi_c1 = 0.3_dp
    real(dp), parameter :: aphi_c2 = -0.4_dp, aphi_c3 = 0.15_dp
    real(dp), parameter :: torflux_val = 1.0_dp
    ! The reader derives rmajor as the (theta,zeta)-average of sqrt(x^2+y^2)
    ! on the innermost rho surface (cm -> m). For the circular-torus geometry
    ! below that average is exactly r0_cm/1e2 on the uniform theta grid.
    real(dp), parameter :: rmajor_val = 7.5_dp
    real(dp), parameter :: r0_cm = rmajor_val*1.0e2_dp, a_cm = 1.0e2_dp
    integer, parameter :: n_rho = 65, n_theta = 8, n_zeta = 8
    integer, parameter :: nfp = 2
    character(len=*), parameter :: fname = 'test_aphi_abscissa.nc'

    real(dp) :: rho(n_rho)
    real(dp) :: s_grid(n_rho)
    integer :: i, n_fail
    real(dp) :: rho_min, s, rho_eval
    real(dp) :: A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3
    real(dp) :: B_vth, dB_vth, d2B_vth, B_vph, dB_vph, d2B_vph
    real(dp) :: Bmod, dBmod(3), d2Bmod(6), B_r, dB_r(3), d2B_r(6)
    real(dp) :: aphi_ref, daphi_ds_ref, d2_ref, d3_ref, err_val, err_der
    real(dp) :: err_d2, err_d3, rr
    real(dp), parameter :: tol_val = 1.0e-13_dp, tol_der = 1.0e-13_dp
    real(dp), parameter :: tol_d2 = 5.0e-10_dp, tol_d3 = 2.0e-9_dp

    rho_min = 0.05_dp
    do i = 1, n_rho
        rho(i) = rho_min + real(i - 1, dp)*(1.0_dp - rho_min)/real(n_rho - 1, dp)
        s_grid(i) = rho_min**2 + real(i - 1, dp)*(1.0_dp - rho_min**2) &
                    /real(n_rho - 1, dp)
    end do

    call write_synthetic_chartmap()
    call check_reader_contract(n_fail)

    call reset_boozer_batch_splines
    call load_boozer_from_chartmap(fname)

    do i = 8, n_rho - 8, 6
        rho_eval = 0.5_dp*(rho(i) + rho(i + 1))   ! off-node
        s = rho_eval**2
        call splint_boozer_coord(s, 0.0_dp, 0.0_dp, 1, &
                                 A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, &
                                 B_vth, dB_vth, d2B_vth, &
                                 B_vph, dB_vph, d2B_vph, &
                                 Bmod, dBmod, d2Bmod, B_r, dB_r, d2B_r)

        rr = rho_eval
        aphi_ref = aphi_profile(s)
        daphi_ds_ref = aphi_c1 + 2.0_dp*aphi_c2*s + 3.0_dp*aphi_c3*s**2
        d2_ref = 2.0_dp*aphi_c2 + 6.0_dp*aphi_c3*s
        d3_ref = 6.0_dp*aphi_c3

        err_val = abs(A_phi - aphi_ref)
        err_der = abs(dA_phi_dr - daphi_ds_ref)/max(abs(daphi_ds_ref), 1.0_dp)
        err_d2 = abs(d2A_phi_dr2 - d2_ref)/max(abs(d2_ref), 1.0_dp)
        err_d3 = abs(d3A_phi_dr3 - d3_ref)/max(abs(d3_ref), 1.0_dp)

        write (*, '(a,f6.3,4(a,es11.3))') &
            'rho=', rr, '  Aphi err=', err_val, '  d/ds=', err_der, &
            '  d2/ds2=', err_d2, '  d3/ds3=', err_d3
        if (err_val > tol_val) then
            write (*, '(a,es12.4,a,es12.4)') '  FAIL A_phi value: got ', A_phi, &
                ' expected ', aphi_ref
            n_fail = n_fail + 1
        end if
        if (err_der > tol_der) then
            write (*, '(a,es12.4,a,es12.4)') '  FAIL dA_phi/ds: got ', dA_phi_dr, &
                ' expected ', daphi_ds_ref
            n_fail = n_fail + 1
        end if
        if (err_d2 > tol_d2) then
            write (*, '(a,es12.4,a,es12.4)') '  FAIL d2A_phi/ds2: got ', d2A_phi_dr2, &
                ' expected ', d2_ref
            n_fail = n_fail + 1
        end if
        if (err_d3 > tol_d3) then
            write (*, '(a,es12.4,a,es12.4)') '  FAIL d3A_phi/ds3: got ', d3A_phi_dr3, &
                ' expected ', d3_ref
            n_fail = n_fail + 1
        end if
    end do

    call reset_boozer_batch_splines

    if (n_fail > 0) then
        write (*, '(a,i0,a)') 'test_chartmap_aphi_abscissa: ', n_fail, ' FAILURES'
        error stop 1
    end if
    print *, 'test_chartmap_aphi_abscissa: PASSED'

contains

    pure function aphi_profile(s) result(a)
        real(dp), intent(in) :: s
        real(dp) :: a
        a = aphi_c0 + aphi_c1*s + aphi_c2*s**2 + aphi_c3*s**3
    end function aphi_profile

    pure function bmod_model(ir, it, iz) result(b)
        integer, intent(in) :: ir, it, iz
        real(dp) :: b
        b = 1.0_dp + 0.1_dp*ir + 0.01_dp*it + 0.001_dp*iz
    end function bmod_model

    subroutine check_reader_contract(n_fail)
        integer, intent(inout) :: n_fail
        type(boozer_chartmap_data_t) :: d

        n_fail = 0
        call read_boozer_chartmap(fname, d)
        call expect_int(d%n_rho, n_rho, 'n_rho', n_fail)
        call expect_int(d%n_theta, n_theta + 1, 'n_theta spline grid', n_fail)
        call expect_int(d%n_phi, n_zeta + 1, 'n_phi spline grid', n_fail)
        call expect_int(d%nfp, nfp, 'nfp', n_fail)
        call expect_real(d%torflux, torflux_val, 'torflux', n_fail)
        call expect_real(d%rmajor, rmajor_val, 'rmajor (geometry-derived)', n_fail)
        call expect_int(d%n_s, n_rho, 'n_s', n_fail)
        call expect_real(d%s(1), rho_min**2, 's min', n_fail)
        call expect_real(d%s(n_rho), 1.0_dp, 's max', n_fail)
        call expect_real(real(d%n_theta - 1, dp)*d%h_theta, twopi, &
                         'poloidal period', n_fail)
        call expect_real(real(d%n_phi - 1, dp)*d%h_phi, twopi/real(nfp, dp), &
                         'toroidal period', n_fail)
        call expect_real(d%Bmod(2, 3, 4), bmod_model(2, 3, 4), 'Bmod(2,3,4)', n_fail)
        call expect_real(d%Bmod(2, n_theta + 1, 4), bmod_model(2, 1, 4), &
                         'Bmod theta endpoint copy', n_fail)
        call expect_real(d%Bmod(2, 3, n_zeta + 1), bmod_model(2, 3, 1), &
                         'Bmod zeta endpoint copy', n_fail)
        call expect_real(d%Bmod(2, n_theta + 1, n_zeta + 1), bmod_model(2, 1, 1), &
                         'Bmod corner endpoint copy', n_fail)
    end subroutine check_reader_contract

    subroutine expect_int(got, want, label, n_fail)
        integer, intent(in) :: got, want
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (got /= want) then
            print *, 'FAIL: ', label, ' got ', got, ' want ', want
            n_fail = n_fail + 1
        end if
    end subroutine expect_int

    subroutine expect_real(got, want, label, n_fail)
        real(dp), intent(in) :: got, want
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (abs(got - want) > 1.0e-12_dp*max(1.0_dp, abs(want))) then
            print *, 'FAIL: ', label, ' got ', got, ' want ', want
            n_fail = n_fail + 1
        end if
    end subroutine expect_real

    subroutine write_synthetic_chartmap()
        integer :: ncid, did_rho, did_s, did_th, did_ze
        integer :: vid_rho, vid_s, vid_th, vid_ze, vid_aphi, vid_bth, vid_bph
        integer :: vid_bmod, vid_nfp, vid_x, vid_y, vid_z
        real(dp) :: theta(n_theta), zeta(n_zeta)
        real(dp) :: a_phi_arr(n_rho), b_theta_arr(n_rho), b_phi_arr(n_rho)
        real(dp) :: bmod_arr(n_rho, n_theta, n_zeta)
        real(dp) :: x_arr(n_rho, n_theta, n_zeta), y_arr(n_rho, n_theta, n_zeta)
        real(dp) :: z_arr(n_rho, n_theta, n_zeta), radius
        integer :: ir, it, iz, j

        do j = 1, n_theta
            theta(j) = real(j - 1, dp)*twopi/real(n_theta, dp)
        end do
        do j = 1, n_zeta
            zeta(j) = real(j - 1, dp)*twopi/real(nfp*n_zeta, dp)
        end do
        do j = 1, n_rho
            a_phi_arr(j) = aphi_profile(s_grid(j))
            b_theta_arr(j) = 0.5_dp*rho(j)
            b_phi_arr(j) = 2.0_dp
        end do
        do iz = 1, n_zeta
            do it = 1, n_theta
                do ir = 1, n_rho
                    bmod_arr(ir, it, iz) = bmod_model(ir, it, iz)
                end do
            end do
        end do
        do iz = 1, n_zeta
            do it = 1, n_theta
                do ir = 1, n_rho
                    radius = r0_cm + a_cm*rho(ir)*cos(theta(it))
                    x_arr(ir, it, iz) = radius*cos(zeta(iz))
                    y_arr(ir, it, iz) = radius*sin(zeta(iz))
                    z_arr(ir, it, iz) = a_cm*rho(ir)*sin(theta(it))
                end do
            end do
        end do

        call nc(nf90_create(fname, nf90_clobber, ncid), 'create')
        call nc(nf90_def_dim(ncid, 'rho', n_rho, did_rho), 'dim rho')
        call nc(nf90_def_dim(ncid, 's', n_rho, did_s), 'dim s')
        call nc(nf90_def_dim(ncid, 'theta', n_theta, did_th), 'dim theta')
        call nc(nf90_def_dim(ncid, 'zeta', n_zeta, did_ze), 'dim zeta')

        call nc(nf90_def_var(ncid, 'rho', nf90_double, [did_rho], vid_rho), 'var rho')
        call nc(nf90_def_var(ncid, 's', nf90_double, [did_s], vid_s), 'var s')
        call nc(nf90_def_var(ncid, 'theta', nf90_double, [did_th], vid_th), 'var th')
        call nc(nf90_def_var(ncid, 'zeta', nf90_double, [did_ze], vid_ze), 'var ze')
        call nc(nf90_def_var(ncid, 'A_phi', nf90_double, [did_s], &
                             vid_aphi), 'var aphi')
        call nc(nf90_put_att(ncid, vid_aphi, 'radial_abscissa', 's'), 'att aphi')
        call nc(nf90_def_var(ncid, 'B_theta', nf90_double, [did_rho], vid_bth), 'var bth')
        call nc(nf90_def_var(ncid, 'B_phi', nf90_double, [did_rho], vid_bph), 'var bph')
        call nc(nf90_def_var(ncid, 'x', nf90_double, &
                             [did_rho, did_th, did_ze], vid_x), 'var x')
        call nc(nf90_def_var(ncid, 'y', nf90_double, &
                             [did_rho, did_th, did_ze], vid_y), 'var y')
        call nc(nf90_def_var(ncid, 'z', nf90_double, &
                             [did_rho, did_th, did_ze], vid_z), 'var z')
        call nc(nf90_def_var(ncid, 'Bmod', nf90_double, &
                             [did_rho, did_th, did_ze], vid_bmod), 'var bmod')
        call nc(nf90_def_var(ncid, 'num_field_periods', nf90_int, vid_nfp), 'var nfp')

        call nc(nf90_put_att(ncid, nf90_global, 'torflux', torflux_val), 'att torflux')
        call nc(nf90_put_att(ncid, nf90_global, 'boozer_field', 1), 'att boozer')
        call nc(nf90_enddef(ncid), 'enddef')

        call nc(nf90_put_var(ncid, vid_rho, rho), 'put rho')
        call nc(nf90_put_var(ncid, vid_s, s_grid), 'put s')
        call nc(nf90_put_var(ncid, vid_th, theta), 'put th')
        call nc(nf90_put_var(ncid, vid_ze, zeta), 'put ze')
        call nc(nf90_put_var(ncid, vid_x, x_arr), 'put x')
        call nc(nf90_put_var(ncid, vid_y, y_arr), 'put y')
        call nc(nf90_put_var(ncid, vid_z, z_arr), 'put z')
        call nc(nf90_put_var(ncid, vid_aphi, a_phi_arr), 'put aphi')
        call nc(nf90_put_var(ncid, vid_bth, b_theta_arr), 'put bth')
        call nc(nf90_put_var(ncid, vid_bph, b_phi_arr), 'put bph')
        call nc(nf90_put_var(ncid, vid_bmod, bmod_arr), 'put bmod')
        call nc(nf90_put_var(ncid, vid_nfp, nfp), 'put nfp')
        call nc(nf90_close(ncid), 'close')
    end subroutine write_synthetic_chartmap

    subroutine nc(status, loc)
        integer, intent(in) :: status
        character(len=*), intent(in) :: loc
        if (status /= nf90_noerr) then
            print *, 'NetCDF error at ', trim(loc), ': ', trim(nf90_strerror(status))
            error stop 2
        end if
    end subroutine nc

end program test_chartmap_aphi_abscissa
