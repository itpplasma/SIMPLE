program test_chartmap_aphi_abscissa
    !> Unit test for the chartmap radial abscissa contract: every 1D profile in a
    !> Boozer chartmap file (A_phi, B_theta, B_phi) is a function of the file's
    !> uniform rho grid. A_phi must therefore be interpolated on rho and
    !> chain-ruled to s, exactly like B_theta/B_phi/Bmod.
    !>
    !> Regression guard for the bug where A_phi was splined on a synthetic
    !> uniform-s grid: that misplaces A_phi radially (rho is uniform, so s=rho^2
    !> is not), which corrupts iota = -dA_phi_ds/dA_theta_ds and hence the
    !> trapped/passing split. We write a synthetic chartmap with a known A_phi(rho)
    !> that is non-polynomial in s, load it, and check A_phi(rho) and dA_phi/ds
    !> against the analytic profile at interior points.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use boozer_sub, only: load_boozer_from_chartmap, splint_boozer_coord, &
        reset_boozer_batch_splines
    implicit none

    real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: aphi_amp = 0.3_dp, aphi_freq = 1.7_dp
    real(dp), parameter :: torflux_val = 1.0_dp
    integer, parameter :: n_rho = 65, n_theta = 7, n_zeta = 7
    integer, parameter :: n_theta_field = 9, n_zeta_field = 9
    character(len=*), parameter :: fname = 'test_aphi_abscissa.nc'

    real(dp) :: rho(n_rho)
    integer :: i, n_fail
    real(dp) :: rho_min, s, rho_eval
    real(dp) :: A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3
    real(dp) :: B_vth, dB_vth, d2B_vth, B_vph, dB_vph, d2B_vph
    real(dp) :: Bmod, dBmod(3), d2Bmod(6), B_r, dB_r(3), d2B_r(6)
    real(dp) :: aphi_ref, daphi_ds_ref, err_val, err_der
    real(dp), parameter :: tol_val = 1.0e-4_dp, tol_der = 5.0e-3_dp

    rho_min = 0.05_dp
    do i = 1, n_rho
        rho(i) = rho_min + real(i - 1, dp)*(1.0_dp - rho_min)/real(n_rho - 1, dp)
    end do

    call write_synthetic_chartmap()

    call reset_boozer_batch_splines
    call load_boozer_from_chartmap(fname)

    ! Evaluate at interior rho points away from the grid boundary and axis.
    n_fail = 0
    do i = 8, n_rho - 8, 6
        rho_eval = 0.5_dp*(rho(i) + rho(i + 1))   ! off-node
        s = rho_eval**2
        call splint_boozer_coord(s, 0.0_dp, 0.0_dp, 1, &
                                 A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, &
                                 B_vth, dB_vth, d2B_vth, &
                                 B_vph, dB_vph, d2B_vph, &
                                 Bmod, dBmod, d2Bmod, B_r, dB_r, d2B_r)

        aphi_ref = aphi_profile(rho_eval)
        ! dA_phi/ds = (dA_phi/drho)/(2 rho)
        daphi_ds_ref = aphi_amp*aphi_freq*cos(aphi_freq*rho_eval)/(2.0_dp*rho_eval)

        err_val = abs(A_phi - aphi_ref)
        err_der = abs(dA_phi_dr - daphi_ds_ref)/max(abs(daphi_ds_ref), 1.0_dp)

        write (*, '(a,f6.3,a,es12.4,a,es12.4)') &
            'rho=', rho_eval, '  |dA_phi|err=', err_val, '  rel d/ds err=', err_der
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
    end do

    call reset_boozer_batch_splines

    if (n_fail > 0) then
        write (*, '(a,i0,a)') 'test_chartmap_aphi_abscissa: ', n_fail, ' FAILURES'
        error stop 1
    end if
    print *, 'test_chartmap_aphi_abscissa: PASSED'

contains

    pure function aphi_profile(r) result(a)
        real(dp), intent(in) :: r
        real(dp) :: a
        a = aphi_amp*sin(aphi_freq*r)
    end function aphi_profile

    subroutine write_synthetic_chartmap()
        integer :: ncid, did_rho, did_th, did_ze, did_thf, did_zef
        integer :: vid_rho, vid_th, vid_ze, vid_aphi, vid_bth, vid_bph
        integer :: vid_bmod, vid_nfp
        real(dp) :: theta(n_theta), zeta(n_zeta)
        real(dp) :: a_phi_arr(n_rho), b_theta_arr(n_rho), b_phi_arr(n_rho)
        real(dp) :: bmod_arr(n_rho, n_theta_field, n_zeta_field)
        integer :: j

        do j = 1, n_theta
            theta(j) = real(j - 1, dp)*twopi/real(n_theta, dp)
        end do
        do j = 1, n_zeta
            zeta(j) = real(j - 1, dp)*twopi/real(n_zeta, dp)
        end do
        do j = 1, n_rho
            a_phi_arr(j) = aphi_profile(rho(j))
            b_theta_arr(j) = 0.5_dp*rho(j)
            b_phi_arr(j) = 2.0_dp
        end do
        bmod_arr = 1.0_dp

        call nc(nf90_create(fname, nf90_clobber, ncid), 'create')
        call nc(nf90_def_dim(ncid, 'rho', n_rho, did_rho), 'dim rho')
        call nc(nf90_def_dim(ncid, 'theta', n_theta, did_th), 'dim theta')
        call nc(nf90_def_dim(ncid, 'zeta', n_zeta, did_ze), 'dim zeta')
        call nc(nf90_def_dim(ncid, 'theta_field', n_theta_field, did_thf), 'dim thf')
        call nc(nf90_def_dim(ncid, 'zeta_field', n_zeta_field, did_zef), 'dim zef')

        call nc(nf90_def_var(ncid, 'rho', nf90_double, [did_rho], vid_rho), 'var rho')
        call nc(nf90_def_var(ncid, 'theta', nf90_double, [did_th], vid_th), 'var th')
        call nc(nf90_def_var(ncid, 'zeta', nf90_double, [did_ze], vid_ze), 'var ze')
        call nc(nf90_def_var(ncid, 'A_phi', nf90_double, [did_rho], vid_aphi), 'var aphi')
        call nc(nf90_def_var(ncid, 'B_theta', nf90_double, [did_rho], vid_bth), 'var bth')
        call nc(nf90_def_var(ncid, 'B_phi', nf90_double, [did_rho], vid_bph), 'var bph')
        call nc(nf90_def_var(ncid, 'Bmod', nf90_double, &
                             [did_rho, did_thf, did_zef], vid_bmod), 'var bmod')
        call nc(nf90_def_var(ncid, 'num_field_periods', nf90_int, vid_nfp), 'var nfp')

        call nc(nf90_put_att(ncid, nf90_global, 'torflux', torflux_val), 'att torflux')
        call nc(nf90_put_att(ncid, nf90_global, 'boozer_field', 1), 'att boozer')
        call nc(nf90_enddef(ncid), 'enddef')

        call nc(nf90_put_var(ncid, vid_rho, rho), 'put rho')
        call nc(nf90_put_var(ncid, vid_th, theta), 'put th')
        call nc(nf90_put_var(ncid, vid_ze, zeta), 'put ze')
        call nc(nf90_put_var(ncid, vid_aphi, a_phi_arr), 'put aphi')
        call nc(nf90_put_var(ncid, vid_bth, b_theta_arr), 'put bth')
        call nc(nf90_put_var(ncid, vid_bph, b_phi_arr), 'put bph')
        call nc(nf90_put_var(ncid, vid_bmod, bmod_arr), 'put bmod')
        call nc(nf90_put_var(ncid, vid_nfp, 1), 'put nfp')
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
