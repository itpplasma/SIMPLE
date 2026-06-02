program test_boozer_vmec_bfield_match
    !> Regression test for the physical equivalence of |B| between the raw VMEC
    !> field and the Boozer reconstruction used by the chartmap path.
    !>
    !> bmin/bmax over a flux surface are coordinate invariant, so the two
    !> representations must agree to interpolation accuracy. This guards against
    !> regressions in the Boozer transform / spline grid that would shift the
    !> trapped/passing boundary (the symptom investigated for the chartmap runs).
    !> Both the true surface extrema (dense scan) and the field-line-trace
    !> empirical extrema (as used by init_starting_surf) are checked.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: netcdffile, multharm, ns_A, ns_s, ns_tp
    use parmot_mod, only: rmu
    use velo_mod, only: isw_field_type
    use boozer_coordinates_mod, only: use_B_r
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, &
        reset_boozer_batch_splines
    use spline_vmec_sub, only: spline_vmec_data, vmec_field
    use vmecin_sub, only: stevvo
    use magfie_sub, only: init_magfie, VMEC, BOOZER
    use alpha_lifetime_sub, only: integrate_mfl_can

    implicit none

    real(dp), parameter :: twopi = 8.0_dp * atan(1.0_dp)
    integer, parameter :: nth = 256, nph = 256
    real(dp), parameter :: tol_surface = 1.0e-3_dp
    real(dp), parameter :: tol_fieldline = 1.0e-3_dp
    real(dp), parameter :: s0 = 0.3_dp

    character(len=256) :: wout_file
    integer :: nargs, nfail
    real(dp) :: bmaxV, bminV, bmaxB, bminB, relmax, relmin

    nargs = command_argument_count()
    if (nargs >= 1) then
        call get_command_argument(1, wout_file)
    else
        wout_file = 'wout.nc'
    end if
    nfail = 0

    isw_field_type = 2
    rmu = 1.0e8_dp
    ns_A = 5
    ns_s = 5
    ns_tp = 5
    use_B_r = .false.
    multharm = 5

    netcdffile = wout_file
    call reset_boozer_batch_splines
    call spline_vmec_data
    call get_boozer_coordinates

    call surface_extrema(s0, bmaxV, bminV, bmaxB, bminB)
    relmax = abs(bmaxB - bmaxV) / abs(bmaxV)
    relmin = abs(bminB - bminV) / abs(bminV)
    print '(a,3es16.8)', 'surface bmax (VMEC, Boozer, rel): ', bmaxV, bmaxB, relmax
    print '(a,3es16.8)', 'surface bmin (VMEC, Boozer, rel): ', bminV, bminB, relmin
    if (relmax > tol_surface) then
        print *, 'FAIL: surface bmax mismatch ', relmax, ' > ', tol_surface
        nfail = nfail + 1
    end if
    if (relmin > tol_surface) then
        print *, 'FAIL: surface bmin mismatch ', relmin, ' > ', tol_surface
        nfail = nfail + 1
    end if

    call fieldline_extrema(s0, bmaxV, bminV, bmaxB, bminB)
    relmax = abs(bmaxB - bmaxV) / abs(bmaxV)
    relmin = abs(bminB - bminV) / abs(bminV)
    print '(a,3es16.8)', 'fieldline bmax (VMEC, Boozer, rel): ', bmaxV, bmaxB, relmax
    print '(a,3es16.8)', 'fieldline bmin (VMEC, Boozer, rel): ', bminV, bminB, relmin
    if (relmax > tol_fieldline) then
        print *, 'FAIL: field-line bmax mismatch ', relmax, ' > ', tol_fieldline
        nfail = nfail + 1
    end if
    if (relmin > tol_fieldline) then
        print *, 'FAIL: field-line bmin mismatch ', relmin, ' > ', tol_fieldline
        nfail = nfail + 1
    end if

    if (nfail /= 0) then
        print *, nfail, ' B-field match tests failed'
        error stop 'test_boozer_vmec_bfield_match failed'
    end if
    print *, 'PASS: Boozer |B| matches raw VMEC on the surface and along a field line'

contains

    subroutine surface_extrema(s, bmaxV, bminV, bmaxB, bminB)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: bmaxV, bminV, bmaxB, bminB
        integer :: i, j
        real(dp) :: th, ph, hth, hph, b

        hth = twopi / real(nth - 1, dp)
        hph = twopi / real(nph - 1, dp)
        bmaxV = -1.0e30_dp; bminV = 1.0e30_dp
        bmaxB = -1.0e30_dp; bminB = 1.0e30_dp
        do i = 1, nth
            th = real(i - 1, dp) * hth
            do j = 1, nph
                ph = real(j - 1, dp) * hph
                b = bmod_vmec(s, th, ph)
                bmaxV = max(bmaxV, b); bminV = min(bminV, b)
                b = bmod_boozer(s, th, ph)
                bmaxB = max(bmaxB, b); bminB = min(bminB, b)
            end do
        end do
    end subroutine surface_extrema

    subroutine fieldline_extrema(s, bmaxV, bminV, bmaxB, bminB)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: bmaxV, bminV, bmaxB, bminB
        integer, parameter :: npoiper = 100, nfl = 1000
        integer :: npoi, ierr, L1i
        real(dp) :: dphi, bmod00, RT0, R0i, cbfi, bz0i, bf0
        real(dp), allocatable :: xstart(:, :), bstart(:), volstart(:)

        call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0)
        npoi = nfl * npoiper
        dphi = twopi / real(L1i * npoiper, dp)
        allocate (xstart(3, npoi), bstart(npoi), volstart(npoi))

        call init_magfie(VMEC)
        ierr = 0
        call integrate_mfl_can(npoi, dphi, s, 0.0_dp, 0.0_dp, &
                               xstart, bstart, volstart, bmod00, ierr)
        bmaxV = maxval(bstart); bminV = minval(bstart)

        call init_magfie(BOOZER)
        ierr = 0
        call integrate_mfl_can(npoi, dphi, s, 0.0_dp, 0.0_dp, &
                               xstart, bstart, volstart, bmod00, ierr)
        bmaxB = maxval(bstart); bminB = minval(bstart)
        deallocate (xstart, bstart, volstart)
    end subroutine fieldline_extrema

    function bmod_vmec(s, theta, varphi) result(b)
        real(dp), intent(in) :: s, theta, varphi
        real(dp) :: b
        real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota
        real(dp) :: sqg, alam, dl_ds, dl_dt, dl_dp
        real(dp) :: Bctrvr_vartheta, Bctrvr_varphi
        real(dp) :: Bcovar_s, Bcovar_vartheta, Bcovar_varphi

        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, &
            dA_phi_ds, aiota, sqg, alam, dl_ds, dl_dt, dl_dp, &
            Bctrvr_vartheta, Bctrvr_varphi, Bcovar_s, Bcovar_vartheta, &
            Bcovar_varphi)
        b = sqrt(Bctrvr_vartheta * Bcovar_vartheta + &
                 Bctrvr_varphi * Bcovar_varphi)
    end function bmod_vmec

    function bmod_boozer(s, vartheta_B, varphi_B) result(b)
        real(dp), intent(in) :: s, vartheta_B, varphi_B
        real(dp) :: b
        real(dp) :: A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3
        real(dp) :: B_vartheta_B, dB_vartheta_B, d2B_vartheta_B
        real(dp) :: B_varphi_B, dB_varphi_B, d2B_varphi_B
        real(dp) :: Bmod_B, B_r
        real(dp), dimension(3) :: dBmod_B, dB_r
        real(dp), dimension(6) :: d2Bmod_B, d2B_r

        call splint_boozer_coord(s, vartheta_B, varphi_B, 0, &
            A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3, &
            B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
            B_varphi_B, dB_varphi_B, d2B_varphi_B, &
            Bmod_B, dBmod_B, d2Bmod_B, B_r, dB_r, d2B_r)
        b = Bmod_B
    end function bmod_boozer

end program test_boozer_vmec_bfield_match
