module beam_stopping
    !> Neutral hydrogen beam stopping cross section of S. Suzuki et al.,
    !> Plasma Phys. Control. Fusion 40 (1998) 2097, doi:10.1088/0741-3335/40/12/009,
    !> Eqs. (24), (26) and (28) with Tables 2 and 3 (B = 5 T columns).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use finite_check, only: is_finite
    implicit none
    private
    public :: suzuki_stopping, suzuki_impurity_supported

    ! A(1:10, beam isotope H, D, T); Table 3(a), 10 <= E < 100 keV/amu.
    real(dp), parameter :: a_low(10, 3) = reshape([ &
                                   -52.9_dp, -1.36_dp, 0.0719_dp, 0.0137_dp, 0.454_dp, &
                              0.403_dp, -0.220_dp, 0.0666_dp, -0.0677_dp, -1.48e-3_dp, &
                                   -67.9_dp, -1.22_dp, 0.0814_dp, 0.0139_dp, 0.454_dp, &
                              0.465_dp, -0.273_dp, 0.0751_dp, -0.0630_dp, -5.08e-4_dp, &
                                   -74.2_dp, -1.18_dp, 0.0843_dp, 0.0139_dp, 0.453_dp, &
                      0.491_dp, -0.294_dp, 0.0788_dp, -0.0612_dp, -1.85e-4_dp], [10, 3])
    real(dp), parameter, public :: e_min_kev_amu = 9.0_dp

    ! Table 2(a), 100 <= E <= 10000 keV/amu.
    real(dp), parameter :: a_high(10, 3) = reshape([ &
                                      12.7_dp, 1.25_dp, 0.452_dp, 0.0105_dp, 0.547_dp, &
                              -0.102_dp, 0.360_dp, -0.0298_dp, -0.0959_dp, 4.21e-3_dp, &
                                      14.1_dp, 1.11_dp, 0.408_dp, 0.0105_dp, 0.547_dp, &
                             -0.0403_dp, 0.345_dp, -0.0288_dp, -0.0971_dp, 4.74e-3_dp, &
                                      12.7_dp, 1.26_dp, 0.449_dp, 0.0105_dp, 0.547_dp, &
                    -0.00577_dp, 0.336_dp, -0.0282_dp, -0.0974_dp, 4.87e-3_dp], [10, 3])

    ! Impurity charges with tables: He and B from Tables 2(b) and 3(b),
    ! C, N, O and Fe from Tables 2(c) and 3(c) (1 < Zeff < 5).
    integer, parameter :: n_imp = 6
    integer, parameter :: z_imp(n_imp) = [2, 5, 6, 7, 8, 26]
    ! B(i,j,k) stored in the order B111, B112, B121, B122, B211, ..., B322.
    real(dp), parameter :: b_low(12, n_imp) = reshape([ &
            -7.92e-1_dp, 4.20e-2_dp, 5.30e-2_dp, -1.39e-2_dp, 3.01e-1_dp, -2.64e-2_dp, &
             -2.99e-2_dp, 6.07e-3_dp, 2.72e-4_dp, 6.11e-3_dp, 3.47e-3_dp, -9.19e-4_dp, &
           1.22e-1_dp, 5.27e-2_dp, -4.30e-4_dp, -3.18e-3_dp, -1.51e-1_dp, -3.64e-2_dp, &
             3.43e-3_dp, 1.51e-3_dp, 4.20e-2_dp, 6.92e-3_dp, -1.41e-3_dp, -2.90e-4_dp, &
           1.61e-1_dp, 5.98e-2_dp, -3.36e-3_dp, -4.26e-3_dp, -1.57e-1_dp, -3.96e-2_dp, &
             4.60e-3_dp, 2.19e-3_dp, 3.91e-2_dp, 7.11e-3_dp, -1.44e-3_dp, -3.85e-4_dp, &
           1.39e-1_dp, 6.06e-2_dp, -3.06e-3_dp, -4.55e-3_dp, -1.33e-1_dp, -3.94e-2_dp, &
             3.99e-3_dp, 2.36e-3_dp, 3.35e-2_dp, 6.90e-3_dp, -1.24e-3_dp, -4.05e-4_dp, &
           1.11e-1_dp, 5.41e-2_dp, -3.46e-4_dp, -3.68e-3_dp, -1.08e-1_dp, -3.47e-2_dp, &
             1.93e-3_dp, 1.81e-3_dp, 2.80e-2_dp, 6.04e-3_dp, -8.41e-4_dp, -3.17e-4_dp, &
           -1.10e-2_dp, 2.02e-2_dp, 9.46e-4_dp, -4.09e-3_dp, -6.66e-3_dp, -1.17e-2_dp, &
           -2.36e-4_dp, 2.02e-3_dp, 4.08e-3_dp, 1.85e-3_dp, -6.48e-5_dp, -3.13e-4_dp], &
                                                      [12, n_imp])
    real(dp), parameter :: b_high(12, n_imp) = reshape([ &
            2.31e-1_dp, 3.43e-1_dp, -1.85e-1_dp, -1.62e-2_dp, 1.05e-1_dp, -7.03e-2_dp, &
            5.31e-2_dp, 3.42e-3_dp, -8.38e-3_dp, 4.15e-3_dp, -3.35e-3_dp, -2.21e-4_dp, &
            -7.32e-1_dp, 1.83e-2_dp, -1.55e-1_dp, -1.72e-2_dp, 3.21e-1_dp, 9.46e-3_dp, &
           3.97e-2_dp, 4.20e-3_dp, -2.04e-2_dp, -6.19e-4_dp, -2.24e-3_dp, -2.54e-4_dp, &
            -1.01e0_dp, -8.65e-3_dp, -1.24e-1_dp, -1.45e-2_dp, 3.91e-1_dp, 1.61e-2_dp, &
           2.98e-2_dp, 3.32e-3_dp, -2.48e-2_dp, -1.04e-3_dp, -1.52e-3_dp, -1.89e-4_dp, &
            -1.02e0_dp, -1.39e-2_dp, -9.79e-2_dp, -1.17e-2_dp, 3.75e-1_dp, 1.56e-2_dp, &
           2.24e-2_dp, 2.54e-3_dp, -2.26e-2_dp, -8.89e-4_dp, -1.04e-3_dp, -1.39e-4_dp, &
            -1.02e0_dp, -1.48e-2_dp, -6.74e-2_dp, -9.17e-3_dp, 3.59e-1_dp, 1.43e-2_dp, &
           1.39e-2_dp, 1.84e-3_dp, -2.09e-2_dp, -7.32e-4_dp, -5.02e-4_dp, -9.49e-5_dp, &
             -8.20e-1_dp, -6.36e-3_dp, 5.42e-2_dp, 3.95e-3_dp, 2.02e-1_dp, 8.06e-4_dp, &
           -2.00e-2_dp, -1.78e-3_dp, -6.10e-3_dp, 6.51e-4_dp, 1.75e-3_dp, 1.46e-4_dp], &
                                                       [12, n_imp])

contains

    pure logical function suzuki_impurity_supported(z) result(ok)
        integer, intent(in) :: z
        ok = any(z_imp == z)
    end function suzuki_impurity_supported

    function suzuki_stopping(e_kev_amu, ne_m3, te_ev, beam_mass, zeff, z_impurity) &
        result(sigma)
        !> Beam stopping cross section [m^2] for a hydrogen beam of mass number
        !> beam_mass (1, 2 or 3) at energy e_kev_amu [keV/amu] in a hydrogen
        !> plasma with one impurity of charge z_impurity (0 for none).
        !> Electron density and temperature are clamped to the fit ranges of the
        !> paper, 1e18 <= ne <= 1e21 m^-3 and E/100 <= Te <= E/2 (low energies)
        !> or 1 <= Te <= 50 keV (high energies), instead of extrapolating.
        real(dp), intent(in) :: e_kev_amu, ne_m3, te_ev, zeff
        integer, intent(in) :: beam_mass, z_impurity
        real(dp) :: sigma
        real(dp) :: eps, ln_n, u, n19, te_kev, s_z, a(10), b(12)
        integer :: k
        logical :: low

        if (.not. all(is_finite([e_kev_amu, ne_m3, te_ev, zeff]))) &
            error stop 'beam stopping: nonfinite input'
        if (beam_mass < 1 .or. beam_mass > 3) &
            error stop 'beam stopping: beam mass number must be 1, 2 or 3'
        ! The low-energy fit holds from 10 keV/amu. It is extrapolated down to
        ! 9 keV/amu, as in ASCOT5, for the E/3 component of 55-60 keV
        ! deuterium beams.
        if (e_kev_amu < e_min_kev_amu .or. e_kev_amu > 1.0e4_dp) &
            error stop 'beam stopping: energy outside [9,10000] keV/amu'
        if (ne_m3 <= 0.0_dp .or. te_ev <= 0.0_dp) &
            error stop 'beam stopping: nonpositive density or temperature'
        if (zeff < 1.0_dp) error stop 'beam stopping: Zeff below one'

        low = e_kev_amu < 100.0_dp
        te_kev = te_ev*1.0e-3_dp
        if (low) then
            te_kev = min(max(te_kev, e_kev_amu/100.0_dp), e_kev_amu/2.0_dp)
            a = a_low(:, beam_mass)
        else
            te_kev = min(max(te_kev, 1.0_dp), 50.0_dp)
            a = a_high(:, beam_mass)
        end if
        n19 = min(max(ne_m3, 1.0e18_dp), 1.0e21_dp)*1.0e-19_dp
        eps = log(e_kev_amu)
        ln_n = log(n19)
        u = log(te_kev)

        ! Eq. (28) in cm^2.
        sigma = a(1)*1.0e-16_dp/e_kev_amu*(1.0_dp + a(2)*eps + a(3)*eps**2) &
           *(1.0_dp + (1.0_dp - exp(-a(4)*n19))**a(5)*(a(6) + a(7)*eps + a(8)*eps**2)) &
                *(1.0_dp + a(9)*u + a(10)*u**2)

        s_z = 0.0_dp
        if (z_impurity /= 0 .and. zeff > 1.0_dp) then
            k = findloc(z_imp, z_impurity, dim=1)
            if (k == 0) error stop 'beam stopping: impurity charge without Suzuki table'
            if (low) then
                b = b_low(:, k)
            else
                b = b_high(:, k)
            end if
            ! Eq. (26) with i, j, k the powers of eps, ln N and U.
            s_z = b(1) + b(2)*u + b(3)*ln_n + b(4)*ln_n*u &
                  + eps*(b(5) + b(6)*u + b(7)*ln_n + b(8)*ln_n*u) &
                  + eps**2*(b(9) + b(10)*u + b(11)*ln_n + b(12)*ln_n*u)
        end if
        ! Eq. (24), cm^2 to m^2.
        sigma = sigma*(1.0_dp + (zeff - 1.0_dp)*s_z)*1.0e-4_dp
        if (.not. is_finite(sigma) .or. sigma <= 0.0_dp) &
            error stop 'beam stopping: nonpositive cross section'
    end function suzuki_stopping

end module beam_stopping
