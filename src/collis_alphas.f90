module collis_alp
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none

    integer, parameter :: nsorts = 3
    integer, parameter :: N_S_GRID = 101

    real(wp), dimension(nsorts) :: efcolf, velrat, enrat
    real(wp), dimension(nsorts, N_S_GRID) :: efcolf_grid, velrat_grid, enrat_grid

contains

    subroutine coleff(p, dpp, dhh, fpeff)
        real(wp), intent(in) :: p
        real(wp), intent(out) :: dpp, dhh, fpeff

        call coleff_local(p, efcolf, velrat, enrat, dpp, dhh, fpeff)
    end subroutine coleff

    subroutine coleff_local(p, efcolf_loc, velrat_loc, enrat_loc, dpp, dhh, fpeff)
        real(wp), intent(in) :: p
        real(wp), intent(in) :: efcolf_loc(nsorts), velrat_loc(nsorts), enrat_loc(nsorts)
        real(wp), intent(out) :: dpp, dhh, fpeff

        integer :: i
        real(wp) :: plim, xbeta, dp, dh, dpd

        plim = max(p, 1.d-8)

        dpp = 0.0d0
        dhh = 0.0d0
        fpeff = 0.0d0

        do i = 1, nsorts
            xbeta = p*velrat_loc(i)
            call onseff(xbeta, dp, dh, dpd)
            dpp = dpp + dp*efcolf_loc(i)
            dhh = dhh + dh*efcolf_loc(i)
            fpeff = fpeff + (dpd/plim - 2.0d0*dp*p*enrat_loc(i))*efcolf_loc(i)
        end do

        dhh = dhh/plim**2
    end subroutine coleff_local

    subroutine onseff(v, dp, dh, dpd)
        real(wp), intent(in) :: v
        real(wp), intent(out) :: dp, dh, dpd

        real(wp), parameter :: sqp = 1.7724538d0
        real(wp), parameter :: cons = 0.75225278d0
        real(wp) :: v2, v3, ex, er

        v2 = v**2
        v3 = v2*v
        if (v < 0.01d0) then
            dp = cons*(1.d0 - 0.6d0*v2)
            dh = cons*(1.d0 - 0.2d0*v2)
            dpd = 2.d0*cons*(1.d0 - 1.2d0*v2)
        else if (v > 6.d0) then
            dp = 1.d0/v3
            dh = (1.d0 - 0.5d0/v2)/v
            dpd = -1.d0/v3
        else
            ex = exp(-v2)/sqp
            er = erf(v)
            dp = er/v3 - 2.d0*ex/v2
            dh = er*(1.d0 - 0.5d0/v2)/v + ex/v2
            dpd = 4.d0*ex - dp
        end if
    end subroutine onseff

    subroutine loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, &
                            tempe, ealpha, v0, dchichi, slowrate, &
                            dchichi_norm, slowrate_norm)
        real(wp), intent(in) :: am1, am2, Z1, Z2, densi1, densi2
        real(wp), intent(in) :: tempi1, tempi2, tempe, ealpha
        real(wp), intent(out) :: v0, dchichi, slowrate, dchichi_norm, slowrate_norm

        call compute_collision_coeffs(am1, am2, Z1, Z2, densi1, densi2, &
                                      tempi1, tempi2, tempe, ealpha, v0, &
                                      efcolf, velrat, enrat)

        call compute_rates(v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
    end subroutine loacol_alpha

    subroutine compute_collision_coeffs(am1, am2, Z1, Z2, densi1, densi2, &
                                        tempi1, tempi2, tempe, ealpha, v0, &
                                        efcolf_out, velrat_out, enrat_out)
        real(wp), intent(in) :: am1, am2, Z1, Z2, densi1, densi2
        real(wp), intent(in) :: tempi1, tempi2, tempe, ealpha
        real(wp), intent(out) :: v0
        real(wp), intent(out) :: efcolf_out(nsorts), velrat_out(nsorts), enrat_out(nsorts)

        real(wp), parameter :: pi = 3.14159265358979d0
        real(wp), parameter :: pmass = 1.6726d-24
        real(wp), parameter :: emass = 9.1094d-28
        real(wp), parameter :: e_cgs = 4.8032d-10
        real(wp), parameter :: ev_cgs = 1.6022d-12

        real(wp) :: dense, vti1, vti2, vte
        real(wp) :: alami1, alami2, alame, frecol_base

        enrat_out(1) = ealpha/tempi1
        enrat_out(2) = ealpha/tempi2
        enrat_out(3) = ealpha/tempe

        v0 = sqrt(2.d0*ealpha*ev_cgs/(4.d0*pmass))
        vti1 = sqrt(2.d0*tempi1*ev_cgs/(pmass*am1))
        vti2 = sqrt(2.d0*tempi2*ev_cgs/(pmass*am2))
        vte = sqrt(2.d0*tempe*ev_cgs/emass)

        velrat_out(1) = v0/vti1
        velrat_out(2) = v0/vti2
        velrat_out(3) = v0/vte

        dense = densi1*Z1 + densi2*Z2
        alami1 = 23.d0 - log(max(epsilon(1.d0), &
                                 sqrt(densi1*Z1**2/tempi1)*2.d0*Z1*(4.d0 + am1)/(4.d0*tempi1 + am1*ealpha)))
        alami2 = 23.d0 - log(max(epsilon(1.d0), &
                                 sqrt(densi2*Z2**2/tempi2)*2.d0*Z2*(4.d0 + am2)/(4.d0*tempi2 + am2*ealpha)))
        alame = 24.d0 - log(sqrt(dense)/tempe)

        frecol_base = 2.d0*pi*dense*e_cgs**4*2.d0**2/((4.d0*pmass)**2*v0**3)
        frecol_base = frecol_base/v0

        efcolf_out(1) = frecol_base*Z1**2*alami1*densi1/dense
        efcolf_out(2) = frecol_base*Z2**2*alami2*densi2/dense
        efcolf_out(3) = frecol_base*alame

        efcolf_out = efcolf_out*velrat_out
    end subroutine compute_collision_coeffs

    subroutine compute_rates(v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        real(wp), intent(in) :: v0
        real(wp), intent(out) :: dchichi, slowrate, dchichi_norm, slowrate_norm

        real(wp) :: p, dpp, dhh, fpeff

        p = 1.d0
        call coleff(p, dpp, dhh, fpeff)

        dchichi = dhh*v0
        slowrate = abs(fpeff)*v0
        dchichi_norm = dhh
        slowrate_norm = abs(fpeff)
    end subroutine compute_rates

    subroutine init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)
        use simple_profiles, only: get_plasma_params

        real(wp), intent(in) :: am1, am2, Z1, Z2, ealpha, v0
        real(wp) :: s, Te, Ti1, Ti2, ni1, ni2
        real(wp) :: densi1_loc, densi2_loc, tempi1_loc, tempi2_loc, tempe_loc
        real(wp) :: v0_loc
        integer :: i

        do i = 1, N_S_GRID
            s = dble(i - 1)/dble(N_S_GRID - 1)

            call get_plasma_params(s, Te, Ti1, Ti2, ni1, ni2)

            tempe_loc = max(Te, 1.0d0)
            tempi1_loc = max(Ti1, 1.0d0)
            tempi2_loc = max(Ti2, 1.0d0)
            densi1_loc = max(ni1*1.0d-6, 1.0d0)
            densi2_loc = max(ni2*1.0d-6, 1.0d0)

            call compute_collision_coeffs(am1, am2, Z1, Z2, &
                                          densi1_loc, densi2_loc, &
                                          tempi1_loc, tempi2_loc, tempe_loc, &
                                          ealpha, v0_loc, &
                                          efcolf_grid(:, i), velrat_grid(:, i), enrat_grid(:, i))
        end do
    end subroutine init_collision_profiles

    subroutine get_local_coeffs(s, efcolf_loc, velrat_loc, enrat_loc)
        real(wp), intent(in) :: s
        real(wp), intent(out) :: efcolf_loc(nsorts), velrat_loc(nsorts), enrat_loc(nsorts)

        real(wp) :: s_clamped, s_norm, weight
        integer :: idx_lo, idx_hi

        s_clamped = max(0.0d0, min(1.0d0, s))
        s_norm = s_clamped*dble(N_S_GRID - 1)
        idx_lo = max(1, min(N_S_GRID - 1, int(s_norm) + 1))
        idx_hi = idx_lo + 1
        weight = s_norm - dble(idx_lo - 1)

        efcolf_loc = (1.0d0 - weight)*efcolf_grid(:, idx_lo) + weight*efcolf_grid(:, idx_hi)
        velrat_loc = (1.0d0 - weight)*velrat_grid(:, idx_lo) + weight*velrat_grid(:, idx_hi)
        enrat_loc = (1.0d0 - weight)*enrat_grid(:, idx_lo) + weight*enrat_grid(:, idx_hi)
    end subroutine get_local_coeffs

    subroutine stost(z, dtauc, iswmode, ierr)
        integer, intent(in) :: iswmode
        real(wp), intent(in) :: dtauc
        real(wp), intent(inout) :: z(5)
        integer, intent(out) :: ierr

        real(wp), parameter :: pmin = 1.d-8
        real(wp) :: p, dpp, dhh, fpeff, alam, dalam, coala
        real(wp) :: efcolf_loc(nsorts), velrat_loc(nsorts), enrat_loc(nsorts)
        real :: ur

        p = z(4)

        call get_local_coeffs(z(1), efcolf_loc, velrat_loc, enrat_loc)
        call coleff_local(p, efcolf_loc, velrat_loc, enrat_loc, dpp, dhh, fpeff)

        ierr = 0

        if (iswmode == 1 .or. iswmode == 4) then
            alam = z(5)
            coala = 1.d0 - alam**2

            if (coala < 0.d0) then
                ierr = 1
                return
            end if

            call getran(1, ur)

            dalam = sqrt(2.d0*dhh*coala*dtauc)*dble(ur) - 2.d0*alam*dhh*dtauc

            if (abs(dalam) > 1.d0) then
                ierr = 2
                call random_number(ur)
                alam = 2.d0*(dble(ur) - 0.5d0)
            else
                alam = alam + dalam
                if (alam > 1.d0) then
                    ierr = 3
                    alam = 2.d0 - alam
                else if (alam < -1.d0) then
                    ierr = 3
                    alam = -2.d0 - alam
                end if
            end if

            z(5) = alam
            if (iswmode == 4) return
        end if

        if (iswmode < 3) then
            call getran(0, ur)
            z(4) = z(4) + sqrt(abs(2.d0*dpp*dtauc))*dble(ur) + fpeff*dtauc
        else
            z(4) = z(4) + fpeff*dtauc
        end if

        if (z(4) < pmin) then
            ierr = ierr + 10
            z(4) = pmin + abs(pmin - z(4))
        end if
    end subroutine stost

    subroutine getran(irand, ur)
        integer, intent(in) :: irand
        real, intent(out) :: ur

        call random_number(ur)

        if (irand == 0) then
            ur = 3.464102*(ur - 0.5)
        else
            if (ur > 0.5) then
                ur = 1.0
            else
                ur = -1.0
            end if
        end if
    end subroutine getran

end module collis_alp
