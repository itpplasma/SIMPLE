program test_dd_beam_rate
    !> D(d,n)3He thermal reactivity against Bosch and Hale (1992) Table VIII and
    !> the beam-Maxwellian rate against its limits and an independent quadrature.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fusion_source, only: dd_neutron_reactivity, dd_beam_maxwellian_rate, &
                             dd_neutron_cross_section
    implicit none
    real(dp), parameter :: ev = 1.602176634e-19_dp, m = 2.0_dp*1.67262192369e-27_dp
    real(dp), parameter :: table_t(4) = [1.0_dp, 5.0_dp, 10.0_dp, 50.0_dp]
    real(dp), parameter :: table_rate(4) = [9.933e-23_dp, 9.128e-20_dp, 6.023e-19_dp, &
                                            1.133e-17_dp]
    real(dp) :: v, rate, w
    integer :: i

    ! Table VIII, four significant digits, cm^3/s.
    do i = 1, 4
        rate = dd_neutron_reactivity(1.0e3_dp*table_t(i))*1.0e6_dp
        if (abs(rate/table_rate(i) - 1.0_dp) > 1.0e-3_dp) then
            print *, table_t(i), rate, table_rate(i)
            error stop 'D(d,n) reactivity differs from Bosch and Hale Table VIII'
        end if
    end do

    ! A deuteron at rest in deuterons at 2 T sees the relative-speed
    ! distribution of a thermal D-D plasma at T. The cross-section fit
    ! (Table IV) deviates by up to 2.5 % from the data behind Table VII.
    do i = 1, 3
        rate = dd_beam_maxwellian_rate(0.0_dp, 2.0e3_dp*table_t(i), m)
     if (abs(rate/dd_neutron_reactivity(1.0e3_dp*table_t(i)) - 1.0_dp) > 2.5e-2_dp) then
            print *, table_t(i), rate, dd_neutron_reactivity(1.0e3_dp*table_t(i))
            error stop 'rest-frame rate differs from the thermal reactivity at T/2'
        end if
    end do

    ! Cold target: sigma(E/2) v.
    v = sqrt(2.0_dp*60.0e3_dp*ev/m)
   if (abs(dd_beam_maxwellian_rate(v, 0.0_dp, m)/(v*dd_neutron_cross_section(30.0_dp)) &
            - 1.0_dp) > 1.0e-12_dp) &
        error stop 'cold-target rate is not sigma(E/2) v'

    ! Gauss-Hermite branch (v >= 2 w) against the Simpson branch evaluated
    ! for the same integral just below the switch, and against a direct sum.
    w = sqrt(2.0_dp*5.0e3_dp*ev/m)
    rate = dd_beam_maxwellian_rate(2.0_dp*w, 5.0e3_dp, m)
    if (abs(dd_beam_maxwellian_rate(2.0_dp*w*(1.0_dp - 1.0e-9_dp), 5.0e3_dp, m)/rate &
            - 1.0_dp) > 1.0e-4_dp) error stop 'quadrature branches disagree at v = 2 w'
    rate = dd_beam_maxwellian_rate(v, 5.0e3_dp, m)
    if (abs(rate/direct(v, w) - 1.0_dp) > 1.0e-4_dp) then
        print *, rate, direct(v, w)
        error stop 'beam rate differs from direct quadrature'
    end if
    ! A warm target raises the 60 keV rate (sigma grows with E below 1 MeV).
    if (rate <= v*dd_neutron_cross_section(30.0_dp)) error stop 'warm target lowers the rate'
    print *, 'dd beam rate OK, 60 keV on 5 keV / cold:', rate/(v*dd_neutron_cross_section(30.0_dp))

contains

    function direct(v, w) result(r)
        ! Midpoint sum of the speed integral with both exponentials.
        real(dp), intent(in) :: v, w
        real(dp) :: r, u, du
        integer :: k
        integer, parameter :: n = 200000
        du = (v + 12.0_dp*w)/n
        r = 0.0_dp
        do k = 1, n
            u = (k - 0.5_dp)*du
            r = r + u**2*(exp(-((u - v)/w)**2) - exp(-((u + v)/w)**2))* &
                dd_neutron_cross_section(0.25_dp*m*u**2/ev*1.0e-3_dp)
        end do
        r = r*du/(sqrt(acos(-1.0_dp))*v*w)
    end function direct

end program test_dd_beam_rate
