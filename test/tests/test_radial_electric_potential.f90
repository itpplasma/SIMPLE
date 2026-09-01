program test_radial_electric_potential
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: field_can_t, apply_radial_electric_potential, &
        set_radial_electric_potential_slope
    use field_can_mod, only: get_derivatives, get_derivatives2
    implicit none

    type(field_can_t) :: field
    real(dp), parameter :: slope = 3.25_dp, s = 0.4_dp

    field%Aph = 0.0_dp; field%Ath = 0.0_dp
    field%hph = 1.0_dp; field%hth = 0.0_dp
    field%Bmod = 1.0_dp; field%mu = 0.0_dp; field%ro0 = 1.0_dp
    field%dAph = 0.0_dp; field%dAth = 0.0_dp
    field%dhph = 0.0_dp; field%dhth = 0.0_dp; field%dBmod = 0.0_dp
    field%d2Aph = 0.0_dp; field%d2Ath = 0.0_dp
    field%d2hph = 0.0_dp; field%d2hth = 0.0_dp; field%d2Bmod = 0.0_dp

    call set_radial_electric_potential_slope(0.0_dp)
    call apply_radial_electric_potential(field, s)
    call get_derivatives2(field, 0.0_dp)
    if (field%H /= 0.0_dp .or. maxval(abs(field%dH)) /= 0.0_dp) &
        error stop 'zero potential changed the magnetic Hamiltonian'

    call set_radial_electric_potential_slope(slope)
    call apply_radial_electric_potential(field, s)
    call get_derivatives(field, 0.0_dp)
    if (abs(field%H - slope*s) > 1.0e-14_dp) error stop 'potential energy mismatch'
    if (abs(field%dH(1) - slope) > 1.0e-14_dp) error stop 'potential force mismatch'
    if (maxval(abs(field%dH(2:3))) > 1.0e-14_dp) error stop 'spurious angular force'

    call get_derivatives2(field, 0.0_dp)
    if (maxval(abs(field%d2Phie)) > 1.0e-14_dp) error stop 'linear potential curvature'
end program test_radial_electric_potential
