module field_can_test

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_base, only: FieldCan, n_field_evaluations

implicit none

contains

subroutine evaluate(f, r, th_c, ph_c, mode_secders)
    type(FieldCan), intent(inout) :: f

    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    call eval_field_test(f, r, th_c, ph_c, mode_secders)

    n_field_evaluations = n_field_evaluations + 1
end subroutine evaluate


! for testing -> circular tokamak
subroutine eval_field_test(f, r, th, ph, mode_secders)
    type(FieldCan), intent(inout) :: f
    double precision, intent(in) :: r, th, ph
    integer, intent(in) :: mode_secders

    double precision :: B0, iota0, a, R0, cth, sth

    B0 = 1.0    ! magnetic field modulus normalization
    iota0 = 1.0 ! constant part of rotational transform
    a = 0.5     ! (equivalent) minor radius
    R0 = 1.0    ! (equivalent) major radius

    cth = cos(th)
    sth = sin(th)

    f%Ath = B0*(r**2/2d0 - r**3/(3d0*R0)*cth)
    f%Aph     = -B0*iota0*(r**2/2d0-r**4/(4d0*a**2))
    f%hth     = iota0*(1d0-r**2/a**2)*r**2/R0
    f%hph     = R0 + r*cth
    f%Bmod     = B0*(1d0 - r/R0*cth)

    f%dAth(1) = B0*(r - r**2/R0*cth)
    f%dAth(2) = B0*r**3*sth/(3.0*R0)
    f%dAth(3) = 0d0

    f%dAph(1) = -B0*iota0*(r-r**3/a**2)
    f%dAph(2) = 0d0
    f%dAph(3) = 0d0

    f%dhth(1) = 2d0*iota0*r*(a**2-2d0*r**2)/(a**2*R0)
    f%dhth(2) = 0d0
    f%dhth(3) = 0d0

    f%dhph(1) = cth
    f%dhph(2) = -r*sth
    f%dhph(3) = 0d0

    f%dBmod(1) = -B0/R0*cth
    f%dBmod(2) = B0*r/R0*sth
    f%dBmod(3) = 0d0

    if (mode_secders <= 0) return

    f%d2Ath(1) = B0*(1d0 - 2d0*r/R0*cth)
    f%d2Ath(4) = B0*r**3*cth/(3d0*R0)
    f%d2Ath(3) = 0d0
    f%d2Ath(2) = B0*r**2/R0*sth
    f%d2Ath(5) = 0d0
    f%d2Ath(6) = 0d0

    f%d2Aph(1) = -B0*iota0*(1d0-3d0*r**2/a**2)
    f%d2Aph(4) = 0d0
    f%d2Aph(3) = 0d0
    f%d2Aph(2) = 0d0
    f%d2Aph(5) = 0d0
    f%d2Aph(6) = 0d0

    f%d2hth = 0d0
    f%d2hth(1) = 2d0*iota0*(a**2-6d0*r**2)/(a**2*R0)

    f%d2hph(1) = 0d0
    f%d2hph(4) = -r*cth
    f%d2hph(3) = 0d0
    f%d2hph(2) = -sth
    f%d2hph(5) = 0d0
    f%d2hph(6) = 0d0

    f%d2Bmod(1) = 0d0
    f%d2Bmod(4) = B0*r/R0*cth
    f%d2Bmod(3) = 0d0
    f%d2Bmod(2) = B0/R0*sth
    f%d2Bmod(5) = 0d0
    f%d2Bmod(6) = 0d0

end subroutine eval_field_test

end module field_can_test
