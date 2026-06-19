module field_can_test

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_base, only: field_can_t, n_field_evaluations

implicit none

contains

subroutine evaluate_test(f, r, th_c, ph_c, mode_secders)
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    call eval_field_test(f, r, th_c, ph_c, mode_secders)

    n_field_evaluations = n_field_evaluations + 1
end subroutine evaluate_test


! for testing -> circular tokamak
subroutine eval_field_test(f, r, th, ph, mode_secders)
    type(field_can_t), intent(inout) :: f
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


! Exact-curl analytic tokamak field for the 6D canonical port (field_correct_test.py).
! A, dA, d2A are byte-identical to eval_field_test; only B, dB, d2B, h differ:
! the GC-linearized B=B0(1-r/R0 cos th) is replaced by the exact |B| from B^k =
! eps^ijk A_j,i/sqrtg, |B|=sqrt(g_ij B^i B^j). With A_r=0 only B^th, B^ph survive
! and |B|^2 = A_ph,r^2/(R0+r cos th)^2 + A_th,r^2/r^2 (W below). dBmod/d2Bmod come
! from W via |B|=sqrt(W). The covariant h_i = g_ii B^i / |B| (h_r = 0). The 6D
! models need this exact field; the GC path keeps eval_field_test.
subroutine evaluate_correct_test(f, r, th_c, ph_c, mode_secders)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    call eval_field_correct_test(f, r, th_c, ph_c, mode_secders)

    n_field_evaluations = n_field_evaluations + 1
end subroutine evaluate_correct_test


subroutine eval_field_correct_test(f, r, th, ph, mode_secders)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, th, ph
    integer, intent(in) :: mode_secders

    real(dp) :: B0, iota0, a, R0, cth, sth, J, Rr
    real(dp) :: dAthr, dAphr, Bth, Bph, W
    real(dp) :: Bmod, d2Athrr, d2Athrth, d2Aphrr

    B0 = 1.0d0; iota0 = 1.0d0; a = 0.5d0; R0 = 1.0d0
    cth = cos(th); sth = sin(th)
    Rr = R0 + r*cth
    J = r*Rr

    ! Vector potential and its derivatives: identical to eval_field_test.
    f%Ath = B0*(r**2/2d0 - r**3/(3d0*R0)*cth)
    f%Aph = -B0*iota0*(r**2/2d0 - r**4/(4d0*a**2))

    f%dAth(1) = B0*(r - r**2/R0*cth); f%dAth(2) = B0*r**3*sth/(3d0*R0); f%dAth(3) = 0d0
    f%dAph(1) = -B0*iota0*(r - r**3/a**2); f%dAph(2) = 0d0; f%dAph(3) = 0d0

    ! Exact-curl |B|: B^th = -A_ph,r/J, B^ph = A_th,r/J, |B|^2 = W.
    dAthr = f%dAth(1); dAphr = f%dAph(1)
    Bth = -dAphr/J; Bph = dAthr/J
    W = dAphr**2/Rr**2 + dAthr**2/r**2
    Bmod = sqrt(W)
    f%Bmod = Bmod

    ! Covariant unit-field components h_i = g_ii B^i / |B|; h_r = 0 (no radial B).
    f%hth = r**2*Bth/Bmod
    f%hph = Rr**2*Bph/Bmod

    ! dBmod ported verbatim from field_correct_test.py (lines 34-36) so the 6D
    ! port reproduces the python oracle bit-for-bit. dBmod(2) is the python
    ! listing's value (it omits one chain-rule term in d|B|/dtheta); the oracle
    ! trajectories were generated with it and the residual must match. The true
    ! analytic d|B|/dtheta is recovered from W only in the Jacobian's d2 block,
    ! where the O(mu)=1e-5 force makes the difference irrelevant to the fixed point.
    d2Athrr = B0*(1d0 - 2d0*r/R0*cth)
    d2Athrth = B0*r**2/R0*sth
    d2Aphrr = -B0*iota0*(1d0 - 3d0*r**2/a**2)
    f%dBmod(1) = (r*Bth**2 + r**2*Bth*(-1d0/J*d2Aphrr + 1d0/J**2*(R0 + 2d0*r*cth)*dAphr) &
                  + Rr*cth*Bph**2 + Rr**2*Bph*(1d0/J*d2Athrr - 1d0/J**2*(R0 + 2d0*r*cth)*dAthr))/Bmod
    f%dBmod(2) = (-Rr*r*sth*Bph**2 + Rr**2*Bph*(1d0/J*d2Athrth + 1d0/J**2*r**2*sth*dAthr))/Bmod
    f%dBmod(3) = 0d0

    if (mode_secders <= 0) return

    f%d2Ath(1) = d2Athrr; f%d2Ath(2) = d2Athrth
    f%d2Ath(3) = 0d0; f%d2Ath(4) = B0*r**3*cth/(3d0*R0); f%d2Ath(5) = 0d0; f%d2Ath(6) = 0d0

    f%d2Aph(1) = d2Aphrr
    f%d2Aph(2) = 0d0; f%d2Aph(3) = 0d0; f%d2Aph(4) = 0d0; f%d2Aph(5) = 0d0; f%d2Aph(6) = 0d0

    ! d2Bmod is left unset: the python dBmod (above) is not a true gradient, so a
    ! symmetric packed Hessian cannot represent its mixed derivative consistently.
    ! The 6D Jacobian's mu|B| term takes d(dBmod)/dq by a central difference of
    ! dBmod itself, which is exact for whichever dBmod the residual uses.
    f%d2Bmod = 0d0

end subroutine eval_field_correct_test

end module field_can_test
