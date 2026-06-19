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
! and |B|^2 = A_ph,r^2/(R0+r cos th)^2 + A_th,r^2/r^2 (W below). dBmod and the
! Hessian d2Bmod are the exact closed-form derivatives of |B|=sqrt(W):
! d_k|B| = dW_k/(2|B|), d2_kl|B| = d2W_kl/(2|B|) - dW_k dW_l/(4|B|^3), FD-verified
! to ~1e-9 against the field. The covariant h_i = g_ii B^i / |B| (h_r = 0). The 6D
! models need this exact field; the GC path keeps eval_field_test.
subroutine evaluate_correct_test(f, r, th_c, ph_c, mode_secders)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    call eval_field_correct_test(f, r, th_c, ph_c, mode_secders)

    n_field_evaluations = n_field_evaluations + 1
end subroutine evaluate_correct_test


subroutine eval_field_correct_test(f, r, th, ph, mode_secders)
    !$acc routine seq
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, th, ph
    integer, intent(in) :: mode_secders

    real(dp) :: B0, iota0, a, R0, cth, sth, J, Rr
    real(dp) :: dAthr, dAphr, Bth, Bph, W, Bmod
    real(dp) :: d2Athrr, d2Athrth, d2Aphrr, d3Athrrr, d3Athrrth, d3Athrthth, d3Aphrrr
    real(dp) :: dWdr, dWdth, d2Wrr, d2Wrth, d2Wthth

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

    ! A_i derivatives entering W and its derivatives (A_r=0; A_th,A_ph are phi-free).
    d2Athrr = B0*(R0 - 2d0*r*cth)/R0       ! A_th,rr
    d2Athrth = B0*r**2*sth/R0              ! A_th,rth
    d2Aphrr = -B0*iota0*(a**2 - 3d0*r**2)/a**2  ! A_ph,rr

    ! Exact closed-form d|B| = dW/(2|B|). dW keeps both the Rr-dependence and
    ! A_th,r's theta-dependence (A_th,rth) -- the latter is the chain-rule term the
    ! python listing dropped in d|B|/dtheta. FD-verified to ~1e-9 against |B|.
    dWdr = 2d0*dAphr*d2Aphrr/Rr**2 - 2d0*dAphr**2*cth/Rr**3 &
           + 2d0*dAthr*d2Athrr/r**2 - 2d0*dAthr**2/r**3
    dWdth = 2d0*r*sth*dAphr**2/Rr**3 + 2d0*dAthr*d2Athrth/r**2
    f%dBmod(1) = dWdr/(2d0*Bmod)
    f%dBmod(2) = dWdth/(2d0*Bmod)
    f%dBmod(3) = 0d0

    if (mode_secders <= 0) return

    f%d2Ath(1) = d2Athrr; f%d2Ath(2) = d2Athrth
    f%d2Ath(3) = 0d0; f%d2Ath(4) = B0*r**3*cth/(3d0*R0); f%d2Ath(5) = 0d0; f%d2Ath(6) = 0d0

    f%d2Aph(1) = d2Aphrr
    f%d2Aph(2) = 0d0; f%d2Aph(3) = 0d0; f%d2Aph(4) = 0d0; f%d2Aph(5) = 0d0; f%d2Aph(6) = 0d0

    ! True analytic Hessian of the corrected |B|: d2|B|_kl = d2W_kl/(2|B|) -
    ! dW_k dW_l/(4|B|^3). Third A_i derivatives entering d2W (only rr,rth,thth couple):
    d3Athrrr = -2d0*B0*cth/R0; d3Athrrth = 2d0*B0*r*sth/R0; d3Athrthth = B0*r**2*cth/R0
    d3Aphrrr = 6d0*B0*iota0*r/a**2
    d2Wrr = 2d0*d2Aphrr**2/Rr**2 + 2d0*dAphr*d3Aphrrr/Rr**2 - 8d0*dAphr*d2Aphrr*cth/Rr**3 &
            + 6d0*dAphr**2*cth**2/Rr**4 &
            + 2d0*d2Athrr**2/r**2 + 2d0*dAthr*d3Athrrr/r**2 - 8d0*dAthr*d2Athrr/r**3 &
            + 6d0*dAthr**2/r**4
    d2Wrth = 4d0*dAphr*d2Aphrr*r*sth/Rr**3 - 2d0*dAphr**2*(-sth/Rr**3 + 3d0*r*sth*cth/Rr**4) &
             + 2d0*(d2Athrth*d2Athrr + dAthr*d3Athrrth)/r**2 - 4d0*dAthr*d2Athrth/r**3
    d2Wthth = 2d0*r*dAphr**2*(cth/Rr**3 + 3d0*r*sth**2/Rr**4) &
              + 2d0*(d2Athrth**2 + dAthr*d3Athrthth)/r**2
    f%d2Bmod = 0d0
    f%d2Bmod(1) = d2Wrr/(2d0*Bmod)   - dWdr*dWdr/(4d0*Bmod**3)     ! d2|B|/drr
    f%d2Bmod(2) = d2Wrth/(2d0*Bmod)  - dWdr*dWdth/(4d0*Bmod**3)    ! d2|B|/drth
    f%d2Bmod(4) = d2Wthth/(2d0*Bmod) - dWdth*dWdth/(4d0*Bmod**3)   ! d2|B|/dthth

end subroutine eval_field_correct_test

end module field_can_test
