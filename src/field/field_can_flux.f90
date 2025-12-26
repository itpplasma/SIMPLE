module field_can_flux

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_base, only: field_can_t, n_field_evaluations, twopi

implicit none
private

public :: evaluate_flux, integ_to_ref_flux, ref_to_integ_flux, eval_field_can
public :: twopi

contains

subroutine evaluate_flux(f, r, th_c, ph_c, mode_secders)
    ! Note: !$acc routine seq removed - n_field_evaluations is threadprivate
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    call eval_field_can(f, r, th_c, ph_c, mode_secders)

    ! Note: !$acc atomic removed - n_field_evaluations is threadprivate
    n_field_evaluations = n_field_evaluations + 1
end subroutine evaluate_flux


subroutine integ_to_ref_flux(xinteg, xref)
    use get_can_sub, only: can_to_vmec
    real(dp), intent(in) :: xinteg(3)
    real(dp), intent(out) :: xref(3)

    xref(1) = xinteg(1)
    call can_to_vmec(xinteg(1), xinteg(2), xinteg(3), xref(2), xref(3))
    xref(2) = mod(xref(2), twopi)
    xref(3) = mod(xref(3), twopi)
end subroutine integ_to_ref_flux


subroutine ref_to_integ_flux(xref, xinteg)
    use get_can_sub, only: vmec_to_can
    real(dp), intent(in) :: xref(3)
    real(dp), intent(out) :: xinteg(3)

    xinteg(1) = xref(1)
    call vmec_to_can(xref(1), xref(2), xref(3), xinteg(2), xinteg(3))
    xinteg(2) = mod(xinteg(2), twopi)
    xinteg(3) = mod(xinteg(3), twopi)
end subroutine ref_to_integ_flux


!> Evaluate magnetic field in canonical coordinates (r, th_c, ph_c)
!> Works for A_th linear in r (toroidal flux as radial variable)
!> mode_secders = 0: no second derivatives
!> mode_secders = 1: second derivatives only in d/dr^2
!> mode_secders = 2: all second derivatives, including mixed
subroutine eval_field_can(f, r, th_c, ph_c, mode_secders)
    ! Note: !$acc routine seq removed - splint_can_coord uses threadprivate variables
    use get_can_sub, only: splint_can_coord

    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    real(dp) :: Bctr_vartheta, Bctr_varphi, bmod2, sqg, dsqg(3), &
        d2sqg(6), d3Aphdr3, dummy, Bth, Bph, dBth(3), dBph(3), d2Bth(6), &
        d2Bph(6), twobmod, dbmod2(3)

    ! Initialize to zero - no angular derivatives for straight field line Ath(r) Aph(r)
    f%dAth = 0.0_dp
    f%dAph = 0.0_dp

    ! Initialize all 2nd derivatives to zero; mode_secders controls which are computed
    f%d2Ath = 0.0_dp
    f%d2Aph = 0.0_dp
    f%d2hth = 0.0_dp
    f%d2hph = 0.0_dp
    f%d2Bmod = 0.0_dp

    call splint_can_coord(.false., mode_secders, r, th_c, ph_c, &
        f%Ath, f%Aph, f%dAth(1), f%dAph(1), f%d2Aph(1), d3Aphdr3, &
        sqg, dsqg(1), dsqg(2), dsqg(3), &
        Bth, dBth(1), dBth(2), dBth(3), &
        Bph, dBph(1), dBph(2), dBph(3), &
        d2sqg(1), d2sqg(2), d2sqg(3), d2sqg(4), d2sqg(5), d2sqg(6), &
        d2Bth(1), d2Bth(2), d2Bth(3), d2Bth(4), d2Bth(5), d2Bth(6), &
        d2Bph(1), d2Bph(2), d2Bph(3), d2Bph(4), d2Bph(5), d2Bph(6), dummy)

    Bctr_vartheta = -f%dAph(1)/sqg
    Bctr_varphi = f%dAth(1)/sqg

    bmod2 = Bctr_vartheta*Bth + Bctr_varphi*Bph
    f%Bmod = sqrt(abs(bmod2))
    twobmod = 2.0_dp*f%Bmod

    dbmod2(1) = (f%dAth(1)*dBph(1) - f%dAph(1)*dBth(1) &
        - f%d2Aph(1)*Bth - bmod2*dsqg(1))/sqg
    dbmod2(2) = (f%dAth(1)*dBph(2) - f%dAph(1)*dBth(2) - bmod2*dsqg(2))/sqg
    dbmod2(3) = (f%dAth(1)*dBph(3) - f%dAph(1)*dBth(3) - bmod2*dsqg(3))/sqg

    f%dBmod = dbmod2/twobmod

    f%hth = Bth/f%Bmod
    f%hph = Bph/f%Bmod
    f%dhth = dBth/f%Bmod - Bth*f%dBmod/bmod2
    f%dhph = dBph/f%Bmod - Bph*f%dBmod/bmod2

    if (mode_secders > 0) then
        ! d2dr2
        f%d2Bmod(1) = (d2Bph(1)*f%dAth(1) - d2Bth(1)*f%dAph(1) &
            - 2.0_dp*dBth(1)*f%d2Aph(1) - f%Bmod*f%hth*d3Aphdr3 &
            - 2.0_dp*dsqg(1)*dbmod2(1) - bmod2*d2sqg(1))/sqg

        f%d2Bmod(1) = f%d2Bmod(1)/twobmod - f%dBmod(1)**2/f%Bmod

        f%d2hth(1) = d2Bth(1)/f%Bmod - 2.0_dp*dBth(1)*f%dBmod(1)/bmod2 &
            + Bth/bmod2*(2.0_dp*f%dBmod(1)**2/f%Bmod - f%d2Bmod(1))
        f%d2hph(1) = d2Bph(1)/f%Bmod - 2.0_dp*dBph(1)*f%dBmod(1)/bmod2 &
            + Bph/bmod2*(2.0_dp*f%dBmod(1)**2/f%Bmod - f%d2Bmod(1))
    end if

    if (mode_secders == 2) then
        ! d2dth2, d2dph2
        f%d2Bmod(4) = (d2Bph(4)*f%dAth(1) - d2Bth(4)*f%dAph(1) &
            - 2.0_dp*dsqg(2)*dbmod2(2) - bmod2*d2sqg(4))/sqg
        f%d2Bmod(6) = (d2Bph(6)*f%dAth(1) - d2Bth(6)*f%dAph(1) &
            - 2.0_dp*dsqg(3)*dbmod2(3) - bmod2*d2sqg(6))/sqg

        f%d2Bmod(4) = f%d2Bmod(4)/twobmod - f%dBmod(2)**2/f%Bmod
        f%d2Bmod(6) = f%d2Bmod(6)/twobmod - f%dBmod(3)**2/f%Bmod

        f%d2hth((/4,6/)) = d2Bth((/4,6/))/f%Bmod &
            - 2.0_dp*dBth((/2,3/))*f%dBmod((/2,3/))/bmod2 &
            + Bth/bmod2*(2.0_dp*f%dBmod((/2,3/))**2/f%Bmod - f%d2Bmod((/4,6/)))
        f%d2hph((/4,6/)) = d2Bph((/4,6/))/f%Bmod &
            - 2.0_dp*dBph((/2,3/))*f%dBmod((/2,3/))/bmod2 &
            + Bph/bmod2*(2.0_dp*f%dBmod((/2,3/))**2/f%Bmod - f%d2Bmod((/4,6/)))

        ! d2drdth, d2drdph, d2dthdph
        f%d2Bmod(2) = (d2Bph(2)*f%dAth(1) - d2Bth(2)*f%dAph(1) &
            - dBth(2)*f%d2Aph(1) - dsqg(1)*dbmod2(2) &
            - dsqg(2)*dbmod2(1) - bmod2*d2sqg(2))/sqg
        f%d2Bmod(3) = (d2Bph(3)*f%dAth(1) - d2Bth(3)*f%dAph(1) &
            - dBth(3)*f%d2Aph(1) - dsqg(1)*dbmod2(3) &
            - dsqg(3)*dbmod2(1) - bmod2*d2sqg(3))/sqg
        f%d2Bmod(5) = (d2Bph(5)*f%dAth(1) - d2Bth(5)*f%dAph(1) &
            - dsqg(2)*dbmod2(3) - dsqg(3)*dbmod2(2) - bmod2*d2sqg(5))/sqg

        f%d2Bmod(2) = f%d2Bmod(2)/twobmod - f%dBmod(1)*f%dBmod(2)/f%Bmod
        f%d2Bmod(3) = f%d2Bmod(3)/twobmod - f%dBmod(1)*f%dBmod(3)/f%Bmod
        f%d2Bmod(5) = f%d2Bmod(5)/twobmod - f%dBmod(2)*f%dBmod(3)/f%Bmod

        f%d2hth(2) = d2Bth(2)/f%Bmod &
            - (dBth(1)*f%dBmod(2) + dBth(2)*f%dBmod(1))/bmod2 &
            + Bth/bmod2*(2.0_dp*f%dBmod(1)*f%dBmod(2)/f%Bmod - f%d2Bmod(2))
        f%d2hph(2) = d2Bph(2)/f%Bmod &
            - (dBph(1)*f%dBmod(2) + dBph(2)*f%dBmod(1))/bmod2 &
            + Bph/bmod2*(2.0_dp*f%dBmod(1)*f%dBmod(2)/f%Bmod - f%d2Bmod(2))

        f%d2hth(3) = d2Bth(3)/f%Bmod &
            - (dBth(1)*f%dBmod(3) + dBth(3)*f%dBmod(1))/bmod2 &
            + Bth/bmod2*(2.0_dp*f%dBmod(1)*f%dBmod(3)/f%Bmod - f%d2Bmod(3))
        f%d2hph(3) = d2Bph(3)/f%Bmod &
            - (dBph(1)*f%dBmod(3) + dBph(3)*f%dBmod(1))/bmod2 &
            + Bph/bmod2*(2.0_dp*f%dBmod(1)*f%dBmod(3)/f%Bmod - f%d2Bmod(3))

        f%d2hth(5) = d2Bth(5)/f%Bmod &
            - (dBth(2)*f%dBmod(3) + dBth(3)*f%dBmod(2))/bmod2 &
            + Bth/bmod2*(2.0_dp*f%dBmod(2)*f%dBmod(3)/f%Bmod - f%d2Bmod(5))
        f%d2hph(5) = d2Bph(5)/f%Bmod &
            - (dBph(2)*f%dBmod(3) + dBph(3)*f%dBmod(2))/bmod2 &
            + Bph/bmod2*(2.0_dp*f%dBmod(2)*f%dBmod(3)/f%Bmod - f%d2Bmod(5))
    end if
end subroutine eval_field_can

end module field_can_flux
