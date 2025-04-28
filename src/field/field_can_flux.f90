module field_can_flux

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_base, only: FieldCan, n_field_evaluations, twopi

implicit none

contains

subroutine evaluate_flux(f, r, th_c, ph_c, mode_secders)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    call eval_field_can(f, r, th_c, ph_c, mode_secders)

    n_field_evaluations = n_field_evaluations + 1
end subroutine evaluate_flux


subroutine can_to_ref_flux(xcan, xref)
    use get_can_sub, only: can_to_vmec
    real(dp), intent(in) :: xcan(3)
    real(dp), intent(out) :: xref(3)

    xref(1) = xcan(1)
    call can_to_vmec(xcan(1), xcan(2), xcan(3), xref(2), xref(3))
    xref(2) = mod(xref(2), twopi)
    xref(3) = mod(xref(3), twopi)
end subroutine can_to_ref_flux


subroutine ref_to_can_flux(xref, xcan)
    use get_can_sub, only: vmec_to_can
    real(dp), intent(in) :: xref(3)
    real(dp), intent(out) :: xcan(3)

    xcan(1) = xref(1)
    call vmec_to_can(xref(1), xref(2), xref(3), xcan(2), xcan(3))
    xcan(2) = mod(xcan(2), twopi)
    xcan(3) = mod(xcan(3), twopi)
end subroutine ref_to_can_flux


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine eval_field_can(f, r, th_c, ph_c, mode_secders)
!
! Evaluates magnetic field in canonical coordinates (r, th_c, ph_c)
! and stores results in variable f
! Works for A_th linear in r (toroidal flux as radial variable)
!
! mode_secders = 0: no second derivatives
! mode_secders = 1: second derivatives only in d/dr^2
! mode_secders = 2: all second derivatives, including mixed
!
! tested in test_magfie.f90, 2018-10-23, C. Albert <albert@alumni.tugraz.at>
!

    use get_can_sub, only : splint_can_coord

    type(FieldCan), intent(inout) :: f
    double precision, intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    double precision :: Bctr_vartheta, Bctr_varphi, bmod2, sqg, dsqg(3), &
    d2sqg(6), d3Aphdr3, dummy, Bth, Bph, dBth(3), dBph(3), d2Bth(6), &
    d2Bph(6), twobmod, dbmod2(3)

    ! Count evaluations for profiling

    ! initialize to zero - no angular derivatives will be set due to straight field line Ath(r) Aph(r)
    f%dAth = 0d0
    f%dAph = 0d0

    ! initialize all 2nd derivatives to zero, as the mode decides, which one to use
    f%d2Ath = 0d0
    f%d2Aph = 0d0
    f%d2hth = 0d0
    f%d2hph = 0d0
    f%d2Bmod = 0d0

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
    !if (bmod2<0) print *, r, th_c, ph_c, bmod2
    f%Bmod = sqrt(abs(bmod2))
    twobmod = 2.d0*f%Bmod

    dbmod2(1) = (f%dAth(1)*dBph(1)-f%dAph(1)*dBth(1)-f%d2Aph(1)*Bth-bmod2*dsqg(1))/sqg
    dbmod2(2) = (f%dAth(1)*dBph(2)-f%dAph(1)*dBth(2)-bmod2*dsqg(2))/sqg
    dbmod2(3) = (f%dAth(1)*dBph(3)-f%dAph(1)*dBth(3)-bmod2*dsqg(3))/sqg

    f%dBmod = dbmod2/twobmod

    f%hth = Bth/f%Bmod
    f%hph = Bph/f%Bmod
    f%dhth = dBth/f%Bmod - Bth*f%dBmod/bmod2
    f%dhph = dBph/f%Bmod - Bph*f%dBmod/bmod2

    if(mode_secders > 0) then
    ! d2dr2
    f%d2Bmod(1)=(d2Bph(1)*f%dAth(1)-d2Bth(1)*f%dAph(1)-2.d0*dBth(1)*f%d2Aph(1)-f%Bmod*f%hth*d3Aphdr3 &
                - 2.d0*dsqg(1)*dbmod2(1)-bmod2*d2sqg(1))/sqg

    f%d2Bmod(1)=f%d2Bmod(1)/twobmod-f%dBmod(1)**2/f%Bmod

    f%d2hth(1) = d2Bth(1)/f%Bmod - 2d0*dBth(1)*f%dBmod(1)/bmod2 + Bth/bmod2*(2d0*f%dBmod(1)**2/f%Bmod - f%d2Bmod(1))
    f%d2hph(1) = d2Bph(1)/f%Bmod - 2d0*dBph(1)*f%dBmod(1)/bmod2 + Bph/bmod2*(2d0*f%dBmod(1)**2/f%Bmod - f%d2Bmod(1))
    endif

    if(mode_secders.eq.2) then
    ! d2dth2, d2dph2
    f%d2Bmod(4)=(d2Bph(4)*f%dAth(1)-d2Bth(4)*f%dAph(1)-2.d0*dsqg(2)*dbmod2(2)-bmod2*d2sqg(4))/sqg
    f%d2Bmod(6)=(d2Bph(6)*f%dAth(1)-d2Bth(6)*f%dAph(1)-2.d0*dsqg(3)*dbmod2(3)-bmod2*d2sqg(6))/sqg

    f%d2Bmod(4) = f%d2Bmod(4)/twobmod - f%dBmod(2)**2/f%Bmod
    f%d2Bmod(6) = f%d2Bmod(6)/twobmod - f%dBmod(3)**2/f%Bmod

    f%d2hth((/4,6/)) = d2Bth((/4,6/))/f%Bmod - 2d0*dBth((/2,3/))*f%dBmod((/2,3/))/bmod2 &
    + Bth/bmod2*(2d0*f%dBmod((/2,3/))**2/f%Bmod - f%d2Bmod((/4,6/)))
    f%d2hph((/4,6/)) = d2Bph((/4,6/))/f%Bmod - 2d0*dBph((/2,3/))*f%dBmod((/2,3/))/bmod2 &
    + Bph/bmod2*(2d0*f%dBmod((/2,3/))**2/f%Bmod - f%d2Bmod((/4,6/)))

    ! d2drdth, d2drdph, d2dthdph
    f%d2Bmod(2)=(d2Bph(2)*f%dAth(1)-d2Bth(2)*f%dAph(1)-dBth(2)*f%d2Aph(1) &
        - dsqg(1)*dbmod2(2)-dsqg(2)*dbmod2(1)-bmod2*d2sqg(2))/sqg
    f%d2Bmod(3)=(d2Bph(3)*f%dAth(1)-d2Bth(3)*f%dAph(1)-dBth(3)*f%d2Aph(1) &
        - dsqg(1)*dbmod2(3)-dsqg(3)*dbmod2(1)-bmod2*d2sqg(3))/sqg
    f%d2Bmod(5)=(d2Bph(5)*f%dAth(1)-d2Bth(5)*f%dAph(1)-dsqg(2)*dbmod2(3) &
        - dsqg(3)*dbmod2(2)-bmod2*d2sqg(5))/sqg

    f%d2Bmod(2) = f%d2Bmod(2)/twobmod - f%dBmod(1)*f%dBmod(2)/f%Bmod
    f%d2Bmod(3) = f%d2Bmod(3)/twobmod - f%dBmod(1)*f%dBmod(3)/f%Bmod
    f%d2Bmod(5) = f%d2Bmod(5)/twobmod - f%dBmod(2)*f%dBmod(3)/f%Bmod

    f%d2hth(2) = d2Bth(2)/f%Bmod - (dBth(1)*f%dBmod(2) + dBth(2)*f%dBmod(1))/bmod2 &
        + Bth/bmod2*(2d0*f%dBmod(1)*f%dBmod(2)/f%Bmod - f%d2Bmod(2))
    f%d2hph(2) = d2Bph(2)/f%Bmod - (dBph(1)*f%dBmod(2) + dBph(2)*f%dBmod(1))/bmod2 &
        + Bph/bmod2*(2d0*f%dBmod(1)*f%dBmod(2)/f%Bmod - f%d2Bmod(2))

    f%d2hth(3) = d2Bth(3)/f%Bmod - (dBth(1)*f%dBmod(3) + dBth(3)*f%dBmod(1))/bmod2 &
        + Bth/bmod2*(2d0*f%dBmod(1)*f%dBmod(3)/f%Bmod - f%d2Bmod(3))
    f%d2hph(3) = d2Bph(3)/f%Bmod - (dBph(1)*f%dBmod(3) + dBph(3)*f%dBmod(1))/bmod2 &
        + Bph/bmod2*(2d0*f%dBmod(1)*f%dBmod(3)/f%Bmod - f%d2Bmod(3))

    f%d2hth(5) = d2Bth(5)/f%Bmod - (dBth(2)*f%dBmod(3) + dBth(3)*f%dBmod(2))/bmod2 &
        + Bth/bmod2*(2d0*f%dBmod(2)*f%dBmod(3)/f%Bmod - f%d2Bmod(5))
    f%d2hph(5) = d2Bph(5)/f%Bmod - (dBph(2)*f%dBmod(3) + dBph(3)*f%dBmod(2))/bmod2 &
        + Bph/bmod2*(2d0*f%dBmod(2)*f%dBmod(3)/f%Bmod - f%d2Bmod(5))
    endif

end subroutine eval_field_can

end module field_can_flux
