module field_can_mod

use diag_mod, only : icounter

implicit none

type :: FieldCan
  double precision :: Ath, Aph
  double precision :: hth, hph
  double precision :: Bmod

  double precision, dimension(3) :: dAth, dAph
  double precision, dimension(3) :: dhth, dhph
  double precision, dimension(3) :: dBmod

  ! second derivatives: drdr, drdhth, drdph, dthdth, dthdph, dphdph
  double precision, dimension(6) :: d2Ath, d2Aph
  double precision, dimension(6) :: d2hth, d2hph
  double precision, dimension(6) :: d2Bmod

  double precision :: H, pth, vpar
  double precision, dimension(4) :: dvpar, dH, dpth
  ! order of second derivatives:
  ! d2dr2, d2drdth, d2drph, d2dth2, d2dthdph, d2dph2,
  ! d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
  double precision, dimension(10) :: d2vpar, d2H, d2pth

  double precision :: mu, ro0
end type FieldCan

interface eval_field
  module procedure eval_field_can
end interface

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine get_val(f, pphi)
!
! computes values of H, pth and vpar at z=(r, th, ph, pphi)
!
!
  type(FieldCan), intent(inout) :: f
  double precision, intent(in) :: pphi

  f%vpar = (pphi - f%Aph/f%ro0)/f%hph
  f%H = f%vpar**2/2d0 + f%mu*f%Bmod
  f%pth = f%hth*f%vpar + f%Ath/f%ro0

end subroutine get_val


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine get_derivatives(f, pphi)
!
! computes H, pth and vpar at z=(r, th, ph, pphi) and their derivatives
!
!
  type(FieldCan), intent(inout) :: f
  double precision, intent(in) :: pphi

  call get_val(f, pphi)

  f%dvpar(1:3) = -(f%dAph/f%ro0 + f%dhph*f%vpar)/f%hph
  f%dvpar(4)   = 1d0/f%hph

  f%dH(1:3) = f%vpar*f%dvpar(1:3) + f%mu*f%dBmod
  f%dH(4)   = f%vpar/f%hph

  f%dpth(1:3) = f%dvpar(1:3)*f%hth + f%vpar*f%dhth + f%dAth/f%ro0

  f%dpth(4) = f%hth/f%hph

end subroutine get_derivatives

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine get_derivatives2(f, pphi)
!
! computes H, pth and vpar at z=(r, th, ph, pphi) up to 2nd derivatives
! order of second derivatives:
! d2dr2, d2drdth, d2drph, d2dth2, d2dthdph, d2dph2,
! d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
!
  type(FieldCan), intent(inout) :: f
  double precision, intent(in) :: pphi

  call get_derivatives(f, pphi)

  f%d2vpar = 0d0
  f%d2H    = 0d0
  f%d2pth  = 0d0

  f%d2vpar(1:6) = -f%d2Aph/f%ro0 - f%d2hph*f%vpar
  f%d2vpar((/1,4,6/)) = f%d2vpar((/1,4,6/)) - 2d0*f%dhph*f%dvpar(1:3)
  f%d2vpar(2) = f%d2vpar(2) - (f%dhph(1)*f%dvpar(2) + f%dhph(2)*f%dvpar(1))
  f%d2vpar(3) = f%d2vpar(3) - (f%dhph(1)*f%dvpar(3) + f%dhph(3)*f%dvpar(1))
  f%d2vpar(5) = f%d2vpar(5) - (f%dhph(2)*f%dvpar(3) + f%dhph(3)*f%dvpar(2))
  f%d2vpar(1:6) = f%d2vpar(1:6)/f%hph

  f%d2H(1:6) = f%vpar*f%d2vpar(1:6) + f%mu*f%d2Bmod ! + qi*d2Phie
  f%d2H((/1,4,6/)) = f%d2H((/1,4,6/)) + f%dvpar(1:3)**2
  f%d2H(2) = f%d2H(2) + f%dvpar(1)*f%dvpar(2)
  f%d2H(3) = f%d2H(3) + f%dvpar(1)*f%dvpar(3)
  f%d2H(5) = f%d2H(5) + f%dvpar(2)*f%dvpar(3)

  f%d2pth(1:6) = f%d2vpar(1:6)*f%hth + f%vpar*f%d2hth + f%d2Ath/f%ro0
  f%d2pth((/1,4,6/)) = f%d2pth((/1,4,6/)) + 2d0*f%dvpar(1:3)*f%dhth
  f%d2pth(2) = f%d2pth(2) + f%dvpar(1)*f%dhth(2) + f%dvpar(2)*f%dhth(1)
  f%d2pth(3) = f%d2pth(3) + f%dvpar(1)*f%dhth(3) + f%dvpar(3)*f%dhth(1)
  f%d2pth(5) = f%d2pth(5) + f%dvpar(2)*f%dhth(3) + f%dvpar(3)*f%dhth(2)

  f%d2vpar(7:9) = -f%dhph/f%hph**2
  f%d2H(7:9) = f%dvpar(1:3)/f%hph + f%vpar*f%d2vpar(7:9)
  f%d2pth(7:9) = f%dhth/f%hph + f%hth*f%d2vpar(7:9)
end subroutine get_derivatives2


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

  ! TODO check if this is working:
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

! for testing -> circular tokamak
! subroutine eval_field_test(r, th, ph, mode_secders)
! !

!  double precision, intent(in) :: r, th, ph
!  integer, intent(in) :: mode_secders

!  double precision :: B0th, B0ph, cth, sth

!    ! Count evaluations for profiling
!   icounter = icounter + 1

!   B0th = .99d0
!   B0ph = sqrt(1d0-B0th**2)

!   cth = cos(th)
!   sth = sin(th)

!   f%Ath      = B0ph*(r**2/2d0 - r**3/3d0*cth)
!   f%dAth(1) = B0ph*(r - r**2*cth)
!   f%dAth(2) = B0ph*r**3/3d0*sth
!   f%dAth(3) = 0d0

!   f%Aph     = -B0th*r
!   f%dAph(1) = -B0th
!   f%dAph(2) = 0d0
!   f%dAph(3) = 0d0

!   f%Bth      = B0th*r*(1d0 - r*cth)
!   f%dBth(1) = B0th*(1d0 - 2d0*r*cth)
!   f%dBth(2) = B0th*r**2*sth
!   f%dBth(3) = 0d0

!   f%Bph      = B0ph*(1d0 - (r*cth)**2)
!   Bph(1) = -2d0*B0ph*r*cth**2
!   Bph(2) = 2d0*B0ph*r**2*cth*sth
!   Bph(3) = 0d0

!   f%Bmod   = 1d0 - r*cth
!   f%dBmod(1) = -cth
!   f%dBmod(2) = r*sth
!   f%dBmod(3) = 0d0

! !  TODO: second derivatives
!   f%d2Ath = 0d0
!   f%d2Aph = 0d0
!   f%d2Bth = 0d0
!   f%d2Bph = 0d0
!   f%d2Bmod = 0d0

! end subroutine eval_field_test


end module field_can_mod
