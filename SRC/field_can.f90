module field_can_mod

type :: field_can
    double precision :: Ath, Aph
    double precision :: Bth, Bph
    double precision :: Bmod
end type field_can

type :: d_field_can
    double precision, dimension(3) :: dAth, dAph
    double precision, dimension(3) :: dBth, dBph
    double precision, dimension(3) :: dBmod
end type d_field_can

! second derivatives: drdr, drdt, drdp, dtdt, dtdp, dpdp
type :: d2_field_can 
    double precision, dimension(6) :: d2Ath, d2Aph  
    double precision, dimension(6) :: d2Bth, d2Bph  
    double precision, dimension(6) :: d2Bmod
end type d2_field_can

contains

! subroutine eval_field(r, th_c, ph_c, mode_secders, f, df, d2f)

!   use orbit_symplectic, only: field_can, d_field_can, d2_field_can

!   implicit none

!   double precision, intent(in) :: r, th_c, ph_c     
!   integer, intent(in) :: mode_secders          
                             
!   type(field_can), intent(out) :: f
!   type(d_field_can), intent(out) :: df
!   type(d2_field_can), intent(out) :: d2f

!   double precision :: Bctr_vartheta, Bctr_varphi, bmod2, sqg, dsqg(3), d2sqg(6), d3Aphdr3, dummy

!   ! initialize to zero - no angular derivatives will be set due to straight field line Ath(r) Aph(r)
!   df%dAth = 0d0
!   df%dAph = 0d0
!   d2f%d2Ath = 0d0
!   d2f%d2Aph = 0d0

!   call splint_can_coord(.false., mode_secders, r, th_c, ph_c,                           &
!     f%Ath, f%Aph, df%dAth(1), df%dAph(1), d2f%d2Aph(1), d3Aphdr3,                       &
!     sqg, dsqg(1), dsqg(2), dsqg(3),                                                     &
!     f%Bth, df%dBth(1), df%dBth(2), df%dBth(3),                                          &
!     f%Bph, df%dBph(1), df%dBph(2), df%dBph(3),                                          &
!     d2sqg(1), d2sqg(2), d2sqg(3), d2sqg(4), d2sqg(5), d2sqg(6),                         &
!     d2f%d2Bth(1), d2f%d2Bth(2), d2f%d2Bth(3), d2f%d2Bth(4), d2f%d2Bth(5), d2f%d2Bth(6), &
!     d2f%d2Bph(1), d2f%d2Bph(2), d2f%d2Bph(3), d2f%d2Bph(4), d2f%d2Bph(5), d2f%d2Bph(6), dummy)
 
!   Bctr_vartheta = -df%dAph(1)/sqg
!   Bctr_varphi = df%dAth(1)/sqg
  
!   bmod2 = Bctr_vartheta*f%Bth + Bctr_varphi*f%Bph
!   f%Bmod = sqrt(bmod2)
  
!   df%dBmod(1) = 0.5d0*((df%dAth(1)*df%dBph(1)-df%dAph(1)*df%dBth(1)-d2f%d2Aph(1)*f%Bth)/f%Bmod-dsqg(1)*f%Bmod)/sqg
!   df%dBmod(2) = 0.5d0*((df%dAth(1)*df%dBph(2)-df%dAth(1)*df%dBth(2))/f%Bmod-dsqg(2)*f%Bmod)/sqg
!   df%dBmod(3) = 0.5d0*((df%dAth(1)*df%dBph(3)-df%dAth(1)*df%dBth(3))/f%Bmod-dsqg(3)*f%Bmod)/sqg

! end subroutine eval_field

! for testing -> circular tokamak
subroutine eval_field(r, th, ph, mode_secders, f, df, d2f)
  implicit none

  double precision, intent(in) :: r, th, ph   
  integer, intent(in) :: mode_secders        
                             
  type(field_can), intent(out) :: f
  type(d_field_can), intent(out) :: df
  type(d2_field_can), intent(out) :: d2f

  double precision :: B0th, B0ph, cth, sth 
  B0th = .99d0
  B0ph = sqrt(1d0-B0th**2)

  cth = cos(th)
  sth = sin(th)
  
  f%Ath      = B0ph*(r**2/2d0 - r**3/3d0*cth)
  df%dAth(1) = B0ph*(r - r**2*cth)
  df%dAth(2) = B0ph*r**3/3d0*sth
  df%dAth(3) = 0d0

  f%Aph     = -B0th*r
  df%dAph(1) = -B0th
  df%dAph(2) = 0d0
  df%dAph(3) = 0d0

  f%Bth      = B0th*r*(1d0 - r*cth)
  df%dBth(1) = B0th*(1d0 - 2d0*r*cth)
  df%dBth(2) = B0th*r**2*sth
  df%dBth(3) = 0d0
  
  f%Bph      = B0ph*(1d0 - (r*cth)**2)
  df%dBph(1) = -2d0*B0ph*r*cth**2
  df%dBph(2) = 2d0*B0ph*r**2*cth*sth
  df%dBph(3) = 0d0

  f%Bmod   = 1d0 - r*cth
  df%dBmod(1) = -cth
  df%dBmod(2) = r*sth
  df%dBmod(3) = 0d0

  ! TODO: second derivatives
  d2f%d2Ath = 0d0
  d2f%d2Aph = 0d0
  d2f%d2Bth = 0d0
  d2f%d2Bph = 0d0
  d2f%d2Bmod = 0d0

end subroutine eval_field


end module field_can_mod