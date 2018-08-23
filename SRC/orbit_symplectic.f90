module orbit_symplectic

use parmot_mod, only : ro0_parmot => ro0

implicit none

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

double precision, dimension(6) :: yold
double precision :: rmumag
double precision :: dt, ro0

double precision :: coala
double precision :: derphi(3)
double precision :: alambd, pabs

type(field_can) :: f
type(d_field_can) :: df
type(d2_field_can) :: d2f

contains

subroutine orbit_sympl_init(z)

  implicit none
  
  double precision, intent(in) :: z(5)

  call eval_field(z(1), z(2), z(3), 0, f, df, d2f)

  pabs=z(4)
  alambd=z(5)

  rmumag = .5d0*pabs**2*(1.d0-alambd**2)/f%Bmod*2d0 ! rmumag by factor 2 from other modules
  ro0 = ro0_parmot/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules

  yold(1:3) = z(1:3) ! r, theta, varphi
  yold(4) = pabs*alambd*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules
  yold(5) = yold(4)*f%Bth/f%Bmod + f%Ath/ro0
  yold(6) = yold(4)*f%Bph/f%Bmod + f%Aph/ro0

end subroutine orbit_sympl_init

subroutine f_sympl_euler(n, y, fvec, iflag)

  implicit none

  integer, parameter :: mode = 2

  integer, intent(in) :: n
  double precision, intent(in) :: y(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(out) :: iflag

  double precision, dimension(2) :: q, p, w
  double precision, dimension(2) :: dqdt, dpdt, dHdq, dHdw, pqw
  double precision, dimension(2,2) :: dwdq, dwdp

  if (mode==1) then
    q = y(2:3)
    p = yold(5:6)
    w = y((/1,4/))
  else if (mode==2) then ! TODO: debug mode 2
    q = yold(2:3)
    p = y(5:6)
    w = y((/1,4/))
  end if

  dqdt = (y(2:3)-yold(2:3))/dt
  dpdt = (y(5:6)-yold(5:6))/dt

  call get_derivatives(w(1), q(1), q(2), w(2), dHdq, dHdw, pqw, dwdq, dwdp)
  call step_sympl(q, p, w, dqdt, dpdt, dHdq, dHdw, pqw, dwdq, dwdp, fvec)
!  write(5001,*) y
!  write(5002,*) fvec
  !fvec = fvec + exp(-w(1)*1e2) ! to remove negative radius
end

subroutine step_sympl(q, p, w, dqdt, dpdt, dHdq, dHdw, pqw, dwdq, dwdp, ret)

  implicit none

  double precision, intent(out) :: ret(:)
  double precision, intent(in) :: q(:), p(:), w(:)
  double precision, intent(in) :: dqdt(:), dpdt(:)
  double precision, intent(in) :: dHdq(:), dHdw(:)
  double precision, intent(in) :: pqw(:)
  double precision, intent(in) :: dwdq(:,:), dwdp(:,:)

  integer :: k, nq

  nq = size(q)

  ret = 0d0
  do k = 1, nq
    ret(k) = sum(dHdw*dwdp(:,k)) - dqdt(k)
    ret(k + nq) = -sum(dHdw*dwdq(:,k)) - dHdq(k) - dpdt(k)
    ret(k + 2*nq) = p(k) - pqw(k)
  end do

end subroutine step_sympl

subroutine get_derivatives(r, th, ph, vpar, dHdq, dHdw, pqw, dwdq, dwdp)

  implicit none

  double precision, intent(in) :: r, th, ph, vpar
  double precision, dimension(2), intent(out) :: dHdq, dHdw, pqw
  double precision, dimension(2,2), intent(out) :: dwdq, dwdp
  double precision, dimension(3) :: x
  double precision, dimension(2,2) :: dpdw, dpdq

  call eval_field(r, th, ph, 0, f, df, d2f)

  x(1) = r
  x(2) = th
  x(3) = ph
  call elefie_can(x, derphi)

  dHdq(1) = rmumag*df%dBmod(2) + derphi(2)
  dHdq(2) = rmumag*df%dBmod(3) + derphi(3)
  dHdw(1) = rmumag*df%dBmod(1) + derphi(1)
  dHdw(2) = vpar

  pqw(1) = vpar*f%Bth/f%Bmod + f%Ath/ro0
  pqw(2) = vpar*f%Bph/f%Bmod + f%Aph/ro0

  dpdq(1,1) = vpar*(df%dBth(2) - f%Bth*df%dBmod(2)/f%Bmod)/f%Bmod + df%dAth(2)/ro0
  dpdq(1,2) = vpar*(df%dBth(3) - f%Bth*df%dBmod(3)/f%Bmod)/f%Bmod + df%dAth(3)/ro0
  dpdq(2,1) = vpar*(df%dBph(2) - f%Bph*df%dBmod(2)/f%Bmod)/f%Bmod + df%dAph(2)/ro0
  dpdq(2,2) = vpar*(df%dBph(3) - f%Bph*df%dBmod(3)/f%Bmod)/f%Bmod + df%dAph(3)/ro0

  dpdw(1,1) = vpar*(df%dBth(1) - f%Bth*df%dBmod(1)/f%Bmod)/f%Bmod + df%dAth(1)/ro0
  dpdw(1,2) = f%Bth/f%Bmod
  dpdw(2,1) = vpar*(df%dBph(1) - f%Bth*df%dBmod(1)/f%Bmod)/f%Bmod + df%dAph(1)/ro0
  dpdw(2,2) = f%Bph/f%Bmod

  ! inverse matrix
  dwdp(1,1) = dpdw(2,2)
  dwdp(1,2) = -dpdw(1,2)
  dwdp(2,1) = -dpdw(2,1)
  dwdp(2,2) = dpdw(1,1)
  dwdp = dwdp/(dpdw(1,1)*dpdw(2,2)-dpdw(1,2)*dpdw(2,1))

  ! matrix multiplication
  dwdq = -matmul(dwdp,dpdq)

end subroutine get_derivatives

subroutine step_forward(y)
  double precision, dimension(6), intent(inout) :: y
  double precision, dimension(2) :: dqdt, dpdt, dHdq, dHdw, pqw
  double precision, dimension(2,2) :: dwdq, dwdp

  double precision, dimension(2) :: dHdq_can, dHdp_can

  !call eval_field(y(1), y(2), y(3), 0, f, df, d2f)
  call get_derivatives(y(1), y(2), y(3), y(4), dHdq, dHdw, pqw, dwdq, dwdp)

  ! dH(q,p)/dq
  dHdq_can(1) = dHdq(1) + sum(dHdw*dwdq(:,1)) 
  dHdq_can(2) = dHdq(2) + sum(dHdw*dwdq(:,2))
  dHdp_can(1) = sum(dHdw*dwdp(:,1)) 
  dHdp_can(2) = sum(dHdw*dwdp(:,2))

  y(1) = y(1) + sum(dwdq(1,:)*dHdp_can - dwdp(1,:)*dHdq_can)*dt
  y(4) = y(4) + sum(dwdq(2,:)*dHdp_can - dwdp(2,:)*dHdq_can)*dt

  ! for using new r and vpar for angle timestep (pseudo-symplectic Euler)
  !call eval_field(y(1), y(2), y(3), 0, f, df, d2f)
  call get_derivatives(y(1), y(2), y(3), y(4), dHdq, dHdw, pqw, dwdq, dwdp)
  
  y(2) = y(2) + sum(dHdw*dwdp(:,1))*dt
  y(3) = y(3) + sum(dHdw*dwdp(:,2))*dt
  
end subroutine step_forward

subroutine orbit_timestep_sympl(z, dtau, dtaumin, ierr)

  implicit none
  
  integer, parameter :: ndim = 5
  integer, parameter :: n = 6
  double precision :: tol

  integer, intent(out) :: ierr
  double precision, dimension(ndim), intent(inout) :: z
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  double precision, dimension(6) :: y, fvec

  double precision tau2
  
  ierr = 0

  ! for nleq1

  ! INTEGER NMAX, LIOPT, LIWK, LRWK, LUPRT
  ! PARAMETER ( NMAX=50 )
  ! PARAMETER ( LIOPT=50, LIWK=NMAX+50, LRWK=(NMAX+13)*NMAX+60 )
  ! PARAMETER ( LUPRT=6 )
  ! EXTERNAL FCN, NLEQ1
  ! DOUBLE PRECISION XSCAL(NMAX)
  ! INTEGER IOPT(LIOPT), IWK(LIWK)
  ! DOUBLE PRECISION RWK(LRWK)
  ! INTEGER I, NIW, NRW, DUMMY
  ! CHARACTER CHGDAT*20, PRODCT*8

  ! end for nleq1
  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0.0
  do while(tau2.lt.dtau)
    y = yold
    call step_forward(y)
    tol = 1d-10
    call hybrd1 (f_sympl_euler, n, y, fvec, tol, ierr)
    !IOPT = 0
    !iopt(11) = 1 ! print error messages
    !IWK = 0
    !iwk(31) = 1000 ! maximum iterations, default: 50
    !RWK = 0.0D0
    !XSCAL = 0.0D0
    !CALL nleq1(n,f_sympl_euler,dummy,y,XSCAL,tol,IOPT,ierr,LIWK,IWK,LRWK,RWK)
    !write(5001,*) '---------------------------------------------'
    !write(5002,*) '---------------------------------------------'
    if(ierr > 1) then
      print *, 'orbit_timestep_sympl: warning ', ierr, ', step at half accuracy'
      tol = 1d-5
      ierr = 0
      !y = yold
      !call step_forward(y)
      call hybrd1 (f_sympl_euler, n, y, fvec, tol, ierr)
      if(ierr > 1) stop
      !if(ierr /= 4 .and. ierr /= 5) stop
    end if
    yold = y
    tau2=tau2+dtaumin
  enddo
  z = 0.0
  z(1:4) = y(1:4)
  z(5) = y(4)**2/2d0 + rmumag*f%Bmod
end

end module orbit_symplectic

subroutine eval_field(r, th_c, ph_c, mode_secders, f, df, d2f)

  use orbit_symplectic, only: field_can, d_field_can, d2_field_can

  implicit none

  double precision, intent(in) :: r, th_c, ph_c     
  integer, intent(in) :: mode_secders          
                             
  type(field_can), intent(out) :: f
  type(d_field_can), intent(out) :: df
  type(d2_field_can), intent(out) :: d2f

  double precision :: Bctr_vartheta, Bctr_varphi, bmod2, sqg, dsqg(3), d2sqg(6), d3Aphdr3, dummy

  ! initialize to zero - no angular derivatives will be set due to straight field line Ath(r) Aph(r)
  df%dAth = 0d0
  df%dAph = 0d0
  d2f%d2Ath = 0d0
  d2f%d2Aph = 0d0

  call splint_can_coord(.false., mode_secders, r, th_c, ph_c,                             &
    f%Ath, f%Aph, df%dAth(1), df%dAph(1), d2f%d2Aph(1), d3Aphdr3,                       &
    sqg, dsqg(1), dsqg(2), dsqg(3),                                                     &
    f%Bth, df%dBth(1), df%dBth(2), df%dBth(3),                                          &
    f%Bph, df%dBph(1), df%dBph(2), df%dBph(3),                                          &
    d2sqg(1), d2sqg(2), d2sqg(3), d2sqg(4), d2sqg(5), d2sqg(6), &
    d2f%d2Bth(1), d2f%d2Bth(2), d2f%d2Bth(3), d2f%d2Bth(4), d2f%d2Bth(5), d2f%d2Bth(6), &
    d2f%d2Bph(1), d2f%d2Bph(2), d2f%d2Bph(3), d2f%d2Bph(4), d2f%d2Bph(5), d2f%d2Bph(6), dummy)
 
  Bctr_vartheta=-df%dAph(1)/sqg
  Bctr_varphi=df%dAth(1)/sqg
!
  bmod2=Bctr_vartheta*f%Bth+Bctr_varphi*f%Bph
  f%Bmod=sqrt(bmod2)
!
  df%dBmod(1) = 0.5d0*((df%dAth(1)*df%dBph(1)-df%dAph(1)*df%dBth(1)-d2f%d2Aph(1)*f%Bth)/f%Bmod-dsqg(1)*f%Bmod)/sqg
  df%dBmod(2) = 0.5d0*((df%dAth(1)*df%dBph(2)-df%dAth(1)*df%dBth(2))/f%Bmod-dsqg(2)*f%Bmod)/sqg
  df%dBmod(3) = 0.5d0*((df%dAth(1)*df%dBph(3)-df%dAth(1)*df%dBth(3))/f%Bmod-dsqg(3)*f%Bmod)/sqg

end subroutine eval_field

! for testing -> circular tokamak
! subroutine eval_field(r, th, ph, mode_secders, f, df, d2f)
!   use orbit_symplectic, only: field_can, d_field_can, d2_field_can
!   implicit none

!   double precision, intent(in) :: r, th, ph   
!   integer, intent(in) :: mode_secders        
                             
!   type(field_can), intent(out) :: f
!   type(d_field_can), intent(out) :: df
!   type(d2_field_can), intent(out) :: d2f

!   double precision :: B0th, B0ph, cth, sth 
!   B0th = .99d0
!   B0ph = sqrt(1d0-B0th**2)

!   cth = cos(th)
!   sth = sin(th)
  
!   f%Ath      = B0ph*(r**2/2d0 - r**3/3d0*cth)
!   df%dAth(1) = B0ph*(r - r**2*cth)
!   df%dAth(2) = B0ph*r**3/3d0*sth
!   df%dAth(3) = 0d0

!   f%Aph     = -B0th*r
!   df%dAph(1) = -B0th
!   df%dAph(2) = 0d0
!   df%dAph(3) = 0d0

!   f%Bth      = B0th*r*(1d0 - r*cth)
!   df%dBth(1) = B0th*(1d0 - 2d0*r*cth)
!   df%dBth(2) = B0th*r**2*sth
!   df%dBth(3) = 0d0
  
!   f%Bph      = B0ph*(1d0 - (r*cth)**2)
!   df%dBph(1) = -2d0*B0ph*r*cth**2
!   df%dBph(2) = 2d0*B0ph*r**2*cth*sth
!   df%dBph(3) = 0d0

!   f%Bmod   = 1d0 - r*cth
!   df%dBmod(1) = -cth
!   df%dBmod(2) = r*sth
!   df%dBmod(3) = 0d0

!   ! TODO: second derivatives
!   d2f%d2Ath = 0d0
!   d2f%d2Aph = 0d0
!   d2f%d2Bth = 0d0
!   d2f%d2Bph = 0d0
!   d2f%d2Bmod = 0d0

! end subroutine eval_field

! TODO: Check with VMEC field
! !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! !
!   subroutine splint_vmec_extra(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
! !
! ! Computes magnetic field module in units of the magnetic code  - bmod,
! ! square root of determinant of the metric tensor               - sqrtg,
! ! derivatives of the logarythm of the magnetic field module
! ! over coordinates                                              - bder,
! ! covariant componets of the unit vector of the magnetic
! ! field direction                                               - hcovar,
! ! contravariant components of this vector                       - hctrvr,
! ! contravariant component of the curl of this vector            - hcurl
! ! Order of coordinates is the following: x(1)=R (big radius),
! ! x(2)=phi (toroidal angle), x(3)=Z (altitude).
! !
! !  Input parameters:
! !            formal:  x(3)             - array of VMEC coordinates
! !  Output parameters:
! !            formal:  bmod
! !                     sqrtg
! !                     bder(3)          - derivatives of $\log(B)$
! !                     hcovar(3)        - covariant components of unit vector $\bh$ along $\bB$
! !                     hctrvr(3)        - contra-variant components of unit vector $\bh$ along $\bB$
! !                     hcurl(3)         - contra-variant components of curl of $\bh$
! !
! !  Called routines: vmec_field
! !
!   implicit none
! !
!   double precision, parameter :: twopi=2.d0*3.14159265358979d0, hs=1.d-3, ht=hs*twopi, hp=ht/5.d0
! !
!   double precision :: bmod,sqrtg
!   double precision :: s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                       sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                       Bcovar_r,Bcovar_vartheta,Bcovar_varphi
!   double precision :: cjac,bcov_s_vmec,bcov_t_vmec,bcov_p_vmec
!   double precision :: dhs_dt,dhs_dp,dht_ds,dht_dp,dhp_ds,dhp_dt
!   double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
! !
! ! Begin derivatives over s
! !
!   theta=x(2)
!   varphi=x(3)
!   s=x(1)+hs
! !
!   call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                   sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                   Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
! !
!   bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
!   bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
!   bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
!   bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
!   bder(1)=bmod
!   dht_ds=bcov_t_vmec/bmod
!   dhp_ds=bcov_p_vmec/bmod
! !
!   s=x(1)-hs
! !
!   call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                   sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                   Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
! !
!   bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
!   bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
!   bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
!   bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
!   bder(1)=(bder(1)-bmod)/(2.d0*hs)
!   dht_ds=(dht_ds-bcov_t_vmec/bmod)/(2.d0*hs)
!   dhp_ds=(dhp_ds-bcov_p_vmec/bmod)/(2.d0*hs)
! !
! ! End derivatives over s
! !
! !-------------------------
! !
! ! Begin derivatives over theta
! !
!   s=x(1)
!   theta=x(2)+ht
! !
!   call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                   sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                   Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
! !
!   bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
!   bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
!   bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
!   bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
!   bder(2)=bmod
!   dhs_dt=bcov_s_vmec/bmod
!   dhp_dt=bcov_p_vmec/bmod
! !
!   theta=x(2)-ht
! !
!   call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                   sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                   Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
! !
!   bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
!   bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
!   bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
!   bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
!   bder(2)=(bder(2)-bmod)/(2.d0*ht)
!   dhs_dt=(dhs_dt-bcov_s_vmec/bmod)/(2.d0*ht)
!   dhp_dt=(dhp_dt-bcov_p_vmec/bmod)/(2.d0*ht)
! !
! ! End derivatives over theta
! !
! !-------------------------
! !
! ! Begin derivatives over varphi
! !
!   theta=x(2)
!   varphi=x(3)+hp
! !
!   call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                   sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                   Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
! !
!   bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
!   bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
!   bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
!   bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
!   bder(3)=bmod
!   dhs_dp=bcov_s_vmec/bmod
!   dht_dp=bcov_t_vmec/bmod
! !
!   varphi=x(3)-hp
! !
!   call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                   sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                   Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
! !
!   bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
!   bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
!   bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
!   bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
!   bder(3)=(bder(3)-bmod)/(2.d0*hp)
!   dhs_dp=(dhs_dp-bcov_s_vmec/bmod)/(2.d0*hp)
!   dht_dp=(dht_dp-bcov_t_vmec/bmod)/(2.d0*hp)
! !
! ! End derivatives over varphi
! !
! !-------------------------
! !
!   varphi=x(3)
! !
!   call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                   sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                   Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
! !
!   bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
!   cjac=1.d0+dl_dt
!   sqrtg=sqg*cjac
!   bder=bder/bmod
!   bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
!   bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
!   bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
!   hcovar(1)=bcov_s_vmec/bmod
!   hcovar(2)=bcov_t_vmec/bmod
!   hcovar(3)=bcov_p_vmec/bmod
!   hctrvr(1)=0.d0
!   hctrvr(2)=(Bctrvr_vartheta-dl_dp*Bctrvr_varphi)/(cjac*bmod)
!   hctrvr(3)=Bctrvr_varphi/bmod
!   hcurl(1)=(dhp_dt-dht_dp)/sqrtg
!   hcurl(2)=(dhs_dp-dhp_ds)/sqrtg
!   hcurl(3)=(dht_ds-dhs_dt)/sqrtg
! !
!   end subroutine splint_vmec_extra

