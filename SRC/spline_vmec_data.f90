!
  subroutine spline_vmec_data
!
  use new_vmec_stuff_mod
  use vector_potentail_mod, only : ns,hs,torflux,sA_phi
!
  implicit none
!
  integer :: i,k,m,n,is,i_theta,i_phi,m_max,n_max,nsize_exp_imt,nsize_exp_inp,iexpt,iexpp
  integer :: iss,ist,isp,nrho,nheal
  double precision :: twopi,cosphase,sinphase
  complex(8)   :: base_exp_imt,base_exp_inp,base_exp_inp_inv,expphase
  double precision, dimension(:,:), allocatable :: splcoe
  double precision, dimension(:,:), allocatable :: almn_rho,rmn_rho,zmn_rho
  complex(8),   dimension(:),   allocatable :: exp_imt,exp_inp
!
  print *,'Splining VMEC data: ns_A = ',ns_A,'  ns_s = ',ns_s,'  ns_tp = ',ns_tp
!
  call new_allocate_vmec_stuff
!
  call vmecin(rmn,zmn,almn,aiota,phi,sps,axm,axn,s,    &
              nsurfm,nstrm,kpar,torflux)
!
  ns=kpar+1
  allocate(splcoe(0:ns_A,ns))
  hs=s(2)-s(1)
!
  nrho=ns
  allocate(almn_rho(nstrm,0:nrho-1),rmn_rho(nstrm,0:nrho-1),zmn_rho(nstrm,0:nrho-1))
!
  do i=1,nstrm
!
    m=nint(abs(axm(i)))
!
    nheal=min(m,10)
!
    call s_to_rho_healaxis(m,ns,nrho,nheal,rmn(i,:),rmn_rho(i,:))
!
    call s_to_rho_healaxis(m,ns,nrho,nheal,zmn(i,:),zmn_rho(i,:))
!
    call s_to_rho_healaxis(m,ns,nrho,nheal,almn(i,:),almn_rho(i,:))
!
  enddo
!
!------------------------------------
! Begin poloidal flux ($A_\varphi$):
!
  splcoe(0,:)=aiota
!
  call spl_reg(ns_A-1,ns,hs,splcoe(0:ns_A-1,:))
!
  do i=ns_A,1,-1
    splcoe(i,:)=splcoe(i-1,:)/dble(i)
  enddo
!
  splcoe(0,1)=0.d0
  do is=1,ns-1
    splcoe(0,is+1)=splcoe(ns_A,is)
    do k=ns_A-1,0,-1
      splcoe(0,is+1)=splcoe(k,is)+hs*splcoe(0,is+1)
    enddo
  enddo
!
  if(allocated(sA_phi)) deallocate(sA_phi)
  allocate(sA_phi(ns_A+1,ns))
  do k=0,ns_A
    sA_phi(k+1,:)=-torflux*splcoe(k,:)
  enddo
!
  deallocate(splcoe)
! End poloidal flux
!------------------------------
!
! Begin angular grid, sin and cos
!
  m_max=nint(maxval(axm))
  n_max=nint(maxval(axn))
!
  print *,'VMEC ns = ',ns,' m_max = ',m_max,' n_max = ',n_max
!
  n_theta = m_max*multharm+1
  n_phi = n_max*multharm+1
  twopi=8.d0*atan2(1.d0,1.d0)
  h_theta=twopi/dble(n_theta-1)
  h_phi=twopi/dble((n_phi-1)*nper)
!
  nsize_exp_imt=(n_theta-1)*m_max
  nsize_exp_inp=(n_phi-1)*n_max
! 
  allocate(exp_imt(0:nsize_exp_imt),exp_inp(-nsize_exp_inp:nsize_exp_inp))
!
  base_exp_imt=exp(dcmplx(0.d0,h_theta))
  base_exp_inp=exp(dcmplx(0.d0,h_phi))
  base_exp_inp_inv=(1.d0,0.d0)/base_exp_inp
  exp_imt(0)=(1.d0,0.d0)
  exp_inp(0)=(1.d0,0.d0)
!
  do i=1,nsize_exp_imt
    exp_imt(i)=exp_imt(i-1)*base_exp_imt
  enddo
!
  do i=1,nsize_exp_inp
    exp_inp(i)=exp_inp(i-1)*base_exp_inp
    exp_inp(-i)=exp_inp(1-i)*base_exp_inp_inv
  enddo
!
  if(allocated(sR)) deallocate(sR)
  allocate(sR(ns_s+1,ns_tp+1,ns_tp+1,ns,n_theta,n_phi))
  if(allocated(sZ)) deallocate(sZ)
  allocate(sZ(ns_s+1,ns_tp+1,ns_tp+1,ns,n_theta,n_phi))
  if(allocated(slam)) deallocate(slam)
  allocate(slam(ns_s+1,ns_tp+1,ns_tp+1,ns,n_theta,n_phi))
!
  sR(1,1,1,:,:,:)=0.d0
  sZ(1,1,1,:,:,:)=0.d0
  slam(1,1,1,:,:,:)=0.d0
!
!$omp parallel private(m, n, i_theta, i_phi, i, is, iexpt, iexpp, &
!$omp&  expphase, cosphase, sinphase, k, splcoe)
!$omp do
  do i_theta=1,n_theta
    do i_phi=1,n_phi
      do i=1,nstrm
        m=nint(axm(i))
        n=nint(axn(i))
        iexpt=m*(i_theta-1)
        iexpp=n*(i_phi-1)
        expphase=exp_imt(iexpt)*exp_inp(-iexpp)
        cosphase=dble(expphase)
        sinphase=dimag(expphase)
        do is=1,ns
          sR(1,1,1,is,i_theta,i_phi) = sR(1,1,1,is,i_theta,i_phi)      &
                                     + rmn_rho(i,is-1)*cosphase
          sZ(1,1,1,is,i_theta,i_phi) = sZ(1,1,1,is,i_theta,i_phi)      &
                                     + zmn_rho(i,is-1)*sinphase
          slam(1,1,1,is,i_theta,i_phi) = slam(1,1,1,is,i_theta,i_phi)  &
                                       + almn_rho(i,is-1)*sinphase
        enddo
      enddo
    enddo
  enddo
!$omp end do
!
!$omp barrier
!$omp single
  call new_deallocate_vmec_stuff
!$omp end single
!
! splining over $\varphi$:
!
  allocate(splcoe(0:ns_tp,n_phi))
!
!$omp do
  do is=1,ns
    do i_theta=1,n_theta
!
      splcoe(0,:)=sR(1,1,1,is,i_theta,:)
!
      call spl_per(ns_tp,n_phi,h_phi,splcoe)
!
      do k=1,ns_tp
        sR(1,1,k+1,is,i_theta,:)=splcoe(k,:)
      enddo
!
      splcoe(0,:)=sZ(1,1,1,is,i_theta,:)
!
      call spl_per(ns_tp,n_phi,h_phi,splcoe)
!
      do k=1,ns_tp
        sZ(1,1,k+1,is,i_theta,:)=splcoe(k,:)
      enddo
!
      splcoe(0,:)=slam(1,1,1,is,i_theta,:)
!
      call spl_per(ns_tp,n_phi,h_phi,splcoe)
!
      do k=1,ns_tp
        slam(1,1,k+1,is,i_theta,:)=splcoe(k,:)
      enddo
!
    enddo
  enddo
!$omp end do
!
  deallocate(splcoe)
!
! splining over $\vartheta$:
!
  allocate(splcoe(0:ns_tp,n_theta))
!
!$omp do
  do is=1,ns
    do i_phi=1,n_phi
      do isp=1,ns_tp+1
!
        splcoe(0,:)=sR(1,1,isp,is,:,i_phi)
!
        call spl_per(ns_tp,n_theta,h_theta,splcoe)
!
        do k=1,ns_tp
          sR(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
        enddo
!
        splcoe(0,:)=sZ(1,1,isp,is,:,i_phi)
!
        call spl_per(ns_tp,n_theta,h_theta,splcoe)
!
        do k=1,ns_tp
          sZ(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
        enddo
!
        splcoe(0,:)=slam(1,1,isp,is,:,i_phi)
!
        call spl_per(ns_tp,n_theta,h_theta,splcoe)
!
        do k=1,ns_tp
          slam(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
        enddo
!
      enddo
    enddo
  enddo
!$omp end do
!
  deallocate(splcoe)
!
! splining over $s$:
!
  allocate(splcoe(0:ns_s,ns))
!
!$omp do
  do i_theta=1,n_theta
    do i_phi=1,n_phi
      do ist=1,ns_tp+1
        do isp=1,ns_tp+1
!
          splcoe(0,:)=sR(1,ist,isp,:,i_theta,i_phi)
!
          call spl_reg(ns_s,ns,hs,splcoe)
!
          do k=1,ns_s
            sR(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
          enddo
!
          splcoe(0,:)=sZ(1,ist,isp,:,i_theta,i_phi)
!
          call spl_reg(ns_s,ns,hs,splcoe)
!
          do k=1,ns_s
            sZ(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
          enddo
!
          splcoe(0,:)=slam(1,ist,isp,:,i_theta,i_phi)
!
          call spl_reg(ns_s,ns,hs,splcoe)
!
          do k=1,ns_s
            slam(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
          enddo
!
        enddo
      enddo
    enddo
  enddo
!$omp end do
!
  deallocate(splcoe)
!$omp end parallel
  deallocate(exp_imt,exp_inp)
!
  end subroutine spline_vmec_data
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine deallocate_vmec_spline(mode)
!
  use new_vmec_stuff_mod
!
  implicit none
!
  integer :: mode
!
  if(mode.eq.0) then
    deallocate(sR,sZ,slam)
  elseif(mode.eq.1) then
    deallocate(sR,sZ)
  elseif(mode.eq.2) then
    deallocate(slam)
  else
    print *,'deallocate_vmec_spline: unknown mode'
  endif
!
  end subroutine deallocate_vmec_spline
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine splint_vmec_data(s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                              R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
  use new_vmec_stuff_mod,   only : n_theta,n_phi,h_theta,h_phi,sR,sZ,slam,nper,ns_A,ns_s,ns_tp
  use vector_potentail_mod, only : ns,hs,torflux,sA_phi
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  integer :: is,i_theta,i_phi,k
  double precision :: ds,dtheta,dphi,rho_tor
  double precision :: s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                      R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
  integer, parameter :: ns_max=6
!
  integer :: nstp
!
  double precision, dimension(ns_max)        :: sp_R,sp_Z,sp_lam
  double precision, dimension(ns_max)        :: dsp_R_ds,dsp_Z_ds,dsp_lam_ds
  double precision, dimension(ns_max)        :: dsp_R_dt,dsp_Z_dt,dsp_lam_dt
  double precision, dimension(ns_max,ns_max) :: stp_R,stp_Z,stp_lam
  double precision, dimension(ns_max,ns_max) :: dstp_R_ds,dstp_Z_ds,dstp_lam_ds
!
  nstp=ns_tp+1
!
  A_theta=torflux*s
  dA_theta_ds=torflux
!
  ds=s/hs
  is=max(0,min(ns-1,int(ds)))
  ds=(ds-dble(is))*hs
  is=is+1
!
  dtheta=modulo(theta,twopi)/h_theta
  i_theta=max(0,min(n_theta-1,int(dtheta)))
  dtheta=(dtheta-dble(i_theta))*h_theta
  i_theta=i_theta+1
!
  dphi=modulo(varphi,twopi/dble(nper))/h_phi
  i_phi=max(0,min(n_phi-1,int(dphi)))
  dphi=(dphi-dble(i_phi))*h_phi
  i_phi=i_phi+1
!
! Begin interpolation over $s$
!
! Vector potential $A_\varphi$ and its derivative:
  A_phi=sA_phi(ns_A+1,is)
  dA_phi_ds=0.d0
!
  do k=ns_A,1,-1
    A_phi=sA_phi(k,is)+ds*A_phi
    dA_phi_ds=sA_phi(k+1,is)*dble(k)+ds*dA_phi_ds
  enddo
!
! R, Z and $\lambda$ and their derivatives over $s$:
!
  rho_tor=sqrt(s)
  ds=rho_tor/hs
  is=max(0,min(ns-2,int(ds)))
  ds=(ds-dble(is))*hs
  is=is+1
!
  stp_R(1:nstp,1:nstp)=sR(ns_s+1,:,:,is,i_theta,i_phi)
  dstp_R_ds(1:nstp,1:nstp)=0.d0
  stp_Z(1:nstp,1:nstp)=sZ(ns_s+1,:,:,is,i_theta,i_phi)
  dstp_Z_ds(1:nstp,1:nstp)=0.d0
  stp_lam(1:nstp,1:nstp)=slam(ns_s+1,:,:,is,i_theta,i_phi)
  dstp_lam_ds(1:nstp,1:nstp)=0.d0
!
  do k=ns_s,1,-1
    stp_R(1:nstp,1:nstp)=sR(k,:,:,is,i_theta,i_phi)+ds*stp_R(1:nstp,1:nstp)
    dstp_R_ds(1:nstp,1:nstp)=sR(k+1,:,:,is,i_theta,i_phi)*dble(k)+ds*dstp_R_ds(1:nstp,1:nstp)
    stp_Z(1:nstp,1:nstp)=sZ(k,:,:,is,i_theta,i_phi)+ds*stp_Z(1:nstp,1:nstp)
    dstp_Z_ds(1:nstp,1:nstp)=sZ(k+1,:,:,is,i_theta,i_phi)*dble(k)+ds*dstp_Z_ds(1:nstp,1:nstp)
    stp_lam(1:nstp,1:nstp)=slam(k,:,:,is,i_theta,i_phi)+ds*stp_lam(1:nstp,1:nstp)
    dstp_lam_ds(1:nstp,1:nstp)=slam(k+1,:,:,is,i_theta,i_phi)*dble(k)+ds*dstp_lam_ds(1:nstp,1:nstp)
  enddo
!
! End interpolation over $s$
!----------------------------
!
  aiota=-dA_phi_ds/dA_theta_ds
!
! Begin interpolation over $\theta$
!
  sp_R(1:nstp)=stp_R(nstp,1:nstp)
  dsp_R_ds(1:nstp)=dstp_R_ds(nstp,1:nstp)
  dsp_R_dt(1:nstp)=0.d0
  sp_Z(1:nstp)=stp_Z(nstp,1:nstp)
  dsp_Z_ds(1:nstp)=dstp_Z_ds(nstp,1:nstp)
  dsp_Z_dt(1:nstp)=0.d0
  sp_lam(1:nstp)=stp_lam(nstp,1:nstp)
  dsp_lam_ds(1:nstp)=dstp_lam_ds(nstp,1:nstp)
  dsp_lam_dt(1:nstp)=0.d0
!
  do k=ns_tp,1,-1
    sp_R(1:nstp)=stp_R(k,1:nstp)+dtheta*sp_R(1:nstp)
    dsp_R_ds(1:nstp)=dstp_R_ds(k,1:nstp)+dtheta*dsp_R_ds(1:nstp)
    dsp_R_dt(1:nstp)=stp_R(k+1,1:nstp)*dble(k)+dtheta*dsp_R_dt(1:nstp)
!
    sp_Z(1:nstp)=stp_Z(k,1:nstp)+dtheta*sp_Z(1:nstp)
    dsp_Z_ds(1:nstp)=dstp_Z_ds(k,1:nstp)+dtheta*dsp_Z_ds(1:nstp)
    dsp_Z_dt(1:nstp)=stp_Z(k+1,1:nstp)*dble(k)+dtheta*dsp_Z_dt(1:nstp)
!
    sp_lam(1:nstp)=stp_lam(k,1:nstp)+dtheta*sp_lam(1:nstp)
    dsp_lam_ds(1:nstp)=dstp_lam_ds(k,1:nstp)+dtheta*dsp_lam_ds(1:nstp)
    dsp_lam_dt(1:nstp)=stp_lam(k+1,1:nstp)*dble(k)+dtheta*dsp_lam_dt(1:nstp)
  enddo
!
! End interpolation over $\theta$
!--------------------------------
!
! Begin interpolation over $\varphi$
!
  R=sp_R(nstp)
  dR_ds=dsp_R_ds(nstp)
  dR_dt=dsp_R_dt(nstp)
  dR_dp=0.d0
  Z=sp_Z(nstp)
  dZ_ds=dsp_Z_ds(nstp)
  dZ_dt=dsp_Z_dt(nstp)
  dZ_dp=0.d0
  alam=sp_lam(nstp)
  dl_ds=dsp_lam_ds(nstp)
  dl_dt=dsp_lam_dt(nstp)
  dl_dp=0.d0
!
  do k=ns_tp,1,-1
    R=sp_R(k)+dphi*R
    dR_ds=dsp_R_ds(k)+dphi*dR_ds
    dR_dt=dsp_R_dt(k)+dphi*dR_dt
    dR_dp=sp_R(k+1)*dble(k)+dphi*dR_dp
!
    Z=sp_Z(k)+dphi*Z
    dZ_ds=dsp_Z_ds(k)+dphi*dZ_ds
    dZ_dt=dsp_Z_dt(k)+dphi*dZ_dt
    dZ_dp=sp_Z(k+1)*dble(k)+dphi*dZ_dp
!
    alam=sp_lam(k)+dphi*alam
    dl_ds=dsp_lam_ds(k)+dphi*dl_ds
    dl_dt=dsp_lam_dt(k)+dphi*dl_dt
    dl_dp=sp_lam(k+1)*dble(k)+dphi*dl_dp
  enddo
!
! End interpolation over $\varphi$
!
  dR_ds=0.5d0*dR_ds/rho_tor
  dZ_ds=0.5d0*dZ_ds/rho_tor
  dl_ds=0.5d0*dl_ds/rho_tor
!
  end subroutine splint_vmec_data
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                        sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                        Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
                        
!
  implicit none
!
  double precision :: s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,        &
                      R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
  double precision :: Bctrvr_vartheta,Bctrvr_varphi,Bcovar_vartheta,Bcovar_varphi,sqg
  double precision :: cjac,sqgV,Bcovar_r
  double precision, dimension(3,3) :: cmat,gV,g
!
!
  call splint_vmec_data(s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,      &
                        R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
  gV(1,1)=dR_ds**2+dZ_ds**2
  gV(1,2)=dR_ds*dR_dt+dZ_ds*dZ_dt
  gV(1,3)=dR_ds*dR_dp+dZ_ds*dZ_dp
  gV(2,1)=gV(1,2)
  gV(2,2)=dR_dt**2+dZ_dt**2
  gV(2,3)=dR_dt*dR_dp+dZ_dt*dZ_dp
  gV(3,1)=gV(1,3)
  gV(3,2)=gV(2,3)
  gV(3,3)=R**2+dR_dp**2+dZ_dp**2
  sqgV=R*(dR_dt*dZ_ds-dR_ds*dZ_dt)
!
  cjac=1.d0/(1.d0+dl_dt)
  sqg=sqgV*cjac
  Bctrvr_vartheta=-dA_phi_ds/sqg
  Bctrvr_varphi=dA_theta_ds/sqg
!
  cmat(1,2:3)=0.d0
  cmat(3,1:2)=0.d0
  cmat(1,1)=1.d0
  cmat(3,3)=1.d0
  cmat(2,1)=-dl_ds*cjac
  cmat(2,2)=cjac
  cmat(2,3)=-dl_dp*cjac
!
  g=matmul(transpose(cmat),matmul(gV,cmat))
!
  Bcovar_r=g(1,2)*Bctrvr_vartheta+g(1,3)*Bctrvr_varphi
  Bcovar_vartheta=g(2,2)*Bctrvr_vartheta+g(2,3)*Bctrvr_varphi
  Bcovar_varphi=g(3,2)*Bctrvr_vartheta+g(3,3)*Bctrvr_varphi
!
  end subroutine vmec_field
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine splint_iota(s,aiota,daiota_ds)
!
  use vector_potentail_mod, only : ns,hs,torflux,sA_phi
  use new_vmec_stuff_mod,   only : ns_A
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  integer :: is,i_theta,i_phi,k
  double precision :: ds,dtheta,dphi
  double precision :: s,dA_phi_ds,dA_theta_ds,d2A_phi_ds2,aiota,daiota_ds
!
  integer, parameter :: ns_max=6
!
  integer :: nstp
!
  dA_theta_ds=torflux
!
  ds=s/hs
  is=max(0,min(ns-1,int(ds)))
  ds=(ds-dble(is))*hs
  is=is+1
!
  dA_phi_ds=0.d0
!
  do k=ns_A,1,-1
    dA_phi_ds=sA_phi(k+1,is)*dble(k)+ds*dA_phi_ds
  enddo
!
  d2A_phi_ds2=0.d0
!
  do k=ns_A,2,-1
    d2A_phi_ds2=sA_phi(k+1,is)*dble(k)*dble(k-1)+ds*d2A_phi_ds2
  enddo
!
  aiota=-dA_phi_ds/dA_theta_ds
  daiota_ds=-d2A_phi_ds2/dA_theta_ds
!
  end subroutine splint_iota
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine splint_lambda(s,theta,varphi,alam,dl_dt)
!
  use new_vmec_stuff_mod,   only : n_theta,n_phi,h_theta,h_phi,slam,nper,ns_s,ns_tp
  use vector_potentail_mod, only : ns,hs
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  integer :: is,i_theta,i_phi,k
  double precision :: ds,dtheta,dphi
  double precision :: s,theta,varphi,alam,dl_dt
!
  integer, parameter :: ns_max=6
!
  integer :: nstp
!
  double precision, dimension(ns_max)        :: sp_lam
  double precision, dimension(ns_max)        :: dsp_lam_dt
  double precision, dimension(ns_max,ns_max) :: stp_lam
!
  nstp=ns_tp+1
!
  ds=sqrt(s)/hs
  is=max(0,min(ns-1,int(ds)))
  ds=(ds-dble(is))*hs
  is=is+1
!
  dtheta=modulo(theta,twopi)/h_theta
  i_theta=max(0,min(n_theta-1,int(dtheta)))
  dtheta=(dtheta-dble(i_theta))*h_theta
  i_theta=i_theta+1
!
  dphi=modulo(varphi,twopi/dble(nper))/h_phi
  i_phi=max(0,min(n_phi-1,int(dphi)))
  dphi=(dphi-dble(i_phi))*h_phi
  i_phi=i_phi+1
!
! Begin interpolation over $s$
!
  stp_lam(1:nstp,1:nstp)=slam(ns_s+1,:,:,is,i_theta,i_phi)
!
  do k=ns_s,1,-1
    stp_lam(1:nstp,1:nstp)=slam(k,:,:,is,i_theta,i_phi)+ds*stp_lam(1:nstp,1:nstp)
  enddo
!
! End interpolation over $s$
!----------------------------
!
! Begin interpolation over $\theta$
!
  sp_lam(1:nstp)=stp_lam(nstp,1:nstp)
  dsp_lam_dt(1:nstp)=0.d0
!
  do k=ns_tp,1,-1
    sp_lam(1:nstp)=stp_lam(k,1:nstp)+dtheta*sp_lam(1:nstp)
    dsp_lam_dt(1:nstp)=stp_lam(k+1,1:nstp)*dble(k)+dtheta*dsp_lam_dt(1:nstp)
  enddo
!
! End interpolation over $\theta$
!--------------------------------
!
! Begin interpolation over $\varphi$
!
  alam=sp_lam(nstp)
  dl_dt=dsp_lam_dt(nstp)
!
  do k=ns_tp,1,-1
    alam=sp_lam(k)+dphi*alam
    dl_dt=dsp_lam_dt(k)+dphi*dl_dt
  enddo
!
  end subroutine splint_lambda
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine s_to_rho_healaxis(m,ns,nrho,nheal,arr_in,arr_out)
!
  use new_vmec_stuff_mod, only : ns_s
!
  implicit none
!
  integer :: m,ns,nrho,nheal,irho,is,k,nhe
!
  double precision :: hs,hrho,s,ds,rho,a,b,c
!
  double precision, dimension(ns)   :: arr_in
  double precision, dimension(nrho) :: arr_out
  double precision, dimension(:,:), allocatable :: splcoe
!
!do is=1,ns
!write(2001,*) is,arr_in(is)
!enddo
!close(2001)
!print *,m
!
  hs=1.d0/dble(ns-1)
  hrho=1.d0/dble(nrho-1)
!
  nhe=max(1,nheal)+1
!
  do is=nhe,ns
    if(m.gt.0) then
      rho=sqrt(hs*dble(is-1))
      arr_out(is)=arr_in(is)/rho**m
    else
      arr_out(is)=arr_in(is)
    endif
  enddo
!
  a=arr_out(nhe)
  b=0.5d0*(4.d0*arr_out(nhe+1)-3.d0*arr_out(nhe)-arr_out(nhe+2))
  c=0.5d0*(arr_out(nhe)+arr_out(nhe+2)-2.d0*arr_out(nhe+1))
!
  do is=1,nhe-1
    arr_out(is)=a+b*dble(is-nhe)+c*dble(is-nhe)**2
  enddo
!do is=1,ns
!write(2002,*) is,arr_out(is)
!enddo
!close(2002)
!
  allocate(splcoe(0:ns_s,ns))
!
  splcoe(0,:)=arr_out
!
  call spl_reg(ns_s,ns,hs,splcoe)
!
  do irho=1,nrho
    rho=hrho*dble(irho-1)
    s=rho**2
!
    ds=s/hs
    is=max(0,min(ns-1,int(ds)))
    ds=(ds-dble(is))*hs
    is=is+1
!
    arr_out(irho)=splcoe(ns_s,is)
!
    do k=ns_s-1,0,-1
      arr_out(irho)=splcoe(k,is)+ds*arr_out(irho)
    enddo
!
    if(m.gt.0) arr_out(irho)=arr_out(irho)*rho**m
  enddo
!
  deallocate(splcoe)
!
  end subroutine s_to_rho_healaxis
