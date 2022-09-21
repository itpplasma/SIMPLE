!
  module exchange_get_cancoord_mod
    implicit none

    logical :: onlytheta
    double precision :: vartheta_c,varphi_c,sqg,aiota,Bcovar_vartheta,&
      Bcovar_varphi,A_theta,A_phi,theta,Bctrvr_vartheta,Bctrvr_varphi
!$omp threadprivate(onlytheta, vartheta_c, varphi_c, sqg, aiota)
!$omp threadprivate(Bcovar_vartheta,Bcovar_varphi,A_theta,A_phi,theta,Bctrvr_vartheta,Bctrvr_varphi)
  end module exchange_get_cancoord_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_canonical_coordinates
!
  use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,           &
                                        hs_c,h_theta_c,h_phi_c,           &
                                        ns_s_c,ns_tp_c,                   &
                                        nh_stencil,G_c,sqg_c,             &
                                        B_vartheta_c,B_varphi_c
  use vector_potentail_mod, only : ns,hs
  use exchange_get_cancoord_mod, only : vartheta_c,varphi_c,sqg,aiota,    &
                                        Bcovar_vartheta,Bcovar_varphi,    &
                                        onlytheta
  use new_vmec_stuff_mod, only : n_theta,n_phi,h_theta,h_phi,ns_s,ns_tp
!
  implicit none
!
  logical :: fullset
  double precision, parameter :: relerr=1d-10
  integer :: i_theta,i_phi,i_sten,ndim,is_beg
  integer,          dimension(:),     allocatable :: ipoi_t,ipoi_p
  double precision, dimension(:),     allocatable :: y,dy
  double precision :: dstencil_theta(-nh_stencil:nh_stencil), &
                      dstencil_phi(-nh_stencil:nh_stencil)
!
  double precision :: r,r1,r2,G_beg,dG_c_dt,dG_c_dp
  integer :: is
  integer :: i_ctr ! for nice counting in parallel
!
  external rhs_cancoord
!
!
  ns_c=ns
  n_theta_c=n_theta
  n_phi_c=n_phi
  h_theta_c=h_theta
  h_phi_c=h_phi
  hs_c=hs
!
  if(nh_stencil.eq.1) then
    dstencil_theta(-1)=-0.5d0
    dstencil_theta(0)=0.0d0
    dstencil_theta(1)=0.5d0
  elseif(nh_stencil.eq.2) then
    dstencil_theta(-2)=1.d0/12.d0
    dstencil_theta(-1)=-2.d0/3.d0
    dstencil_theta(0)=0.0d0
    dstencil_theta(1)=2.d0/3.d0
    dstencil_theta(2)=-1.d0/12.d0
  elseif(nh_stencil.eq.3) then
    dstencil_theta(-3)=-1.d0/60.d0
    dstencil_theta(-2)=0.15d0
    dstencil_theta(-1)=-0.75d0
    dstencil_theta(0)=0.0d0
    dstencil_theta(1)=0.75d0
    dstencil_theta(2)=-0.15d0
    dstencil_theta(3)=1.d0/60.d0
  endif
!
  dstencil_phi=dstencil_theta
  dstencil_theta=dstencil_theta/h_theta_c
  dstencil_phi=dstencil_phi/h_phi_c
!
  allocate(ipoi_t(1-nh_stencil:n_theta_c+nh_stencil))
  allocate(ipoi_p(1-nh_stencil:n_phi_c+nh_stencil))
!
  do i_theta=1,n_theta_c
    ipoi_t(i_theta)=i_theta
  enddo
!
  do i_phi=1,n_phi_c
    ipoi_p(i_phi)=i_phi
  enddo
!
  do i_sten=1,nh_stencil
    ipoi_t(1-i_sten)=ipoi_t(n_theta-i_sten)
    ipoi_t(n_theta_c+i_sten)=ipoi_t(1+i_sten)
    ipoi_p(1-i_sten)=ipoi_p(n_phi_c-i_sten)
    ipoi_p(n_phi_c+i_sten)=ipoi_p(1+i_sten)
  enddo
!
  allocate(G_c(ns_c,n_theta_c,n_phi_c))
  allocate(sqg_c(ns_c,n_theta_c,n_phi_c))
  allocate(B_vartheta_c(ns_c,n_theta_c,n_phi_c))
  allocate(B_varphi_c(ns_c,n_theta_c,n_phi_c))
!
  onlytheta=.false.
  ndim=1
!  is_beg=ns_c/2 !<=OLD
  is_beg=1       !<=NEW
!  G_beg=1.d-5 !<=OLD
  G_beg=1.d-8   !<=NEW

  i_ctr=0
!$omp parallel private(y, dy, i_theta, i_phi, is, r1, r2, r, dG_c_dt, dG_c_dp)
!$omp critical
  allocate(y(ndim),dy(ndim))
!$omp end critical
!
!$omp do
  do i_theta=1,n_theta_c
!$omp critical
    i_ctr = i_ctr + 1
    print *,'integrate ODE: ',i_ctr,' of ',n_theta_c
!$omp end critical
    vartheta_c=h_theta_c*dble(i_theta-1)
    do i_phi=1,n_phi_c
      varphi_c=h_phi_c*dble(i_phi-1)
!
      G_c(is_beg,i_theta,i_phi)=G_beg
      y(1)=G_beg
!
!      do is=is_beg-1,1,-1
      do is=is_beg-1,2,-1
!        r1=(hs_c*dble(is))**2    !<=OLD
!        r2=(hs_c*dble(is-1))**2  !<=OLD
        r1=hs_c*dble(is)          !<=NEW
        r2=hs_c*dble(is-1)        !<=NEW
!        if(is.eq.1) r2=hs_c*1d-5
!        if(is.eq.1) r2=(hs_c*1d-5)**2
!
        call odeint_allroutines(y,ndim,r1,r2,relerr,rhs_cancoord)
!
        G_c(is,i_theta,i_phi)=y(1)
      enddo
!
      y(1)=G_beg
!
      do is=is_beg+1,ns_c
!        r1=(hs_c*dble(is-2))**2  !<=OLD
!        r2=(hs_c*dble(is-1))**2  !<=OLD
        r1=hs_c*dble(is-2)        !<=NEW
        r2=hs_c*dble(is-1)        !<=NEW
        if(is.eq.2) r1=1.d-8
!
        call odeint_allroutines(y,ndim,r1,r2,relerr,rhs_cancoord)
!
        G_c(is,i_theta,i_phi)=y(1)
      enddo
    enddo
  enddo
!$omp end do
!
i_ctr=0
!$omp barrier
!$omp do
  do i_theta=1,n_theta_c
!$omp critical
    i_ctr = i_ctr + 1
    print *,'compute components: ',i_ctr,' of ',n_theta_c
!$omp end critical
    vartheta_c=h_theta_c*dble(i_theta-1)
    do i_phi=1,n_phi_c
      varphi_c=h_phi_c*dble(i_phi-1)
!      do is=1,ns_c
      do is=2,ns_c
!        r=(hs_c*dble(is-1))**2  !<=OLD
        r=hs_c*dble(is-1)        !<=NEW
        y(1)=G_c(is,i_theta,i_phi)
!
        call rhs_cancoord(r,y,dy)
!
        dG_c_dt=sum(dstencil_theta*G_c(is,ipoi_t(i_theta-nh_stencil:i_theta+nh_stencil),i_phi))
        dG_c_dp=sum(dstencil_phi*G_c(is,i_theta,ipoi_p(i_phi-nh_stencil:i_phi+nh_stencil)))
        sqg_c(is,i_theta,i_phi)=sqg*(1.d0+aiota*dG_c_dt+dG_c_dp)
        B_vartheta_c(is,i_theta,i_phi)=Bcovar_vartheta+(aiota*Bcovar_vartheta+Bcovar_varphi)*dG_c_dt
        B_varphi_c(is,i_theta,i_phi)=Bcovar_varphi+(aiota*Bcovar_vartheta+Bcovar_varphi)*dG_c_dp
      enddo
!First point is=1 (on axis) is bad, extrapolate with parabola:
      sqg_c(1,i_theta,i_phi) = 3.d0*(sqg_c(2,i_theta,i_phi)-sqg_c(3,i_theta,i_phi))                      &
                             + sqg_c(4,i_theta,i_phi)
!      B_vartheta_c(1,i_theta,i_phi) = 3.d0*(B_vartheta_c(2,i_theta,i_phi)-B_vartheta_c(3,i_theta,i_phi)) & !<=OLD
!                                    + B_vartheta_c(4,i_theta,i_phi)                                        !<=OLD
      B_vartheta_c(1,i_theta,i_phi) = 0.d0                                                                  !<=NEW
      B_varphi_c(1,i_theta,i_phi) = 3.d0*(B_varphi_c(2,i_theta,i_phi)-B_varphi_c(3,i_theta,i_phi))       &
                                  + B_varphi_c(4,i_theta,i_phi)
    enddo
  enddo
!$omp end do
!$omp critical
deallocate(y,dy)
!$omp end critical
!$omp end parallel
!
  ns_s_c=ns_s
  ns_tp_c=ns_tp
  fullset=.true.
!
!  call deallocate_vmec_spline(1)
!
  onlytheta=.true.
!
!do is=1,ns_c
!write(400,*) G_c(is,:,10)
!write(401,*) B_vartheta_c(is,:,10)
!write(402,*) B_varphi_c(is,:,10)
!enddo
!is=90
!do i_phi=1,n_phi_c
!write(500,*) G_c(is,:,i_phi)
!write(501,*) B_vartheta_c(is,:,i_phi)
!write(502,*) B_varphi_c(is,:,i_phi)
!enddo
!stop
  call spline_can_coord(fullset)
  deallocate(ipoi_t,ipoi_p,sqg_c,B_vartheta_c,B_varphi_c,G_c)
!
  end subroutine get_canonical_coordinates
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs_cancoord(r,y,dy)
!
  use exchange_get_cancoord_mod, only : vartheta_c,varphi_c,sqg,aiota,Bcovar_vartheta,Bcovar_varphi, &
                                        theta,onlytheta
!
  implicit none
!
  double precision, parameter :: epserr=1.d-14
  integer :: iter
  double precision :: s,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,                 &
                      alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,Bcovar_r
!
  double precision :: r,vartheta,daiota_ds,deltheta
  double precision, dimension(1) :: y,dy
!
!  s=r   !<=OLD
  s=r**2 !<=NEW
!
  call splint_iota(s,aiota,daiota_ds)
!
  vartheta=vartheta_c+aiota*y(1)
  varphi=varphi_c+y(1)
!
! Begin Newton iteration to find VMEC theta
!
  theta = vartheta
!
  do iter=1,100
!
    call splint_lambda(s,theta,varphi,alam,dl_dt)
!
    deltheta = (vartheta-theta-alam)/(1.d0+dl_dt)
    theta = theta + deltheta
    if(abs(deltheta).lt.epserr) exit
  enddo
!
! End Newton iteration to find VMEC theta
!
  if(onlytheta) return
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  dy(1)=-(Bcovar_r+daiota_ds*Bcovar_vartheta*y(1))/(aiota*Bcovar_vartheta+Bcovar_varphi)
  dy(1)=2.d0*r*dy(1)  !<=NEW
!
  end subroutine rhs_cancoord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spline_can_coord(fullset)
!
  use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,hs_c,h_theta_c,h_phi_c,    &
                                        ns_s_c,ns_tp_c,G_c,sqg_c,B_vartheta_c,B_varphi_c, &
                                        s_sqg_Bt_Bp,s_G_c,ns_max,derf1,derf2,derf3
!
  implicit none
!
  logical :: fullset
  integer :: k,is,i_theta,i_phi,i_qua
  integer :: iss,ist,isp
  double precision, dimension(:,:), allocatable :: splcoe
!
  if (.not. allocated(s_sqg_Bt_Bp)) &
    allocate(s_sqg_Bt_Bp(3,ns_s_c+1,ns_tp_c+1,ns_tp_c+1,ns_c,n_theta_c,n_phi_c))
  if(fullset .and. (.not. allocated(s_G_c))) &
    allocate(s_G_c(ns_s_c+1,ns_tp_c+1,ns_tp_c+1,ns_c,n_theta_c,n_phi_c))
!
  s_sqg_Bt_Bp(1,1,1,1,:,:,:)=sqg_c
  s_sqg_Bt_Bp(2,1,1,1,:,:,:)=B_vartheta_c
  s_sqg_Bt_Bp(3,1,1,1,:,:,:)=B_varphi_c
  if(fullset) s_G_c(1,1,1,:,:,:)=G_c
!
! splining over $\varphi$:
!
  allocate(splcoe(0:ns_tp_c,n_phi_c))
!
  do is=1,ns_c
    do i_theta=1,n_theta_c
      do i_qua=1,3
!
        splcoe(0,:)=s_sqg_Bt_Bp(i_qua,1,1,1,is,i_theta,:)
!
        call spl_per(ns_tp_c,n_phi_c,h_phi_c,splcoe)
!
        do k=1,ns_tp_c
          s_sqg_Bt_Bp(i_qua,1,1,k+1,is,i_theta,:)=splcoe(k,:)
        enddo
!
      enddo
!
      if(fullset) then
!
        splcoe(0,:)=s_G_c(1,1,1,is,i_theta,:)
!
        call spl_per(ns_tp_c,n_phi_c,h_phi_c,splcoe)
!
        do k=1,ns_tp_c
          s_G_c(1,1,k+1,is,i_theta,:)=splcoe(k,:)
        enddo
!
      endif
    enddo
  enddo
!
  deallocate(splcoe)
!
! splining over $\vartheta$:
!
  allocate(splcoe(0:ns_tp_c,n_theta_c))
!
  do is=1,ns_c
    do i_phi=1,n_phi_c
      do isp=1,ns_tp_c+1
        do i_qua=1,3
!
          splcoe(0,:)=s_sqg_Bt_Bp(i_qua,1,1,isp,is,:,i_phi)
!
          call spl_per(ns_tp_c,n_theta_c,h_theta_c,splcoe)
!
          do k=1,ns_tp_c
            s_sqg_Bt_Bp(i_qua,1,k+1,isp,is,:,i_phi)=splcoe(k,:)
          enddo
!
        enddo
!
        if(fullset) then
!
          splcoe(0,:)=s_G_c(1,1,isp,is,:,i_phi)
!
          call spl_per(ns_tp_c,n_theta_c,h_theta_c,splcoe)
!
          do k=1,ns_tp_c
            s_G_c(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
          enddo
!
        endif
!
      enddo
    enddo
  enddo
!
  deallocate(splcoe)
!
! splining over $s$:
!
  allocate(splcoe(0:ns_s_c,ns_c))
!
  do i_theta=1,n_theta_c
    do i_phi=1,n_phi_c
      do ist=1,ns_tp_c+1
        do isp=1,ns_tp_c+1
          do i_qua=1,3
!
            splcoe(0,:)=s_sqg_Bt_Bp(i_qua,1,ist,isp,:,i_theta,i_phi)
!
            call spl_reg(ns_s_c,ns_c,hs_c,splcoe)
!
            do k=1,ns_s_c
              s_sqg_Bt_Bp(i_qua,k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
            enddo
!
          enddo
!
          if(fullset) then
!
            splcoe(0,:)=s_G_c(1,ist,isp,:,i_theta,i_phi)
!
            call spl_reg(ns_s_c,ns_c,hs_c,splcoe)
!
            do k=1,ns_s_c
              s_G_c(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
            enddo
!
          endif
!
        enddo
      enddo
    enddo
  enddo
!
  deallocate(splcoe)
!
  do k=1,ns_max
    derf1(k)=dble(k-1)
    derf2(k)=dble((k-1)*(k-2))
    derf3(k)=dble((k-1)*(k-2)*(k-3))
  enddo
!
  end subroutine spline_can_coord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine splint_can_coord(fullset,mode_secders,r,vartheta_c,varphi_c,                      &
                              A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                              sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                              B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                              B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                              d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                              d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                              d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp,G_c)
!
  use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,hs_c,h_theta_c,h_phi_c,    &
                                        ns_s_c,ns_tp_c,ns_max,n_qua,derf1,derf2,derf3,    &
                                        s_sqg_Bt_Bp,s_G_c
  use vector_potentail_mod, only : ns,hs,torflux,sA_phi
  use new_vmec_stuff_mod,   only : nper,ns_A
  use chamb_mod,            only : rnegflag
use diag_mod, only : icounter
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  logical :: fullset
!
  integer :: mode_secders,nstp,ns_A_p1,ns_s_p1
  integer :: k,is,i_theta,i_phi
  integer :: iss,ist,isp
!
  double precision :: r,vartheta_c,varphi_c,                                           &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c,     &
                      d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                      d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                      d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp
  double precision :: s,ds,dtheta,dphi,rho_tor,drhods,drhods2,d2rhods2m
!
  double precision, dimension(ns_max)              :: sp_G
  double precision, dimension(ns_max,ns_max)       :: stp_G
!
  double precision, dimension(n_qua)               :: qua,dqua_dr,dqua_dt,dqua_dp
  double precision, dimension(n_qua)               :: d2qua_dr2,d2qua_drdt,d2qua_drdp,d2qua_dt2,d2qua_dtdp,d2qua_dp2
  double precision, dimension(n_qua,ns_max)        :: sp_all,dsp_all_ds,dsp_all_dt
  double precision, dimension(n_qua,ns_max)        :: d2sp_all_ds2,d2sp_all_dsdt,d2sp_all_dt2
  double precision, dimension(n_qua,ns_max,ns_max) :: stp_all,dstp_all_ds,d2stp_all_ds2
!$omp atomic
icounter=icounter+1
  if(r.le.0.d0) then
    rnegflag=.true.
    r=abs(r)
  endif
!
  A_theta=torflux*r
  dA_theta_dr=torflux
!
  dtheta=modulo(vartheta_c,twopi)/h_theta_c
  i_theta=max(0,min(n_theta_c-1,int(dtheta)))
  dtheta=(dtheta-dble(i_theta))*h_theta_c
  i_theta=i_theta+1
!
  dphi=modulo(varphi_c,twopi/dble(nper))/h_phi_c
  i_phi=max(0,min(n_phi_c-1,int(dphi)))
  dphi=(dphi-dble(i_phi))*h_phi_c
  i_phi=i_phi+1
!
! Begin interpolation of vector potentials over $s$
!
  ds=r/hs
  is=max(0,min(ns-1,int(ds)))
  ds=(ds-dble(is))*hs
  is=is+1
!
  ns_A_p1=ns_A+1
  A_phi=sA_phi(ns_A_p1,is)
  dA_phi_dr=A_phi*derf1(ns_A_p1)
  d2A_phi_dr2=A_phi*derf2(ns_A_p1)
!
  do k=ns_A,3,-1
    A_phi=sA_phi(k,is)+ds*A_phi
    dA_phi_dr=sA_phi(k,is)*derf1(k)+ds*dA_phi_dr
    d2A_phi_dr2=sA_phi(k,is)*derf2(k)+ds*d2A_phi_dr2
  enddo
!
  A_phi=sA_phi(1,is)+ds*(sA_phi(2,is)+ds*A_phi)
  dA_phi_dr=sA_phi(2,is)+ds*dA_phi_dr
!
  if(mode_secders.gt.0) then
    d3A_phi_dr3=sA_phi(ns_A_p1,is)*derf3(ns_A_p1)
!
    do k=ns_A,4,-1
      d3A_phi_dr3=sA_phi(k,is)*derf3(k)+ds*d3A_phi_dr3
    enddo
!
  endif
!
! End interpolation of vector potentials over $s$
!
!--------------------------------
!-------------------------------
!
  rho_tor=sqrt(r)
  !hs_c=hs !added by Johanna in alalogy to get_canonical_coordinites to make test_orbits_vmec working
  ds=rho_tor/hs_c
  is=max(0,min(ns_c-1,int(ds)))
  ds=(ds-dble(is))*hs_c
  is=is+1
!
  nstp=ns_tp_c+1
!
  if(fullset) then
!
! Begin interpolation of G over $s$
!
    stp_G(1:nstp,1:nstp)=s_G_c(ns_s_c+1,:,:,is,i_theta,i_phi)
!
    do k=ns_s_c,1,-1
      stp_G(1:nstp,1:nstp)=s_G_c(k,:,:,is,i_theta,i_phi)+ds*stp_G(1:nstp,1:nstp)
    enddo
!
! End interpolation of G over $s$
!----------------------------
! Begin interpolation of G over $\theta$
!
    sp_G(1:nstp)=stp_G(nstp,1:nstp)
!
    do k=ns_tp_c,1,-1
      sp_G(1:nstp)=stp_G(k,1:nstp)+dtheta*sp_G(1:nstp)
    enddo
!
! End interpolation of G over $\theta$
!--------------------------------
! Begin interpolation of G over $\varphi$
!
    G_c=sp_G(nstp)
!
    do k=ns_tp_c,1,-1
      G_c=sp_G(k)+dphi*G_c
    enddo
!
! End interpolation of G over $\varphi$
!
  endif
!
!--------------------------------
  if(mode_secders.eq.2) then
!--------------------------------
!
! Begin interpolation of all over $s$
!
    ns_s_p1=ns_s_c+1
    stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,ns_s_p1,:,:,is,i_theta,i_phi)
    dstp_all_ds(:,1:nstp,1:nstp)=stp_all(:,1:nstp,1:nstp)*derf1(ns_s_p1)
    d2stp_all_ds2(:,1:nstp,1:nstp)=stp_all(:,1:nstp,1:nstp)*derf2(ns_s_p1)
!
    do k=ns_s_c,3,-1
      stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,k,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp)
      dstp_all_ds(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,k,:,:,is,i_theta,i_phi)*derf1(k)+ds*dstp_all_ds(:,1:nstp,1:nstp)
      d2stp_all_ds2(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,k,:,:,is,i_theta,i_phi)*derf2(k)+ds*d2stp_all_ds2(:,1:nstp,1:nstp)
    enddo
!
    stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,1,:,:,is,i_theta,i_phi)                                 &
                            +ds*(s_sqg_Bt_Bp(:,2,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp))
    dstp_all_ds(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,2,:,:,is,i_theta,i_phi)+ds*dstp_all_ds(:,1:nstp,1:nstp)
!
! End interpolation of all over $s$
!-------------------------------
! Begin interpolation of all over $\theta$
!
    sp_all(:,1:nstp)=stp_all(:,nstp,1:nstp)
    dsp_all_ds(:,1:nstp)=dstp_all_ds(:,nstp,1:nstp)
    d2sp_all_ds2(:,1:nstp)=d2stp_all_ds2(:,nstp,1:nstp)
    dsp_all_dt(:,1:nstp)=sp_all(:,1:nstp)*derf1(nstp)
    d2sp_all_dsdt(:,1:nstp)=dsp_all_ds(:,1:nstp)*derf1(nstp)
    d2sp_all_dt2(:,1:nstp)=sp_all(:,1:nstp)*derf2(nstp)
!
    do k=ns_tp_c,3,-1
      sp_all(:,1:nstp)=stp_all(:,k,1:nstp)+dtheta*sp_all(:,1:nstp)
      dsp_all_ds(:,1:nstp)=dstp_all_ds(:,k,1:nstp)+dtheta*dsp_all_ds(:,1:nstp)
      d2sp_all_ds2(:,1:nstp)=d2stp_all_ds2(:,k,1:nstp)+dtheta*d2sp_all_ds2(:,1:nstp)
      dsp_all_dt(:,1:nstp)=stp_all(:,k,1:nstp)*derf1(k)+dtheta*dsp_all_dt(:,1:nstp)
      d2sp_all_dsdt(:,1:nstp)=dstp_all_ds(:,k,1:nstp)*derf1(k)+dtheta*d2sp_all_dsdt(:,1:nstp)
      d2sp_all_dt2(:,1:nstp)=stp_all(:,k,1:nstp)*derf2(k)+dtheta*d2sp_all_dt2(:,1:nstp)
    enddo
!
    sp_all(:,1:nstp)=stp_all(:,1,1:nstp)                                                    &
                    +dtheta*(stp_all(:,2,1:nstp)+dtheta*sp_all(:,1:nstp))
    dsp_all_ds(:,1:nstp)=dstp_all_ds(:,1,1:nstp)                                            &
                        +dtheta*(dstp_all_ds(:,2,1:nstp)+dtheta*dsp_all_ds(:,1:nstp))
    d2sp_all_ds2(:,1:nstp)=d2stp_all_ds2(:,1,1:nstp)                                        &
                          +dtheta*(d2stp_all_ds2(:,2,1:nstp)+dtheta*d2sp_all_ds2(:,1:nstp))
    dsp_all_dt(:,1:nstp)=stp_all(:,2,1:nstp)+dtheta*dsp_all_dt(:,1:nstp)
    d2sp_all_dsdt(:,1:nstp)=dstp_all_ds(:,2,1:nstp)+dtheta*d2sp_all_dsdt(:,1:nstp)
!
! End interpolation of all over $\theta$
!--------------------------------
! Begin interpolation of all over $\varphi$
!
    qua=sp_all(:,nstp)
    dqua_dr=dsp_all_ds(:,nstp)
    dqua_dt=dsp_all_dt(:,nstp)
    dqua_dp=qua*derf1(nstp)
!
    d2qua_dr2=d2sp_all_ds2(:,nstp)
    d2qua_drdt=d2sp_all_dsdt(:,nstp)
    d2qua_drdp=dqua_dr*derf1(nstp)
    d2qua_dt2=d2sp_all_dt2(:,nstp)
    d2qua_dtdp=dqua_dt*derf1(nstp)
    d2qua_dp2=qua*derf2(nstp)
!
    do k=ns_tp_c,3,-1
      qua=sp_all(:,k)+dphi*qua
      dqua_dr=dsp_all_ds(:,k)+dphi*dqua_dr
      dqua_dt=dsp_all_dt(:,k)+dphi*dqua_dt
      dqua_dp=sp_all(:,k)*derf1(k)+dphi*dqua_dp
!
      d2qua_dr2=d2sp_all_ds2(:,k)+dphi*d2qua_dr2
      d2qua_drdt=d2sp_all_dsdt(:,k)+dphi*d2qua_drdt
      d2qua_drdp=dsp_all_ds(:,k)*derf1(k)+dphi*d2qua_drdp
      d2qua_dt2=d2sp_all_dt2(:,k)+dphi*d2qua_dt2
      d2qua_dtdp=dsp_all_dt(:,k)*derf1(k)+dphi*d2qua_dtdp
      d2qua_dp2=sp_all(:,k)*derf2(k)+dphi*d2qua_dp2
    enddo
!
    qua=sp_all(:,1)+dphi*(sp_all(:,2)+dphi*qua)
    dqua_dr=dsp_all_ds(:,1)+dphi*(dsp_all_ds(:,2)+dphi*dqua_dr)
    dqua_dt=dsp_all_dt(:,1)+dphi*(dsp_all_dt(:,2)+dphi*dqua_dt)
!
    d2qua_dr2=d2sp_all_ds2(:,1)+dphi*(d2sp_all_ds2(:,2)+dphi*d2qua_dr2)
    d2qua_drdt=d2sp_all_dsdt(:,1)+dphi*(d2sp_all_dsdt(:,2)+dphi*d2qua_drdt)
    d2qua_dt2=d2sp_all_dt2(:,1)+dphi*(d2sp_all_dt2(:,2)+dphi*d2qua_dt2)
!
    dqua_dp=sp_all(:,2)+dphi*dqua_dp
    d2qua_drdp=dsp_all_ds(:,2)+dphi*d2qua_drdp
    d2qua_dtdp=dsp_all_dt(:,2)+dphi*d2qua_dtdp
!
! End interpolation of all over $\varphi$
!
    drhods=0.5d0/rho_tor
    drhods2=drhods**2
    d2rhods2m=drhods2/rho_tor
!
    d2qua_dr2=d2qua_dr2*drhods2-dqua_dr*d2rhods2m
    dqua_dr=dqua_dr*drhods
    d2qua_drdt=d2qua_drdt*drhods
    d2qua_drdp=d2qua_drdp*drhods
!
    sqg_c=qua(1)
    B_vartheta_c=qua(2)
    B_varphi_c=qua(3)
!
    dsqg_c_dr=dqua_dr(1)
    dB_vartheta_c_dr=dqua_dr(2)
    dB_varphi_c_dr=dqua_dr(3)
!
    dsqg_c_dt=dqua_dt(1)
    dB_vartheta_c_dt=dqua_dt(2)
    dB_varphi_c_dt=dqua_dt(3)
!
    dsqg_c_dp=dqua_dp(1)
    dB_vartheta_c_dp=dqua_dp(2)
    dB_varphi_c_dp=dqua_dp(3)
!
    d2sqg_rr=d2qua_dr2(1)
    d2bth_rr=d2qua_dr2(2)
    d2bph_rr=d2qua_dr2(3)
!
    d2sqg_rt=d2qua_drdt(1)
    d2bth_rt=d2qua_drdt(2)
    d2bph_rt=d2qua_drdt(3)
!
    d2sqg_rp=d2qua_drdp(1)
    d2bth_rp=d2qua_drdp(2)
    d2bph_rp=d2qua_drdp(3)
!
    d2sqg_tt=d2qua_dt2(1)
    d2bth_tt=d2qua_dt2(2)
    d2bph_tt=d2qua_dt2(3)
!
    d2sqg_tp=d2qua_dtdp(1)
    d2bth_tp=d2qua_dtdp(2)
    d2bph_tp=d2qua_dtdp(3)
!
    d2sqg_pp=d2qua_dp2(1)
    d2bth_pp=d2qua_dp2(2)
    d2bph_pp=d2qua_dp2(3)
!
!--------------------------------
  elseif(mode_secders.eq.1) then
!--------------------------------
!
! Begin interpolation of all over $s$
!
    ns_s_p1=ns_s_c+1
    stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,ns_s_p1,:,:,is,i_theta,i_phi)
    dstp_all_ds(:,1:nstp,1:nstp)=stp_all(:,1:nstp,1:nstp)*derf1(ns_s_p1)
    d2stp_all_ds2(:,1:nstp,1:nstp)=stp_all(:,1:nstp,1:nstp)*derf2(ns_s_p1)
!
    do k=ns_s_c,3,-1
      stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,k,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp)
      dstp_all_ds(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,k,:,:,is,i_theta,i_phi)*derf1(k)+ds*dstp_all_ds(:,1:nstp,1:nstp)
      d2stp_all_ds2(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,k,:,:,is,i_theta,i_phi)*derf2(k)+ds*d2stp_all_ds2(:,1:nstp,1:nstp)
    enddo
!
    stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,1,:,:,is,i_theta,i_phi)                                 &
                            +ds*(s_sqg_Bt_Bp(:,2,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp))
    dstp_all_ds(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,2,:,:,is,i_theta,i_phi)+ds*dstp_all_ds(:,1:nstp,1:nstp)
!
! End interpolation of all over $s$
!-------------------------------
! Begin interpolation of all over $\theta$
!
    sp_all(:,1:nstp)=stp_all(:,nstp,1:nstp)
    dsp_all_ds(:,1:nstp)=dstp_all_ds(:,nstp,1:nstp)
    d2sp_all_ds2(:,1:nstp)=d2stp_all_ds2(:,nstp,1:nstp)
    dsp_all_dt(:,1:nstp)=sp_all(:,1:nstp)*derf1(nstp)
!
    do k=ns_tp_c,2,-1
      sp_all(:,1:nstp)=stp_all(:,k,1:nstp)+dtheta*sp_all(:,1:nstp)
      dsp_all_ds(:,1:nstp)=dstp_all_ds(:,k,1:nstp)+dtheta*dsp_all_ds(:,1:nstp)
      d2sp_all_ds2(:,1:nstp)=d2stp_all_ds2(:,k,1:nstp)+dtheta*d2sp_all_ds2(:,1:nstp)
      dsp_all_dt(:,1:nstp)=stp_all(:,k,1:nstp)*derf1(k)+dtheta*dsp_all_dt(:,1:nstp)
    enddo
!
    sp_all(:,1:nstp)=stp_all(:,1,1:nstp)+dtheta*sp_all(:,1:nstp)
    dsp_all_ds(:,1:nstp)=dstp_all_ds(:,1,1:nstp)+dtheta*dsp_all_ds(:,1:nstp)
    d2sp_all_ds2(:,1:nstp)=d2stp_all_ds2(:,1,1:nstp)+dtheta*d2sp_all_ds2(:,1:nstp)
!
! End interpolation of all over $\theta$
!--------------------------------
! Begin interpolation of all over $\varphi$
!
    qua=sp_all(:,nstp)
    dqua_dr=dsp_all_ds(:,nstp)
    dqua_dt=dsp_all_dt(:,nstp)
    dqua_dp=qua*derf1(nstp)
!
    d2qua_dr2=d2sp_all_ds2(:,nstp)
!
    do k=ns_tp_c,2,-1
      qua=sp_all(:,k)+dphi*qua
      dqua_dr=dsp_all_ds(:,k)+dphi*dqua_dr
      dqua_dt=dsp_all_dt(:,k)+dphi*dqua_dt
      dqua_dp=sp_all(:,k)*derf1(k)+dphi*dqua_dp
!
      d2qua_dr2=d2sp_all_ds2(:,k)+dphi*d2qua_dr2
    enddo
!
    qua=sp_all(:,1)+dphi*qua
    dqua_dr=dsp_all_ds(:,1)+dphi*dqua_dr
    dqua_dt=dsp_all_dt(:,1)+dphi*dqua_dt
!
    d2qua_dr2=d2sp_all_ds2(:,1)+dphi*d2qua_dr2
!
! End interpolation of all over $\varphi$
!
    drhods=0.5d0/rho_tor
    drhods2=drhods**2
    d2rhods2m=drhods2/rho_tor
!
    d2qua_dr2=d2qua_dr2*drhods2-dqua_dr*d2rhods2m
    dqua_dr=dqua_dr*drhods
!
    sqg_c=qua(1)
    B_vartheta_c=qua(2)
    B_varphi_c=qua(3)
!
    dsqg_c_dr=dqua_dr(1)
    dB_vartheta_c_dr=dqua_dr(2)
    dB_varphi_c_dr=dqua_dr(3)
!
    dsqg_c_dt=dqua_dt(1)
    dB_vartheta_c_dt=dqua_dt(2)
    dB_varphi_c_dt=dqua_dt(3)
!
    dsqg_c_dp=dqua_dp(1)
    dB_vartheta_c_dp=dqua_dp(2)
    dB_varphi_c_dp=dqua_dp(3)
!
    d2sqg_rr=d2qua_dr2(1)
    d2bth_rr=d2qua_dr2(2)
    d2bph_rr=d2qua_dr2(3)
!
!--------------------------------
  else
!--------------------------------
!
! Begin interpolation of all over $s$
!
    ns_s_p1=ns_s_c+1
    stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,ns_s_p1,:,:,is,i_theta,i_phi)
    dstp_all_ds(:,1:nstp,1:nstp)=stp_all(:,1:nstp,1:nstp)*derf1(ns_s_p1)
!
    do k=ns_s_c,2,-1
      stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,k,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp)
      dstp_all_ds(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,k,:,:,is,i_theta,i_phi)*derf1(k)+ds*dstp_all_ds(:,1:nstp,1:nstp)
    enddo
!
    stp_all(:,1:nstp,1:nstp)=s_sqg_Bt_Bp(:,1,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp)
!
! End interpolation of all over $s$
!-------------------------------
! Begin interpolation of all over $\theta$
!
    sp_all(:,1:nstp)=stp_all(:,nstp,1:nstp)
    dsp_all_ds(:,1:nstp)=dstp_all_ds(:,nstp,1:nstp)
    dsp_all_dt(:,1:nstp)=sp_all(:,1:nstp)*derf1(nstp)
!
    do k=ns_tp_c,2,-1
      sp_all(:,1:nstp)=stp_all(:,k,1:nstp)+dtheta*sp_all(:,1:nstp)
      dsp_all_ds(:,1:nstp)=dstp_all_ds(:,k,1:nstp)+dtheta*dsp_all_ds(:,1:nstp)
      dsp_all_dt(:,1:nstp)=stp_all(:,k,1:nstp)*derf1(k)+dtheta*dsp_all_dt(:,1:nstp)
    enddo
!
    sp_all(:,1:nstp)=stp_all(:,1,1:nstp)+dtheta*sp_all(:,1:nstp)
    dsp_all_ds(:,1:nstp)=dstp_all_ds(:,1,1:nstp)+dtheta*dsp_all_ds(:,1:nstp)
!
! End interpolation of all over $\theta$
!--------------------------------
! Begin interpolation of all over $\varphi$
!
    qua=sp_all(:,nstp)
    dqua_dr=dsp_all_ds(:,nstp)
    dqua_dt=dsp_all_dt(:,nstp)
    dqua_dp=qua*derf1(nstp)
!
    do k=ns_tp_c,2,-1
      qua=sp_all(:,k)+dphi*qua
      dqua_dr=dsp_all_ds(:,k)+dphi*dqua_dr
      dqua_dt=dsp_all_dt(:,k)+dphi*dqua_dt
      dqua_dp=sp_all(:,k)*derf1(k)+dphi*dqua_dp
    enddo
!
    qua=sp_all(:,1)+dphi*qua
    dqua_dr=dsp_all_ds(:,1)+dphi*dqua_dr
    dqua_dt=dsp_all_dt(:,1)+dphi*dqua_dt
!
! End interpolation of all over $\varphi$
!
    drhods=0.5d0/rho_tor
!
    dqua_dr=dqua_dr*drhods
!
    sqg_c=qua(1)
    B_vartheta_c=qua(2)
    B_varphi_c=qua(3)
!
    dsqg_c_dr=dqua_dr(1)
    dB_vartheta_c_dr=dqua_dr(2)
    dB_varphi_c_dr=dqua_dr(3)
!
    dsqg_c_dt=dqua_dt(1)
    dB_vartheta_c_dt=dqua_dt(2)
    dB_varphi_c_dt=dqua_dt(3)
!
    dsqg_c_dp=dqua_dp(1)
    dB_vartheta_c_dp=dqua_dp(2)
    dB_varphi_c_dp=dqua_dp(3)
!
!--------------------------------
  endif
!--------------------------------
!
  end subroutine splint_can_coord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine can_to_vmec(r,vartheta_c_in,varphi_c_in,theta_vmec,varphi_vmec)
!
  use exchange_get_cancoord_mod, only : vartheta_c,varphi_c,theta
!
  implicit none
!
  logical :: fullset
  integer :: mode_secders
  double precision, intent(in) :: r,vartheta_c_in,varphi_c_in
  double precision, intent(out) :: theta_vmec,varphi_vmec
  double precision :: A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c,     &
                      d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                      d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                      d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp
  double precision, dimension(1) :: y,dy
!
  fullset=.true.
  mode_secders=0
!
  call splint_can_coord(fullset,mode_secders,r,vartheta_c_in,varphi_c_in,                &
                        A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                        sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                        B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                        B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                        d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                        d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                        d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp,G_c)
!
  vartheta_c=vartheta_c_in
  varphi_c=varphi_c_in
  y(1)=G_c
!
!  call rhs_cancoord(r,y,dy)      !<=OLD
  call rhs_cancoord(sqrt(r),y,dy) !<=NEW
!
  theta_vmec=theta
  varphi_vmec=varphi_c_in+G_c
!
  end subroutine can_to_vmec
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine deallocate_can_coord
!
  use canonical_coordinates_mod, only : s_sqg_Bt_Bp,s_G_c
!
  implicit none
!
  if(allocated(s_sqg_Bt_Bp)) deallocate(s_sqg_Bt_Bp)
  if(allocated(s_G_c)) deallocate(s_G_c)
!
  end subroutine deallocate_can_coord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine vmec_to_can(r,theta,varphi,vartheta_c,varphi_c)
!
! Input : r,theta,varphi      - VMEC coordinates
! Output: vartheta_c,varphi_c - canonical coordinates
!
  implicit none
!
  double precision, parameter :: epserr=1.d-14
  integer,          parameter :: niter=100
  integer          :: iter
  double precision :: r,theta,varphi
  double precision, intent(out) :: vartheta_c,varphi_c
  double precision :: delthe,delphi,alam,dl_dt,vartheta
!
  call splint_lambda(r,theta,varphi,alam,dl_dt)
!
  vartheta=theta+alam
!
  vartheta_c=vartheta
  varphi_c=varphi
!
  do iter=1,niter
!
    call newt_step
!
    vartheta_c=vartheta_c+delthe
    varphi_c=varphi_c+delphi
    if(abs(delthe)+abs(delphi).lt.epserr) exit
  enddo
!
!------------------------------------------
!
  contains
!
!------------------------------------------
!
  subroutine newt_step
!
  use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,hs_c,h_theta_c,h_phi_c,    &
                                        ns_s_c,ns_tp_c,ns_max,derf1,s_G_c
  use vector_potentail_mod, only : ns,hs,torflux,sA_phi
  use new_vmec_stuff_mod,   only : nper,ns_A
  use chamb_mod,            only : rnegflag
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  integer :: nstp,ns_A_p1,ns_s_p1
  integer :: k,is,i_theta,i_phi
  integer :: iss,ist,isp
!
  double precision :: A_phi,A_theta,dA_phi_dr,dA_theta_dr
  double precision :: s,ds,dtheta,dphi,rho_tor,drhods,drhods2,d2rhods2m
  double precision :: aiota,G_c,dG_c_dt,dG_c_dp
  double precision :: ts,ps,dts_dtc,dts_dpc,dps_dtc,dps_dpc,det
!
  double precision, dimension(ns_max)              :: sp_G,dsp_G_dt
  double precision, dimension(ns_max,ns_max)       :: stp_G
!
  if(r.le.0.d0) then
    rnegflag=.true.
    r=abs(r)
  endif
!
  dA_theta_dr=torflux
!
  dtheta=modulo(vartheta_c,twopi)/h_theta_c
  i_theta=max(0,min(n_theta_c-1,int(dtheta)))
  dtheta=(dtheta-dble(i_theta))*h_theta_c
  i_theta=i_theta+1
!
  dphi=modulo(varphi_c,twopi/dble(nper))/h_phi_c
  i_phi=max(0,min(n_phi_c-1,int(dphi)))
  dphi=(dphi-dble(i_phi))*h_phi_c
  i_phi=i_phi+1
!
! Begin interpolation of vector potentials over $s$
!
  ds=r/hs
  is=max(0,min(ns-1,int(ds)))
  ds=(ds-dble(is))*hs
  is=is+1
!
  ns_A_p1=ns_A+1
  A_phi=sA_phi(ns_A_p1,is)
  dA_phi_dr=A_phi*derf1(ns_A_p1)
!
  do k=ns_A,3,-1
    dA_phi_dr=sA_phi(k,is)*derf1(k)+ds*dA_phi_dr
  enddo
!
  dA_phi_dr=sA_phi(2,is)+ds*dA_phi_dr
!
  aiota=-dA_phi_dr/dA_theta_dr
!
! End interpolation of vector potentials over $s$
!
  rho_tor=sqrt(r)
  !hs_c=hs !added by Johanna in alalogy to get_canonical_coordinites to make test_orbits_vmec working
  ds=rho_tor/hs_c
  is=max(0,min(ns_c-1,int(ds)))
  ds=(ds-dble(is))*hs_c
  is=is+1
!
  nstp=ns_tp_c+1
!
! Begin interpolation of G over $s$
!
  stp_G(1:nstp,1:nstp)=s_G_c(ns_s_c+1,:,:,is,i_theta,i_phi)
!
  do k=ns_s_c,1,-1
    stp_G(1:nstp,1:nstp)=s_G_c(k,:,:,is,i_theta,i_phi)+ds*stp_G(1:nstp,1:nstp)
  enddo
!
! End interpolation of G over $s$
! Begin interpolation of G over $\theta$
!
  sp_G(1:nstp)=stp_G(nstp,1:nstp)
  dsp_G_dt(1:nstp)=sp_G(1:nstp)*derf1(nstp)
!
  do k=ns_tp_c,1,-1
    sp_G(1:nstp)=stp_G(k,1:nstp)+dtheta*sp_G(1:nstp)
    if(k.gt.1) dsp_G_dt(1:nstp)=stp_G(k,1:nstp)*derf1(k)+dtheta*dsp_G_dt(1:nstp)
  enddo
!
! End interpolation of G over $\theta$
! Begin interpolation of G over $\varphi$
!
  G_c=sp_G(nstp)
  dG_c_dt=dsp_G_dt(nstp)
  dG_c_dp=G_c*derf1(nstp)
!
  do k=ns_tp_c,1,-1
    G_c=sp_G(k)+dphi*G_c
    dG_c_dt=dsp_G_dt(k)+dphi*dG_c_dt
    if(k.gt.1) dG_c_dp=sp_G(k)*derf1(k)+dphi*dG_c_dp
  enddo
!
! End interpolation of G over $\varphi$
!
  ts=vartheta_c+aiota*G_c-vartheta
  ps=varphi_c+G_c-varphi
  dts_dtc=1.d0+aiota*dG_c_dt
  dts_dpc=aiota*dG_c_dp
  dps_dtc=dG_c_dt
  dps_dpc=1.d0+dG_c_dp
  det=1.d0+aiota*dG_c_dt+dG_c_dp
!
  delthe=(ps*dts_dpc-ts*dps_dpc)/det
  delphi=(ts*dps_dtc-ps*dts_dtc)/det
!
  end subroutine newt_step
!
!------------------------------------------
!
  end subroutine vmec_to_can
