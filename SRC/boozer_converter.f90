!
  subroutine boozer_converter
!
  use vector_potentail_mod,   only : ns,hs
  use new_vmec_stuff_mod,     only : n_theta,n_phi,h_theta,h_phi,ns_s,ns_tp
  use boozer_coordinates_mod, only : ns_s_B,ns_tp_B,ns_B,n_theta_B,n_phi_B, &
                                     hs_B,h_theta_B,h_phi_B,                &
                                     s_Bcovar_tp_B,                         &
                                     s_Bmod_B,s_Bcovar_r_B,                 &
                                     s_delt_delp_V,s_delt_delp_B,           &
                                     ns_max,derf1,derf2,derf3,              &
                                     use_B_r, use_del_tp_B
!
  implicit none
!
  double precision, parameter :: s_min=1.d-6, rho_min=sqrt(s_min)
!
  integer :: i,k,i_rho,i_theta,i_phi,npoilag,nder,nshift,ibeg,iend,i_qua,nqua,ist,isp
  double precision :: s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                      sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                      Bcovar_r,Bcovar_vartheta,Bcovar_varphi
  double precision :: Bcovar_vartheta_B,Bcovar_varphi_B
  double precision :: denomjac,G00,Gbeg,aper,per_theta,per_phi,gridcellnum
  double precision, dimension(:),       allocatable :: wint_t,wint_p,theta_V,theta_B,   &
                                                       phi_V,phi_B,aiota_arr,rho_tor
  double precision, dimension(:,:),     allocatable :: Bcovar_theta_V,Bcovar_varphi_V,  &
                                                       bmod_Vg,bmod_Bg,alam_2D,         &
                                                       deltheta_BV_Vg,delphi_BV_Vg,     &
                                                       deltheta_BV_Bg,delphi_BV_Bg,     &
                                                       splcoe_r,splcoe_t,splcoe_p,coef, &
                                                       perqua_t,perqua_p
  double precision, dimension(:,:,:),   allocatable :: perqua_2D,Gfunc
  double precision, dimension(:,:,:,:), allocatable :: Bcovar_symfl
!
  nqua=6
!
  ns_s_B=ns_s
  ns_tp_B=ns_tp
  ns_B=ns
  n_theta_B=n_theta
  n_phi_B=n_phi
!
  gridcellnum=dble((n_theta_B-1)*(n_phi_B-1))
!
  npoilag=ns_tp_B+1
  nder=0
  nshift=npoilag/2
!
  hs_B=hs*dble(ns-1)/dble(ns_B-1)
  h_theta_B=h_theta*dble(n_theta-1)/dble(n_theta_B-1)
  h_phi_B=h_phi*dble(n_phi-1)/dble(n_phi_B-1)
!
!
!-----------------------------------------------------------------------
!
!
! Compute Boozer data
!
  print *,'Transforming to Boozer coordinates'
!
  if(use_B_r) then
    print *,'B_r is computed'
  else
    print *,'B_r is not computed'
  endif
!
  G00=0.d0
!
  allocate(rho_tor(ns_B))
  if(use_B_r) then
    allocate(aiota_arr(ns_B))
    allocate(Gfunc(ns_B,n_theta_B,n_phi_B), Bcovar_symfl(3,ns_B,n_theta_B,n_phi_B))
  endif
!
  allocate(Bcovar_theta_V(n_theta_B,n_phi_B),Bcovar_varphi_V(n_theta_B,n_phi_B), &
                          bmod_Vg(n_theta_B,n_phi_B),alam_2D(n_theta_B,n_phi_B))
  allocate(deltheta_BV_Vg(n_theta_B,n_phi_B),delphi_BV_Vg(n_theta_B,n_phi_B),    &
           deltheta_BV_Bg(n_theta_B,n_phi_B),delphi_BV_Bg(n_theta_B,n_phi_B))
  allocate(wint_t(0:ns_tp_B),wint_p(0:ns_tp_B))
  allocate(coef(0:nder,npoilag))
  allocate(theta_V(2-n_theta_B:2*n_theta_B-1),theta_B(2-n_theta_B:2*n_theta_B-1))
  allocate(phi_V(2-n_phi_B:2*n_phi_B-1),phi_B(2-n_phi_B:2*n_phi_B-1))
  allocate(perqua_t(nqua,2-n_theta_B:2*n_theta_B-1),perqua_p(nqua,2-n_phi_B:2*n_phi_B-1))
  allocate(perqua_2D(nqua,n_theta_B,n_phi_B))
!
  allocate(splcoe_t(0:ns_tp_B,n_theta_B),splcoe_p(0:ns_tp_B,n_phi_B))
!
! allocate spline coefficients for Boozer data:
  if(.not.allocated(s_Bcovar_tp_B)) &
          allocate(s_Bcovar_tp_B(2,ns_s_B+1,ns_B))
  if(.not.allocated(s_Bmod_B)) &
          allocate(s_Bmod_B(ns_s_B+1,ns_tp_B+1,ns_tp_B+1,ns_B,n_theta_B,n_phi_B))
  if(use_B_r.and..not.allocated(s_Bcovar_r_B)) &
          allocate(s_Bcovar_r_B(ns_s_B+1,ns_tp_B+1,ns_tp_B+1,ns_B,n_theta_B,n_phi_B))
  if(.not.allocated(s_delt_delp_V)) &
          allocate(s_delt_delp_V(2,ns_s_B+1,ns_tp_B+1,ns_tp_B+1,ns_B,n_theta_B,n_phi_B))
  if(use_del_tp_B.and..not.allocated(s_delt_delp_B)) &
          allocate(s_delt_delp_B(2,ns_s_B+1,ns_tp_B+1,ns_tp_B+1,ns_B,n_theta_B,n_phi_B))
!
  do i=0,ns_tp_B
    wint_t(i)=h_theta_B**(i+1)/dble(i+1)
    wint_p(i)=h_phi_B**(i+1)/dble(i+1)
  enddo

  ! Set theta_V and phi_V linear, with value 0 at index 1 and stepsize h.
  ! Then expand this in both directions beyond 1:n_theta_B.
  do i_theta=1,n_theta_B
    theta_V(i_theta)=dble(i_theta-1)*h_theta_B
  enddo
  per_theta=dble(n_theta_B-1)*h_theta_B
  theta_V(2-n_theta_B:0)=theta_V(1:n_theta_B-1)-per_theta
  theta_V(n_theta_B+1:2*n_theta_B-1)=theta_V(2:n_theta_B)+per_theta

  do i_phi=1,n_phi_B
    phi_V(i_phi)=dble(i_phi-1)*h_phi_B
  enddo
  per_phi=dble(n_phi_B-1)*h_phi_B
  phi_V(2-n_phi_B:0)=phi_V(1:n_phi_B-1)-per_phi
  phi_V(n_phi_B+1:2*n_phi_B-1)=phi_V(2:n_phi_B)+per_phi

  do i_rho=1,ns_B
    rho_tor(i_rho)=max(dble(i_rho-1)*hs_B,rho_min)
    s=rho_tor(i_rho)**2
!
    do i_theta=1,n_theta_B
      theta=dble(i_theta-1)*h_theta_B
      do i_phi=1,n_phi_B
        varphi=dble(i_phi-1)*h_phi_B
!
        call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                        sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                        Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
        alam_2D(i_theta,i_phi)=alam
        bmod_Vg(i_theta,i_phi)=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
        Bcovar_theta_V(i_theta,i_phi)=Bcovar_vartheta*(1.d0+dl_dt)
        Bcovar_varphi_V(i_theta,i_phi)=Bcovar_varphi+Bcovar_vartheta*dl_dp
        perqua_2D(4,i_theta,i_phi)=Bcovar_r
        perqua_2D(5,i_theta,i_phi)=Bcovar_vartheta
        perqua_2D(6,i_theta,i_phi)=Bcovar_varphi
      enddo
    enddo
!
! covariant components $B_\vartheta$ and $B_\varphi$ of Boozer coordinates:
    Bcovar_vartheta_B=sum(Bcovar_theta_V(2:n_theta_B,2:n_phi_B))/gridcellnum
    Bcovar_varphi_B=sum(Bcovar_varphi_V(2:n_theta_B,2:n_phi_B))/gridcellnum
    s_Bcovar_tp_B(1,1,i_rho)=Bcovar_vartheta_B
    s_Bcovar_tp_B(2,1,i_rho)=Bcovar_varphi_B
!
    denomjac=1.d0/(aiota*Bcovar_vartheta_B+Bcovar_varphi_B)
    Gbeg=G00+Bcovar_vartheta_B*denomjac*alam_2D(1,1)
!
    splcoe_t(0,:)=Bcovar_theta_V(:,1)
!
    call spl_per(ns_tp_B,n_theta_B,h_theta_B,splcoe_t)
!
    delphi_BV_Vg(1,1)=0.d0
    do i_theta=1,n_theta_B-1
      delphi_BV_Vg(i_theta+1,1)=delphi_BV_Vg(i_theta,1)+sum(wint_t*splcoe_t(:,i_theta))
    enddo
    ! Remove linear increasing component from delphi_BV_Vg
    aper=(delphi_BV_Vg(n_theta_B,1)-delphi_BV_Vg(1,1))/dble(n_theta_B-1)
    do i_theta=2,n_theta_B
      delphi_BV_Vg(i_theta,1)=delphi_BV_Vg(i_theta,1)-aper*dble(i_theta-1)
    enddo
!
    do i_theta=1,n_theta_B
      splcoe_p(0,:)=Bcovar_varphi_V(i_theta,:)
!
      call spl_per(ns_tp_B,n_phi_B,h_phi_B,splcoe_p)
!
      do i_phi=1,n_phi_B-1
        delphi_BV_Vg(i_theta,i_phi+1)=delphi_BV_Vg(i_theta,i_phi)+sum(wint_p*splcoe_p(:,i_phi))
      enddo
      aper=(delphi_BV_Vg(i_theta,n_phi_B)-delphi_BV_Vg(i_theta,1))/dble(n_phi_B-1)
      do i_phi=2,n_phi_B
        delphi_BV_Vg(i_theta,i_phi)=delphi_BV_Vg(i_theta,i_phi)-aper*dble(i_phi-1)
      enddo
    enddo
!
! difference between Boozer and VMEC toroidal angle, $\Delta \varphi_{BV}=\varphi_B-\varphi=G$:
    delphi_BV_Vg=denomjac*delphi_BV_Vg+Gbeg
! difference between Boozer and VMEC poloidal angle, $\Delta \vartheta_{BV}=\vartheta_B-\theta$:
    deltheta_BV_Vg=aiota*delphi_BV_Vg+alam_2D
!
    s_delt_delp_V(1,1,1,1,i_rho,:,:)=deltheta_BV_Vg
    s_delt_delp_V(2,1,1,1,i_rho,:,:)=delphi_BV_Vg
!
! At this point, all quantities are specified on equidistant grid in VMEC angles $(\theta,\varphi)$
!
! Re-interpolate to equidistant grid in $(\vartheta_B,\varphi)$:
!
    do i_phi=1,n_phi_B
      perqua_t(1,1:n_theta_B)=deltheta_BV_Vg(:,i_phi)
      perqua_t(2,1:n_theta_B)=delphi_BV_Vg(:,i_phi)
      perqua_t(3,1:n_theta_B)=bmod_Vg(:,i_phi)
      perqua_t(4:6,1:n_theta_B)=perqua_2D(4:6,:,i_phi)
      ! Extend range of theta values
      perqua_t(:,2-n_theta_B:0)=perqua_t(:,1:n_theta_B-1)
      perqua_t(:,n_theta_B+1:2*n_theta_B-1)=perqua_t(:,2:n_theta_B)
      theta_B=theta_V+perqua_t(1,:)
      do i_theta=1,n_theta_B
!
        call binsrc(theta_B,2-n_theta_B,2*n_theta_B-1,theta_V(i_theta),i)
!
        ibeg=i-nshift
        iend=ibeg+ns_tp_B
!
        call plag_coeff(npoilag,nder,theta_V(i_theta),theta_B(ibeg:iend),coef)
!
        perqua_2D(:,i_theta,i_phi)=matmul(perqua_t(:,ibeg:iend),coef(0,:))
      enddo
    enddo
!
! End re-interpolate to equidistant grid in $(\vartheta_B,\varphi)$
!
! Re-interpolate to equidistant grid in $(\vartheta_B,\varphi_B)$:
!
    do i_theta=1,n_theta_B
      perqua_p(:,1:n_phi_B)=perqua_2D(:,i_theta,:)
      perqua_p(:,2-n_phi_B:0)=perqua_p(:,1:n_phi_B-1)
      ! Extend range of phi values
      perqua_p(:,n_phi_B+1:2*n_phi_B-1)=perqua_p(:,2:n_phi_B)
      phi_B=phi_V+perqua_p(2,:)
      do i_phi=1,n_phi_B
!
        call binsrc(phi_B,2-n_phi_B,2*n_phi_B-1,phi_V(i_phi),i)
!
        ibeg=i-nshift
        iend=ibeg+ns_tp_B
!
        call plag_coeff(npoilag,nder,phi_V(i_phi),phi_B(ibeg:iend),coef)
!
        perqua_2D(:,i_theta,i_phi)=matmul(perqua_p(:,ibeg:iend),coef(0,:))
      enddo
    enddo
!
    if(use_del_tp_B) s_delt_delp_B(:,1,1,1,i_rho,:,:)=perqua_2D(1:2,:,:)
    s_Bmod_B(1,1,1,i_rho,:,:)=perqua_2D(3,:,:)
!
! End re-interpolate to equidistant grid in $(\vartheta_B,\varphi_B)$
!
    if(use_B_r) then
      aiota_arr(i_rho)=aiota
      Gfunc(i_rho,:,:)=perqua_2D(2,:,:)
! covariant components $B_k$ in symmetry flux coordinates on equidistant grid of Boozer coordinates:
      Bcovar_symfl(:,i_rho,:,:)=perqua_2D(4:6,:,:)
    endif
!
  enddo
!
! Compute radial covariant magnetic field component in Boozer coordinates
!
  if(use_B_r) then
    if(allocated(coef)) deallocate(coef)
    nder=1
    npoilag=5
    nshift=npoilag/2
    allocate(coef(0:nder,npoilag))
!
    do i_rho=1,ns_B
      ibeg=i_rho-nshift
      iend=ibeg+npoilag-1
      if(ibeg.lt.1) then
        ibeg=1
        iend=ibeg+npoilag-1
      elseif(iend.gt.ns_B) then
        iend=ns_B
        ibeg=iend-npoilag+1
      endif
!
      call plag_coeff(npoilag,nder,rho_tor(i_rho),rho_tor(ibeg:iend),coef)
!
! We spline covariant component $B_\rho$ instead of $B_s$:
      do i_phi=1,n_phi_B
        s_Bcovar_r_B(1,1,1,i_rho,:,i_phi)=2.d0*rho_tor(i_rho)*Bcovar_symfl(1,i_rho,:,i_phi) &
                           -matmul(coef(1,:)*aiota_arr(ibeg:iend),Gfunc(ibeg:iend,:,i_phi)) &
                           *Bcovar_symfl(2,i_rho,:,i_phi) &
                           -matmul(coef(1,:),Gfunc(ibeg:iend,:,i_phi)) &
                           *Bcovar_symfl(3,i_rho,:,i_phi)
      enddo
!
    enddo
    deallocate(aiota_arr,Gfunc,Bcovar_symfl)
  endif
!
! End compute radial covariant magnetic field component in Boozer coordinates
!
!
  deallocate(Bcovar_theta_V,Bcovar_varphi_V,bmod_Vg,alam_2D,          &
             deltheta_BV_Vg,delphi_BV_Vg,deltheta_BV_Bg,delphi_BV_Bg, &
             wint_t,wint_p,coef,theta_V,theta_B,phi_V,phi_B,          &
             perqua_t,perqua_p,perqua_2D)
!
  print *,'done'
!
! End compute Boozer data
!
!-----------------------------------------------------------------------
!
! Spline Boozer data
!
  print *,'Splining Boozer data'
!
  if(use_del_tp_B) then
    print *,'Delta theta and Delta phi in Boozer coordinates are splined'
  else
    print *,'Delta theta and Delta phi in Boozer coordinates are not splined'
  endif
!
! splining over $\varphi$:
!
  do i_rho=1,ns_B
    do i_theta=1,n_theta_B
!
      do i_qua=1,2
        splcoe_p(0,:)=s_delt_delp_V(i_qua,1,1,1,i_rho,i_theta,:)
!
        call spl_per(ns_tp_B,n_phi_B,h_phi_B,splcoe_p)
!
        do k=1,ns_tp_B
          s_delt_delp_V(i_qua,1,1,k+1,i_rho,i_theta,:)=splcoe_p(k,:)
        enddo
!
        if(use_del_tp_B) then
          splcoe_p(0,:)=s_delt_delp_B(i_qua,1,1,1,i_rho,i_theta,:)
!
          call spl_per(ns_tp_B,n_phi_B,h_phi_B,splcoe_p)
!
          do k=1,ns_tp_B
            s_delt_delp_B(i_qua,1,1,k+1,i_rho,i_theta,:)=splcoe_p(k,:)
          enddo
        endif
      enddo
!
      splcoe_p(0,:)=s_Bmod_B(1,1,1,i_rho,i_theta,:)
!
      call spl_per(ns_tp_B,n_phi_B,h_phi_B,splcoe_p)
!
      do k=1,ns_tp_B
        s_Bmod_B(1,1,k+1,i_rho,i_theta,:)=splcoe_p(k,:)
      enddo
!
      if(use_B_r) then
        splcoe_p(0,:)=s_Bcovar_r_B(1,1,1,i_rho,i_theta,:)
!
        call spl_per(ns_tp_B,n_phi_B,h_phi_B,splcoe_p)
!
        do k=1,ns_tp_B
          s_Bcovar_r_B(1,1,k+1,i_rho,i_theta,:)=splcoe_p(k,:)
        enddo
      endif
!
    enddo
  enddo
!
! splining over $\vartheta$:
!
  do i_rho=1,ns_B
    do i_phi=1,n_phi_B
      do isp=1,ns_tp_B+1
!
        do i_qua=1,2
          splcoe_t(0,:)=s_delt_delp_V(i_qua,1,1,isp,i_rho,:,i_phi)
!
          call spl_per(ns_tp_B,n_theta_B,h_theta_B,splcoe_t)
!
          do k=1,ns_tp_B
            s_delt_delp_V(i_qua,1,k+1,isp,i_rho,:,i_phi)=splcoe_t(k,:)
          enddo
!
          if(use_del_tp_B) then
            splcoe_t(0,:)=s_delt_delp_B(i_qua,1,1,isp,i_rho,:,i_phi)
!
            call spl_per(ns_tp_B,n_theta_B,h_theta_B,splcoe_t)
!
            do k=1,ns_tp_B
              s_delt_delp_B(i_qua,1,k+1,isp,i_rho,:,i_phi)=splcoe_t(k,:)
            enddo
          endif
        enddo
!
        splcoe_t(0,:)=s_Bmod_B(1,1,isp,i_rho,:,i_phi)
!
        call spl_per(ns_tp_B,n_theta_B,h_theta_B,splcoe_t)
!
        do k=1,ns_tp_B
          s_Bmod_B(1,k+1,isp,i_rho,:,i_phi)=splcoe_t(k,:)
        enddo
!
        if(use_B_r) then
          splcoe_t(0,:)=s_Bcovar_r_B(1,1,isp,i_rho,:,i_phi)
!
          call spl_per(ns_tp_B,n_theta_B,h_theta_B,splcoe_t)
!
          do k=1,ns_tp_B
            s_Bcovar_r_B(1,k+1,isp,i_rho,:,i_phi)=splcoe_t(k,:)
          enddo
        endif
!
      enddo
    enddo
  enddo
!
! splining over $\rho$:
!
  allocate(splcoe_r(0:ns_s_B,ns_B))
!
  do i_qua=1,2
    splcoe_r(0,:)=s_Bcovar_tp_B(i_qua,1,:)
!
    call spl_reg(ns_s_B,ns_B,hs_B,splcoe_r)
!
    do k=1,ns_s_B
      s_Bcovar_tp_B(i_qua,k+1,:)=splcoe_r(k,:)
    enddo
  enddo
!
  do i_theta=1,n_theta_B
    do i_phi=1,n_phi_B
      do ist=1,ns_tp_B+1
        do isp=1,ns_tp_B+1
!
          do i_qua=1,2
            splcoe_r(0,:)=s_delt_delp_V(i_qua,1,ist,isp,:,i_theta,i_phi)
!
            call spl_reg(ns_s_B,ns_B,hs_B,splcoe_r)
!
            do k=1,ns_s_B
              s_delt_delp_V(i_qua,k+1,ist,isp,:,i_theta,i_phi)=splcoe_r(k,:)
            enddo
!
            if(use_del_tp_B) then
              splcoe_r(0,:)=s_delt_delp_B(i_qua,1,ist,isp,:,i_theta,i_phi)
!
              call spl_reg(ns_s_B,ns_B,hs_B,splcoe_r)
!
              do k=1,ns_s_B
                s_delt_delp_B(i_qua,k+1,ist,isp,:,i_theta,i_phi)=splcoe_r(k,:)
              enddo
            endif
          enddo
!
          splcoe_r(0,:)=s_Bmod_B(1,ist,isp,:,i_theta,i_phi)
!
          call spl_reg(ns_s_B,ns_B,hs_B,splcoe_r)
!
          do k=1,ns_s_B
            s_Bmod_B(k+1,ist,isp,:,i_theta,i_phi)=splcoe_r(k,:)
          enddo
!
          if(use_B_r) then
            splcoe_r(0,:)=s_Bcovar_r_B(1,ist,isp,:,i_theta,i_phi)
!
            call spl_reg(ns_s_B,ns_B,hs_B,splcoe_r)
!
            do k=1,ns_s_B
              s_Bcovar_r_B(k+1,ist,isp,:,i_theta,i_phi)=splcoe_r(k,:)
            enddo
          endif
!
        enddo
      enddo
    enddo
  enddo
!
  deallocate(splcoe_r,splcoe_t,splcoe_p)
!
  do k=1,ns_max
    derf1(k)=dble(k-1)
    derf2(k)=dble((k-1)*(k-2))
    derf3(k)=dble((k-1)*(k-2)*(k-3))
  enddo
!
  print *,'done'
!
! End spline Boozer data
!
  end subroutine boozer_converter
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine splint_boozer_coord(r,vartheta_B,varphi_B,                                       &
                                 A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                                 B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                                 B_varphi_B,dB_varphi_B,d2B_varphi_B,                         &
                                 Bmod_B,dBmod_B,d2Bmod_B,                                     &
                                 B_r,dB_r,d2B_r)
!
  use boozer_coordinates_mod, only : ns_s_B,ns_tp_B,ns_B,n_theta_B,n_phi_B, &
                                     hs_B,h_theta_B,h_phi_B,                &
                                     s_Bcovar_tp_B,                         &
                                     s_Bmod_B,s_Bcovar_r_B,                 &
                                     ns_max,derf1,derf2,derf3,              &
                                     use_B_r
  use vector_potentail_mod, only : ns,hs,torflux,sA_phi
  use new_vmec_stuff_mod,   only : nper,ns_A
  use chamb_mod,            only : rnegflag
use diag_mod, only : icounter
!
  implicit none
!
  integer, parameter          :: mode_secders=1
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  integer :: nstp,ns_A_p1,ns_s_p1
  integer :: k,is,i_theta,i_phi
  integer :: iss,ist,isp
!
  double precision :: r,vartheta_B,varphi_B,                                       &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3, &
                      B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                      B_varphi_B,dB_varphi_B,d2B_varphi_B,Bmod_B,B_r
  double precision, dimension(3) :: dBmod_B,dB_r
  double precision, dimension(6) :: d2Bmod_B,d2B_r
!
  double precision :: s,ds,dtheta,dphi,rho_tor,drhods,drhods2,d2rhods2m
  double precision :: qua,dqua_dr,dqua_dt,dqua_dp
  double precision :: d2qua_dr2,d2qua_drdt,d2qua_drdp,d2qua_dt2,d2qua_dtdp,d2qua_dp2
  double precision, dimension(ns_max)        :: sp_all,dsp_all_ds,dsp_all_dt
  double precision, dimension(ns_max)        :: d2sp_all_ds2,d2sp_all_dsdt,d2sp_all_dt2
  double precision, dimension(ns_max,ns_max) :: stp_all,dstp_all_ds,d2stp_all_ds2
!
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
  dtheta=modulo(vartheta_B,twopi)/h_theta_B
  i_theta=max(0,min(n_theta_B-1,int(dtheta)))
  dtheta=(dtheta-dble(i_theta))*h_theta_B
  i_theta=i_theta+1
!
  dphi=modulo(varphi_B,twopi/dble(nper))/h_phi_B
  i_phi=max(0,min(n_phi_B-1,int(dphi)))
  dphi=(dphi-dble(i_phi))*h_phi_B
  i_phi=i_phi+1
!
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
!
  rho_tor=sqrt(r)
  ds=rho_tor/hs_B
  is=max(0,min(ns_B-1,int(ds)))
  ds=(ds-dble(is))*hs_B
  is=is+1
!
  nstp=ns_tp_B+1
  ns_s_p1=ns_s_B+1
!
!--------------------------------
! Interpolation of mod-B:
!--------------------------------
!
! Begin interpolation of mod-B over $rho$
!
  stp_all(1:nstp,1:nstp)=s_Bmod_B(ns_s_p1,:,:,is,i_theta,i_phi)
  dstp_all_ds(1:nstp,1:nstp)=stp_all(1:nstp,1:nstp)*derf1(ns_s_p1)
  d2stp_all_ds2(1:nstp,1:nstp)=stp_all(1:nstp,1:nstp)*derf2(ns_s_p1)
!
  do k=ns_s_B,3,-1
    stp_all(1:nstp,1:nstp)=s_Bmod_B(k,:,:,is,i_theta,i_phi)+ds*stp_all(1:nstp,1:nstp)
    dstp_all_ds(1:nstp,1:nstp)=s_Bmod_B(k,:,:,is,i_theta,i_phi)*derf1(k)+ds*dstp_all_ds(1:nstp,1:nstp)
    d2stp_all_ds2(1:nstp,1:nstp)=s_Bmod_B(k,:,:,is,i_theta,i_phi)*derf2(k)+ds*d2stp_all_ds2(1:nstp,1:nstp)
  enddo
!
  stp_all(1:nstp,1:nstp)=s_Bmod_B(1,:,:,is,i_theta,i_phi)                                 &
                        +ds*(s_Bmod_B(2,:,:,is,i_theta,i_phi)+ds*stp_all(1:nstp,1:nstp))
  dstp_all_ds(1:nstp,1:nstp)=s_Bmod_B(2,:,:,is,i_theta,i_phi)+ds*dstp_all_ds(1:nstp,1:nstp)
!
! End interpolation of mod-B over $rho$
!-------------------------------
! Begin interpolation of mod-B over $\theta$
!
  sp_all(1:nstp)=stp_all(nstp,1:nstp)
  dsp_all_ds(1:nstp)=dstp_all_ds(nstp,1:nstp)
  d2sp_all_ds2(1:nstp)=d2stp_all_ds2(nstp,1:nstp)
  dsp_all_dt(1:nstp)=sp_all(1:nstp)*derf1(nstp)
  d2sp_all_dsdt(1:nstp)=dsp_all_ds(1:nstp)*derf1(nstp)
  d2sp_all_dt2(1:nstp)=sp_all(1:nstp)*derf2(nstp)
!
  do k=ns_tp_B,3,-1
    sp_all(1:nstp)=stp_all(k,1:nstp)+dtheta*sp_all(1:nstp)
    dsp_all_ds(1:nstp)=dstp_all_ds(k,1:nstp)+dtheta*dsp_all_ds(1:nstp)
    d2sp_all_ds2(1:nstp)=d2stp_all_ds2(k,1:nstp)+dtheta*d2sp_all_ds2(1:nstp)
    dsp_all_dt(1:nstp)=stp_all(k,1:nstp)*derf1(k)+dtheta*dsp_all_dt(1:nstp)
    d2sp_all_dsdt(1:nstp)=dstp_all_ds(k,1:nstp)*derf1(k)+dtheta*d2sp_all_dsdt(1:nstp)
    d2sp_all_dt2(1:nstp)=stp_all(k,1:nstp)*derf2(k)+dtheta*d2sp_all_dt2(1:nstp)
  enddo
!
  sp_all(1:nstp)=stp_all(1,1:nstp)                                                    &
                +dtheta*(stp_all(2,1:nstp)+dtheta*sp_all(1:nstp))
  dsp_all_ds(1:nstp)=dstp_all_ds(1,1:nstp)                                            &
                    +dtheta*(dstp_all_ds(2,1:nstp)+dtheta*dsp_all_ds(1:nstp))
  d2sp_all_ds2(1:nstp)=d2stp_all_ds2(1,1:nstp)                                        &
                      +dtheta*(d2stp_all_ds2(2,1:nstp)+dtheta*d2sp_all_ds2(1:nstp))
  dsp_all_dt(1:nstp)=stp_all(2,1:nstp)+dtheta*dsp_all_dt(1:nstp)
  d2sp_all_dsdt(1:nstp)=dstp_all_ds(2,1:nstp)+dtheta*d2sp_all_dsdt(1:nstp)
!
! End interpolation of mod-B over $\theta$
!--------------------------------
! Begin interpolation of mod-B over $\varphi$
!
  qua=sp_all(nstp)
  dqua_dr=dsp_all_ds(nstp)
  dqua_dt=dsp_all_dt(nstp)
  dqua_dp=qua*derf1(nstp)
!
  d2qua_dr2=d2sp_all_ds2(nstp)
  d2qua_drdt=d2sp_all_dsdt(nstp)
  d2qua_drdp=dqua_dr*derf1(nstp)
  d2qua_dt2=d2sp_all_dt2(nstp)
  d2qua_dtdp=dqua_dt*derf1(nstp)
  d2qua_dp2=qua*derf2(nstp)
!
  do k=ns_tp_B,3,-1
    qua=sp_all(k)+dphi*qua
    dqua_dr=dsp_all_ds(k)+dphi*dqua_dr
    dqua_dt=dsp_all_dt(k)+dphi*dqua_dt
    dqua_dp=sp_all(k)*derf1(k)+dphi*dqua_dp
!
    d2qua_dr2=d2sp_all_ds2(k)+dphi*d2qua_dr2
    d2qua_drdt=d2sp_all_dsdt(k)+dphi*d2qua_drdt
    d2qua_drdp=dsp_all_ds(k)*derf1(k)+dphi*d2qua_drdp
    d2qua_dt2=d2sp_all_dt2(k)+dphi*d2qua_dt2
    d2qua_dtdp=dsp_all_dt(k)*derf1(k)+dphi*d2qua_dtdp
    d2qua_dp2=sp_all(k)*derf2(k)+dphi*d2qua_dp2
  enddo
!
  qua=sp_all(1)+dphi*(sp_all(2)+dphi*qua)
  dqua_dr=dsp_all_ds(1)+dphi*(dsp_all_ds(2)+dphi*dqua_dr)
  dqua_dt=dsp_all_dt(1)+dphi*(dsp_all_dt(2)+dphi*dqua_dt)
!
  d2qua_dr2=d2sp_all_ds2(1)+dphi*(d2sp_all_ds2(2)+dphi*d2qua_dr2)
  d2qua_drdt=d2sp_all_dsdt(1)+dphi*(d2sp_all_dsdt(2)+dphi*d2qua_drdt)
  d2qua_dt2=d2sp_all_dt2(1)+dphi*(d2sp_all_dt2(2)+dphi*d2qua_dt2)
!
  dqua_dp=sp_all(2)+dphi*dqua_dp
  d2qua_drdp=dsp_all_ds(2)+dphi*d2qua_drdp
  d2qua_dtdp=dsp_all_dt(2)+dphi*d2qua_dtdp
!
! End interpolation of mod-B over $\varphi$
!
! Coversion coefficients for derivatives over s
  drhods=0.5d0/rho_tor
  drhods2=drhods**2
  d2rhods2m=drhods2/rho_tor    !$-\dr^2 \rho / \rd s^2$ (second derivative with minus sign)
!
  d2qua_dr2=d2qua_dr2*drhods2-dqua_dr*d2rhods2m
  dqua_dr=dqua_dr*drhods
  d2qua_drdt=d2qua_drdt*drhods
  d2qua_drdp=d2qua_drdp*drhods
!
  Bmod_B=qua
!
  dBmod_B(1)=dqua_dr
  dBmod_B(2)=dqua_dt
  dBmod_B(3)=dqua_dp
!
  d2Bmod_B(1)=d2qua_dr2
  d2Bmod_B(2)=d2qua_drdt
  d2Bmod_B(3)=d2qua_drdp
  d2Bmod_B(4)=d2qua_dt2
  d2Bmod_B(5)=d2qua_dtdp
  d2Bmod_B(6)=d2qua_dp2
!
!--------------------------------
! End Interpolation of mod-B
!--------------------------------
! Interpolation of B_\vartheta and B_\varphi:
!--------------------------------
!
  B_vartheta_B=s_Bcovar_tp_B(1,ns_s_p1,is)
  dB_vartheta_B=B_vartheta_B*derf1(ns_s_p1)
  d2B_vartheta_B=B_vartheta_B*derf2(ns_s_p1)
  B_varphi_B=s_Bcovar_tp_B(2,ns_s_p1,is)
  dB_varphi_B=B_varphi_B*derf1(ns_s_p1)
  d2B_varphi_B=B_varphi_B*derf2(ns_s_p1)
!
  do k=ns_s_B,3,-1
    B_vartheta_B=s_Bcovar_tp_B(1,k,is)+ds*B_vartheta_B
    dB_vartheta_B=s_Bcovar_tp_B(1,k,is)*derf1(k)+ds*dB_vartheta_B
    d2B_vartheta_B=s_Bcovar_tp_B(1,k,is)*derf2(k)+ds*d2B_vartheta_B
    B_varphi_B=s_Bcovar_tp_B(2,k,is)+ds*B_varphi_B
    dB_varphi_B=s_Bcovar_tp_B(2,k,is)*derf1(k)+ds*dB_varphi_B
    d2B_varphi_B=s_Bcovar_tp_B(2,k,is)*derf2(k)+ds*d2B_varphi_B
  enddo
!
  B_vartheta_B=s_Bcovar_tp_B(1,1,is)+ds*(s_Bcovar_tp_B(1,2,is)+ds*B_vartheta_B)
  dB_vartheta_B=s_Bcovar_tp_B(1,2,is)+ds*dB_vartheta_B
  B_varphi_B=s_Bcovar_tp_B(2,1,is)+ds*(s_Bcovar_tp_B(2,2,is)+ds*B_varphi_B)
  dB_varphi_B=s_Bcovar_tp_B(2,2,is)+ds*dB_varphi_B
!
  d2B_vartheta_B=d2B_vartheta_B*drhods2-dB_vartheta_B*d2rhods2m
  d2B_varphi_B=d2B_varphi_B*drhods2-dB_varphi_B*d2rhods2m
  dB_vartheta_B=dB_vartheta_B*drhods
  dB_varphi_B=dB_varphi_B*drhods
!
!--------------------------------
! End interpolation of B_\vartheta and B_\varphi
!--------------------------------
! Interpolation of B_r:
!--------------------------------
! Note that splined quantity is $B_\rho$, not $B_s$
!
  if(use_B_r) then
!
! Begin interpolation of B_rho over $rho$
!
    stp_all(1:nstp,1:nstp)=s_Bcovar_r_B(ns_s_p1,:,:,is,i_theta,i_phi)
    dstp_all_ds(1:nstp,1:nstp)=stp_all(1:nstp,1:nstp)*derf1(ns_s_p1)
    d2stp_all_ds2(1:nstp,1:nstp)=stp_all(1:nstp,1:nstp)*derf2(ns_s_p1)
!
    do k=ns_s_B,3,-1
      stp_all(1:nstp,1:nstp)=s_Bcovar_r_B(k,:,:,is,i_theta,i_phi)+ds*stp_all(1:nstp,1:nstp)
      dstp_all_ds(1:nstp,1:nstp)=s_Bcovar_r_B(k,:,:,is,i_theta,i_phi)*derf1(k)+ds*dstp_all_ds(1:nstp,1:nstp)
      d2stp_all_ds2(1:nstp,1:nstp)=s_Bcovar_r_B(k,:,:,is,i_theta,i_phi)*derf2(k)+ds*d2stp_all_ds2(1:nstp,1:nstp)
    enddo
!
    stp_all(1:nstp,1:nstp)=s_Bcovar_r_B(1,:,:,is,i_theta,i_phi)                                 &
                          +ds*(s_Bcovar_r_B(2,:,:,is,i_theta,i_phi)+ds*stp_all(1:nstp,1:nstp))
    dstp_all_ds(1:nstp,1:nstp)=s_Bcovar_r_B(2,:,:,is,i_theta,i_phi)+ds*dstp_all_ds(1:nstp,1:nstp)
!
! End interpolation of B_rho over $rho$
!-------------------------------
! Begin interpolation of B_rho over $\theta$
!
    sp_all(1:nstp)=stp_all(nstp,1:nstp)
    dsp_all_ds(1:nstp)=dstp_all_ds(nstp,1:nstp)
    d2sp_all_ds2(1:nstp)=d2stp_all_ds2(nstp,1:nstp)
    dsp_all_dt(1:nstp)=sp_all(1:nstp)*derf1(nstp)
    d2sp_all_dsdt(1:nstp)=dsp_all_ds(1:nstp)*derf1(nstp)
    d2sp_all_dt2(1:nstp)=sp_all(1:nstp)*derf2(nstp)
!
    do k=ns_tp_B,3,-1
      sp_all(1:nstp)=stp_all(k,1:nstp)+dtheta*sp_all(1:nstp)
      dsp_all_ds(1:nstp)=dstp_all_ds(k,1:nstp)+dtheta*dsp_all_ds(1:nstp)
      d2sp_all_ds2(1:nstp)=d2stp_all_ds2(k,1:nstp)+dtheta*d2sp_all_ds2(1:nstp)
      dsp_all_dt(1:nstp)=stp_all(k,1:nstp)*derf1(k)+dtheta*dsp_all_dt(1:nstp)
      d2sp_all_dsdt(1:nstp)=dstp_all_ds(k,1:nstp)*derf1(k)+dtheta*d2sp_all_dsdt(1:nstp)
      d2sp_all_dt2(1:nstp)=stp_all(k,1:nstp)*derf2(k)+dtheta*d2sp_all_dt2(1:nstp)
    enddo
!
    sp_all(1:nstp)=stp_all(1,1:nstp)                                                    &
                  +dtheta*(stp_all(2,1:nstp)+dtheta*sp_all(1:nstp))
    dsp_all_ds(1:nstp)=dstp_all_ds(1,1:nstp)                                            &
                      +dtheta*(dstp_all_ds(2,1:nstp)+dtheta*dsp_all_ds(1:nstp))
    d2sp_all_ds2(1:nstp)=d2stp_all_ds2(1,1:nstp)                                        &
                        +dtheta*(d2stp_all_ds2(2,1:nstp)+dtheta*d2sp_all_ds2(1:nstp))
    dsp_all_dt(1:nstp)=stp_all(2,1:nstp)+dtheta*dsp_all_dt(1:nstp)
    d2sp_all_dsdt(1:nstp)=dstp_all_ds(2,1:nstp)+dtheta*d2sp_all_dsdt(1:nstp)
!
! End interpolation of B_rho over $\theta$
!--------------------------------
! Begin interpolation of B_rho over $\varphi$
!
    qua=sp_all(nstp)
    dqua_dr=dsp_all_ds(nstp)
    dqua_dt=dsp_all_dt(nstp)
    dqua_dp=qua*derf1(nstp)
!
    d2qua_dr2=d2sp_all_ds2(nstp)
    d2qua_drdt=d2sp_all_dsdt(nstp)
    d2qua_drdp=dqua_dr*derf1(nstp)
    d2qua_dt2=d2sp_all_dt2(nstp)
    d2qua_dtdp=dqua_dt*derf1(nstp)
    d2qua_dp2=qua*derf2(nstp)
!
    do k=ns_tp_B,3,-1
      qua=sp_all(k)+dphi*qua
      dqua_dr=dsp_all_ds(k)+dphi*dqua_dr
      dqua_dt=dsp_all_dt(k)+dphi*dqua_dt
      dqua_dp=sp_all(k)*derf1(k)+dphi*dqua_dp
!
      d2qua_dr2=d2sp_all_ds2(k)+dphi*d2qua_dr2
      d2qua_drdt=d2sp_all_dsdt(k)+dphi*d2qua_drdt
      d2qua_drdp=dsp_all_ds(k)*derf1(k)+dphi*d2qua_drdp
      d2qua_dt2=d2sp_all_dt2(k)+dphi*d2qua_dt2
      d2qua_dtdp=dsp_all_dt(k)*derf1(k)+dphi*d2qua_dtdp
      d2qua_dp2=sp_all(k)*derf2(k)+dphi*d2qua_dp2
    enddo
!
    qua=sp_all(1)+dphi*(sp_all(2)+dphi*qua)
    dqua_dr=dsp_all_ds(1)+dphi*(dsp_all_ds(2)+dphi*dqua_dr)
    dqua_dt=dsp_all_dt(1)+dphi*(dsp_all_dt(2)+dphi*dqua_dt)
!
    d2qua_dr2=d2sp_all_ds2(1)+dphi*(d2sp_all_ds2(2)+dphi*d2qua_dr2)
    d2qua_drdt=d2sp_all_dsdt(1)+dphi*(d2sp_all_dsdt(2)+dphi*d2qua_drdt)
    d2qua_dt2=d2sp_all_dt2(1)+dphi*(d2sp_all_dt2(2)+dphi*d2qua_dt2)
!
    dqua_dp=sp_all(2)+dphi*dqua_dp
    d2qua_drdp=dsp_all_ds(2)+dphi*d2qua_drdp
    d2qua_dtdp=dsp_all_dt(2)+dphi*d2qua_dtdp
!
! End interpolation of B_rho over $\varphi$
!
! Coversion coefficients for derivatives over s
    drhods=0.5d0/rho_tor
    drhods2=drhods**2
    d2rhods2m=drhods2/rho_tor    !$-\dr^2 \rho / \rd s^2$ (second derivative with minus sign)
!
    d2qua_dr2=d2qua_dr2*drhods2-dqua_dr*d2rhods2m
    dqua_dr=dqua_dr*drhods
    d2qua_drdt=d2qua_drdt*drhods
    d2qua_drdp=d2qua_drdp*drhods
!
    B_r=qua*drhods
!
    dB_r(1)=dqua_dr*drhods-qua*d2rhods2m
    dB_r(2)=dqua_dt*drhods
    dB_r(3)=dqua_dp*drhods
!
    d2B_r(1)=d2qua_dr2*drhods-2.d0*dqua_dr*d2rhods2m+qua*drhods*3.d0/r**2
    d2B_r(2)=d2qua_drdt*drhods-dqua_dt*d2rhods2m
    d2B_r(3)=d2qua_drdp*drhods-dqua_dp*d2rhods2m
    d2B_r(4)=d2qua_dt2*drhods
    d2B_r(5)=d2qua_dtdp*drhods
    d2B_r(6)=d2qua_dp2*drhods
!
  else
    B_r=0.d0
    dB_r=0.d0
    d2B_r=0.d0
  endif
!--------------------------------
! End interpolation of B_r
!--------------------------------
  end subroutine splint_boozer_coord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine delthe_delphi_BV(isw,r,vartheta,varphi,deltheta_BV,delphi_BV, &
                              ddeltheta_BV,ddelphi_BV)
!
! Computes $\Delta \vartheta = \vartheta_B - \theta_V$ and $\Delta \varphi = \varphi_B - \varphi_V$
! and their first derivatives over angles for two cases:
! isw=0 - if they are given as functions of VMEC coordinates (r,vartheta,varphi)
! isw=1 - if they are given as functions of Boozer coordinates (r,vartheta,varphi)
!
  use boozer_coordinates_mod, only : ns_s_B,ns_tp_B,ns_B,n_theta_B,n_phi_B, &
                                     hs_B,h_theta_B,h_phi_B,                &
                                     s_delt_delp_V,s_delt_delp_B,           &
                                     ns_max,derf1,derf2,derf3,              &
                                     use_del_tp_B
  use new_vmec_stuff_mod,   only : nper
  use chamb_mod,            only : rnegflag
!
  implicit none
!
  integer :: isw
  double precision :: r,vartheta,varphi,deltheta_BV,delphi_BV
  double precision, dimension(2) :: ddeltheta_BV,ddelphi_BV
!
  integer, parameter :: n_qua=2
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
  integer :: nstp,ns_A_p1,ns_s_p1
  integer :: k,is,i_theta,i_phi
  integer :: iss,ist,isp
!
  double precision :: s,ds,dtheta,dphi,rho_tor,drhods,drhods2,d2rhods2m
!
!
  double precision, dimension(n_qua)               :: qua,dqua_dt,dqua_dp
  double precision, dimension(n_qua,ns_max)        :: sp_all,dsp_all_dt
  double precision, dimension(n_qua,ns_max,ns_max) :: stp_all
!
  if(r.le.0.d0) then
    rnegflag=.true.
    r=abs(r)
  endif
!
  dtheta=modulo(vartheta,twopi)/h_theta_B
  i_theta=max(0,min(n_theta_B-1,int(dtheta)))
  dtheta=(dtheta-dble(i_theta))*h_theta_B
  i_theta=i_theta+1
!
  dphi=modulo(varphi,twopi/dble(nper))/h_phi_B
  i_phi=max(0,min(n_phi_B-1,int(dphi)))
  dphi=(dphi-dble(i_phi))*h_phi_B
  i_phi=i_phi+1
!
!
!--------------------------------
!
  rho_tor=sqrt(r)
  ds=rho_tor/hs_B
  is=max(0,min(ns_B-1,int(ds)))
  ds=(ds-dble(is))*hs_B
  is=is+1
!
  nstp=ns_tp_B+1
  ns_s_p1=ns_s_B+1
!
!--------------------------------
!
! Begin interpolation of all over $rho$
!
  if(isw.eq.0) then
    stp_all(:,1:nstp,1:nstp)=s_delt_delp_V(:,ns_s_p1,:,:,is,i_theta,i_phi)
!
    do k=ns_s_B,1,-1
      stp_all(:,1:nstp,1:nstp)=s_delt_delp_V(:,k,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp)
    enddo
  elseif(isw.eq.1) then
    if(.not.use_del_tp_B) then
      print *,'delthe_delphi_BV : Boozer data is not loaded'
      return
    endif
    stp_all(:,1:nstp,1:nstp)=s_delt_delp_B(:,ns_s_p1,:,:,is,i_theta,i_phi)
!
    do k=ns_s_B,1,-1
      stp_all(:,1:nstp,1:nstp)=s_delt_delp_B(:,k,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp)
    enddo
  else
    print *,'delthe_delphi_BV : unknown value of switch isw'
    return
  endif
!
! End interpolation of all over $rho$
!-------------------------------
! Begin interpolation of all over $\theta$
!
  sp_all(:,1:nstp)=stp_all(:,nstp,1:nstp)
  dsp_all_dt(:,1:nstp)=sp_all(:,1:nstp)*derf1(nstp)
!
  do k=ns_tp_B,2,-1
    sp_all(:,1:nstp)=stp_all(:,k,1:nstp)+dtheta*sp_all(:,1:nstp)
    dsp_all_dt(:,1:nstp)=stp_all(:,k,1:nstp)*derf1(k)+dtheta*dsp_all_dt(:,1:nstp)
  enddo
!
  sp_all(:,1:nstp)=stp_all(:,1,1:nstp)+dtheta*sp_all(:,1:nstp)
!
! End interpolation of all over $\theta$
!--------------------------------
! Begin interpolation of all over $\varphi$
!
  qua=sp_all(:,nstp)
  dqua_dt=dsp_all_dt(:,nstp)
  dqua_dp=qua*derf1(nstp)
!
  do k=ns_tp_B,2,-1
    qua=sp_all(:,k)+dphi*qua
    dqua_dt=dsp_all_dt(:,k)+dphi*dqua_dt
    dqua_dp=sp_all(:,k)*derf1(k)+dphi*dqua_dp
  enddo
!
  qua=sp_all(:,1)+dphi*qua
  dqua_dt=dsp_all_dt(:,1)+dphi*dqua_dt
!
! End interpolation of all over $\varphi$
!
  deltheta_BV=qua(1)
  delphi_BV=qua(2)
!
  ddeltheta_BV(1)=dqua_dt(1)
  ddelphi_BV(1)=dqua_dt(2)
!
  ddeltheta_BV(2)=dqua_dp(1)
  ddelphi_BV(2)=dqua_dp(2)
!
  end subroutine delthe_delphi_BV
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine vmec_to_boozer(r,theta,varphi,vartheta_B,varphi_B)
!
! Input : r,theta,varphi      - VMEC coordinates
! Output: vartheta_B,varphi_B - Boozer coordinates
!
  use new_vmec_stuff_mod,   only : nper
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  double precision :: r,theta,varphi,vartheta_B,varphi_B
  double precision :: deltheta_BV,delphi_BV
  double precision, dimension(2) :: ddeltheta_BV,ddelphi_BV
!
  call delthe_delphi_BV(0,r,theta,varphi,deltheta_BV,delphi_BV, &
                        ddeltheta_BV,ddelphi_BV)
!
  vartheta_B=modulo(theta+deltheta_BV,twopi)
  varphi_B=modulo(varphi+delphi_BV,twopi/dble(nper))
!
  end subroutine vmec_to_boozer
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine boozer_to_vmec(r,vartheta_B,varphi_B,theta,varphi)
!
! Input : r,vartheta_B,varphi_B - Boozer coordinates
! Output: theta,varphi          - VMEC coordinates
!
  use new_vmec_stuff_mod,   only : nper
  use boozer_coordinates_mod, only : use_del_tp_B
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
  double precision, parameter :: epserr=1.d-14
  integer,          parameter :: niter=100
!
  integer          :: iter
  double precision :: r,theta,varphi,vartheta_B,varphi_B
  double precision :: deltheta_BV,delphi_BV
  double precision :: f1,f2,f11,f12,f21,f22,delthe,delphi,det
  double precision, dimension(2) :: ddeltheta_BV,ddelphi_BV
!
  if(use_del_tp_B) then
!
    call delthe_delphi_BV(1,r,vartheta_B,varphi_B,deltheta_BV,delphi_BV, &
                          ddeltheta_BV,ddelphi_BV)
!
    theta=vartheta_B-deltheta_BV
    varphi=varphi_B-delphi_BV
  else
    theta=vartheta_B
    varphi=varphi_B
  endif
!
! Newton method:
!
  do iter=1,niter
!
    call delthe_delphi_BV(0,r,theta,varphi,deltheta_BV,delphi_BV, &
                          ddeltheta_BV,ddelphi_BV)
!
    f1=theta+deltheta_BV-vartheta_B
    f2=varphi+delphi_BV-varphi_B
    f11=1.d0+ddeltheta_BV(1)
    f12=ddeltheta_BV(2)
    f21=ddelphi_BV(1)
    f22=1.d0+ddelphi_BV(2)
!
    det=f11*f22-f12*f21
    delthe=(f2*f12-f1*f22)/det
    delphi=(f1*f21-f2*f11)/det
!
    theta=theta+delthe
    varphi=varphi+delphi
    if(abs(delthe)+abs(delphi).lt.epserr) exit
  enddo
!
!  theta=modulo(theta,twopi)
!  varphi=modulo(varphi,twopi/dble(nper))
!
  end subroutine boozer_to_vmec
!
