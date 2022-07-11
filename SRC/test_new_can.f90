
module new_can

  implicit none

  character(*), parameter :: filename = 'wout.nc'

  contains

  subroutine init_own
    use new_vmec_stuff_mod, only : netcdffile, multharm

    netcdffile = filename
    multharm = 3

    call spline_vmec_data
  end subroutine init_own

  subroutine get_canonical_coordinates_new
    use canonical_coordinates_new_mod, only : ns_c,n_theta_c,n_phi_c,       &
                                          hs_c,h_theta_c,h_phi_c,           &
                                          ns_s_c,ns_tp_c,                   &
                                          nh_stencil,                       &
                                          h_vartheta_c,h_varphi_c,          &
                                          A_vartheta_c, A_varphi_c,         &
                                          Delta_varphi_c, chi_gauge, Bmod_c, n_qua
    use vector_potentail_mod, only : ns,hs
    use exchange_get_cancoord_mod, only : vartheta_c,varphi_c, &
                                          Bcovar_vartheta,Bcovar_varphi, &
                                          A_theta,A_phi,onlytheta, &
                                          Bctrvr_vartheta, Bctrvr_varphi
    use new_vmec_stuff_mod, only : n_theta,n_phi,h_theta,h_phi,ns_s,ns_tp

    implicit none

    logical :: fullset
    double precision, parameter :: relerr=1d-10
    integer :: i_theta,i_phi,i_sten,ndim,is_beg
    integer,          dimension(:),     allocatable :: ipoi_t,ipoi_p
    double precision, dimension(:),     allocatable :: y,dy
    double precision :: dstencil_theta(-nh_stencil:nh_stencil), &
                        dstencil_phi(-nh_stencil:nh_stencil)

    double precision :: r,r1,r2,Delta_varphi_c_beg,dchi_dt,dchi_dp, &
      dDelphi_dt, dDelphi_dp
    integer :: is
    integer :: i_ctr ! for nice counting in parallel

    ns_c=ns
    n_theta_c=n_theta
    n_phi_c=n_phi
    h_theta_c=h_theta
    h_phi_c=h_phi
    hs_c=hs

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

    dstencil_phi=dstencil_theta
    dstencil_theta=dstencil_theta/h_theta_c
    dstencil_phi=dstencil_phi/h_phi_c

    allocate(ipoi_t(1-nh_stencil:n_theta_c+nh_stencil))
    allocate(ipoi_p(1-nh_stencil:n_phi_c+nh_stencil))

    do i_theta=1,n_theta_c
      ipoi_t(i_theta)=i_theta
    enddo

    do i_phi=1,n_phi_c
      ipoi_p(i_phi)=i_phi
    enddo

    do i_sten=1,nh_stencil
      ipoi_t(1-i_sten)=ipoi_t(n_theta-i_sten)
      ipoi_t(n_theta_c+i_sten)=ipoi_t(1+i_sten)
      ipoi_p(1-i_sten)=ipoi_p(n_phi_c-i_sten)
      ipoi_p(n_phi_c+i_sten)=ipoi_p(1+i_sten)
    enddo

    allocate(Bmod_c(ns_c,n_theta_c,n_phi_c))
    allocate(h_vartheta_c(ns_c,n_theta_c,n_phi_c))
    allocate(h_varphi_c(ns_c,n_theta_c,n_phi_c))
    allocate(A_vartheta_c(ns_c,n_theta_c,n_phi_c))
    allocate(A_varphi_c(ns_c,n_theta_c,n_phi_c))
    allocate(Delta_varphi_c(ns_c,n_theta_c,n_phi_c))
    allocate(chi_gauge(ns_c,n_theta_c,n_phi_c))

    onlytheta=.false.
    ndim=2
    is_beg=1
    Delta_varphi_c_beg=1d-8

    i_ctr=0
    !$omp parallel private(y, dy, i_theta, i_phi, is, r1, r2, r, &
    !$omp&  dchi_dt, dchi_dp, dDelphi_dt, dDelphi_dp)
    !$omp critical
    allocate(y(ndim),dy(ndim))
    !$omp end critical

    !$omp do
    do i_theta=1,n_theta_c
    !$omp critical
      i_ctr = i_ctr + 1
      print *,'integrate ODE: ',i_ctr,' of ',n_theta_c
    !$omp end critical
      vartheta_c=h_theta_c*dble(i_theta-1)
      do i_phi=1,n_phi_c
        varphi_c=h_phi_c*dble(i_phi-1)

        Delta_varphi_c(is_beg,i_theta,i_phi)=Delta_varphi_c_beg
        y(1)=Delta_varphi_c_beg
        y(2)=0d0

        do is=is_beg-1,2,-1
          r1=hs_c*dble(is)
          r2=hs_c*dble(is-1)

          call odeint_allroutines(y,ndim,r1,r2,relerr,rhs_cancoord_new)

          Delta_varphi_c(is,i_theta,i_phi)=y(1)
          chi_gauge(is,i_theta,i_phi)=y(2)
        enddo

        y(1)=Delta_varphi_c_beg
        y(2)=0d0

        do is=is_beg+1,ns_c
          r1=hs_c*dble(is-2)
          r2=hs_c*dble(is-1)
          if(is.eq.2) r1=1.d-8

          call odeint_allroutines(y,ndim,r1,r2,relerr,rhs_cancoord_new)

          Delta_varphi_c(is,i_theta,i_phi)=y(1)
          chi_gauge(is,i_theta,i_phi)=y(2)
        enddo
      enddo
    enddo
    !$omp end do

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
          do is=2,ns_c
            r=hs_c*dble(is-1)
            y(1)=Delta_varphi_c(is,i_theta,i_phi)

            call rhs_cancoord_new(r,y,dy)

            Bmod_c(is,i_theta,i_phi) = sqrt( &
              abs(Bcovar_vartheta*Bctrvr_vartheta+Bcovar_varphi*Bctrvr_varphi))
            dDelphi_dt=sum(dstencil_theta*chi_gauge( &
              is,ipoi_t(i_theta-nh_stencil:i_theta+nh_stencil),i_phi))
            dDelphi_dp=sum(dstencil_phi*chi_gauge( &
              is,i_theta,ipoi_p(i_phi-nh_stencil:i_phi+nh_stencil)))
            dchi_dt=sum(dstencil_theta*chi_gauge( &
              is,ipoi_t(i_theta-nh_stencil:i_theta+nh_stencil),i_phi))
            dchi_dp=sum(dstencil_phi*chi_gauge( &
              is,i_theta,ipoi_p(i_phi-nh_stencil:i_phi+nh_stencil)))
            h_vartheta_c(is,i_theta,i_phi)=(Bcovar_vartheta &
              + Bcovar_varphi*dDelphi_dt)/Bmod_c(is,i_theta,i_phi)
            h_varphi_c(is,i_theta,i_phi)=(Bcovar_varphi &
              + Bcovar_varphi*dDelphi_dp)/Bmod_c(is,i_theta,i_phi)
            A_vartheta_c(is,i_theta,i_phi)=A_theta &
              + A_phi*dDelphi_dt + dchi_dt
            A_varphi_c(is,i_theta,i_phi)=A_phi &
              + A_phi*dDelphi_dp + dchi_dp
          enddo
          !First point is=1 (on axis) is bad, extrapolate with parabola:
          Bmod_c(1,i_theta,i_phi) = 3.d0*(Bmod_c(2,i_theta,i_phi) &
            -Bmod_c(3,i_theta,i_phi)) + Bmod_c(4,i_theta,i_phi)
          h_vartheta_c(1,i_theta,i_phi) = 0.d0
          h_varphi_c(1,i_theta,i_phi) = 3.d0*(h_varphi_c(2,i_theta,i_phi) &
            -h_varphi_c(3,i_theta,i_phi)) + h_varphi_c(4,i_theta,i_phi)
          A_vartheta_c(1,i_theta,i_phi) = 0.d0
          A_varphi_c(1,i_theta,i_phi) = 3.d0*(A_varphi_c(2,i_theta,i_phi) &
            -A_varphi_c(3,i_theta,i_phi)) + A_varphi_c(4,i_theta,i_phi)
        enddo
      enddo
    !$omp end do
    !$omp critical
    deallocate(y,dy)
    !$omp end critical
    !$omp end parallel

    ns_s_c=ns_s
    ns_tp_c=ns_tp
    fullset=.true.

    onlytheta=.true.

    call spline_can_coord_new(fullset)
    deallocate(ipoi_t,ipoi_p,Bmod_c,h_vartheta_c,h_varphi_c,&
      A_vartheta_c,A_varphi_c,Delta_varphi_c,chi_gauge)

  end subroutine get_canonical_coordinates_new


  subroutine spline_can_coord_new(fullset)
    use canonical_coordinates_new_mod, only : ns_c,n_theta_c,n_phi_c,hs_c,h_theta_c,h_phi_c,    &
                                          ns_s_c,ns_tp_c,h_vartheta_c,h_varphi_c, &
                                          ns_max,derf1,derf2,derf3, s_Bmod_B_A, s_Delta_varphi_c, &
                                          A_vartheta_c, A_varphi_c, Delta_varphi_c, Bmod_c, n_qua

    implicit none

    logical :: fullset
    integer :: k,is,i_theta,i_phi,i_qua
    integer :: ist,isp
    double precision, dimension(:,:), allocatable :: splcoe

    if (.not. allocated(s_Bmod_B_A)) &
      allocate(s_Bmod_B_A(n_qua,ns_s_c+1,ns_tp_c+1,ns_tp_c+1,ns_c,n_theta_c,n_phi_c))
    if(fullset .and. (.not. allocated(s_Delta_varphi_c))) &
      allocate(s_Delta_varphi_c(ns_s_c+1,ns_tp_c+1,ns_tp_c+1,ns_c,n_theta_c,n_phi_c))

    s_Bmod_B_A(1,1,1,1,:,:,:)=Bmod_c
    s_Bmod_B_A(2,1,1,1,:,:,:)=h_vartheta_c
    s_Bmod_B_A(3,1,1,1,:,:,:)=h_varphi_c
    s_Bmod_B_A(4,1,1,1,:,:,:)=A_vartheta_c
    s_Bmod_B_A(5,1,1,1,:,:,:)=A_varphi_c
    if(fullset) s_Delta_varphi_c(1,1,1,:,:,:)=Delta_varphi_c

    !
    ! splining over $\varphi$:
    !
    allocate(splcoe(0:ns_tp_c,n_phi_c))

    do is=1,ns_c
      do i_theta=1,n_theta_c
        do i_qua=1,n_qua

          splcoe(0,:)=s_Bmod_B_A(i_qua,1,1,1,is,i_theta,:)

          call spl_per(ns_tp_c,n_phi_c,h_phi_c,splcoe)

          do k=1,ns_tp_c
            s_Bmod_B_A(i_qua,1,1,k+1,is,i_theta,:)=splcoe(k,:)
          enddo

        enddo

        if(fullset) then

          splcoe(0,:)=s_Delta_varphi_c(1,1,1,is,i_theta,:)

          call spl_per(ns_tp_c,n_phi_c,h_phi_c,splcoe)

          do k=1,ns_tp_c
            s_Delta_varphi_c(1,1,k+1,is,i_theta,:)=splcoe(k,:)
          enddo

        endif
      enddo
    enddo

    deallocate(splcoe)
    !
    ! splining over $\vartheta$:
    !
    allocate(splcoe(0:ns_tp_c,n_theta_c))

    do is=1,ns_c
      do i_phi=1,n_phi_c
        do isp=1,ns_tp_c+1
          do i_qua=1,n_qua

            splcoe(0,:)=s_Bmod_B_A(i_qua,1,1,isp,is,:,i_phi)

            call spl_per(ns_tp_c,n_theta_c,h_theta_c,splcoe)

            do k=1,ns_tp_c
              s_Bmod_B_A(i_qua,1,k+1,isp,is,:,i_phi)=splcoe(k,:)
            enddo

          enddo

          if(fullset) then

            splcoe(0,:)=s_Delta_varphi_c(1,1,isp,is,:,i_phi)

            call spl_per(ns_tp_c,n_theta_c,h_theta_c,splcoe)

            do k=1,ns_tp_c
              s_Delta_varphi_c(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
            enddo

          endif

        enddo
      enddo
    enddo

    deallocate(splcoe)
    !
    ! splining over $s$:
    !
    allocate(splcoe(0:ns_s_c,ns_c))

    do i_theta=1,n_theta_c
      do i_phi=1,n_phi_c
        do ist=1,ns_tp_c+1
          do isp=1,ns_tp_c+1
            do i_qua=1,n_qua

              splcoe(0,:)=s_Bmod_B_A(i_qua,1,ist,isp,:,i_theta,i_phi)

              call spl_reg(ns_s_c,ns_c,hs_c,splcoe)

              do k=1,ns_s_c
                s_Bmod_B_A(i_qua,k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
              enddo

            enddo

            if(fullset) then

              splcoe(0,:)=s_Delta_varphi_c(1,ist,isp,:,i_theta,i_phi)

              call spl_reg(ns_s_c,ns_c,hs_c,splcoe)

              do k=1,ns_s_c
                s_Delta_varphi_c(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
              enddo

            endif

          enddo
        enddo
      enddo
    enddo

    deallocate(splcoe)

    do k=1,ns_max
      derf1(k)=dble(k-1)
      derf2(k)=dble((k-1)*(k-2))
      derf3(k)=dble((k-1)*(k-2)*(k-3))
    enddo

  end subroutine spline_can_coord_new


  subroutine splint_can_coord_new(fullset,mode_secders,r,vartheta_c,varphi_c,f, &
      Delta_varphi_c)

    use canonical_coordinates_new_mod, only : ns_c,n_theta_c,n_phi_c,hs_c, &
                  h_theta_c,h_phi_c,ns_s_c,ns_tp_c,ns_max,n_qua, &
                  derf1,derf2,derf3,s_Bmod_B_A,s_Delta_varphi_c
    use vector_potentail_mod, only : ns,hs,torflux,sA_phi
    use new_vmec_stuff_mod,   only : nper,ns_A
    use chamb_mod,            only : rnegflag
    use diag_mod, only : icounter
    use field_can_mod, only : FieldCan

    implicit none

    double precision, intent(inout) :: r
    double precision, intent(in) :: vartheta_c,varphi_c
    type(FieldCan), intent(inout) :: f
    double precision, intent(out) :: Delta_varphi_c

    double precision, parameter :: twopi=2.d0*3.14159265358979d0

    logical :: fullset

    integer :: mode_secders,nstp,ns_A_p1,ns_s_p1
    integer :: k,is,i_theta,i_phi
    integer :: iss,ist,isp

    double precision :: s,ds,dtheta,dphi,rho_tor,drhods,drhods2,d2rhods2m

    double precision, dimension(ns_max)              :: sp_G
    double precision, dimension(ns_max,ns_max)       :: stp_G

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

    dtheta=modulo(vartheta_c,twopi)/h_theta_c
    i_theta=max(0,min(n_theta_c-1,int(dtheta)))
    dtheta=(dtheta-dble(i_theta))*h_theta_c
    i_theta=i_theta+1

    dphi=modulo(varphi_c,twopi/dble(nper))/h_phi_c
    i_phi=max(0,min(n_phi_c-1,int(dphi)))
    dphi=(dphi-dble(i_phi))*h_phi_c
    i_phi=i_phi+1

    rho_tor=sqrt(r)
    ds=rho_tor/hs_c
    is=max(0,min(ns_c-1,int(ds)))
    ds=(ds-dble(is))*hs_c
    is=is+1
    nstp=ns_tp_c+1

    if(fullset) then
    !
    ! Begin interpolation of Delta_phi over $s$
    !
    stp_G(1:nstp,1:nstp)=s_Delta_varphi_c(ns_s_c+1,:,:,is,i_theta,i_phi)
    !
    do k=ns_s_c,1,-1
      stp_G(1:nstp,1:nstp)=s_Delta_varphi_c(k,:,:,is,i_theta,i_phi)+ds*stp_G(1:nstp,1:nstp)
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
      Delta_varphi_c=sp_G(nstp)
    !
    do k=ns_tp_c,1,-1
      Delta_varphi_c=sp_G(k)+dphi*Delta_varphi_c
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
      stp_all(:,1:nstp,1:nstp)=s_Bmod_B_A(:,ns_s_p1,:,:,is,i_theta,i_phi)
      dstp_all_ds(:,1:nstp,1:nstp)=stp_all(:,1:nstp,1:nstp)*derf1(ns_s_p1)
      d2stp_all_ds2(:,1:nstp,1:nstp)=stp_all(:,1:nstp,1:nstp)*derf2(ns_s_p1)
      !
      do k=ns_s_c,3,-1
        stp_all(:,1:nstp,1:nstp)=s_Bmod_B_A(:,k,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp)
        dstp_all_ds(:,1:nstp,1:nstp)=s_Bmod_B_A(:,k,:,:,is,i_theta,i_phi)*derf1(k)+ds*dstp_all_ds(:,1:nstp,1:nstp)
        d2stp_all_ds2(:,1:nstp,1:nstp)=s_Bmod_B_A(:,k,:,:,is,i_theta,i_phi)*derf2(k)+ds*d2stp_all_ds2(:,1:nstp,1:nstp)
      enddo
      !
      stp_all(:,1:nstp,1:nstp)=s_Bmod_B_A(:,1,:,:,is,i_theta,i_phi)                                 &
        +ds*(s_Bmod_B_A(:,2,:,:,is,i_theta,i_phi)+ds*stp_all(:,1:nstp,1:nstp))
      dstp_all_ds(:,1:nstp,1:nstp)=s_Bmod_B_A(:,2,:,:,is,i_theta,i_phi)+ds*dstp_all_ds(:,1:nstp,1:nstp)
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

      d2qua_dr2=d2qua_dr2*drhods2-dqua_dr*d2rhods2m
      dqua_dr=dqua_dr*drhods
      d2qua_drdt=d2qua_drdt*drhods
      d2qua_drdp=d2qua_drdp*drhods

      f%Bmod=qua(1)
      f%hth=qua(2)
      f%hph=qua(3)
      f%Ath=qua(4)
      f%Aph=qua(5)

      f%dBmod(1)=dqua_dr(1)
      f%dhth(1)=dqua_dr(2)
      f%dhph(1)=dqua_dr(3)
      f%dAth(1)=dqua_dr(4)
      f%dAph(1)=dqua_dr(5)

      f%dBmod(2)=dqua_dt(1)
      f%dhth(2)=dqua_dt(2)
      f%dhph(2)=dqua_dt(3)
      f%dAth(2)=dqua_dt(4)
      f%dAph(2)=dqua_dt(5)

      f%dBmod(3)=dqua_dp(1)
      f%dhth(3)=dqua_dp(2)
      f%dhph(3)=dqua_dp(3)
      f%dAth(3)=dqua_dp(4)
      f%dAph(3)=dqua_dp(5)

      f%d2Bmod(1)=d2qua_dr2(1)
      f%d2hth(1)=d2qua_dr2(2)
      f%d2hph(1)=d2qua_dr2(3)
      f%d2Ath(1)=d2qua_dr2(4)
      f%d2Aph(1)=d2qua_dr2(5)

      f%d2Bmod(2)=d2qua_drdt(1)
      f%d2hth(2)=d2qua_drdt(2)
      f%d2hph(2)=d2qua_drdt(3)
      f%d2Ath(2)=d2qua_drdt(4)
      f%d2Aph(2)=d2qua_drdt(5)

      f%d2Bmod(3)=d2qua_drdp(1)
      f%d2hth(3)=d2qua_drdp(2)
      f%d2hph(3)=d2qua_drdp(3)
      f%d2Ath(3)=d2qua_drdp(4)
      f%d2Aph(3)=d2qua_drdp(5)

      f%d2Bmod(4)=d2qua_dt2(1)
      f%d2hth(4)=d2qua_dt2(2)
      f%d2hph(4)=d2qua_dt2(3)
      f%d2Ath(4)=d2qua_dt2(4)
      f%d2Aph(4)=d2qua_dt2(5)

      f%d2Bmod(5)=d2qua_dtdp(1)
      f%d2hth(5)=d2qua_dtdp(2)
      f%d2hph(5)=d2qua_dtdp(3)
      f%d2Ath(5)=d2qua_dtdp(4)
      f%d2Aph(5)=d2qua_dtdp(5)

      f%d2Bmod(6)=d2qua_dp2(1)
      f%d2hth(6)=d2qua_dp2(2)
      f%d2hph(6)=d2qua_dp2(3)
      f%d2Ath(6)=d2qua_dp2(4)
      f%d2Aph(6)=d2qua_dp2(5)
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
  end subroutine splint_can_coord_new


  subroutine convert_can
    use new_vmec_stuff_mod, only : n_theta,n_phi,h_theta,h_phi
    use vector_potentail_mod, only : ns,hs
    use canonical_coordinates_new_mod, only : ns_c,n_theta_c,n_phi_c,           &
                                          hs_c,h_theta_c,h_phi_c,           &
                                          nh_stencil
    use exchange_get_cancoord_mod, only : vartheta_c,varphi_c
    use util, only : twopi
    implicit none

    integer, parameter :: ntheta = 100, nvarphi = 100
    double precision, parameter :: relerr=1d-10

    integer :: i_theta
    double precision :: r1, r2
    double precision :: y(2)

    call init_own

    ns_c=ns
    n_theta_c=n_theta
    n_phi_c=n_phi
    h_theta_c=h_theta
    h_phi_c=h_phi
    hs_c=hs

    varphi_c = 0.2d0

    do i_theta = 1,ntheta
      vartheta_c = i_theta*twopi/ntheta
      y = 0d0

      r1 = 0.1d0
      r2 = 0.9d0
      call odeint_allroutines(y,2,r1,r2,relerr,rhs_cancoord_new)
      print *, vartheta_c, y(2)
    end do
  end subroutine convert_can

  subroutine rhs_cancoord_new(r, y, dy)
    ! Returns
    ! dy_1/ds = d\Delta\varphi/ds   ... toroidal angle shift
    ! dy_2/ds = d\chi/ds            ... gauge potential

    use exchange_get_cancoord_mod, only : vartheta_c,varphi_c,sqg,aiota,       &
      Bcovar_vartheta,Bcovar_varphi,theta,A_theta,A_phi,Bctrvr_vartheta,Bctrvr_varphi, &
      onlytheta

    implicit none

    double precision, parameter :: epserr=1.d-14
    integer :: iter
    double precision :: s,varphi,vartheta,deltheta,A_s,                        &
                        dA_theta_ds,dA_phi_ds,alam,dl_ds,dl_dt,dl_dp,          &
                        Bcovar_r

    double precision :: r
    double precision, dimension(2), intent(in) :: y
    double precision, dimension(2), intent(out) :: dy

    s=r**2

    vartheta = vartheta_c
    varphi = varphi_c - y(1)


  ! Begin Newton iteration to find VMEC theta

    theta = vartheta

    do iter=1,100

      call splint_lambda(s,theta,varphi,alam,dl_dt)

      deltheta = (vartheta-theta-alam)/(1.d0+dl_dt)
      theta = theta + deltheta
      if(abs(deltheta).lt.epserr) exit
    enddo

    if(onlytheta) return

  ! End Newton iteration to find VMEC theta

    call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,  &
                    sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,  &
                    Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
    A_s = 0d0  ! No covariant radial component for VMEC

    dy(1) = -Bcovar_r/Bcovar_varphi
    dy(2) = -0d0 - dy(1)*A_phi

  end subroutine rhs_cancoord_new

end module new_can

program test_new_can
  use new_can
  implicit none

  ! call convert_can
  call init_own
  call get_canonical_coordinates_new
end program test_new_can
