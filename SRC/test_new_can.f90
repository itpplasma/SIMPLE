
module new_can
  implicit none

  character(*), parameter :: filename = 'wout.nc'

  contains

  subroutine init_own
    use new_vmec_stuff_mod, only : netcdffile, multharm

    netcdffile = filename
    multharm = 7

    call spline_vmec_data
  end subroutine init_own

  subroutine get_canonical_coordinates
    use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,           &
                                          hs_c,h_theta_c,h_phi_c,           &
                                          ns_s_c,ns_tp_c,                   &
                                          nh_stencil,sqg_c,             &
                                          B_vartheta_c,B_varphi_c,          &
                                          A_vartheta_c, A_varphi_c,         &
                                          Delta_varphi_c, chi_gauge
    use vector_potentail_mod, only : ns,hs
    use exchange_get_cancoord_mod, only : vartheta_c,varphi_c,sqg, &
                                          Bcovar_vartheta,Bcovar_varphi, &
                                          A_theta,A_phi,onlytheta
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

    allocate(sqg_c(ns_c,n_theta_c,n_phi_c))
    allocate(B_vartheta_c(ns_c,n_theta_c,n_phi_c))
    allocate(B_varphi_c(ns_c,n_theta_c,n_phi_c))
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

            dDelphi_dt=sum(dstencil_theta*chi_gauge( &
              is,ipoi_t(i_theta-nh_stencil:i_theta+nh_stencil),i_phi))
            dDelphi_dp=sum(dstencil_phi*chi_gauge( &
              is,i_theta,ipoi_p(i_phi-nh_stencil:i_phi+nh_stencil)))
            dchi_dt=sum(dstencil_theta*chi_gauge( &
              is,ipoi_t(i_theta-nh_stencil:i_theta+nh_stencil),i_phi))
            dchi_dp=sum(dstencil_phi*chi_gauge( &
              is,i_theta,ipoi_p(i_phi-nh_stencil:i_phi+nh_stencil)))
            sqg_c(is,i_theta,i_phi)=sqg*(1.d0+dDelphi_dp)
            B_vartheta_c(is,i_theta,i_phi)=Bcovar_vartheta &
              + Bcovar_varphi*dDelphi_dt
            B_varphi_c(is,i_theta,i_phi)=Bcovar_varphi &
              + Bcovar_varphi*dDelphi_dp
            A_vartheta_c(is,i_theta,i_phi)=A_theta &
              + A_phi*dDelphi_dt + dchi_dt
            A_varphi_c(is,i_theta,i_phi)=A_phi &
              + A_phi*dDelphi_dp + dchi_dp
          enddo
          !First point is=1 (on axis) is bad, extrapolate with parabola:
          sqg_c(1,i_theta,i_phi) = 3.d0*(sqg_c(2,i_theta,i_phi) &
            -sqg_c(3,i_theta,i_phi)) + sqg_c(4,i_theta,i_phi)
          B_vartheta_c(1,i_theta,i_phi) = 0.d0
          B_varphi_c(1,i_theta,i_phi) = 3.d0*(B_varphi_c(2,i_theta,i_phi) &
            -B_varphi_c(3,i_theta,i_phi)) + B_varphi_c(4,i_theta,i_phi)
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

    !call spline_can_coord(fullset)
    deallocate(ipoi_t,ipoi_p,sqg_c,B_vartheta_c,B_varphi_c,&
      A_vartheta_c,A_varphi_c,Delta_varphi_c,chi_gauge)

  end subroutine get_canonical_coordinates

  subroutine convert_can
    use new_vmec_stuff_mod, only : n_theta,n_phi,h_theta,h_phi
    use vector_potentail_mod, only : ns,hs
    use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,           &
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
      Bcovar_vartheta,Bcovar_varphi,theta,A_theta,A_phi

    implicit none

    double precision, parameter :: epserr=1.d-14
    integer :: iter
    double precision :: s,varphi,vartheta,deltheta,A_s,                        &
                        dA_theta_ds,dA_phi_ds,alam,dl_ds,dl_dt,dl_dp,          &
                        Bctrvr_vartheta,Bctrvr_varphi,Bcovar_r

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
  call spline_vmec_data
  call get_canonical_coordinates
end program test_new_can
