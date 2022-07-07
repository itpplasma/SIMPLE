
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

  subroutine convert_can
    use new_vmec_stuff_mod, only : n_theta,n_phi,h_theta,h_phi
    use vector_potentail_mod, only : ns,hs
    use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,           &
                                          hs_c,h_theta_c,h_phi_c,           &
                                          ns_s_c,ns_tp_c,                   &
                                          nh_stencil,G_c,sqg_c,             &
                                          B_vartheta_c,B_varphi_c
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
      Bcovar_vartheta,Bcovar_varphi,theta

    implicit none

    double precision, parameter :: epserr=1.d-14
    integer :: iter
    double precision :: s,varphi,vartheta,deltheta,A_s,A_theta,A_phi,          &
                        dA_theta_ds,dA_phi_ds,alam,dl_ds,dl_dt,dl_dp,          &
                        Bctrvr_vartheta,Bctrvr_varphi,Bcovar_r

    double precision :: r
    double precision, dimension(2), intent(in) :: y
    double precision, dimension(2), intent(out) :: dy

    s=r**2

    vartheta = vartheta_c
    varphi = varphi_c - y(1)

  !
  ! Begin Newton iteration to find VMEC theta
  !
    theta = vartheta

    do iter=1,100

      call splint_lambda(s,theta,varphi,alam,dl_dt)

      deltheta = (vartheta-theta-alam)/(1.d0+dl_dt)
      theta = theta + deltheta
      if(abs(deltheta).lt.epserr) exit
    enddo
  !
  ! End Newton iteration to find VMEC theta
  !
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

  call convert_can
end program test_new_can
