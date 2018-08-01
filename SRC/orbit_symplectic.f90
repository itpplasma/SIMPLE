module orbit_symplectic

use parmot_mod, only : rmu, ro0

implicit none

double precision, dimension(6) :: yold
double precision :: rmumag
double precision :: dt

contains

subroutine orbit_sympl_init(z)

  implicit none
  
  double precision, intent(in) :: z(5)
  double precision :: r,vartheta_c,varphi_c,                                           &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,                 &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c

  double precision :: Bctr_vartheta, Bctr_varphi, bmod, bmod2, coala
  double precision :: derphi(3)
  double precision :: p, alambd

  call splint_can_coord(z(1),z(2),z(3),                                           &
                        A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                        sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                        B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                        B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                        .false., G_c)

  Bctr_vartheta=-dA_phi_dr/sqg_c
  Bctr_varphi=dA_theta_dr/sqg_c
  bmod2=Bctr_vartheta*B_vartheta_c+Bctr_varphi*B_varphi_c
  bmod=sqrt(bmod2)

  rmumag=.5d0*p**2*(1.d0-alambd**2)/bmod

  p=z(4)
  alambd=z(5)
  
  yold(1:3) = z(1:3) ! r, theta, varphi
  yold(4) = p*alambd/dsqrt(p**2*2.d0/rmu+1.d0) ! vpar_bar = vpar/sqrt(2T/m)
  yold(5) = yold(4)*B_vartheta_c/bmod + A_theta/ro0
  yold(6) = yold(4)*B_varphi_c/bmod + A_phi/ro0

  print *, 'Symplectic old variables: '
  print *, yold

end subroutine orbit_sympl_init

subroutine f_sympl_euler(n, y, fvec, iflag)

  implicit none

  integer, parameter :: mode = 1

  integer, intent(in) :: n
  double precision, intent(in) :: y(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(out) :: iflag

  double precision, dimension(2) :: q, p, w
  double precision, dimension(2) :: dqdt, dpdt, dHdq, dHdw, pqw
  double precision, dimension(2,2) :: dpdq, dpdw, dwdq, dwdp

  double precision :: r,vartheta_c,varphi_c,                                           &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,                 &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c

  double precision :: Bctr_vartheta, Bctr_varphi, bmod2, bmod, coala
  double precision :: bder(3), derphi(3)

  if (mode==1) then
    q = y(2:3)
    p = yold(5:6)
    w = y((/1,4/))
  else if (mode==2) then
    q = y(2:3)
    p = yold(5:6)
    w = y((/1,4/))
  end if

  ! factor 2.0 due to phase space volume factor by normalised variables
  dqdt = (y(2:3)-yold(2:3))/dt*2.0
  dpdt = (y(5:6)-yold(5:6))/dt*2.0

  call splint_can_coord(y(1),y(2),y(3),                                           &
                        A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                        sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                        B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                        B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                        .false., G_c)

  Bctr_vartheta=-dA_phi_dr/sqg_c
  Bctr_varphi=dA_theta_dr/sqg_c
  bmod2=Bctr_vartheta*B_vartheta_c+Bctr_varphi*B_varphi_c
  bmod=sqrt(bmod2)
  bder(1)=0.5d0*((dA_theta_dr*dB_varphi_c_dr-dA_phi_dr*dB_vartheta_c_dr-d2A_phi_dr2*B_vartheta_c) &
         /bmod2-dsqg_c_dr)/sqg_c
  bder(2)=0.5d0*((dA_theta_dr*dB_varphi_c_dt-dA_phi_dr*dB_vartheta_c_dt)/bmod2-dsqg_c_dt)/sqg_c
  bder(3)=0.5d0*((dA_theta_dr*dB_varphi_c_dp-dA_phi_dr*dB_vartheta_c_dp)/bmod2-dsqg_c_dp)/sqg_c
  
  call elefie_can(y(1:3),derphi)

  dHdq = rmumag*bmod*bder(2:3) + derphi(2:3)
  dHdw(1) = rmumag*bmod*bder(1) + derphi(1)
  dHdw(2) = 2*w(2) 
  pqw(1) = w(2)*B_vartheta_c/bmod + A_theta/ro0
  pqw(2) = w(2)*B_varphi_c/bmod + A_phi/ro0

  dpdq(1,1) = w(2)*(dB_vartheta_c_dt - B_vartheta_c*bder(2))/bmod
  dpdq(1,2) = w(2)*(dB_vartheta_c_dp - B_vartheta_c*bder(3))/bmod
  dpdq(2,1) = w(2)*(dB_varphi_c_dt - B_varphi_c*bder(2))/bmod
  dpdq(2,2) = w(2)*(dB_varphi_c_dp - B_varphi_c*bder(3))/bmod

  dpdw(1,1) = w(2)*(dB_vartheta_c_dr - B_vartheta_c*bder(1))/bmod + dA_theta_dr/ro0
  dpdw(1,2) = B_vartheta_c/bmod
  dpdw(2,1) = w(2)*(dB_varphi_c_dr- B_varphi_c*bder(1))/bmod + dA_phi_dr/ro0
  dpdw(2,2) = B_varphi_c/bmod

  ! inverse matrix
  dwdp(1,1) = dpdw(2,2)
  dwdp(1,2) = -dpdw(1,2)
  dwdp(2,1) = -dpdw(2,1)
  dwdp(2,2) = dpdw(1,1)
  dwdp = dwdp/(dpdw(1,1)*dpdw(2,2)-dpdw(1,2)*dpdw(2,1))

  ! matrix multiplication
  dwdq = -matmul(dwdp,dpdq)

  call step_sympl(q, p, w, dqdt, dpdt, dHdq, dHdw, pqw, dwdq, dwdp, fvec)
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

  ret = 0.0
  do k = 1, nq
    ret(k) = sum(dHdw*dwdp(:,k)) - dqdt(k)
    ret(k + nq) = -sum(dHdw*dwdq(:,k)) - dHdq(k) - dpdt(k)
    ret(k + 2*nq) = p(k) - pqw(k)
  end do

end subroutine step_sympl

subroutine orbit_timestep_sympl(z, dtau, dtaumin, ierr)
  implicit none
  
  integer, parameter :: ndim = 5
  integer, parameter :: n = 6
  double precision, parameter :: tol = 1e-10

  integer, intent(out) :: ierr
  double precision, dimension(ndim), intent(inout) :: z
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  double precision, dimension(6) :: y, fvec
  integer :: info

  double precision tau2

  dt = dtau
  tau2 = 0.0
  do while(tau2.lt.dtau)
    y = yold
    call hybrd1 (f_sympl_euler, n, y, fvec, tol, info)
  !  print *, info
  !  print *, y
  !  print *, fvec
    yold = y
    z = 0.0
    z(1:3) = y(1:3)
    tau2=tau2+dtaumin
  enddo
end


end module orbit_symplectic
