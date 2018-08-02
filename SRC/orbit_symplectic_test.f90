program orbit_symplectic_test

use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl, f_sympl_euler, yold,&
  rmumag, ro0

implicit none

integer, parameter :: n = 6

double precision :: dt0
double precision, dimension(6) :: z, fvec
integer :: info

double precision :: r,vartheta_c,varphi_c,                                           &
                    A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,                 &
                    sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                    B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                    B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c,     &
                    B_mod, dlnB_mod_dr, dlnB_mod_dt, dlnB_mod_dp, dA_phi_dt

integer :: k

rmumag = 0.1d0
ro0 = 1d0
dt0 = 1.0d0*dsqrt(2d0)

yold = 0d0
yold(1) = 0.3d0
yold(2) = 1.5d0
yold(4) = 0.1d0

call splint_can_extra(yold(1),yold(2),yold(3),                                           &
                      A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                      .false., G_c, B_mod, dlnB_mod_dr, dlnB_mod_dt, dlnB_mod_dp, dA_phi_dt)

yold(5) = yold(4)*B_vartheta_c/B_mod + A_theta/ro0
yold(6) = yold(4)*B_varphi_c/B_mod + A_phi/ro0

! call f_sympl_euler(n,yold,fvec,info)
! print *, yold
! print *, fvec

z = 0d0
z(1:3) = yold(1:3)
z(4) = yold(4)**2/2d0 + rmumag*B_mod
write(4001,*) z

do k = 1, 10000
  call orbit_timestep_sympl(z, 10*dt0, dt0, info)
  write(4001,*) z
end do

end program orbit_symplectic_test