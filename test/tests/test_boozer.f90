!
  program test_boozer
!
  use simple, only : Tracer
  use simple_main, only : init_field
  use boozer_coordinates_mod, only : use_B_r, use_del_tp_B
  use new_vmec_stuff_mod,   only : nper
  use boozer_sub, only : splint_boozer_coord, boozer_converter, &
    delthe_delphi_bv, vmec_to_boozer, boozer_to_vmec
!
  implicit none
!
  type(Tracer) :: norb
!
  double precision :: r,vartheta_B,varphi_B,                                       &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3, &
                      B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                      B_varphi_B,dB_varphi_B,d2B_varphi_B,Bmod_B,B_r
  double precision, dimension(3) :: dBmod_B,dB_r
  double precision, dimension(6) :: d2Bmod_B,d2B_r
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
  integer, parameter :: n=100
  integer :: it,ip,isw,i
  double precision, dimension(0:n) :: arr,arr1,arr2,arr3
  double precision :: deltheta_BV,delphi_BV,theta,varphi,eps,xi
  double precision, dimension(2) :: ddeltheta_BV,ddelphi_BV
!
  double precision :: dA_phi_dr_num,d2A_phi_dr2_num,        &
                      dB_vartheta_B_num,d2B_vartheta_B_num, &
                      dB_varphi_B_num,d2B_varphi_B_num
  double precision :: A_phi_p,A_phi_m
  double precision :: B_vartheta_B_p,B_vartheta_B_m
  double precision :: B_varphi_B_p,B_varphi_B_m
  double precision :: Bmod_B_p,Bmod_B_m,Bmod_B_pp,Bmod_B_mm,Bmod_B_pm,Bmod_B_mp
  double precision :: B_r_p,B_r_m,B_r_pp,B_r_mm,B_r_pm,B_r_mp
  double precision, dimension(3) :: dBmod_B_num,dB_r_num
  double precision, dimension(6) :: d2Bmod_B_num,d2B_r_num
!
!
!
  call init_field(norb, 'wout.nc', 5, 5, 5, -1)
!
!  use_B_r=.false.
  use_B_r=.true.        !important for plotting, not needed for canonical coordinates
  use_del_tp_B=.false.  !not very important (saves 1 iteration in coord. transform), don't use to save memory
!  use_del_tp_B=.true.
!
  call get_boozer_coordinates
!
!---------------------------------------------------------
!
! Plot mod-B, B_r, theta_B-theta_V and phi_B-phi_V on a given flux surface:
!

  print *, "enter radius s"
  read(*,*) r

  print *, "r = ", r

  open(1234,file='angle.dat')
  open(12345,file='Bmod.dat')
  open(123456,file='B_r.dat')
  open(21345,file='deltheta.dat')
  open(213456,file='delphi.dat')
!
  do it=0,n
    vartheta_B=dble(it)*twopi/dble(n)
    do ip=0,n
      varphi_B=dble(ip)*twopi/dble(n*nper)
!
      call splint_boozer_coord(r,vartheta_B,varphi_B,                                       &
                               A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                               B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                               B_varphi_B,dB_varphi_B,d2B_varphi_B,                         &
                               Bmod_B,dBmod_B,d2Bmod_B,                                     &
                               B_r,dB_r,d2B_r)
!
      arr(ip)=Bmod_B
      arr1(ip)=B_r
!
      isw=0
      call delthe_delphi_BV(isw,r,vartheta_B,varphi_B,deltheta_BV,delphi_BV, &
                              ddeltheta_BV,ddelphi_BV)
!
      arr2(ip)=deltheta_BV
      arr3(ip)=delphi_BV
    enddo
    write(1234,*) vartheta_B
    write(12345,*) arr
    write(123456,*) arr1
    write(21345,*) arr2
    write(213456,*) arr3
  enddo
!
  close(1234)
  close(12345)
  close(123456)
  close(21345)
  close(213456)
!
!---------------------------------------------------------
!
! Plot B_r, B_theta, B_phi in (s,theta) plane at half period
!
  open(1234,file='B_theta_phi_vs_r.dat')
  open(12345,file='B_r_st_plane.dat')
  varphi_B=0.5d0*twopi/dble(nper)
  do i=1,300
    r=dble(i)/300.d0
    do it=1,n
      vartheta_B=dble(it)*twopi/dble(n)
!
      call splint_boozer_coord(r,vartheta_B,varphi_B,                                       &
                               A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                               B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                               B_varphi_B,dB_varphi_B,d2B_varphi_B,                         &
                               Bmod_B,dBmod_B,d2Bmod_B,                                     &
                               B_r,dB_r,d2B_r)
!
      arr(it)=B_r
    enddo
    write(1234,*) r,B_vartheta_B,B_varphi_B
    write(12345,*) arr
  enddo
  close(1234)
  close(12345)
!
!---------------------------------------------------------
!
! Test forward-backward convertion of VMEC to Boozer angles:
!
  r=0.5d0
  theta=3.2602972560067962d0
  varphi=0.62918037430259344d0
!
  open(1234,file='forw_backw.dat')
  write(1234,*) theta,varphi
!
  do i=1,100
!
    call vmec_to_boozer(r,theta,varphi,vartheta_B,varphi_B)
!
    write(1234,*) vartheta_B,varphi_B
!
    call boozer_to_vmec(r,vartheta_B,varphi_B,theta,varphi)
!
    write(1234,*) theta,varphi
  enddo
!
  close(1234)
!
!
!---------------------------------------------------------
!
! Test derivatives
!
  open(1234,file='derivatives_spline.dat')
  open(12345,file='derivatives_findif.dat')
!
  eps=1.d-3
!
  do i=0,1000
    xi=1.d-3*dble(i)
    r=0.4d0+0.2d0*xi
    vartheta_B=0.1d0+twopi*xi
    varphi_B=0.1d0-twopi*xi
!
    call splint_boozer_coord(r,vartheta_B,varphi_B,                                       &
                             A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B,dBmod_B,d2Bmod_B,                                     &
                             B_r,dB_r,d2B_r)
!
    write(1234,*) xi,dA_phi_dr,d2A_phi_dr2,     &
                  dB_vartheta_B,d2B_vartheta_B, &
                  dB_varphi_B,d2B_varphi_B,     &
                  dBmod_B,d2Bmod_B,dB_r,d2B_r
!
    call splint_boozer_coord(r+eps,vartheta_B,varphi_B,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_p,dBmod_B,d2Bmod_B,                                     &
                             B_r_p,dB_r,d2B_r)
!
    call splint_boozer_coord(r-eps,vartheta_B,varphi_B,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_m,dBmod_B,d2Bmod_B,                                     &
                             B_r_m,dB_r,d2B_r)
!
    dA_phi_dr_num=(A_phi_p-A_phi_m)/(2d0*eps)
    d2A_phi_dr2_num=(A_phi_p+A_phi_m-2.d0*A_phi)/eps**2
    dB_vartheta_B_num=(B_vartheta_B_p-B_vartheta_B_m)/(2d0*eps)
    d2B_vartheta_B_num=(B_vartheta_B_p+B_vartheta_B_m-2.d0*B_vartheta_B)/eps**2
    dB_varphi_B_num=(B_varphi_B_p-B_varphi_B_m)/(2d0*eps)
    d2B_varphi_B_num=(B_varphi_B_p+B_varphi_B_m-2.d0*B_varphi_B)/eps**2
    dBmod_B_num(1)=(Bmod_B_p-Bmod_B_m)/(2d0*eps)
    d2Bmod_B_num(1)=(Bmod_B_p+Bmod_B_m-2.d0*Bmod_B)/eps**2
    dB_r_num(1)=(B_r_p-B_r_m)/(2d0*eps)
    d2B_r_num(1)=(B_r_p+B_r_m-2.d0*B_r)/eps**2
!
!
    call splint_boozer_coord(r,vartheta_B+eps,varphi_B,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_p,dBmod_B,d2Bmod_B,                                     &
                             B_r_p,dB_r,d2B_r)
!
    call splint_boozer_coord(r,vartheta_B-eps,varphi_B,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_m,dBmod_B,d2Bmod_B,                                     &
                             B_r_m,dB_r,d2B_r)
!
    dBmod_B_num(2)=(Bmod_B_p-Bmod_B_m)/(2d0*eps)
    d2Bmod_B_num(4)=(Bmod_B_p+Bmod_B_m-2.d0*Bmod_B)/eps**2
    dB_r_num(2)=(B_r_p-B_r_m)/(2d0*eps)
    d2B_r_num(4)=(B_r_p+B_r_m-2.d0*B_r)/eps**2
!
!
    call splint_boozer_coord(r,vartheta_B,varphi_B+eps,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_p,dBmod_B,d2Bmod_B,                                     &
                             B_r_p,dB_r,d2B_r)
!
    call splint_boozer_coord(r,vartheta_B,varphi_B-eps,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_m,dBmod_B,d2Bmod_B,                                     &
                             B_r_m,dB_r,d2B_r)
!
    dBmod_B_num(3)=(Bmod_B_p-Bmod_B_m)/(2d0*eps)
    d2Bmod_B_num(6)=(Bmod_B_p+Bmod_B_m-2.d0*Bmod_B)/eps**2
    dB_r_num(3)=(B_r_p-B_r_m)/(2d0*eps)
    d2B_r_num(6)=(B_r_p+B_r_m-2.d0*B_r)/eps**2
!
!
    call splint_boozer_coord(r+eps,vartheta_B+eps,varphi_B,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_pp,dBmod_B,d2Bmod_B,                                     &
                             B_r_pp,dB_r,d2B_r)
!
    call splint_boozer_coord(r-eps,vartheta_B-eps,varphi_B,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_mm,dBmod_B,d2Bmod_B,                                     &
                             B_r_mm,dB_r,d2B_r)
!
    call splint_boozer_coord(r+eps,vartheta_B-eps,varphi_B,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_pm,dBmod_B,d2Bmod_B,                                     &
                             B_r_pm,dB_r,d2B_r)
!
    call splint_boozer_coord(r-eps,vartheta_B+eps,varphi_B,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_mp,dBmod_B,d2Bmod_B,                                     &
                             B_r_mp,dB_r,d2B_r)
!
    d2Bmod_B_num(2)=(Bmod_B_pp+Bmod_B_mm-Bmod_B_pm-Bmod_B_mp)/(2.d0*eps)**2
    d2B_r_num(2)=(B_r_pp+B_r_mm-B_r_pm-B_r_mp)/(2d0*eps)**2
!
!
    call splint_boozer_coord(r+eps,vartheta_B,varphi_B+eps,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_pp,dBmod_B,d2Bmod_B,                                     &
                             B_r_pp,dB_r,d2B_r)
!
    call splint_boozer_coord(r-eps,vartheta_B,varphi_B-eps,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_mm,dBmod_B,d2Bmod_B,                                     &
                             B_r_mm,dB_r,d2B_r)
!
    call splint_boozer_coord(r+eps,vartheta_B,varphi_B-eps,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_pm,dBmod_B,d2Bmod_B,                                     &
                             B_r_pm,dB_r,d2B_r)
!
    call splint_boozer_coord(r-eps,vartheta_B,varphi_B+eps,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_mp,dBmod_B,d2Bmod_B,                                     &
                             B_r_mp,dB_r,d2B_r)
!
    d2Bmod_B_num(3)=(Bmod_B_pp+Bmod_B_mm-Bmod_B_pm-Bmod_B_mp)/(2.d0*eps)**2
    d2B_r_num(3)=(B_r_pp+B_r_mm-B_r_pm-B_r_mp)/(2d0*eps)**2
!
!
    call splint_boozer_coord(r,vartheta_B+eps,varphi_B+eps,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_pp,dBmod_B,d2Bmod_B,                                     &
                             B_r_pp,dB_r,d2B_r)
!
    call splint_boozer_coord(r,vartheta_B-eps,varphi_B-eps,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_mm,dBmod_B,d2Bmod_B,                                     &
                             B_r_mm,dB_r,d2B_r)
!
    call splint_boozer_coord(r,vartheta_B-eps,varphi_B+eps,                                       &
                             A_theta,A_phi_p,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_p,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_p,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_pm,dBmod_B,d2Bmod_B,                                     &
                             B_r_pm,dB_r,d2B_r)
!
    call splint_boozer_coord(r,vartheta_B+eps,varphi_B-eps,                                       &
                             A_theta,A_phi_m,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B_m,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B_m,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B_mp,dBmod_B,d2Bmod_B,                                     &
                             B_r_mp,dB_r,d2B_r)
!
    d2Bmod_B_num(5)=(Bmod_B_pp+Bmod_B_mm-Bmod_B_pm-Bmod_B_mp)/(2.d0*eps)**2
    d2B_r_num(5)=(B_r_pp+B_r_mm-B_r_pm-B_r_mp)/(2d0*eps)**2
!
!
    write(12345,*) xi,dA_phi_dr_num,d2A_phi_dr2_num,     &
                   dB_vartheta_B_num,d2B_vartheta_B_num, &
                   dB_varphi_B_num,d2B_varphi_B_num,     &
                   dBmod_B_num,d2Bmod_B_num,dB_r_num,d2B_r_num
  enddo
!
  close(1234)
  close(12345)
!
!
  end program test_boozer
