program test_magfie

use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
use parmot_mod, only : rmu,ro0
use velo_mod,   only : isw_field_type
use orbit_symplectic, only : orbit_sympl_init, orbit_timestep_sympl, z, &
  f_sympl_euler1, f_sympl_euler2, jac_sympl_euler1, jac_sympl_euler2
use field_can_mod, only : field_can, d_field_can, d2_field_can, f, df, d2f,&
        vpar, H, pth, dvpar, dH, dpth, d2vpar, d2H, d2pth, eval_field_can, &
        get_val, get_derivatives, get_derivatives2 

implicit none

double precision, parameter :: pi=3.14159265358979d0
double precision, parameter :: c=2.9979d10
double precision, parameter :: e_charge=4.8032d-10
double precision, parameter :: e_mass=9.1094d-28
double precision, parameter :: p_mass=1.6726d-24
double precision, parameter :: ev=1.6022d-12

integer          :: L1i,nper,npoiper,ntimstep,ntestpart
integer          :: notrace_passing,loopskip

double precision :: phibeg,bmod00,rlarm
double precision :: tau,dtau,v0,bmod_ref,E_alpha,trace_time
double precision :: RT0,R0i,cbfi,bz0i,bf0
double precision :: sbeg,thetabeg
double precision :: z0(5)
integer          :: npoiper2
double precision :: contr_pp
double precision :: facE_al
integer          :: n_e,n_d
double precision :: r,vartheta_c,varphi_c,theta_vmec,varphi_vmec,alam0

open(1,file='alpha_lifetime_m.inp', recl=1024)
read (1,*) notrace_passing   !skip tracing passing prts if notrace_passing=1
read (1,*) nper              !number of periods for initial field line
read (1,*) npoiper           !number of points per period on this field line
read (1,*) ntimstep          !number of time steps per slowing down time
read (1,*) ntestpart         !number of test particles
read (1,*) bmod_ref          !reference field, G, for Boozer $B_{00}$
read (1,*) trace_time        !slowing down time, s
read (1,*) sbeg              !starting s for field line                       !<=2017
read (1,*) phibeg            !starting phi for field line                     !<=2017
read (1,*) thetabeg          !starting theta for field line                   !<=2017
read (1,*) loopskip          !how many loops to skip to shift random numbers
read (1,*) contr_pp          !control of passing particle fraction
read (1,*) facE_al           !facE_al test particle energy reduction factor
read (1,*) npoiper2          !additional integration step split factor
read (1,*) n_e               !test particle charge number (the same as Z)
read (1,*) n_d               !test particle mass number (the same as A)
read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
close(1)

! inverse relativistic temperature
  rmu=1d8
!
! alpha particle energy, eV:
  E_alpha=3.5d6/facE_al
! alpha particle velocity, cm/s
  v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
! 14.04.2013 end
!
! Larmor radius:
  rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)
! normalized slowing down time:
  tau=trace_time*v0
! normalized time step:
  dtau=tau/dfloat(ntimstep-1)
!
bmod00=281679.46317784750d0
! Larmor raidus corresponds to the field stregth egual to $B_{00}$ harmonic
! in Boozer coordinates:
! 14.11.2011  bmod00=bmod_ref  !<=deactivated, use value from the 'alpha_lifetime.inp'
  ro0=rlarm*bmod00  ! 23.09.2013
!
  multharm=3 !7
  ns_A=5
  ns_s=5
  ns_tp=5

call spline_vmec_data
call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0)
call get_canonical_coordinates

!r=0.7d0
!vartheta_c=0.5d0
!varphi_c=0.5d0
!alam0=0.3d0 !0.5d0

r=0.5d0
vartheta_c=0.0d0
varphi_c=0.314d0
alam0=0.3d0

call can_to_vmec(r, vartheta_c, varphi_c, theta_vmec, varphi_vmec)

print *,'can : ', r, vartheta_c, varphi_c
print *,'VMEC: ', r, theta_vmec, varphi_vmec

isw_field_type=0
z0(1)=r
z0(2)=vartheta_c
z0(3)=varphi_c
z0(4)=1.d0
z0(5)=alam0

call orbit_sympl_init(z0, 1.0d0, 1.0d0, 0)
call do_test

contains

function relerr(a, b)
    double precision :: relerr
    double precision, intent(in) :: a, b
    relerr = merge(0d0, (a - b)/b, b == 0d0)
end function relerr

subroutine der2(x0, i, j)
    double precision, intent(in) :: x0(3)
    integer, intent(in) :: i, j
    double precision hi, hj
    type(field_can) :: f00, f01, f10, f11
    type(d2_field_can) :: d2fnum
    double precision :: x(3), dxi(3), dxj(3)
    double precision, dimension(10) ::  d2vparnum, d2Hnum, d2pthnum
    double precision :: vpar00, vpar11, vpar10, vpar01, &
        H00, H11, H10, H01, pth00, pth11, pth10, pth01
    integer :: k

    hi = 1d-4
    hj = 1d-4

    dxi = 0d0
    dxj = 0d0
    dxi(i) = .5d0*hi
    dxj(j) = .5d0*hj

    x = x0 - dxi - dxj
    call eval_field_can(x(1), x(2), x(3), 0)
    call get_val(z(4))
    f00 = f
    vpar00 = vpar
    H00 = H
    pth00 = pth
    x = x0 - dxi + dxj
    call eval_field_can(x(1), x(2), x(3), 0)
    call get_val(z(4))
    f01 = f
    vpar01 = vpar
    H10 = H
    pth01 = pth
    x = x0 + dxi - dxj
    call eval_field_can(x(1), x(2), x(3), 0)
    call get_val(z(4))
    f10 = f
    vpar10 = vpar
    H01 = H
    pth10 = pth
    x = x0 + dxi + dxj
    call eval_field_can(x(1), x(2), x(3), 0)
    call get_val(z(4))
    f11 = f
    vpar11 = vpar
    H11 = H
    pth11 = pth

    call eval_field_can(x0(1), x0(2), x0(3), 2)
    call get_derivatives2(z(4))

    if (i==1 .and. j==1) k=1
    if ((i==1 .and. j==2) .or. (i==2 .and. j==1)) k=2
    if ((i==1 .and. j==3) .or. (i==3 .and. j==1)) k=3
    
    if (i==2 .and. j==2) k=4
    if ((i==2 .and. j==3) .or. (i==3 .and. j==2)) k=5
    
    if (i==3 .and. j==3) k=6

    if(i==j) then
        d2fnum%d2Ath(k) = (f11%Ath - 2d0*f%Ath + f00%Ath)/(hi*hj)
        d2fnum%d2Aph(k) = (f11%Aph - 2d0*f%Aph + f00%Aph)/(hi*hj)
        d2fnum%d2hth(k) = (f11%hth - 2d0*f%hth + f00%hth)/(hi*hj)
        d2fnum%d2hph(k) = (f11%hph - 2d0*f%hph + f00%hph)/(hi*hj)
        d2fnum%d2Bmod(k) = (f11%Bmod - 2d0*f%Bmod + f00%Bmod)/(hi*hj) 
        d2vparnum(k) = (vpar11 - 2d0*vpar + vpar00)/(hi*hj) 
        d2Hnum(k) = (H11 - 2d0*H + H00)/(hi*hj) 
        d2pthnum(k) = (pth11 - 2d0*pth + pth00)/(hi*hj) 
    else
        d2fnum%d2Ath(k) = (f11%Ath - f10%Ath - f01%Ath + f00%Ath)/(hi*hj)
        d2fnum%d2Aph(k) = (f11%Aph - f10%Aph - f01%Aph + f00%Aph)/(hi*hj)
        d2fnum%d2hth(k) = (f11%hth - f10%hth - f01%hth + f00%hth)/(hi*hj)
        d2fnum%d2hph(k) = (f11%hph - f10%hph - f01%hph + f00%hph)/(hi*hj)
        d2fnum%d2Bmod(k) = (f11%Bmod - f10%Bmod - f01%Bmod + f00%Bmod)/(hi*hj)
        d2vparnum(k) = (vpar11 - vpar10 - vpar01 + vpar00)/(hi*hj)
        d2Hnum(k) = (H11 - H10 - H01 + H00)/(hi*hj)
        d2pthnum(k) = (pth11 - pth10 - pth01 + pth00)/(hi*hj)
    end if

    print *, 'd2Ath (',i,j,')', d2f%d2Ath(k), d2fnum%d2Ath(k), relerr(d2fnum%d2Ath(k), d2f%d2Ath(k))
    print *, 'd2Aph (',i,j,')', d2f%d2Aph(k), d2fnum%d2Aph(k), relerr(d2fnum%d2Aph(k), d2f%d2Aph(k))
    print *, 'd2hth (',i,j,')', d2f%d2hth(k), d2fnum%d2hth(k), relerr(d2fnum%d2hth(k), d2f%d2hth(k))
    print *, 'd2hph (',i,j,')', d2f%d2hph(k), d2fnum%d2hph(k), relerr(d2fnum%d2hph(k), d2f%d2hph(k))
    print *, 'd2Bmod(',i,j,')', d2f%d2Bmod(k), d2fnum%d2Bmod(k), relerr(d2fnum%d2Bmod(k), d2f%d2Bmod(k))
    print *, 'd2vpar(',i,j,')', d2vpar(k), d2vparnum(k), relerr(d2vparnum(k), d2vpar(k))
    print *, 'd2H(',i,j,')', d2H(k), d2Hnum(k), relerr(d2Hnum(k), d2H(k))
    print *, 'd2pth(',i,j,')', d2pth(k), d2pthnum(k), relerr(d2pthnum(k), d2pth(k))
end subroutine der2

subroutine test_jac  
  double precision :: x1(2), x2(3), dx1(2), dx2(3), jac1(2,2), jac2(3,3), x10(2), x20(3), &
    h1(2), h2(3), jac1num(2,2), jac2num(3,3), fvec1(2), fvec2(3)
  integer :: k

  h1(1) = 1d-8
  h1(2) = z(4)*1d-8

  h2 = 1d-8

  do k = 1,2
    dx1 = 0d0
    dx1(k) = h1(k)*0.5d0
    x10 = z((/1,4/)) + (/1d-4, 1d-2/)

    x1 = x10 + dx1
    call f_sympl_euler1(2, x1, fvec1, 0)
    jac1num(:, k) = fvec1

    x1 = x10 - dx1
    call f_sympl_euler1(2, x1, fvec1, 0)
    jac1num(:, k) = (jac1num(:, k) - fvec1)/h1(k)

    x1 = x10
    call f_sympl_euler1(2, x1, fvec1, 0)
    call jac_sympl_euler1(x1, jac1)

  end do

  print *, 'jac_sympl_euler1(1,1)', jac1(1,1), jac1num(1,1), relerr(jac1(1,1), jac1num(1,1))
  print *, 'jac_sympl_euler1(1,2)', jac1(1,2), jac1num(1,2), relerr(jac1(1,2), jac1num(1,2))
  print *, 'jac_sympl_euler1(2,1)', jac1(2,1), jac1num(2,1), relerr(jac1(2,1), jac1num(2,1))
  print *, 'jac_sympl_euler1(2,2)', jac1(2,2), jac1num(2,2), relerr(jac1(2,2), jac1num(2,2))


  do k = 1,3
    dx2 = 0d0
    dx2(k) = h2(k)*0.5d0
    x20 = z(1:3) + 1d-4

    x2 = x20 + dx2
    call f_sympl_euler2(3, x2, fvec2, 0)
    jac2num(:, k) = fvec2

    x2 = x20 - dx2
    call f_sympl_euler2(3, x2, fvec2, 0)
    jac2num(:, k) = (jac2num(:, k) - fvec2)/h2(k)

    x2 = x20
    call f_sympl_euler2(3, x2, fvec2, 0)
    call jac_sympl_euler2(x2, jac2)

  end do


  print *, 'jac_sympl_euler2(1,1)', jac2(1,1), jac2num(1,1), relerr(jac2(1,1), jac2num(1,1))
  print *, 'jac_sympl_euler2(1,2)', jac2(1,2), jac2num(1,2), relerr(jac2(1,2), jac2num(1,2))
  print *, 'jac_sympl_euler2(1,3)', jac2(1,3), jac2num(1,3), relerr(jac2(1,3), jac2num(1,3))
  print *, 'jac_sympl_euler2(2,1)', jac2(2,1), jac2num(2,1), relerr(jac2(2,1), jac2num(2,1))
  print *, 'jac_sympl_euler2(2,2)', jac2(2,2), jac2num(2,2), relerr(jac2(2,2), jac2num(2,2))
  print *, 'jac_sympl_euler2(2,3)', jac2(2,3), jac2num(2,3), relerr(jac2(2,3), jac2num(2,3))
  print *, 'jac_sympl_euler2(3,1)', jac2(3,1), jac2num(3,1), relerr(jac2(3,1), jac2num(3,1))
  print *, 'jac_sympl_euler2(3,2)', jac2(3,2), jac2num(3,2), relerr(jac2(3,2), jac2num(3,2))
  print *, 'jac_sympl_euler2(3,3)', jac2(3,3), jac2num(3,3), relerr(jac2(3,3), jac2num(3,3))


end subroutine test_jac

subroutine test_newton
  integer, parameter :: n = 2
  double precision :: x(n), fvec(n), fjac(n,n), ijac(n,n)
  integer :: k

  x = z((/1,4/)) + (/1d-4, 1d-2/)
 
  do k=1,10
    call f_sympl_euler1(n, x, fvec, 1)
    call jac_sympl_euler1(x, fjac)
    print *, x, fvec
    ijac(1,1) = fjac(2,2)
    ijac(1,2) = -fjac(1,2)
    ijac(2,1) = -fjac(2,1)
    ijac(2,2) = fjac(1,1)
    ijac = ijac/(fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1))
    x = x - matmul(ijac, fvec)
  enddo

  call f_sympl_euler1(n, x, fvec, 1)
  print *, x, fvec
  print *, fjac

end subroutine


subroutine test_newton2
  integer, parameter :: n = 3
  double precision :: x(n)
  double precision :: fvec(n), fjac(n,n)
  integer :: pivot(n), info
  integer :: k

  x = z(1:3) + 1d-4

  do k=1,10
    call f_sympl_euler2(n, x, fvec, 1)
    call jac_sympl_euler2(x, fjac)
    print *, x, fvec

    call dgesv(n, 1, fjac, n, pivot, fvec, 3, info) 
    ! after solution: fvec = (xold-xnew)_Newton

    x = x - fvec
  enddo
  call f_sympl_euler2(n, x, fvec, 1)
  print *, x, fvec
  print *, fjac
end subroutine

subroutine do_test

    double precision :: dz(4)
    integer :: i, j, k
    double precision :: dx
    type(d_field_can) :: dfnum
    double precision :: dvparnum(4), dHnum(4), dpthnum(4)

    print *, 'f\t', 'derivative\t', 'numerical derivative\t', 'relative error'

    ! quantities to test: Ath, Aph, hth, hph, Bmod, vpar, H, pth

    do k = 1,3
        dz = 0d0
        dx = 1d-8
        dz(k) = .5d0*dx
        call eval_field_can(z(1) + dz(1), z(2) + dz(2), z(3) + dz(3), 0)
        call get_val(z(4))
        dfnum%dAth(k) = f%Ath
        dfnum%dAph(k) = f%Aph
        dfnum%dhth(k) = f%hth
        dfnum%dhph(k) = f%hph
        dfnum%dBmod(k) = f%Bmod
        dvparnum(k) = vpar
        dHnum(k) = H
        dpthnum(k) = pth
        call eval_field_can(z(1) - dz(1), z(2) - dz(2), z(3) - dz(3), 0)
        call get_val(z(4))
        dfnum%dAth(k) = (dfnum%dAth(k) - f%Ath)/dx
        dfnum%dAph(k) = (dfnum%dAph(k) - f%Aph)/dx
        dfnum%dhth(k) = (dfnum%dhth(k) - f%hth)/dx
        dfnum%dhph(k) = (dfnum%dhph(k) - f%hph)/dx
        dfnum%dBmod(k) = (dfnum%dBmod(k) - f%Bmod)/dx
        dvparnum(k) = (dvparnum(k) - vpar)/dx
        dHnum(k) = (dHnum(k) - H)/dx
        dpthnum(k) = (dpthnum(k) - pth)/dx
        call eval_field_can(z(1), z(2), z(3), 0)
        call get_derivatives(z(4))

        print *, 'dAth (',k,')', df%dAth(k), dfnum%dAth(k), relerr(dfnum%dAth(k), df%dAth(k))
        print *, 'dAph (',k,')', df%dAph(k), dfnum%dAph(k), relerr(dfnum%dAph(k), df%dAph(k))
        print *, 'dhth (',k,')', df%dhth(k), dfnum%dhth(k), relerr(dfnum%dhth(k), df%dhth(k))
        print *, 'dhph (',k,')', df%dhph(k), dfnum%dhph(k), relerr(dfnum%dhph(k), df%dhph(k))
        print *, 'dBmod(',k,')', df%dBmod(k), dfnum%dBmod(k), relerr(dfnum%dBmod(k), df%dBmod(k))
    enddo

    dx = 1d-8*z(4)
    call get_val(z(4) + .5d0*dx)
    dvparnum(4) = vpar
    dHnum(4) = H
    dpthnum(4) = pth
    call get_val(z(4) - .5d0*dx)
    dvparnum(4) = (dvparnum(4) - vpar)/dx
    dHnum(4) = (dHnum(4) - H)/dx
    dpthnum(4) = (dpthnum(4) - pth)/dx
    call get_derivatives(z(4))
    
    do k=1,3
        print *, 'dvpar(',k,')', dvpar(k), dvparnum(k), relerr(dvparnum(k), dvpar(k))
        print *, 'dH   (',k,')', dH(k), dHnum(k), relerr(dHnum(k), dH(k))
        print *, 'dpth (',k,')', dpth(k), dpthnum(k), relerr(dpthnum(k), dpth(k))
    enddo
    print *, 'dvpardpph', dvpar(4), dvparnum(4), relerr(dvparnum(4), dvpar(4))
    print *, 'dHdpph   ', dH(4), dHnum(4), relerr(dHnum(4), dH(4))
    print *, 'dpthdpph ', dpth(4), dpthnum(4), relerr(dpthnum(4), dpth(4))

    do i = 1,3
        do j = 1,3
            call der2(z(1:3), i, j)
        enddo
    enddo

    ! TODO: second ders in pphi and mixed

    call test_jac
    call test_newton
    call test_newton2
end subroutine do_test

end program test_magfie
