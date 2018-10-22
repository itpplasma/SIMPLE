program test_magfie

use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
use chamb_mod,  only : rbig,rcham2
use parmot_mod, only : rmu,ro0,eeff
use velo_mod,   only : isw_field_type, neval_rk
use orbit_symplectic, only : orbit_sympl_init, orbit_timestep_sympl, z
use field_can_mod, only : field, neval, f, df, d2f, eval_field_can, get_derivatives2

implicit none

double precision, parameter :: pi=3.14159265358979d0
double precision, parameter :: c=2.9979d10
double precision, parameter :: e_charge=4.8032d-10
double precision, parameter :: e_mass=9.1094d-28
double precision, parameter :: p_mass=1.6726d-24
double precision, parameter :: ev=1.6022d-12

integer          :: npoi,ierr,L1i,nper,npoiper,i,ntimstep,ntestpart
integer          :: ipart,notrace_passing,loopskip,iskip,ilost
real             :: zzg
double precision :: dphi,rbeg,phibeg,zbeg,bmod00,rcham,rlarm,bmax,bmin
double precision :: tau,dtau,dtaumin,xi,v0,bmod_ref,E_alpha,trace_time
double precision :: RT0,R0i,cbfi,bz0i,bf0,trap_par
double precision :: sbeg,thetabeg
double precision :: z0(5)
integer          :: npoiper2
double precision :: contr_pp
double precision :: facE_al
integer          :: ibins
integer          :: n_e,n_d,n_b
double precision :: r,vartheta_c,varphi_c,theta_vmec,varphi_vmec,alam0
double precision :: alam,alam_prev,par_inv

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
call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)
call get_canonical_coordinates

r=0.7d0
vartheta_c=0.5d0
varphi_c=0.5d0
alam0=0.3d0 !0.5d0

call can_to_vmec(r,vartheta_c,varphi_c,theta_vmec,varphi_vmec)

print *,'can : ',r,vartheta_c,varphi_c
print *,'VMEC: ',r,theta_vmec,varphi_vmec

isw_field_type=0
z0(1)=r
z0(2)=vartheta_c
z0(3)=varphi_c
z0(4)=1.d0
z0(5)=alam0

call orbit_sympl_init(z0)
call do_test(z)

contains

subroutine do_test(z)

    double precision :: z(4), dz(4)
    integer :: k
    type(field) :: dfnum
    double precision, parameter :: dx = 1d-8

    ! quantities to test: f, df, d2f containing Ath, Aph, hth, hph, Bmod

    do k = 1,3
        dz = 0d0
        dz(k) = .5d0*dx
        call eval_field_can(z + dz)
        dfnum%dAth(k) = f%Ath
        dfnum%dAph(k) = f%Aph
        dfnum%dhth(k) = f%hth
        dfnum%dhph(k) = f%hph
        dfnum%dBmod(k) = f%Bmod
        call eval_field_can(z - dz)
        dfnum%dAth(k) = (dfnum%dAth(k) - f%Ath)/dx
        dfnum%dAph(k) = (dfnum%dAph(k) - f%Aph)/dx
        dfnum%dhth(k) = (dfnum%dhth(k) - f%hth)/dx
        dfnum%dhph(k) = (dfnum%dhph(k) - f%hph)/dx
        dfnum%dBmod(k) = (dfnum%dBmod(k) - f%Bmod(k))/dx
        call eval_field_can(z)
        print *, dfnum
        print *, df
    end do


    call get_derivatives2(z(4))

    ! quantities to test: H, pth, vpar
end subroutine do_test

end program test_magfie
