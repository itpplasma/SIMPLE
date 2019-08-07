program test_magfie

! use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
! use parmot_mod, only : rmu,ro0
! use velo_mod,   only : isw_field_type
use orbit_symplectic
use field_can_mod

implicit none

! double precision, parameter :: pi=3.14159265358979d0
! double precision, parameter :: c=2.9979d10
! double precision, parameter :: e_charge=4.8032d-10
! double precision, parameter :: e_mass=9.1094d-28
! double precision, parameter :: p_mass=1.6726d-24
! double precision, parameter :: ev=1.6022d-12

! integer          :: L1i,nper,npoiper,ntimstep,ntestpart
! integer          :: notrace_passing,loopskip

! double precision :: phibeg,bmod00,rlarm
! double precision :: tau,dtau,v0,bmod_ref,E_alpha,trace_time
! double precision :: RT0,R0i,cbfi,bz0i,bf0
! double precision :: sbeg,thetabeg
! double precision :: z0(5)
! integer          :: npoiper2
! double precision :: contr_pp
! double precision :: facE_al
! integer          :: n_e,n_d
! double precision :: r,vartheta_c,varphi_c,theta_vmec,varphi_vmec,alam0

! open(1,file='alpha_lifetime_m.inp', recl=1024)
! read (1,*) notrace_passing   !skip tracing passing prts if notrace_passing=1
! read (1,*) nper              !number of periods for initial field line
! read (1,*) npoiper           !number of points per period on this field line
! read (1,*) ntimstep          !number of time steps per slowing down time
! read (1,*) ntestpart         !number of test particles
! read (1,*) bmod_ref          !reference field, G, for Boozer $B_{00}$
! read (1,*) trace_time        !slowing down time, s
! read (1,*) sbeg              !starting s for field line                       !<=2017
! read (1,*) phibeg            !starting phi for field line                     !<=2017
! read (1,*) thetabeg          !starting theta for field line                   !<=2017
! read (1,*) loopskip          !how many loops to skip to shift random numbers
! read (1,*) contr_pp          !control of passing particle fraction
! read (1,*) facE_al           !facE_al test particle energy reduction factor
! read (1,*) npoiper2          !additional integration step split factor
! read (1,*) n_e               !test particle charge number (the same as Z)
! read (1,*) n_d               !test particle mass number (the same as A)
! read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
! close(1)

! ! inverse relativistic temperature
!   rmu=1d8
! !
! ! alpha particle energy, eV:
!   E_alpha=3.5d6/facE_al
! ! alpha particle velocity, cm/s
!   v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
! ! 14.04.2013 end
! !
! ! Larmor radius:
!   rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)
! ! normalized slowing down time:
!   tau=trace_time*v0
! ! normalized time step:
!   dtau=tau/dfloat(ntimstep-1)
! !
! bmod00=281679.46317784750d0
! ! Larmor raidus corresponds to the field stregth egual to $B_{00}$ harmonic
! ! in Boozer coordinates:
! ! 14.11.2011  bmod00=bmod_ref  !<=deactivated, use value from the 'alpha_lifetime.inp'
!   ro0=rlarm*bmod00  ! 23.09.2013
! !
!   multharm=3 !7
!   ns_A=5
!   ns_s=5
!   ns_tp=5

! call spline_vmec_data
! call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0)
! call get_canonical_coordinates

! !r=0.7d0
! !vartheta_c=0.5d0
! !varphi_c=0.5d0
! !alam0=0.3d0 !0.5d0

! r=0.5d0
! vartheta_c=0.0d0
! varphi_c=0.314d0
! alam0=0.3d0

! call can_to_vmec(r, vartheta_c, varphi_c, theta_vmec, varphi_vmec)

! print *,'can : ', r, vartheta_c, varphi_c
! print *,'VMEC: ', r, theta_vmec, varphi_vmec

! isw_field_type=0
! z0(1)=r
! z0(2)=vartheta_c
! z0(3)=varphi_c
! z0(4)=1.d0
! z0(5)=alam0

! call orbit_sympl_init(si, f, z0, 1.0d0, 1, 1d-12, 0, 0)

double precision :: z0(4), vpar0
type(FieldCan) :: f

! Initial conditions
z0(1) = 0.7d0  ! r
z0(2) = 0.5d0  ! theta
z0(3) = 0.3d0  ! phi
vpar0 = 0.1d0  ! parallel velocity

! Compute toroidal momentum from initial conditions
f%ro0 = 1.0d0
f%mu = 1.0d-5
call eval_field(f, z0(1), z0(2), z0(3), 0)
z0(4) = vpar0*f%hph + f%Aph/f%ro0  ! p_phi
call do_test

contains

function relerr(a, b)
    double precision :: relerr
    double precision, intent(in) :: a, b
    relerr = merge(0d0, (a - b)/b, b == 0d0)
end function relerr

subroutine der2(f, x0, pphi, i, j)
    double precision, intent(in) :: x0(3)
    integer, intent(in) :: i, j
    double precision hi, hj
    type(FieldCan) :: f, f00, f01, f10, f11
    type(FieldCan) :: d2fnum
    double precision :: pphi, x(3), dxi(3), dxj(3)
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
    call eval_field(f, x(1), x(2), x(3), 0)
    call get_val(f, pphi)
    f00 = f
    vpar00 = f%vpar
    H00 = f%H
    pth00 = f%pth
    x = x0 - dxi + dxj
    call eval_field(f, x(1), x(2), x(3), 0)
    call get_val(f, pphi)
    f01 = f
    vpar01 = f%vpar
    H10 = f%H
    pth01 = f%pth
    x = x0 + dxi - dxj
    call eval_field(f, x(1), x(2), x(3), 0)
    call get_val(f, pphi)
    f10 = f
    vpar10 = f%vpar
    H01 = f%H
    pth10 = f%pth
    x = x0 + dxi + dxj
    call eval_field(f, x(1), x(2), x(3), 0)
    call get_val(f, pphi)
    f11 = f
    vpar11 = f%vpar
    H11 = f%H
    pth11 = f%pth

    call eval_field(f, x0(1), x0(2), x0(3), 2)
    call get_derivatives2(f, pphi)

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
        d2vparnum(k) = (vpar11 - 2d0*f%vpar + vpar00)/(hi*hj) 
        d2Hnum(k) = (H11 - 2d0*f%H + H00)/(hi*hj) 
        d2pthnum(k) = (pth11 - 2d0*f%pth + pth00)/(hi*hj) 
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

    print *, 'd2Ath (',i,j,')', f%d2Ath(k), d2fnum%d2Ath(k), relerr(d2fnum%d2Ath(k), f%d2Ath(k))
    print *, 'd2Aph (',i,j,')', f%d2Aph(k), d2fnum%d2Aph(k), relerr(d2fnum%d2Aph(k), f%d2Aph(k))
    print *, 'd2hth (',i,j,')', f%d2hth(k), d2fnum%d2hth(k), relerr(d2fnum%d2hth(k), f%d2hth(k))
    print *, 'd2hph (',i,j,')', f%d2hph(k), d2fnum%d2hph(k), relerr(d2fnum%d2hph(k), f%d2hph(k))
    print *, 'd2Bmod(',i,j,')', f%d2Bmod(k), d2fnum%d2Bmod(k), relerr(d2fnum%d2Bmod(k), f%d2Bmod(k))
    print *, 'd2vpar(',i,j,')', f%d2vpar(k), d2vparnum(k), relerr(d2vparnum(k), f%d2vpar(k))
    print *, 'd2H(',i,j,')', f%d2H(k), d2Hnum(k), relerr(d2Hnum(k), f%d2H(k))
    print *, 'd2pth(',i,j,')', f%d2pth(k), d2pthnum(k), relerr(d2pthnum(k), f%d2pth(k))
end subroutine der2

subroutine test_jac1(si)  
  type(SymplecticIntegrator) :: si
  double precision :: x1(2), x2(3), dx1(2), dx2(3), jac1(2,2), jac2(3,3), x10(2), x20(3), &
    h1(2), h2(3), jac1num(2,2), jac2num(3,3), fvec1(2), fvec2(3)
  integer :: k

  h1(1) = 1d-8
  h1(2) = si%z(4)*1d-8

  do k = 1,2
    dx1 = 0d0
    dx1(k) = h1(k)*0.5d0
    x10 = si%z((/1,4/)) + (/1d-4, 1d-2/)

    x1 = x10 + dx1
    call f_sympl_euler1(si, f, 2, x1, fvec1, 0)
    jac1num(:, k) = fvec1

    x1 = x10 - dx1
    call f_sympl_euler1(si, f, 2, x1, fvec1, 0)
    jac1num(:, k) = (jac1num(:, k) - fvec1)/h1(k)

    x1 = x10
    call f_sympl_euler1(si, f, 2, x1, fvec1, 0)
    call jac_sympl_euler1(si, f, x1, jac1)

  end do

  print *, 'jac_sympl_euler1(1,1)', jac1(1,1), jac1num(1,1), relerr(jac1(1,1), jac1num(1,1))
  print *, 'jac_sympl_euler1(1,2)', jac1(1,2), jac1num(1,2), relerr(jac1(1,2), jac1num(1,2))
  print *, 'jac_sympl_euler1(2,1)', jac1(2,1), jac1num(2,1), relerr(jac1(2,1), jac1num(2,1))
  print *, 'jac_sympl_euler1(2,2)', jac1(2,2), jac1num(2,2), relerr(jac1(2,2), jac1num(2,2))
end subroutine test_jac1

subroutine test_jac2(si)  
  type(SymplecticIntegrator) :: si
  double precision :: x1(2), x2(3), dx1(2), dx2(3), jac1(2,2), jac2(3,3), x10(2), x20(3), &
    h1(2), h2(3), jac1num(2,2), jac2num(3,3), fvec1(2), fvec2(3)
  integer :: k

  h2 = 1d-8

  do k = 1,3
    dx2 = 0d0
    dx2(k) = h2(k)*0.5d0
    x20 = si%z(1:3) + 1d-4

    x2 = x20 + dx2
    call f_sympl_euler2(si, f, 3, x2, fvec2, 0)
    jac2num(:, k) = fvec2

    x2 = x20 - dx2
    call f_sympl_euler2(si, f, 3, x2, fvec2, 0)
    jac2num(:, k) = (jac2num(:, k) - fvec2)/h2(k)

    x2 = x20
    call f_sympl_euler2(si, f, 3, x2, fvec2, 0)
    call jac_sympl_euler2(si, f, x2, jac2)

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

end subroutine test_jac2

subroutine test_jac_midpoint(si)
  type(SymplecticIntegrator) :: si
  type(FieldCan) :: fmid
  double precision :: x2(5), dx2(5), jac2(5,5), x20(5), h2(5), jac2num(5,5), fvec2(5)
  integer :: k

  h2 = 1d-6

  do k = 1,5
    dx2 = 0d0
    dx2(k) = h2(k)*0.5d0
    x20(1:4) = si%z(1:4) + 1d-4
    x20(5) = si%z(1) + 0.5d-4

    x2 = x20 + dx2
    call f_midpoint_part1(si, f, 5, x2, fvec2, 0)
    call f_midpoint_part2(si, f, 5, x2, fvec2, 0)
    jac2num(:, k) = fvec2

    x2 = x20 - dx2
    call f_midpoint_part1(si, f, 5, x2, fvec2, 0)
    call f_midpoint_part2(si, f, 5, x2, fvec2, 0)
    jac2num(:, k) = (jac2num(:, k) - fvec2)/h2(k)

    x2 = x20
    call f_midpoint_part1(si, f, 5, x2, fvec2, 0)
    call jac_midpoint_part1(si, f, x2, jac2)
    fmid = f
    call f_midpoint_part2(si, f, 5, x2, fvec2, 0)
    call jac_midpoint_part2(si, f, fmid, x2, jac2)

  end do


  print *, 'jac_midpoint(1,1)', jac2(1,1), jac2num(1,1), relerr(jac2(1,1), jac2num(1,1))
  print *, 'jac_midpoint(1,2)', jac2(1,2), jac2num(1,2), relerr(jac2(1,2), jac2num(1,2))
  print *, 'jac_midpoint(1,3)', jac2(1,3), jac2num(1,3), relerr(jac2(1,3), jac2num(1,3))
  print *, 'jac_midpoint(1,4)', jac2(1,4), jac2num(1,4), relerr(jac2(1,4), jac2num(1,4))
  print *, 'jac_midpoint(1,5)', jac2(1,5), jac2num(1,5), relerr(jac2(1,5), jac2num(1,5))
  print *, 'jac_midpoint(2,1)', jac2(2,1), jac2num(2,1), relerr(jac2(2,1), jac2num(2,1))
  print *, 'jac_midpoint(2,2)', jac2(2,2), jac2num(2,2), relerr(jac2(2,2), jac2num(2,2))
  print *, 'jac_midpoint(2,3)', jac2(2,3), jac2num(2,3), relerr(jac2(2,3), jac2num(2,3))
  print *, 'jac_midpoint(2,4)', jac2(2,4), jac2num(2,4), relerr(jac2(2,4), jac2num(2,4))
  print *, 'jac_midpoint(2,5)', jac2(2,5), jac2num(2,5), relerr(jac2(2,5), jac2num(2,5))
  print *, 'jac_midpoint(3,1)', jac2(3,1), jac2num(3,1), relerr(jac2(3,1), jac2num(3,1))
  print *, 'jac_midpoint(3,2)', jac2(3,2), jac2num(3,2), relerr(jac2(3,2), jac2num(3,2))
  print *, 'jac_midpoint(3,3)', jac2(3,3), jac2num(3,3), relerr(jac2(3,3), jac2num(3,3))
  print *, 'jac_midpoint(3,4)', jac2(3,4), jac2num(3,4), relerr(jac2(3,4), jac2num(3,4))
  print *, 'jac_midpoint(3,5)', jac2(3,5), jac2num(3,5), relerr(jac2(3,5), jac2num(3,5))
  print *, 'jac_midpoint(4,1)', jac2(4,1), jac2num(4,1), relerr(jac2(4,1), jac2num(4,1))
  print *, 'jac_midpoint(4,2)', jac2(4,2), jac2num(4,2), relerr(jac2(4,2), jac2num(4,2))
  print *, 'jac_midpoint(4,3)', jac2(4,3), jac2num(4,3), relerr(jac2(4,3), jac2num(4,3))
  print *, 'jac_midpoint(4,4)', jac2(4,4), jac2num(4,4), relerr(jac2(4,4), jac2num(4,4))
  print *, 'jac_midpoint(4,5)', jac2(4,5), jac2num(4,5), relerr(jac2(4,5), jac2num(4,5))
  print *, 'jac_midpoint(5,1)', jac2(5,1), jac2num(5,1), relerr(jac2(5,1), jac2num(5,1))
  print *, 'jac_midpoint(5,2)', jac2(5,2), jac2num(5,2), relerr(jac2(5,2), jac2num(5,2))
  print *, 'jac_midpoint(5,3)', jac2(5,3), jac2num(5,3), relerr(jac2(5,3), jac2num(5,3))
  print *, 'jac_midpoint(5,4)', jac2(5,4), jac2num(5,4), relerr(jac2(5,4), jac2num(5,4))
  print *, 'jac_midpoint(5,5)', jac2(5,5), jac2num(5,5), relerr(jac2(5,5), jac2num(5,5))

end subroutine test_jac_midpoint

subroutine test_newton(si, f)
  type(SymplecticIntegrator) :: si
  type(FieldCan) :: f
  integer, parameter :: n = 2
  double precision :: x(n), fvec(n), fjac(n,n), ijac(n,n)
  integer :: k

  x = si%z((/1,4/)) + (/1d-4, 1d-2/)
 
  do k=1,10
    call f_sympl_euler1(si, f, n, x, fvec, 1)
    call jac_sympl_euler1(si, f, x, fjac)
    print *, x, fvec
    ijac(1,1) = fjac(2,2)
    ijac(1,2) = -fjac(1,2)
    ijac(2,1) = -fjac(2,1)
    ijac(2,2) = fjac(1,1)
    ijac = ijac/(fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1))
    x = x - matmul(ijac, fvec)
  enddo

  call f_sympl_euler1(si, f, n, x, fvec, 1)
  print *, x, fvec
  print *, fjac

end subroutine


subroutine test_newton2(si, f)
  type(SymplecticIntegrator) :: si
  type(FieldCan) :: f
  integer, parameter :: n = 3
  double precision :: x(n)
  double precision :: fvec(n), fjac(n,n)
  integer :: pivot(n), info
  integer :: k

  x = si%z(1:3) + 1d-4

  do k=1,10
    call f_sympl_euler2(si, f, n, x, fvec, 1)
    call jac_sympl_euler2(si, f, x, fjac)
    print *, x, fvec

    call dgesv(n, 1, fjac, n, pivot, fvec, 3, info) 
    ! after solution: fvec = (xold-xnew)_Newton

    x = x - fvec
  enddo
  call f_sympl_euler2(si, f, n, x, fvec, 1)
  print *, x, fvec
  print *, fjac
end subroutine

subroutine do_test()

    type(SymplecticIntegrator) :: euler1, euler2, midpoint

    double precision :: dz(4)
    integer :: i, j, k
    double precision :: dx
    type(FieldCan) :: dfnum
    double precision :: dvparnum(4), dHnum(4), dpthnum(4)

    print *, 'f\t', 'derivative\t', 'numerical derivative\t', 'relative error'

    ! quantities to test: Ath, Aph, hth, hph, Bmod, vpar, H, pth

    do k = 1,3
        dz = 0d0
        dx = 1d-8
        dz(k) = .5d0*dx
        call eval_field(f, z0(1) + dz(1), z0(2) + dz(2), z0(3) + dz(3), 0)
        print *, f
        call get_val(f, z0(4))
        dfnum%dAth(k) = f%Ath
        dfnum%dAph(k) = f%Aph
        dfnum%dhth(k) = f%hth
        dfnum%dhph(k) = f%hph
        dfnum%dBmod(k) = f%Bmod
        dvparnum(k) = f%vpar
        dHnum(k) = f%H
        dpthnum(k) = f%pth
        call eval_field(f, z0(1) - dz(1), z0(2) - dz(2), z0(3) - dz(3), 0)
        call get_val(f, z0(4))
        dfnum%dAth(k) = (dfnum%dAth(k) - f%Ath)/dx
        dfnum%dAph(k) = (dfnum%dAph(k) - f%Aph)/dx
        dfnum%dhth(k) = (dfnum%dhth(k) - f%hth)/dx
        dfnum%dhph(k) = (dfnum%dhph(k) - f%hph)/dx
        dfnum%dBmod(k) = (dfnum%dBmod(k) - f%Bmod)/dx
        dvparnum(k) = (dvparnum(k) - f%vpar)/dx
        dHnum(k) = (dHnum(k) - f%H)/dx
        dpthnum(k) = (dpthnum(k) - f%pth)/dx
        call eval_field(f, z0(1), z0(2), z0(3), 0)
        call get_derivatives(f, z0(4))

        print *, 'dAth (',k,')', f%dAth(k), dfnum%dAth(k), relerr(dfnum%dAth(k), f%dAth(k))
        print *, 'dAph (',k,')', f%dAph(k), dfnum%dAph(k), relerr(dfnum%dAph(k), f%dAph(k))
        print *, 'dhth (',k,')', f%dhth(k), dfnum%dhth(k), relerr(dfnum%dhth(k), f%dhth(k))
        print *, 'dhph (',k,')', f%dhph(k), dfnum%dhph(k), relerr(dfnum%dhph(k), f%dhph(k))
        print *, 'dBmod(',k,')', f%dBmod(k), dfnum%dBmod(k), relerr(dfnum%dBmod(k), f%dBmod(k))
    enddo

    dx = 1d-8*z0(4)
    call get_val(f, z0(4) + .5d0*dx)
    dvparnum(4) = f%vpar
    dHnum(4) = f%H
    dpthnum(4) = f%pth
    call get_val(f, z0(4) - .5d0*dx)
    dvparnum(4) = (dvparnum(4) - f%vpar)/dx
    dHnum(4) = (dHnum(4) - f%H)/dx
    dpthnum(4) = (dpthnum(4) - f%pth)/dx
    call get_derivatives(f, z0(4))
    
    do k=1,3
        print *, 'dvpar(',k,')', f%dvpar(k), dvparnum(k), relerr(dvparnum(k), f%dvpar(k))
        print *, 'dH   (',k,')', f%dH(k), dHnum(k), relerr(dHnum(k), f%dH(k))
        print *, 'dpth (',k,')', f%dpth(k), dpthnum(k), relerr(dpthnum(k), f%dpth(k))
    enddo
    print *, 'dvpardpph', f%dvpar(4), dvparnum(4), relerr(dvparnum(4), f%dvpar(4))
    print *, 'dHdpph   ', f%dH(4), dHnum(4), relerr(dHnum(4), f%dH(4))
    print *, 'dpthdpph ', f%dpth(4), dpthnum(4), relerr(dpthnum(4), f%dpth(4))

    do i = 1,3
        do j = 1,3
            call der2(f, z0(1:3), z0(4), i, j)
        enddo
    enddo

    ! TODO: second ders in pphi and mixed


    call orbit_sympl_init(euler1, f, z0, 1.0d0, 1, 1d-12, 0, 0)
    call test_jac1(euler1)
    call test_newton(euler1, f)

    call orbit_sympl_init(euler2, f, z0, 1.0d0, 1, 1d-12, 0, 0)
    call test_jac2(euler2)
    call test_newton2(euler2, f)

    call orbit_sympl_init(midpoint, f, z0, 1.0d0, 1, 1d-12, 0, 3)
    call test_jac_midpoint(midpoint)
end subroutine do_test

end program test_magfie
