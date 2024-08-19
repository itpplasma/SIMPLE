!
  module bminmax_mod
!
  logical :: prop=.true.
!
  integer, parameter :: nsbmnx=100
  double precision                      :: hsbmnx
  double precision, dimension(0:nsbmnx) :: bmin_arr,bmax_arr
!
  end module bminmax_mod
!
module find_bminmax_sub

implicit none

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_bminmax(s,bmin,bmax)
!
  use new_vmec_stuff_mod, only : nper
!
  implicit none
!
  integer,          parameter :: np=100, nt=100
!
  integer :: it,ip,iter
  double precision :: s,bmin,bmax,bmod,hp,ht,twopi
  double precision :: tmin,tmax,pmin,pmax
  double precision, dimension(3) :: x,bder
!
  twopi=8.d0*atan(1.d0)
  ht=twopi/dble(nt)
  hp=twopi/dble(nt*nper)
  x(1)=s
  x(2:3)=0.d0
!
  call getbmod
!
  tmin=0.d0
  tmax=0.d0
  pmin=0.d0
  pmax=0.d0
  bmin=bmod
  bmax=bmod
!
  do it=1,nt
    x(2)=ht*dble(it)
    do ip=1,np
      x(3)=hp*dble(ip)
!
      call getbmod
!
      if(bmod.gt.bmax) then
        bmax=bmod
        tmax=x(2)
        pmax=x(3)
      elseif(bmod.lt.bmin) then
        bmin=bmod
        tmin=x(2)
        pmin=x(3)
      endif
!
    enddo
  enddo
!
  call newtextr(tmin,pmin,bmin)
!
  call newtextr(tmax,pmax,bmax)
!
!-----------
  contains
!-----------
!
  subroutine getbmod
  use magfie_sub, only : magfie
!
  implicit none
!
  double precision :: sqrtg
  double precision, dimension(3) :: hcovar, hctrvr, hcurl
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  end subroutine getbmod
!
!-----------
!
  subroutine newtextr(theta,phi,bextr)
!
  implicit none
!
  integer,          parameter :: niter=10
  double precision, parameter :: hdt=1.d-3, hdp=1.d-3
  double precision :: bt,bp,btt,btp,bpt,bpp,theta,phi,bextr
  double precision :: det,dtheta,dphi,w
!
  do iter=1,niter
    x(3)=phi
    x(2)=theta
!
    call getbmod
!
    bt=bder(2)*bmod
    bp=bder(3)*bmod
!
    x(3)=phi
    x(2)=theta+hdt
!
    call getbmod
!
    btt=bder(2)*bmod
    bpt=bder(3)*bmod
    x(2)=theta-hdt
!
    call getbmod
!
    btt=(btt-bder(2)*bmod)/(2.d0*hdt)
    bpt=(bpt-bder(3)*bmod)/(2.d0*hdt)
    btp=bpt
    x(2)=theta
    x(3)=phi+hdp
!
    call getbmod
!
    bpp=bder(3)*bmod
!
    x(3)=phi-hdp
!
    call getbmod
!
    bpp=(bpp-bder(3)*bmod)/(2.d0*hdp)
    det=btt*bpp-btp*bpt
    dtheta=(bt*bpp-bp*btp)/det
    dphi=(bp*btt-bt*bpt)/det
    w=min(ht/max(abs(dtheta),ht),hp/max(abs(dphi),hp))
    theta=theta-w*dtheta
    phi=phi-w*dphi
  enddo
!
  bextr=bmod
!
  end subroutine newtextr
!
!-----------
!
  end subroutine find_bminmax
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_bminmax(s,bmin,bmax)
!
  use bminmax_mod, only : prop,nsbmnx,hsbmnx,bmin_arr,bmax_arr
!
  implicit none
!
  integer :: k
  double precision :: s,bmin,bmax,ws,s0
!
  if(prop) then
    prop=.false.
    hsbmnx=1.d0/dble(nsbmnx)
    do k=0,nsbmnx
      s0=max(1.d-8,hsbmnx*dble(k))
!
      call find_bminmax(s0,bmin_arr(k),bmax_arr(k))
!
    enddo
  endif
!
  ws=s/hsbmnx
  k=min(nsbmnx-1,max(0,int(ws)))
  ws=ws-dble(k)
!
  bmin=bmin_arr(k)*(1.d0-ws)+bmin_arr(k+1)*ws
  bmax=bmax_arr(k)*(1.d0-ws)+bmax_arr(k+1)*ws
!
  end subroutine get_bminmax

end module find_bminmax_sub
