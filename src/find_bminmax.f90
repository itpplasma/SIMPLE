  !
  module bminmax_mod
  implicit none

  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)
  !
  logical :: prop=.true.
  !
  integer, parameter :: nsbmnx=100
  real(dp)                      :: hsbmnx
  real(dp), dimension(0:nsbmnx) :: bmin_arr,bmax_arr
    real(dp) :: extremum_angles(2, 2, 0:nsbmnx)
  !
  end module bminmax_mod
  !
module find_bminmax_sub
implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)

contains

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine find_bminmax(s, bmin, bmax, angles)
  !
  use new_vmec_stuff_mod, only : nper
  !
  implicit none
  !
  integer,          parameter :: np=100, nt=100
  !
        integer :: it, ip
        real(dp), intent(in) :: s
        real(dp), intent(out) :: bmin, bmax
        real(dp), intent(out), optional :: angles(2, 2)
        real(dp) :: bmod, hp, ht, twopi
  real(dp) :: tmin,tmax,pmin,pmax
  real(dp), dimension(3) :: x,bder
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
        call refine_bextremum(s, tmin, pmin, ht, hp, bmin)
  !
        call refine_bextremum(s, tmax, pmax, ht, hp, bmax)
        if (present(angles)) then
            angles(:, 1) = [tmin, pmin]
            angles(:, 2) = [tmax, pmax]
        end if
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
  real(dp) :: sqrtg
  real(dp), dimension(3) :: hcovar, hctrvr, hcurl
  !
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  !
  end subroutine getbmod
  !
  !-----------
  !
    end subroutine find_bminmax

    subroutine refine_bextremum(s, theta, phi, ht, hp, bextr)
        use magfie_sub, only: magfie
  !
  implicit none
  !
  integer,          parameter :: niter=10
        integer :: iter
        real(dp), intent(in) :: s, ht, hp
        real(dp), intent(inout) :: theta, phi
        real(dp), intent(out) :: bextr
        real(dp) :: x(3), bmod, bder(3), sqrtg, hcovar(3), hctrvr(3), hcurl(3)
  real(dp), parameter :: hdt=1.d-3, hdp=1.d-3
        real(dp) :: bt, bp, btt, btp, bpt, bpp
  real(dp) :: det,dtheta,dphi,w
  !
        x(1) = s
  do iter=1,niter
    x(3)=phi
    x(2)=theta
  !
            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  !
    bt=bder(2)*bmod
    bp=bder(3)*bmod
            if (max(abs(bt*ht), abs(bp*hp)) <= &
                16.d0*epsilon(bmod)*abs(bmod)) exit
  !
    x(3)=phi
    x(2)=theta+hdt
  !
            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  !
    btt=bder(2)*bmod
    bpt=bder(3)*bmod
    x(2)=theta-hdt
  !
            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  !
    btt=(btt-bder(2)*bmod)/(2.d0*hdt)
    bpt=(bpt-bder(3)*bmod)/(2.d0*hdt)
    btp=bpt
    x(2)=theta
    x(3)=phi+hdp
  !
            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  !
    bpp=bder(3)*bmod
  !
    x(3)=phi-hdp
  !
            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  !
    bpp=(bpp-bder(3)*bmod)/(2.d0*hdp)
    det=btt*bpp-btp*bpt
            ! Axisymmetric fields have a singular two-angle Hessian. Refine the
            ! remaining angular direction instead of dividing by a zero determinant.
            if (abs(det) > epsilon(det)*max(abs(btt*bpp), abs(btp*bpt))) then
                dtheta = (bt*bpp - bp*btp)/det
                dphi = (bp*btt - bt*bpt)/det
            else if (abs(btt) > abs(bpp)) then
                dtheta = bt/btt
                dphi = 0.d0
            else if (abs(bpp) > 0.d0) then
                dtheta = 0.d0
                dphi = bp/bpp
            else
                exit
            end if
    w=min(ht/max(abs(dtheta),ht),hp/max(abs(dphi),hp))
    theta=theta-w*dtheta
    phi=phi-w*dphi
  enddo
  !
        ! Evaluate the final stationary point, not the last phi-hdp stencil point.
        x(2:3) = [theta, phi]
        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  bextr=bmod
  !
    end subroutine refine_bextremum
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine init_bminmax_arrays
        use bminmax_mod, only: prop, nsbmnx, hsbmnx, bmin_arr, bmax_arr, &
                               extremum_angles

    integer :: k
    real(dp) :: s0

    if (.not. prop) return
    prop = .false.
    hsbmnx = 1.d0/dble(nsbmnx)
    ! The cache is defined in the coordinate system of the active magfie
    ! backend. The driver initializes it serially after selecting the tracing
    ! backend, before parallel orbit tracing reads it.
    do k = 0, nsbmnx
        s0 = max(1.d-8, hsbmnx*dble(k))
            call find_bminmax(s0, bmin_arr(k), bmax_arr(k), extremum_angles(:, :, k))
    end do
  end subroutine init_bminmax_arrays

  subroutine get_bminmax(s, bmin, bmax)
        use new_vmec_stuff_mod, only: nper
        use bminmax_mod, only: prop, nsbmnx, hsbmnx, extremum_angles

    real(dp), intent(in) :: s
    real(dp), intent(out) :: bmin, bmax
    integer :: k
        real(dp) :: angles(2, 2), ht, hp

    if (prop) call init_bminmax_arrays

        ! Interpolating extrema values can undershoot the actual maximum between
        ! radial knots. Cache locations instead, then refine on the requested surface
        ! using the same active field representation as the orbit and boundary.
        k = min(nsbmnx, max(0, nint(s/hsbmnx)))
        angles = extremum_angles(:, :, k)
        ht = (8.d0*atan(1.d0))/100.d0
        hp = ht/dble(nper)
        call refine_bextremum(s, angles(1, 1), angles(2, 1), ht, hp, bmin)
        call refine_bextremum(s, angles(1, 2), angles(2, 2), ht, hp, bmax)
  end subroutine get_bminmax

end module find_bminmax_sub
