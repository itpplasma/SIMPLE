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
  !
  end module bminmax_mod
  !
module find_bminmax_sub
implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)

contains

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Scan grid for B field extrema
  subroutine scan_grid_for_extrema(s, nt, np, bmin, bmax, tmin, tmax, pmin, pmax)
    use new_vmec_stuff_mod, only : nper
    use magfie_sub, only : magfie
    
    implicit none
    
    real(dp), intent(in) :: s
    integer, intent(in) :: nt, np
    real(dp), intent(out) :: bmin, bmax, tmin, tmax, pmin, pmax
    
    integer :: it, ip
    real(dp) :: bmod, sqrtg, hp, ht, twopi
    real(dp), dimension(3) :: x, bder, hcovar, hctrvr, hcurl
    
    twopi = 8.d0 * atan(1.d0)
    ht = twopi / dble(nt)
    hp = twopi / dble(nt * nper)
    
    x(1) = s
    x(2:3) = 0.d0
    
    ! Initialize with first point
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    
    tmin = 0.d0
    tmax = 0.d0
    pmin = 0.d0
    pmax = 0.d0
    bmin = bmod
    bmax = bmod
    
    ! Scan grid
    do it = 1, nt
      x(2) = ht * dble(it)
      do ip = 1, np
        x(3) = hp * dble(ip)
        
        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        
        if (bmod > bmax) then
          bmax = bmod
          tmax = x(2)
          pmax = x(3)
        elseif (bmod < bmin) then
          bmin = bmod
          tmin = x(2)
          pmin = x(3)
        endif
      enddo
    enddo
  end subroutine scan_grid_for_extrema
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Newton iteration for extremum refinement
  subroutine newton_refine_extremum(s, theta, phi, bextr)
    use magfie_sub, only : magfie
    
    implicit none
    
    real(dp), intent(in) :: s
    real(dp), intent(inout) :: theta, phi
    real(dp), intent(out) :: bextr
    
    integer, parameter :: niter = 10
    real(dp), parameter :: hdt = 1.d-3, hdp = 1.d-3
    integer :: iter
    real(dp) :: bmod, sqrtg, bt, bp, btt, btp, bpt, bpp
    real(dp) :: det, dtheta, dphi, w, ht, hp, twopi
    real(dp), dimension(3) :: x, bder, hcovar, hctrvr, hcurl
    
    twopi = 8.d0 * atan(1.d0)
    ht = twopi / 100.d0  ! Using fixed grid size from original
    hp = twopi / 100.d0
    
    x(1) = s
    
    do iter = 1, niter
      x(2) = theta
      x(3) = phi
      
      ! Get B field and derivatives at current point
      call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
      bt = bder(2) * bmod
      bp = bder(3) * bmod
      
      ! Get second derivatives using finite differences
      call compute_b_second_derivatives(s, theta, phi, hdt, hdp, &
                                        btt, btp, bpt, bpp)
      
      ! Newton step
      det = btt * bpp - btp * bpt
      dtheta = (bt * bpp - bp * btp) / det
      dphi = (bp * btt - bt * bpt) / det
      
      ! Limited step to stay within reasonable bounds
      w = min(ht / max(abs(dtheta), ht), hp / max(abs(dphi), hp))
      theta = theta - w * dtheta
      phi = phi - w * dphi
    enddo
    
    ! Get final B value
    x(2) = theta
    x(3) = phi
    call magfie(x, bextr, sqrtg, bder, hcovar, hctrvr, hcurl)
  end subroutine newton_refine_extremum
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute second derivatives of B field
  subroutine compute_b_second_derivatives(s, theta, phi, hdt, hdp, &
                                          btt, btp, bpt, bpp)
    use magfie_sub, only : magfie
    
    implicit none
    
    real(dp), intent(in) :: s, theta, phi, hdt, hdp
    real(dp), intent(out) :: btt, btp, bpt, bpp
    
    real(dp) :: bmod, sqrtg, bt_plus, bt_minus, bp_plus, bp_minus
    real(dp), dimension(3) :: x, bder, hcovar, hctrvr, hcurl
    
    x(1) = s
    
    ! Compute d²B/dθ² using central differences
    x(2) = theta + hdt
    x(3) = phi
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    bt_plus = bder(2) * bmod
    
    x(2) = theta - hdt
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    bt_minus = bder(2) * bmod
    
    btt = (bt_plus - bt_minus) / (2.d0 * hdt)
    
    ! Compute d²B/dθdφ using central differences
    x(2) = theta + hdt
    x(3) = phi
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    bp_plus = bder(3) * bmod
    
    x(2) = theta - hdt
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    bp_minus = bder(3) * bmod
    
    bpt = (bp_plus - bp_minus) / (2.d0 * hdt)
    btp = bpt  ! Mixed derivative is symmetric
    
    ! Compute d²B/dφ² using central differences
    x(2) = theta
    x(3) = phi + hdp
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    bp_plus = bder(3) * bmod
    
    x(3) = phi - hdp
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    bp_minus = bder(3) * bmod
    
    bpp = (bp_plus - bp_minus) / (2.d0 * hdp)
  end subroutine compute_b_second_derivatives
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Initialize B min/max arrays
  subroutine initialize_bminmax_arrays
    use bminmax_mod, only : prop, nsbmnx, hsbmnx, bmin_arr, bmax_arr
    
    implicit none
    
    integer :: k
    real(dp) :: s0, bmin, bmax
    
    if (.not. prop) return  ! Already initialized
    
    prop = .false.
    hsbmnx = 1.d0 / dble(nsbmnx)
    
    do k = 0, nsbmnx
      s0 = max(1.d-8, hsbmnx * dble(k))
      call find_bminmax(s0, bmin, bmax)
      bmin_arr(k) = bmin
      bmax_arr(k) = bmax
    enddo
  end subroutine initialize_bminmax_arrays
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Interpolate B min/max from arrays
  subroutine interpolate_bminmax(s, bmin, bmax)
    use bminmax_mod, only : nsbmnx, hsbmnx, bmin_arr, bmax_arr
    
    implicit none
    
    real(dp), intent(in) :: s
    real(dp), intent(out) :: bmin, bmax
    
    integer :: k
    real(dp) :: ws
    
    ws = s / hsbmnx
    k = min(nsbmnx - 1, max(0, int(ws)))
    ws = ws - dble(k)
    
    bmin = bmin_arr(k) * (1.d0 - ws) + bmin_arr(k+1) * ws
    bmax = bmax_arr(k) * (1.d0 - ws) + bmax_arr(k+1) * ws
  end subroutine interpolate_bminmax

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine find_bminmax(s, bmin, bmax)
    implicit none
    
    real(dp), intent(in) :: s
    real(dp), intent(out) :: bmin, bmax
    
    integer, parameter :: np = 100, nt = 100
    real(dp) :: tmin, tmax, pmin, pmax
    
    ! Scan grid to find initial extrema
    call scan_grid_for_extrema(s, nt, np, bmin, bmax, tmin, tmax, pmin, pmax)
    
    ! Refine minimum using Newton iteration
    call newton_refine_extremum(s, tmin, pmin, bmin)
    
    ! Refine maximum using Newton iteration
    call newton_refine_extremum(s, tmax, pmax, bmax)
  end subroutine find_bminmax
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine get_bminmax(s, bmin, bmax)
    implicit none
    
    real(dp), intent(in) :: s
    real(dp), intent(out) :: bmin, bmax
    
    ! Initialize arrays on first call
    call initialize_bminmax_arrays
    
    ! Interpolate from precomputed arrays
    call interpolate_bminmax(s, bmin, bmax)
  end subroutine get_bminmax

end module find_bminmax_sub
