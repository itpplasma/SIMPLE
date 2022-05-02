program test_vmec
use read_wout_mod, only: tosuvspace, nfp
use util, only: pi
implicit none

character(*), parameter :: filename = 'wout.nc'
real(8) :: s, theta, varphi

! For testing own routines
real(8) :: gsqrt, bsupu, bsupv, jsupu, jsupv, lam

! For testing own routines
real(8) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota,&
  R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
  real(8) :: Bctrvr_vartheta, Bctrvr_varphi, Bcovar_vartheta, &
  Bcovar_varphi, sqg, Bcovar_r

integer, parameter :: ntheta = 100
integer :: k

call init_own
call init_libstell

s = 0.6d0
varphi = 0.2d0

do k = 1,ntheta
  theta = k*2d0*pi*1.0/ntheta

  ! VMEC coordinate u=nfp*varphi where nfp is the number of field periods
  ! (e.g. nfp=5 in W7-X, built from 5 identical segments toroidally)
  call tosuvspace(s, theta, nfp*varphi, gsqrt, bsupu, bsupv, jsupu, jsupv, lam)
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota, &
    sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi, &
    Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
    write(101, *) lam
    write(102, *) alam

  call get_derivatives() !ToDo
end do



contains

subroutine der2(x0, pphi, i, j)
  double precision, intent(in) :: x0(3)
  integer, intent(in) :: i, j
  double precision hi, hj
  type(FieldCan) :: f00, f01, f10, f11
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

  function relerr(a, b)
    double precision :: relerr
    double precision, intent(in) :: a, b
    relerr = merge(0d0, (a - b)/b, b == 0d0)
  end function relerr

  subroutine init_own
    use new_vmec_stuff_mod, only : netcdffile, multharm

    netcdffile = filename
    multharm = 7

    call spline_vmec_data
  end subroutine init_own


  subroutine init_libstell
    use read_wout_mod, only: read_wout_file

    integer :: ierr

    call read_wout_file(filename, ierr)
  end subroutine init_libstell

end program test_vmec
