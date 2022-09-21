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
end do



contains

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
