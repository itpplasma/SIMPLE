
module new_can
implicit none

contains

  subroutine init_own
    use new_vmec_stuff_mod, only : netcdffile, multharm

    netcdffile = filename
    multharm = 7

    call spline_vmec_data
  end subroutine init_own

  subroutine convert_can
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

    s = 0.6d0
    varphi = 0.2d0

    do k = 1,ntheta
      theta = k*2d0*pi*1.0/ntheta

      call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota, &
        sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi, &
        Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
        write(102, *) alam
    end do
  end subroutine convert_can
end module

program test_new_can
  use new_can
  implicit none

  call convert_can
end program test_new_can
