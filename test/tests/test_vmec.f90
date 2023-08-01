module test_vmec
    use read_wout_mod, only: tosuvspace, nfp
    use util, only: pi
    use funit

    implicit none

double precision :: errtol
character(*), parameter :: filename = 'wout.nc'

contains

    @test
    subroutine compare_libstell()
        real(8) :: s, theta, varphi, ds


        ! For testing libstell routines
        real(8) :: gsqrt, bsupu, bsupv, jsupu, jsupv, lam

        ! For testing own routines
        real(8) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota,&
          R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
        real(8) :: Bctrvr_vartheta, Bctrvr_varphi, Bcovar_vartheta, &
          Bcovar_varphi, sqg, Bcovar_r, Bsqrt_varphi, Bsqrt_vartheta

        integer, parameter :: ntheta = 100
        integer :: k

        write(*,*) "Before own"
        call init_own
        write(*,*) "Before theirs"
        call init_libstell

        errtol = 1e-14
        ds=1e-1
        
        do s = 0d0,1,ds
            do varphi = 0d0, 2d0 * pi
                do k = 1,ntheta
                    theta = k*2d0*pi*1.0/ntheta

                    ! VMEC coordinate v=nfp*varphi where nfp is the number of field periods
                    ! (e.g. nfp=5 in W7-X, built from 5 identical segments toroidally)
                    call tosuvspace(s, theta, nfp*varphi, gsqrt, bsupu, bsupv, jsupu, jsupv, lam)
                    call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota, &
                    sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi, &
                    Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
                    
                    !cccc Assert outputs cccc
                    @assertEqual(lam, alam, tolerance=errtol, message="Lambda mismatch between SIMPLE and STELLOPT.")
                    @assertEqual(Bctrvr_vartheta, bsupu/nfp, tolerance=errtol)
                    @assertEqual(Bctrvr_varphi, bsupv/nfp, tolerance=errtol)
                    !Bctrvr_vartheta == bsupu ausser nfp? see https://princetonuniversity.github.io/STELLOPT/VMEC.html
                    !@assertEqual(varphi, , tolerance=errtol) !v=nfp*varphi
                    !@assertEqual(dA_theta_ds, ,tolerance=errtol)
                    !@assertEqual(gsqrt, , tolerance=errtol)
                    
                    @assertTrue(relerr(dA_theta_ds, A_theta*ds) < errtol, message = "Disparity between da_theta_ds and A_theta ds!")
                    @assertTrue(relerr(dA_phi_ds, A_phi*ds) < errtol, message = "Disparity between da_phi_ds and A_phi ds!")
                    
                    !write(101, *) lam
                    !write(102, *) alam

                    !zeta  == varphi == phi(libstell) == v/nfp
                    !no j
                    !dA_theta_ds comp. A_theta ds
                    
                end do !k
            end do !varphi
        end do !s

    end subroutine compare_libstell

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

end module test_vmec
