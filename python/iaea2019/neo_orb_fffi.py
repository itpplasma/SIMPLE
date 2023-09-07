from fffi import fortran_library, fortran_module

libneo_orb = fortran_library(
    'simple', compiler={'name': 'gfortran', 'version': 9})
            # compiler={'name': 'ifort', 'version': 18})

neo_orb = fortran_module(libneo_orb, 'neo_orb_global')
neo_orb.fdef("""
  subroutine init_field(ans_s, ans_tp, amultharm, aintegmode)
    integer :: ans_s, ans_tp, amultharm, aintegmode
  end

  subroutine init_params(Z_charge, m_mass, E_kin, adtau, adtaumax, arelerr)
    integer :: Z_charge, m_mass
    double precision :: E_kin, adtau, adtaumax
    double precision :: arelerr
  end

  subroutine init_integrator(z0)
    double precision, dimension(:), intent(in) :: z0
  end

  subroutine timestep_z(z, ierr)
    double precision, dimension(:) :: z
    integer :: ierr
  end

  subroutine timestep_sympl_z(z, ierr)
    double precision, dimension(:) :: z
    integer :: ierr
  end

  subroutine spline_vmec(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
    double precision, intent(in) :: s, theta, varphi
    double precision, intent(out) :: A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
  end

  subroutine field_vmec(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds, aiota, sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi, Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
    double precision :: s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota, R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
    double precision :: Bctrvr_vartheta,Bctrvr_varphi,Bcovar_r,Bcovar_vartheta,Bcovar_varphi,sqg
  end
""")

new_vmec_stuff = fortran_module(libneo_orb, 'new_vmec_stuff_mod')
new_vmec_stuff.fdef("""
    double precision :: rmajor, h_theta, h_phi
    integer :: nsurfm,nstrm,nper,kpar
""")


cut_detector = fortran_module(libneo_orb, 'cut_detector_global')
cut_detector.fdef("""
    subroutine init(z)
      double precision, dimension(:), intent(in) :: z
    end

    subroutine trace_to_cut(z, var_cut, cut_type, ierr)
      double precision, dimension(:), intent(inout) :: z
      double precision, dimension(:), intent(inout) :: var_cut
      integer, intent(out) :: cut_type
      integer, intent(out) :: ierr
    end
""")

libneo_orb.fdef("""
  subroutine can_to_vmec(r,vartheta_c_in,varphi_c_in,theta_vmec,varphi_vmec)
    double precision :: theta_vmec,varphi_vmec
    double precision :: r,vartheta_c_in,varphi_c_in
  end
""")
