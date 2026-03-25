module vmec_field_eval
  !> Module providing VMEC field evaluation functions that work with magnetic_field_t classes
  !> This module is always available regardless of GVEC support

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use field_base, only : magnetic_field_t
  use field_vmec, only : vmec_field_t
#ifdef GVEC_AVAILABLE
  use field_gvec, only : gvec_field_t
  use gvec_export_data, only : gvec_profile_a_theta, gvec_profile_a_phi, &
                               gvec_profile_da_theta_ds, gvec_profile_da_phi_ds, &
                               gvec_profile_iota, gvec_geom_r, gvec_geom_z, &
                               gvec_geom_dr_ds, gvec_geom_dr_dt, gvec_geom_dr_dp, &
                               gvec_geom_dz_ds, gvec_geom_dz_dt, gvec_geom_dz_dp, &
                               gvec_geom_lambda, gvec_geom_dlambda_ds, &
                               gvec_geom_dlambda_dt, gvec_geom_dlambda_dp
#endif
  use spline_vmec_sub

  implicit none
  private

  public :: vmec_field_evaluate, vmec_field_evaluate_with_field
  public :: vmec_iota_interpolate, vmec_iota_interpolate_with_field
  public :: vmec_lambda_interpolate, vmec_lambda_interpolate_with_field
  public :: vmec_data_interpolate, vmec_data_interpolate_with_field

contains

  !> Evaluate VMEC field with field object (boozer_converter interface)
  subroutine vmec_field_evaluate_with_field(field, s, theta, varphi, &
                                            A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                                            sqg, alam, dl_ds, dl_dt, dl_dp, &
                                            Bctrvr_vartheta, Bctrvr_varphi, &
                                            Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds
    real(dp), intent(out) :: aiota, sqg, alam
    real(dp), intent(out) :: dl_ds, dl_dt, dl_dp
    real(dp), intent(out) :: Bctrvr_vartheta, Bctrvr_varphi
    real(dp), intent(out) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi

    select type (field)
    type is (vmec_field_t)
      call vmec_field_evaluate(s, theta, varphi, &
                               A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                               sqg, alam, dl_ds, dl_dt, dl_dp, &
                               Bctrvr_vartheta, Bctrvr_varphi, &
                               Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
#ifdef GVEC_AVAILABLE
    type is (gvec_field_t)
      block
        real(dp) :: profiles(5)
        real(dp) :: geometry(18)

        call field%data%evaluate_profiles(s, profiles)
        call field%data%evaluate_geometry([s, theta, varphi], geometry)

        A_theta = profiles(gvec_profile_a_theta)
        A_phi = profiles(gvec_profile_a_phi)
        dA_theta_ds = profiles(gvec_profile_da_theta_ds)
        dA_phi_ds = profiles(gvec_profile_da_phi_ds)
        aiota = profiles(gvec_profile_iota)
        alam = geometry(gvec_geom_lambda)
        dl_ds = geometry(gvec_geom_dlambda_ds)
        dl_dt = geometry(gvec_geom_dlambda_dt)
        dl_dp = geometry(gvec_geom_dlambda_dp)
        call compute_field_components(geometry(gvec_geom_r), &
                                      geometry(gvec_geom_dr_ds), &
                                      geometry(gvec_geom_dr_dt), &
                                      geometry(gvec_geom_dr_dp), &
                                      geometry(gvec_geom_dz_ds), &
                                      geometry(gvec_geom_dz_dt), &
                                      geometry(gvec_geom_dz_dp), dA_theta_ds, &
                                      dA_phi_ds, dl_ds, dl_dt, dl_dp, sqg, &
                                      Bctrvr_vartheta, Bctrvr_varphi, Bcovar_r, &
                                      Bcovar_vartheta, Bcovar_varphi)
      end block
#endif
    class default
      error stop 'vmec_field_evaluate_with_field: Unsupported field type'
    end select
  end subroutine vmec_field_evaluate_with_field

  !> Original VMEC field evaluation using global splines (boozer_converter interface)
  subroutine vmec_field_evaluate(s, theta, varphi, &
                                  A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                                  sqg, alam, dl_ds, dl_dt, dl_dp, &
                                  Bctrvr_vartheta, Bctrvr_varphi, &
                                  Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds
    real(dp), intent(out) :: aiota, sqg, alam
    real(dp), intent(out) :: dl_ds, dl_dt, dl_dp
    real(dp), intent(out) :: Bctrvr_vartheta, Bctrvr_varphi
    real(dp), intent(out) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi
    
    ! Call the existing VMEC routine
    call vmec_field(s, theta, varphi, &
                    A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                    sqg, alam, dl_ds, dl_dt, dl_dp, &
                    Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
  end subroutine vmec_field_evaluate

  !> Interpolate rotational transform with field object
  subroutine vmec_iota_interpolate_with_field(field, s, aiota, daiota_ds)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds

    select type (field)
    type is (vmec_field_t)
      call vmec_iota_interpolate(s, aiota, daiota_ds)
#ifdef GVEC_AVAILABLE
    type is (gvec_field_t)
      block
        real(dp) :: profiles(5)
        real(dp) :: dprofiles(5)

        call field%data%evaluate_profiles(s, profiles, dprofiles)
        aiota = profiles(gvec_profile_iota)
        daiota_ds = dprofiles(gvec_profile_iota)
      end block
#endif
    class default
      error stop 'vmec_iota_interpolate_with_field: Unsupported field type'
    end select
  end subroutine vmec_iota_interpolate_with_field

  !> Original VMEC iota interpolation
  subroutine vmec_iota_interpolate(s, aiota, daiota_ds)
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds

    call splint_iota(s, aiota, daiota_ds)
  end subroutine vmec_iota_interpolate

  !> Interpolate stream function Lambda with field object
  subroutine vmec_lambda_interpolate_with_field(field, s, theta, varphi, alam, dl_dt)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt

    select type (field)
    type is (vmec_field_t)
      call vmec_lambda_interpolate(s, theta, varphi, alam, dl_dt)
#ifdef GVEC_AVAILABLE
    type is (gvec_field_t)
      block
        real(dp) :: geometry(18)

        call field%data%evaluate_geometry([s, theta, varphi], geometry)
        alam = geometry(gvec_geom_lambda)
        dl_dt = geometry(gvec_geom_dlambda_dt)
      end block
#endif
    class default
      error stop 'vmec_lambda_interpolate_with_field: Unsupported field type'
    end select
  end subroutine vmec_lambda_interpolate_with_field

  !> Original VMEC lambda interpolation
  subroutine vmec_lambda_interpolate(s, theta, varphi, alam, dl_dt)
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt

    call splint_lambda(s, theta, varphi, alam, dl_dt)
  end subroutine vmec_lambda_interpolate

  !> Interpolate all VMEC data with field object
  subroutine vmec_data_interpolate_with_field(field, s, theta, varphi, &
                                                A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                                R, Z, alam, &
                                                dR_ds, dR_dt, dR_dp, &
                                                dZ_ds, dZ_dt, dZ_dp, &
                                                dl_ds, dl_dt, dl_dp)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds
    real(dp), intent(out) :: aiota, R, Z, alam
    real(dp), intent(out) :: dR_ds, dR_dt, dR_dp
    real(dp), intent(out) :: dZ_ds, dZ_dt, dZ_dp
    real(dp), intent(out) :: dl_ds, dl_dt, dl_dp

    select type (field)
    type is (vmec_field_t)
      call vmec_data_interpolate(s, theta, varphi, &
                                 A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                 R, Z, alam, &
                                 dR_ds, dR_dt, dR_dp, &
                                 dZ_ds, dZ_dt, dZ_dp, &
                                 dl_ds, dl_dt, dl_dp)
#ifdef GVEC_AVAILABLE
    type is (gvec_field_t)
      block
        real(dp) :: profiles(5)
        real(dp) :: geometry(18)

        call field%data%evaluate_profiles(s, profiles)
        call field%data%evaluate_geometry([s, theta, varphi], geometry)

        A_phi = profiles(gvec_profile_a_phi)
        A_theta = profiles(gvec_profile_a_theta)
        dA_phi_ds = profiles(gvec_profile_da_phi_ds)
        dA_theta_ds = profiles(gvec_profile_da_theta_ds)
        aiota = profiles(gvec_profile_iota)
        R = geometry(gvec_geom_r)
        Z = geometry(gvec_geom_z)
        alam = geometry(gvec_geom_lambda)
        dR_ds = geometry(gvec_geom_dr_ds)
        dR_dt = geometry(gvec_geom_dr_dt)
        dR_dp = geometry(gvec_geom_dr_dp)
        dZ_ds = geometry(gvec_geom_dz_ds)
        dZ_dt = geometry(gvec_geom_dz_dt)
        dZ_dp = geometry(gvec_geom_dz_dp)
        dl_ds = geometry(gvec_geom_dlambda_ds)
        dl_dt = geometry(gvec_geom_dlambda_dt)
        dl_dp = geometry(gvec_geom_dlambda_dp)
      end block
#endif
    class default
      error stop 'vmec_data_interpolate_with_field: Unsupported field type'
    end select
  end subroutine vmec_data_interpolate_with_field

  !> Original VMEC data interpolation
  subroutine vmec_data_interpolate(s, theta, varphi, &
                                    A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                    R, Z, alam, &
                                    dR_ds, dR_dt, dR_dp, &
                                    dZ_ds, dZ_dt, dZ_dp, &
                                    dl_ds, dl_dt, dl_dp)
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds
    real(dp), intent(out) :: aiota, R, Z, alam
    real(dp), intent(out) :: dR_ds, dR_dt, dR_dp
    real(dp), intent(out) :: dZ_ds, dZ_dt, dZ_dp
    real(dp), intent(out) :: dl_ds, dl_dt, dl_dp

    ! Call the existing VMEC routine
    call splint_vmec_data(s, theta, varphi, &
                          A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                          R, Z, alam, &
                          dR_ds, dR_dt, dR_dp, &
                          dZ_ds, dZ_dt, dZ_dp, &
                          dl_ds, dl_dt, dl_dp)
  end subroutine vmec_data_interpolate

end module vmec_field_eval
