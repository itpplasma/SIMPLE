module vmec_field_eval
   !> Module providing VMEC field evaluation functions that work with magnetic_field_t classes
   !> This module is always available regardless of GVEC support

   use, intrinsic :: iso_fortran_env, only: dp => real64
   use field_base, only: magnetic_field_t
   use field_vmec, only: vmec_field_t
#ifdef GVEC_AVAILABLE
   use field_gvec, only: gvec_field_t
   use gvec_export_data, only: gvec_profile_a_theta, gvec_profile_a_phi, &
                               gvec_profile_da_theta_ds, gvec_profile_da_phi_ds, &
                               gvec_profile_iota, gvec_family_logical, &
                               gvec_geom_x, gvec_geom_y, gvec_geom_z, &
                               gvec_geom_dx_ds, gvec_geom_dx_dt, gvec_geom_dx_dp, &
                               gvec_geom_dy_ds, gvec_geom_dy_dt, gvec_geom_dy_dp, &
                               gvec_geom_dz_ds, gvec_geom_dz_dt, gvec_geom_dz_dp, &
                               gvec_geom_lambda, gvec_geom_dlambda_ds, &
                               gvec_geom_dlambda_dt, gvec_geom_dlambda_dp
   use interpolate, only: evaluate_batch_splines_3d_der2
#endif
   use spline_vmec_sub

   implicit none
   private

   public :: vmec_field_evaluate, vmec_field_evaluate_with_field
   public :: vmec_iota_interpolate, vmec_iota_interpolate_with_field
   public :: vmec_lambda_interpolate, vmec_lambda_interpolate_with_field
   public :: vmec_data_interpolate, vmec_data_interpolate_with_field

contains

#ifdef GVEC_AVAILABLE
   subroutine require_logical_gvec(field, caller)
      class(gvec_field_t), intent(in) :: field
      character(*), intent(in) :: caller

      if (field%data%family /= gvec_family_logical) then
         print *, trim(caller), ': GVEC VMEC-compatibility interface only supports'
         print *, 'logical-coordinate exports. Use native GVEC field/reference'
         print *, 'coordinates for straight-fieldline runs.'
         error stop
      end if
   end subroutine require_logical_gvec

   subroutine geometry_to_vmec_primitives(geometry, &
                                          R, Z, dR_ds, dR_dt, dR_dp, &
                                          dZ_ds, dZ_dt, dZ_dp)
      real(dp), intent(in) :: geometry(16)
      real(dp), intent(out) :: R, Z
      real(dp), intent(out) :: dR_ds, dR_dt, dR_dp
      real(dp), intent(out) :: dZ_ds, dZ_dt, dZ_dp

      real(dp) :: x
      real(dp) :: y
      real(dp) :: inv_r

      x = geometry(gvec_geom_x)
      y = geometry(gvec_geom_y)
      Z = geometry(gvec_geom_z)
      R = sqrt(x*x + y*y)

      if (R > 1.0e-14_dp) then
         inv_r = 1.0_dp/R
         dR_ds = (x*geometry(gvec_geom_dx_ds) + y*geometry(gvec_geom_dy_ds))*inv_r
         dR_dt = (x*geometry(gvec_geom_dx_dt) + y*geometry(gvec_geom_dy_dt))*inv_r
         dR_dp = (x*geometry(gvec_geom_dx_dp) + y*geometry(gvec_geom_dy_dp))*inv_r
      else
         dR_ds = 0.0_dp
         dR_dt = 0.0_dp
         dR_dp = 0.0_dp
      end if

      dZ_ds = geometry(gvec_geom_dz_ds)
      dZ_dt = geometry(gvec_geom_dz_dt)
      dZ_dp = geometry(gvec_geom_dz_dp)
   end subroutine geometry_to_vmec_primitives

   subroutine logical_geometry_from_field(field, x, R, Z, alam, dl_ds, dl_dt, dl_dp, &
                                          dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp)
      class(gvec_field_t), intent(in) :: field
      real(dp), intent(in) :: x(3)
      real(dp), intent(out) :: R, Z, alam
      real(dp), intent(out) :: dl_ds, dl_dt, dl_dp
      real(dp), intent(out) :: dR_ds, dR_dt, dR_dp
      real(dp), intent(out) :: dZ_ds, dZ_dt, dZ_dp

      real(dp) :: geometry(16)
      real(dp) :: geom_values(3)
      real(dp) :: geom_derivs(3, 3)
      real(dp) :: geom_der2(6, 3)

      if (field%has_logical_geom) then
         call evaluate_batch_splines_3d_der2(field%logical_geom_spline, x, geom_values, &
                                             geom_derivs, geom_der2)
         R = geom_values(1)
         Z = geom_values(2)
         alam = geom_values(3)
         dR_ds = geom_derivs(1, 1)
         dR_dt = geom_derivs(2, 1)
         dR_dp = geom_derivs(3, 1)
         dZ_ds = geom_derivs(1, 2)
         dZ_dt = geom_derivs(2, 2)
         dZ_dp = geom_derivs(3, 2)
         dl_ds = geom_derivs(1, 3)
         dl_dt = geom_derivs(2, 3)
         dl_dp = geom_derivs(3, 3)
         return
      end if

      call field%data%evaluate_geometry(x, geometry)
      alam = geometry(gvec_geom_lambda)
      dl_ds = geometry(gvec_geom_dlambda_ds)
      dl_dt = geometry(gvec_geom_dlambda_dt)
      dl_dp = geometry(gvec_geom_dlambda_dp)
      call geometry_to_vmec_primitives(geometry, R, Z, dR_ds, dR_dt, dR_dp, &
                                       dZ_ds, dZ_dt, dZ_dp)
   end subroutine logical_geometry_from_field
#endif

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
            real(dp) :: R, Z
            real(dp) :: dR_ds, dR_dt, dR_dp
            real(dp) :: dZ_ds, dZ_dt, dZ_dp

            call require_logical_gvec(field, 'vmec_field_evaluate_with_field')
            call field%data%evaluate_profiles(s, profiles)

            A_theta = profiles(gvec_profile_a_theta)
            A_phi = profiles(gvec_profile_a_phi)
            dA_theta_ds = profiles(gvec_profile_da_theta_ds)
            dA_phi_ds = profiles(gvec_profile_da_phi_ds)
            aiota = profiles(gvec_profile_iota)
            call logical_geometry_from_field(field, [s, theta, varphi], R, Z, alam, &
                                             dl_ds, dl_dt, dl_dp, dR_ds, dR_dt, dR_dp, &
                                             dZ_ds, dZ_dt, dZ_dp)
            call compute_field_components(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, &
                                          dZ_dp, dA_theta_ds, &
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
            real(dp) :: R, Z
            real(dp) :: dl_ds, dl_dp
            real(dp) :: dR_ds, dR_dt, dR_dp
            real(dp) :: dZ_ds, dZ_dt, dZ_dp

            call require_logical_gvec(field, 'vmec_lambda_interpolate_with_field')
            call logical_geometry_from_field(field, [s, theta, varphi], R, Z, alam, &
                                             dl_ds, dl_dt, dl_dp, dR_ds, dR_dt, dR_dp, &
                                             dZ_ds, dZ_dt, dZ_dp)
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

            call require_logical_gvec(field, 'vmec_data_interpolate_with_field')
            call field%data%evaluate_profiles(s, profiles)

            A_phi = profiles(gvec_profile_a_phi)
            A_theta = profiles(gvec_profile_a_theta)
            dA_phi_ds = profiles(gvec_profile_da_phi_ds)
            dA_theta_ds = profiles(gvec_profile_da_theta_ds)
            aiota = profiles(gvec_profile_iota)
            call logical_geometry_from_field(field, [s, theta, varphi], R, Z, alam, &
                                             dl_ds, dl_dt, dl_dp, dR_ds, dR_dt, dR_dp, &
                                             dZ_ds, dZ_dt, dZ_dp)
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
