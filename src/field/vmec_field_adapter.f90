module vmec_field_adapter
  !> Adapter module that provides high-level VMEC field operations
  !> This module wraps VMEC-specific calls to enable future abstraction
  !> while maintaining exact numerical compatibility with existing code
  
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_base, only: MagneticField
  use field, only: VmecField, GvecField
  use spline_vmec_sub
  
  implicit none
  private
  
  public :: vmec_field_evaluate, vmec_field_evaluate_with_field
  public :: vmec_iota_interpolate, vmec_iota_interpolate_with_field
  public :: vmec_lambda_interpolate, vmec_lambda_interpolate_with_field
  public :: vmec_data_interpolate, vmec_data_interpolate_with_field
  
contains
  
  !> Evaluate VMEC magnetic field at given coordinates
  !> This replaces direct calls to vmec_field subroutine
  !> Overloaded version with field object for future abstraction
  subroutine vmec_field_evaluate_with_field(field, s, theta, varphi, &
                                  A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                                  sqg, alam, dl_ds, dl_dt, dl_dp, &
                                  Bctrvr_vartheta, Bctrvr_varphi, &
                                  Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    class(MagneticField), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds
    real(dp), intent(out) :: aiota, sqg, alam
    real(dp), intent(out) :: dl_ds, dl_dt, dl_dp
    real(dp), intent(out) :: Bctrvr_vartheta, Bctrvr_varphi
    real(dp), intent(out) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi
    
    ! Local variables for field evaluation
    real(dp) :: x(3), Acov(3), hcov(3), Bmod, sqgBctr(3)
    real(dp) :: r, ds_dr, Bcovar_s_over_ds_dr
    
    ! Check field type and dispatch accordingly
    select type (field)
    type is (VmecField)
      ! For VMEC, use the existing optimized routine
      call vmec_field(s, theta, varphi, &
                      A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                      sqg, alam, dl_ds, dl_dt, dl_dp, &
                      Bctrvr_vartheta, Bctrvr_varphi, &
                      Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    
    class default
      ! For other field types, use the generic evaluate interface
      ! Convert from (s, theta, phi) to (r, theta, phi)
      r = sqrt(s)
      ds_dr = 2.0_dp * r
      
      ! Prepare coordinates for field evaluation
      x(1) = r
      x(2) = theta
      x(3) = varphi
      
      ! Call the field's evaluate method
      call field%evaluate(x, Acov, hcov, Bmod, sqgBctr)
      
      ! Extract vector potential components
      ! Acov(1) is A_s * ds/dr, so A_theta = 0 (no A_s in VMEC representation)
      A_theta = 0.0_dp
      A_phi = Acov(3)
      
      ! For derivatives, we need to make additional calls
      ! This is a simplified version - full implementation would need numerical derivatives
      dA_theta_ds = 0.0_dp
      dA_phi_ds = 0.0_dp  ! Would need numerical differentiation
      
      ! Extract magnetic field components
      ! Convert from normalized h to B components
      Bcovar_s_over_ds_dr = hcov(1) * Bmod
      Bcovar_r = Bcovar_s_over_ds_dr * ds_dr
      Bcovar_vartheta = hcov(2) * Bmod
      Bcovar_varphi = hcov(3) * Bmod
      
      ! Contravariant components from sqgBctr
      ! sqgBctr(2) = sqrt(g) * B^theta, sqgBctr(3) = sqrt(g) * B^phi
      ! We need to extract sqrt(g) to get B^theta and B^phi
      ! For now, use a simplified approach
      sqg = abs(sqgBctr(2) / (Bcovar_vartheta / Bmod + 1e-10))
      Bctrvr_vartheta = sqgBctr(2) / sqg
      Bctrvr_varphi = sqgBctr(3) / sqg
      
      ! Rotational transform - not directly available from field interface
      ! Would need to be computed from field line following or stored separately
      aiota = 0.0_dp  ! Placeholder
      
      ! Stream function Lambda and derivatives
      alam = 0.0_dp
      dl_ds = 0.0_dp
      dl_dt = 0.0_dp
      dl_dp = 0.0_dp
    end select
  end subroutine vmec_field_evaluate_with_field
  
  !> Evaluate VMEC magnetic field at given coordinates
  !> This replaces direct calls to vmec_field subroutine
  !> Version without field object for backward compatibility
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
    
    ! For now, directly call the existing VMEC routine
    call vmec_field(s, theta, varphi, &
                    A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                    sqg, alam, dl_ds, dl_dt, dl_dp, &
                    Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
  end subroutine vmec_field_evaluate
  
  !> Interpolate rotational transform (with field object)
  !> This replaces direct calls to splint_iota
  subroutine vmec_iota_interpolate_with_field(field, s, aiota, daiota_ds)
    class(MagneticField), intent(in) :: field
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds
    
    select type (field)
    type is (VmecField)
      ! For VMEC, use the existing optimized routine
      call splint_iota(s, aiota, daiota_ds)
    
    type is (GvecField)
      ! For GVEC, iota evaluation would require access to GVEC internals
      ! which are not easily accessible through the field interface
      ! For now, return placeholder values
      aiota = 0.0_dp
      daiota_ds = 0.0_dp
    
    class default
      ! For other field types, iota not directly available
      ! Would need to compute from field line following
      aiota = 0.0_dp
      daiota_ds = 0.0_dp
    end select
  end subroutine vmec_iota_interpolate_with_field
  
  !> Interpolate rotational transform (without field object)
  !> This replaces direct calls to splint_iota
  subroutine vmec_iota_interpolate(s, aiota, daiota_ds)
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds
    
    ! For now, directly call the existing VMEC routine
    call splint_iota(s, aiota, daiota_ds)
  end subroutine vmec_iota_interpolate
  
  !> Interpolate lambda (stream function) with field object
  !> This replaces direct calls to splint_lambda
  subroutine vmec_lambda_interpolate_with_field(field, s, theta, varphi, alam, dl_dt)
    class(MagneticField), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt
    
    select type (field)
    type is (VmecField)
      ! For VMEC, use the existing optimized routine
      call splint_lambda(s, theta, varphi, alam, dl_dt)
    
    type is (GvecField)
      ! For GVEC, Lambda is evaluated as part of the field
      ! The stream function LA is available in GVEC
      ! For now, return zero as this requires more complex evaluation
      alam = 0.0_dp
      dl_dt = 0.0_dp
    
    class default
      ! For other field types, lambda not directly available
      alam = 0.0_dp
      dl_dt = 0.0_dp
    end select
  end subroutine vmec_lambda_interpolate_with_field
  
  !> Interpolate lambda (stream function) without field object
  !> This replaces direct calls to splint_lambda
  subroutine vmec_lambda_interpolate(s, theta, varphi, alam, dl_dt)
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt
    
    ! For now, directly call the existing VMEC routine
    call splint_lambda(s, theta, varphi, alam, dl_dt)
  end subroutine vmec_lambda_interpolate
  
  !> Interpolate complete VMEC data with field object
  !> This replaces direct calls to splint_vmec_data
  subroutine vmec_data_interpolate_with_field(field, s, theta, varphi, &
                                    A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                    R, Z, alam, &
                                    dR_ds, dR_dt, dR_dp, &
                                    dZ_ds, dZ_dt, dZ_dp, &
                                    dl_ds, dl_dt, dl_dp)
    class(MagneticField), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds
    real(dp), intent(out) :: aiota, R, Z, alam
    real(dp), intent(out) :: dR_ds, dR_dt, dR_dp
    real(dp), intent(out) :: dZ_ds, dZ_dt, dZ_dp
    real(dp), intent(out) :: dl_ds, dl_dt, dl_dp
    
    select type (field)
    type is (VmecField)
      ! For VMEC, use the existing optimized routine
      call splint_vmec_data(s, theta, varphi, &
                            A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                            R, Z, alam, &
                            dR_ds, dR_dt, dR_dp, &
                            dZ_ds, dZ_dt, dZ_dp, &
                            dl_ds, dl_dt, dl_dp)
    
    class default
      ! For other field types, we need to compute or approximate these values
      ! This is a placeholder implementation
      A_phi = 0.0_dp
      A_theta = 0.0_dp
      dA_phi_ds = 0.0_dp
      dA_theta_ds = 0.0_dp
      aiota = 0.0_dp
      R = 0.0_dp
      Z = 0.0_dp
      alam = 0.0_dp
      dR_ds = 0.0_dp
      dR_dt = 0.0_dp
      dR_dp = 0.0_dp
      dZ_ds = 0.0_dp
      dZ_dt = 0.0_dp
      dZ_dp = 0.0_dp
      dl_ds = 0.0_dp
      dl_dt = 0.0_dp
      dl_dp = 0.0_dp
    end select
  end subroutine vmec_data_interpolate_with_field
  
  !> Interpolate complete VMEC data without field object
  !> This replaces direct calls to splint_vmec_data
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
    
    ! For now, directly call the existing VMEC routine
    call splint_vmec_data(s, theta, varphi, &
                          A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                          R, Z, alam, &
                          dR_ds, dR_dt, dR_dp, &
                          dZ_ds, dZ_dt, dZ_dp, &
                          dl_ds, dl_dt, dl_dp)
  end subroutine vmec_data_interpolate
  
end module vmec_field_adapter