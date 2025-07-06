module vmec_field_adapter
  !> Adapter module that provides high-level VMEC field operations
  !> This module wraps VMEC-specific calls to enable future abstraction
  !> while maintaining exact numerical compatibility with existing code
  
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_base, only: MagneticField
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
    class(MagneticField), intent(in) :: field  ! For future use when abstracted
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds
    real(dp), intent(out) :: aiota, sqg, alam
    real(dp), intent(out) :: dl_ds, dl_dt, dl_dp
    real(dp), intent(out) :: Bctrvr_vartheta, Bctrvr_varphi
    real(dp), intent(out) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi
    
    ! For now, directly call the existing VMEC routine
    ! In the future, this will dispatch based on field type
    call vmec_field(s, theta, varphi, &
                    A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                    sqg, alam, dl_ds, dl_dt, dl_dp, &
                    Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
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
    class(MagneticField), intent(in) :: field  ! For future use
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds
    
    ! For now, directly call the existing VMEC routine
    call splint_iota(s, aiota, daiota_ds)
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
    class(MagneticField), intent(in) :: field  ! For future use
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt
    
    ! For now, directly call the existing VMEC routine
    call splint_lambda(s, theta, varphi, alam, dl_dt)
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
    class(MagneticField), intent(in) :: field  ! For future use
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