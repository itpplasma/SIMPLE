module vmec_field_adapter
  !> Module providing GVEC field evaluation functions with VMEC-compatible interface
  !> This module extends vmec_field_eval to support GVEC fields

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use field_base, only : magnetic_field_t
  use field, only : vmec_field_t, gvec_field_t
  use spline_vmec_sub
  use MODgvec_ReadState, only: eval_iota_r
  use MODgvec_ReadState_Vars, only: profiles_1d, sbase_prof, &
                                   X1_base_r, X2_base_r, LA_base_r, &
                                   X1_r, X2_r, LA_r

  implicit none
  private

  ! GVEC derivative constants
  integer, parameter :: DERIV_S = 1
  integer, parameter :: DERIV_THET = 2
  integer, parameter :: DERIV_ZETA = 3

  ! Export functions
  public :: vmec_field_evaluate, vmec_field_evaluate_with_field
  public :: vmec_iota_interpolate, vmec_iota_interpolate_with_field
  public :: vmec_lambda_interpolate, vmec_lambda_interpolate_with_field
  public :: vmec_data_interpolate, vmec_data_interpolate_with_field

contains

  ! Compute vector potential derivatives for GVEC field via numerical differentiation
  subroutine compute_vector_potential_derivatives_gvec(field, s, theta, varphi, dA_theta_ds, dA_phi_ds)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: dA_theta_ds, dA_phi_ds
    
    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    real(dp) :: x_pert(3), Acov_pert(3), hcov_pert(3), Bmod_pert
    real(dp), parameter :: ds = 1.0e-7_dp
    
    ! Check if GVEC state is properly initialized before evaluation
    ! This follows the same pattern as field_gvec.f90
    if (.not. allocated(profiles_1d) .or. .not. allocated(sbase_prof) .or. &
        .not. allocated(X1_r) .or. .not. allocated(X2_r) .or. .not. allocated(LA_r)) then
      ! Need to reinitialize GVEC state
      error stop 'GVEC state not initialized in vmec_field_adapter - call field%evaluate first'
    end if
    
    ! Original position
    x = [sqrt(s), theta, varphi]
    call field%evaluate(x, Acov, hcov, Bmod)
    
    ! Perturb in s-direction to get dA_theta_ds
    x_pert = [sqrt(s + ds), theta, varphi]
    call field%evaluate(x_pert, Acov_pert, hcov_pert, Bmod_pert)
    
    dA_theta_ds = (Acov_pert(2) - Acov(2)) / ds
    dA_phi_ds = (Acov_pert(3) - Acov(3)) / ds
  end subroutine compute_vector_potential_derivatives_gvec

  !> VMEC field evaluation wrapper
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

  !> Override evaluate function to handle GVEC fields (boozer_converter interface)
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

    real(dp) :: x(3), Acov(3), hcov(3), Bmod, sqgBctr(3)
    real(dp) :: daiota_ds
    integer :: deriv

    select type (field)
    type is (vmec_field_t)
      ! For VMEC fields, use the existing spline-based evaluation
      call vmec_field_evaluate(s, theta, varphi, &
                              A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                              sqg, alam, dl_ds, dl_dt, dl_dp, &
                              Bctrvr_vartheta, Bctrvr_varphi, &
                              Bcovar_r, Bcovar_vartheta, Bcovar_varphi)

    type is (gvec_field_t)
      ! For GVEC fields, we need to translate to GVEC coordinates and evaluate
      x = [sqrt(s), theta, varphi]
      
      ! Get basic field quantities
      call field%evaluate(x, Acov, hcov, Bmod)
      
      ! Extract vector potential components
      A_theta = Acov(2)
      A_phi = Acov(3)
      
      ! Get derivatives via numerical differentiation
      call compute_vector_potential_derivatives_gvec(field, s, theta, varphi, dA_theta_ds, dA_phi_ds)
      
      ! Get iota
      call vmec_iota_interpolate_with_field(field, s, aiota, daiota_ds)
      
      ! Get Lambda and its derivatives
      call vmec_lambda_interpolate_with_field(field, s, theta, varphi, alam, dl_dt)
      
      ! For dl_ds and dl_dp, we need to compute numerically or from GVEC data
      dl_ds = 0.0_dp  ! Placeholder - needs proper derivative
      dl_dp = 0.0_dp  ! Placeholder - needs proper derivative
      
      ! Magnetic field components
      ! These are covariant/contravariant components that need proper computation
      sqg = sqrt(Bmod)  ! Placeholder - needs proper Jacobian
      Bctrvr_vartheta = hcov(2) * Bmod  ! Placeholder
      Bctrvr_varphi = hcov(3) * Bmod    ! Placeholder
      Bcovar_r = Acov(1)                ! Placeholder
      Bcovar_vartheta = Acov(2)         ! Placeholder
      Bcovar_varphi = Acov(3)           ! Placeholder

    class default
      error stop 'vmec_field_evaluate_with_field: Unsupported field type'
    end select
  end subroutine vmec_field_evaluate_with_field

  !> Override iota interpolation to handle GVEC fields
  subroutine vmec_iota_interpolate_with_field(field, s, aiota, daiota_ds)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds

    select type (field)
    type is (vmec_field_t)
      call vmec_iota_interpolate(s, aiota, daiota_ds)
      
    type is (gvec_field_t)
      ! For GVEC, use eval_iota_r
      ! GVEC eval_iota_r takes only radial coordinate
      aiota = eval_iota_r(sqrt(s))
      ! For derivative, use finite difference
      block
        real(dp) :: rho, drho, iota_plus, iota_minus
        rho = sqrt(s)
        drho = 1.0e-7_dp
        iota_plus = eval_iota_r(rho + drho)
        iota_minus = eval_iota_r(rho - drho)
        daiota_ds = (iota_plus - iota_minus) / (2.0_dp * drho) / (2.0_dp * sqrt(s))  ! Chain rule
      end block
      
    class default
      error stop 'vmec_iota_interpolate_with_field: Unsupported field type'
    end select
  end subroutine vmec_iota_interpolate_with_field

  !> VMEC iota interpolation wrapper
  subroutine vmec_iota_interpolate(s, aiota, daiota_ds)
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds

    call splint_iota(s, aiota, daiota_ds)
  end subroutine vmec_iota_interpolate

  !> Override lambda interpolation to handle GVEC fields
  subroutine vmec_lambda_interpolate_with_field(field, s, theta, varphi, alam, dl_dt)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt
    
    real(dp) :: rho

    select type (field)
    type is (vmec_field_t)
      call vmec_lambda_interpolate(s, theta, varphi, alam, dl_dt)
      
    type is (gvec_field_t)
      ! For GVEC, Lambda is available as LA
      rho = sqrt(s)
      
      ! Use GVEC's Lambda evaluation following field_gvec.f90 pattern
      ! LA_base_r is the base object, LA_r is the DOF array
      block
        real(dp) :: gvec_coords(3)
        integer :: deriv_flags(2)
        
        gvec_coords = [rho, theta, varphi]
        
        ! No derivatives
        deriv_flags = [0, 0]
        alam = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
        
        ! Derivative w.r.t. theta (DERIV_THET = 2)
        deriv_flags = [0, 2]  ! [radial_deriv, angular_deriv]
        dl_dt = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
      end block
      
    class default
      error stop 'vmec_lambda_interpolate_with_field: Unsupported field type'
    end select
  end subroutine vmec_lambda_interpolate_with_field

  !> VMEC lambda interpolation wrapper
  subroutine vmec_lambda_interpolate(s, theta, varphi, alam, dl_dt)
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt

    call splint_lambda(s, theta, varphi, alam, dl_dt)
  end subroutine vmec_lambda_interpolate

  !> Override data interpolation to handle GVEC fields
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

    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    real(dp) :: daiota_ds

    select type (field)
    type is (vmec_field_t)
      call vmec_data_interpolate(s, theta, varphi, &
                                  A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                  R, Z, alam, &
                                  dR_ds, dR_dt, dR_dp, &
                                  dZ_ds, dZ_dt, dZ_dp, &
                                  dl_ds, dl_dt, dl_dp)

    type is (gvec_field_t)
      ! For GVEC fields, we need to extract these from the field evaluation
      x = [sqrt(s), theta, varphi]
      call field%evaluate(x, Acov, hcov, Bmod)
      
      ! Vector potential components
      A_theta = Acov(2)
      A_phi = Acov(3)
      
      ! Get derivatives via numerical differentiation
      call compute_vector_potential_derivatives_gvec(field, s, theta, varphi, dA_theta_ds, dA_phi_ds)
      
      ! Get iota
      call vmec_iota_interpolate_with_field(field, s, aiota, daiota_ds)
      
      ! Get Lambda and its derivatives
      call vmec_lambda_interpolate_with_field(field, s, theta, varphi, alam, dl_dt)
      
      ! For dl_ds and dl_dp, we need to compute numerically or from GVEC data
      dl_ds = 0.0_dp  ! Placeholder - needs proper derivative
      dl_dp = 0.0_dp  ! Placeholder - needs proper derivative
      
      ! Coordinate transformation values - these would need proper GVEC evaluation
      ! For now, using placeholders
      R = 1.0_dp
      Z = 0.0_dp
      dR_ds = 0.0_dp
      dR_dt = 0.0_dp
      dR_dp = 0.0_dp
      dZ_ds = 0.0_dp
      dZ_dt = 0.0_dp
      dZ_dp = 0.0_dp

    class default
      error stop 'vmec_data_interpolate_with_field: Unsupported field type'
    end select
  end subroutine vmec_data_interpolate_with_field

  !> VMEC data interpolation wrapper
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

end module vmec_field_adapter
