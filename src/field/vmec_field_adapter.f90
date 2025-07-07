module vmec_field_adapter

  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_base, only: MagneticField
  use field, only: VmecField, GvecField
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

  public :: vmec_field_evaluate, vmec_field_evaluate_with_field
  public :: vmec_iota_interpolate, vmec_iota_interpolate_with_field
  public :: vmec_lambda_interpolate, vmec_lambda_interpolate_with_field
  public :: vmec_data_interpolate, vmec_data_interpolate_with_field

contains

  ! Compute vector potential derivatives for GVEC field via numerical differentiation
  subroutine compute_vector_potential_derivatives_gvec(field, s, theta, varphi, dA_theta_ds, dA_phi_ds)
    class(MagneticField), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: dA_theta_ds, dA_phi_ds

    ! Local variables for numerical differentiation
    real(dp), parameter :: h_s = 1.0e-6_dp
    real(dp) :: x_eval(3), Acov_plus(3), Acov_minus(3), hcov_dummy(3), Bmod_dummy, sqgBctr_dummy(3)
    real(dp) :: s_plus, s_minus, r_plus, r_minus
    real(dp) :: A_theta_plus, A_phi_plus, A_theta_minus, A_phi_minus
    real(dp) :: gvec_coords(3), LA_val_plus, LA_val_minus, dLA_dt_plus, dLA_dt_minus
    integer :: deriv_flags(2)

    ! Use symmetric difference for derivatives
    if (s > h_s) then
      s_plus = s + h_s
      s_minus = s - h_s
    else
      ! Forward difference near axis
      s_plus = s + h_s
      s_minus = s
    end if

    r_plus = sqrt(s_plus)
    r_minus = sqrt(s_minus)

    ! Evaluate A_theta and A_phi at s_plus
    x_eval = [r_plus, theta, varphi]
    call field%evaluate(x_eval, Acov_plus, hcov_dummy, Bmod_dummy, sqgBctr_dummy)

    ! For GVEC, A_theta depends on Lambda
    select type (field)
    type is (GvecField)
      if (allocated(LA_r)) then
        gvec_coords = [r_plus, theta, -varphi]
        deriv_flags = [0, 0]
        LA_val_plus = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
        deriv_flags = [0, DERIV_THET]
        dLA_dt_plus = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
        A_theta_plus = Acov_plus(2) / (1.0_dp + dLA_dt_plus)
      else
        A_theta_plus = Acov_plus(2)
      end if
      A_phi_plus = Acov_plus(3)
    class default
      A_theta_plus = Acov_plus(2)
      A_phi_plus = Acov_plus(3)
    end select

    ! Evaluate A_theta and A_phi at s_minus
    x_eval = [r_minus, theta, varphi]
    call field%evaluate(x_eval, Acov_minus, hcov_dummy, Bmod_dummy, sqgBctr_dummy)

    select type (field)
    type is (GvecField)
      if (allocated(LA_r)) then
        gvec_coords = [r_minus, theta, -varphi]
        deriv_flags = [0, 0]
        LA_val_minus = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
        deriv_flags = [0, DERIV_THET]
        dLA_dt_minus = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
        A_theta_minus = Acov_minus(2) / (1.0_dp + dLA_dt_minus)
      else
        A_theta_minus = Acov_minus(2)
      end if
      A_phi_minus = Acov_minus(3)
    class default
      A_theta_minus = Acov_minus(2)
      A_phi_minus = Acov_minus(3)
    end select

    ! Compute derivatives
    if (s > h_s) then
      dA_theta_ds = (A_theta_plus - A_theta_minus) / (2.0_dp * h_s)
      dA_phi_ds = (A_phi_plus - A_phi_minus) / (2.0_dp * h_s)
    else
      dA_theta_ds = (A_theta_plus - A_theta_minus) / h_s
      dA_phi_ds = (A_phi_plus - A_phi_minus) / h_s
    end if

  end subroutine compute_vector_potential_derivatives_gvec

  ! Evaluate VMEC magnetic field at given coordinates (overloaded version with field object)
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
    ! Variables for GVEC field evaluation
    real(dp) :: gvec_coords(3), LA_val, dLA_ds, dLA_dt, dLA_dp
    ! Variables removed: x1_val, x2_val, dx1_ds, dx2_ds (unused)
    integer :: deriv_flags(2)

    ! Check field type and dispatch accordingly
    select type (field)
    type is (VmecField)
      ! Use the existing optimized routine
      call vmec_field(s, theta, varphi, &
                      A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                      sqg, alam, dl_ds, dl_dt, dl_dp, &
                      Bctrvr_vartheta, Bctrvr_varphi, &
                      Bcovar_r, Bcovar_vartheta, Bcovar_varphi)

    type is (GvecField)
      ! For GVEC, compute the required quantities
      r = sqrt(s)
      ds_dr = 2.0_dp * r

      ! Prepare coordinates for GVEC evaluation
      gvec_coords(1) = r
      gvec_coords(2) = theta
      gvec_coords(3) = -varphi  ! GVEC uses zeta = -phi

      ! Get Lambda and its derivatives
      if (allocated(LA_r)) then
        deriv_flags = [0, 0]
        LA_val = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

        deriv_flags = [DERIV_S, 0]
        dLA_ds = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r) / ds_dr

        deriv_flags = [0, DERIV_THET]
        dLA_dt = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

        deriv_flags = [0, DERIV_ZETA]
        dLA_dp = -LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)  ! Negative because phi = -zeta
      else
        LA_val = 0.0_dp
        dLA_ds = 0.0_dp
        dLA_dt = 0.0_dp
        dLA_dp = 0.0_dp
      end if

      ! Set Lambda outputs
      alam = LA_val
      dl_ds = dLA_ds
      dl_dt = dLA_dt
      dl_dp = dLA_dp

      ! Get iota
      aiota = eval_iota_r(r)

      ! Call field evaluate to get magnetic field components
      x(1) = r
      x(2) = theta
      x(3) = varphi
      call field%evaluate(x, Acov, hcov, Bmod, sqgBctr)

      ! Extract vector potential - GVEC gives it in different form
      A_theta = Acov(2) / (1.0_dp + dLA_dt)  ! Approximate
      A_phi = Acov(3)

      ! Derivatives of vector potential - compute via numerical differentiation
      call compute_vector_potential_derivatives_gvec(field, s, theta, varphi, dA_theta_ds, dA_phi_ds)

      ! Magnetic field components
      Bcovar_s_over_ds_dr = hcov(1) * Bmod
      Bcovar_r = Bcovar_s_over_ds_dr * ds_dr
      Bcovar_vartheta = hcov(2) * Bmod
      Bcovar_varphi = hcov(3) * Bmod

      ! Contravariant components and sqrt(g)
      ! Use |B|² = B^i * B_i to extract sqrt(g) properly
      ! We have: sqgBctr = sqrt(g) * B^i and Bcovar = B_i
      ! Since B^r = 0 in flux coordinates, |B|² = B^θ*B_θ + B^φ*B_φ
      ! Therefore: |B|² = (sqgBctr(2)/sqrt(g)) * Bcovar_vartheta + (sqgBctr(3)/sqrt(g)) * Bcovar_varphi
      ! Rearranging: sqrt(g) = (sqgBctr(2)*Bcovar_vartheta + sqgBctr(3)*Bcovar_varphi) / |B|²
      sqg = (sqgBctr(2)*Bcovar_vartheta + sqgBctr(3)*Bcovar_varphi) / (Bmod*Bmod)

      ! Now extract the contravariant components
      Bctrvr_vartheta = sqgBctr(2) / sqg
      Bctrvr_varphi = sqgBctr(3) / sqg

    class default
      ! For other field types, use the generic evaluate interface
      r = sqrt(s)
      ds_dr = 2.0_dp * r

      ! Prepare coordinates for field evaluation
      x(1) = r
      x(2) = theta
      x(3) = varphi

      ! Call the field's evaluate method
      call field%evaluate(x, Acov, hcov, Bmod, sqgBctr)

      ! Extract vector potential components
      A_theta = 0.0_dp
      A_phi = Acov(3)

      ! For derivatives, this is a simplified version (full implementation would need numerical derivatives)
      dA_theta_ds = 0.0_dp
      dA_phi_ds = 0.0_dp  ! Would need numerical differentiation

      ! Extract magnetic field components
      Bcovar_s_over_ds_dr = hcov(1) * Bmod
      Bcovar_r = Bcovar_s_over_ds_dr * ds_dr
      Bcovar_vartheta = hcov(2) * Bmod
      Bcovar_varphi = hcov(3) * Bmod

      ! Contravariant components from sqgBctr
      sqg = (sqgBctr(2)*Bcovar_vartheta + sqgBctr(3)*Bcovar_varphi) / (Bmod*Bmod)
      Bctrvr_vartheta = sqgBctr(2) / sqg
      Bctrvr_varphi = sqgBctr(3) / sqg

      ! Rotational transform not directly available from field interface
      aiota = 0.0_dp  ! Placeholder

      ! Stream function Lambda and derivatives
      alam = 0.0_dp
      dl_ds = 0.0_dp
      dl_dt = 0.0_dp
      dl_dp = 0.0_dp
    end select
  end subroutine vmec_field_evaluate_with_field

  ! Evaluate VMEC magnetic field at given coordinates (version without field object for backward compatibility)
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

    call the existing VMEC routine
    call vmec_field(s, theta, varphi, &
                    A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                    sqg, alam, dl_ds, dl_dt, dl_dp, &
                    Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
  end subroutine vmec_field_evaluate

  ! Interpolate rotational transform (with field object)
  subroutine vmec_iota_interpolate_with_field(field, s, aiota, daiota_ds)
    class(MagneticField), intent(in) :: field
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds
    ! Local variables
    real(dp) :: r, h_s

    select type (field)
    type is (VmecField)
      ! Use the existing optimized routine
      call splint_iota(s, aiota, daiota_ds)

    type is (GvecField)
      ! For GVEC, access the iota profile using eval_iota_r
      r = sqrt(s)
      aiota = eval_iota_r(r)

      ! For derivative, use numerical differentiation
      h_s = 1.0e-6_dp
      if (s > h_s) then
        daiota_ds = (eval_iota_r(sqrt(s + h_s)) - eval_iota_r(sqrt(s - h_s))) / (2.0_dp * h_s)
      else
        ! One-sided derivative near axis
        daiota_ds = (eval_iota_r(sqrt(s + h_s)) - aiota) / h_s
      end if

    class default
      ! For other field types, iota not directly available
      aiota = 0.0_dp
      daiota_ds = 0.0_dp
    end select
  end subroutine vmec_iota_interpolate_with_field

  ! Interpolate rotational transform (without field object)
  subroutine vmec_iota_interpolate(s, aiota, daiota_ds)
    real(dp), intent(in) :: s
    real(dp), intent(out) :: aiota, daiota_ds

    call the existing VMEC routine
    call splint_iota(s, aiota, daiota_ds)
  end subroutine vmec_iota_interpolate

  ! Interpolate lambda (stream function) with field object
  subroutine vmec_lambda_interpolate_with_field(field, s, theta, varphi, alam, dl_dt)
    class(MagneticField), intent(in) :: field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt
    ! Local variables
    real(dp) :: r, gvec_coords(3)
    integer :: deriv_flags(2)

    select type (field)
    type is (VmecField)
      ! Use the existing optimized routine
      call splint_lambda(s, theta, varphi, alam, dl_dt)

    type is (GvecField)
      ! For GVEC, Lambda (stream function) is available as LA; evaluate at the given coordinates
      if (.not. allocated(LA_r)) then
        ! GVEC state not loaded - return zeros
        alam = 0.0_dp
        dl_dt = 0.0_dp
        return
      end if

      r = sqrt(s)
      gvec_coords(1) = r
      gvec_coords(2) = theta  ! This is theta_vmec, not theta*
      gvec_coords(3) = -varphi  ! GVEC uses zeta = -phi

      ! Evaluate Lambda
      deriv_flags = [0, 0]
      alam = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

      ! Evaluate dLambda/dtheta
      deriv_flags = [0, DERIV_THET]
      dl_dt = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)

    class default
      ! For other field types, lambda not directly available
      alam = 0.0_dp
      dl_dt = 0.0_dp
    end select
  end subroutine vmec_lambda_interpolate_with_field

  ! Interpolate lambda (stream function) without field object
  subroutine vmec_lambda_interpolate(s, theta, varphi, alam, dl_dt)
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam, dl_dt

    call the existing VMEC routine
    call splint_lambda(s, theta, varphi, alam, dl_dt)
  end subroutine vmec_lambda_interpolate

  ! Interpolate complete VMEC data with field object
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
    ! Local variables
    real(dp) :: r_loc, gvec_coords(3), LA_val, dLA_ds_r, dLA_dt_loc, dLA_dz
    real(dp) :: x1_val, x2_val, dx1_ds, dx1_dt, dx1_dz, dx2_ds, dx2_dt, dx2_dz
    real(dp) :: ds_dr_loc, x_eval(3), Acov_loc(3), hcov_loc(3), Bmod_loc, sqgBctr_loc(3)
    integer :: deriv_flags(2)

    select type (field)
    type is (VmecField)
      ! Use the existing optimized routine
      call splint_vmec_data(s, theta, varphi, &
                            A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                            R, Z, alam, &
                            dR_ds, dR_dt, dR_dp, &
                            dZ_ds, dZ_dt, dZ_dp, &
                            dl_ds, dl_dt, dl_dp)

    type is (GvecField)
      ! For GVEC, get R, Z and other quantities
      r_loc = sqrt(s)
      ds_dr_loc = 2.0_dp * r_loc

      ! Prepare coordinates for GVEC evaluation
      gvec_coords(1) = r_loc
      gvec_coords(2) = theta
      gvec_coords(3) = -varphi  ! GVEC uses zeta = -phi

      ! Check if GVEC state is initialized
      if (.not. allocated(X1_r) .or. .not. allocated(X2_r) .or. .not. allocated(LA_r)) then
        ! Return zero values
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
        return
      end if

      ! Get R (X1) and Z (X2) coordinates
      deriv_flags = [0, 0]
      x1_val = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)
      x2_val = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

      R = x1_val
      Z = x2_val

      ! Get derivatives of R
      deriv_flags = [DERIV_S, 0]
      dx1_ds = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)
      deriv_flags = [0, DERIV_THET]
      dx1_dt = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)
      deriv_flags = [0, DERIV_ZETA]
      dx1_dz = X1_base_r%evalDOF_x(gvec_coords, deriv_flags, X1_r)

      dR_ds = dx1_ds / ds_dr_loc
      dR_dt = dx1_dt
      dR_dp = -dx1_dz  ! phi = -zeta

      ! Get derivatives of Z
      deriv_flags = [DERIV_S, 0]
      dx2_ds = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)
      deriv_flags = [0, DERIV_THET]
      dx2_dt = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)
      deriv_flags = [0, DERIV_ZETA]
      dx2_dz = X2_base_r%evalDOF_x(gvec_coords, deriv_flags, X2_r)

      dZ_ds = dx2_ds / ds_dr_loc
      dZ_dt = dx2_dt
      dZ_dp = -dx2_dz  ! phi = -zeta

      ! Get Lambda and its derivatives
      deriv_flags = [0, 0]
      LA_val = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
      alam = LA_val

      deriv_flags = [DERIV_S, 0]
      dLA_ds_r = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
      dl_ds = dLA_ds_r / ds_dr_loc

      deriv_flags = [0, DERIV_THET]
      dLA_dt_loc = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
      dl_dt = dLA_dt_loc

      deriv_flags = [0, DERIV_ZETA]
      dLA_dz = LA_base_r%evalDOF_x(gvec_coords, deriv_flags, LA_r)
      dl_dp = -dLA_dz  ! phi = -zeta

      ! Get iota
      aiota = eval_iota_r(r_loc)

      ! Get vector potential by calling field evaluate
      x_eval(1) = r_loc
      x_eval(2) = theta
      x_eval(3) = varphi
      call field%evaluate(x_eval, Acov_loc, hcov_loc, Bmod_loc, sqgBctr_loc)

      ! Extract vector potential
      A_theta = Acov_loc(2) / (1.0_dp + dLA_dt_loc)  ! Approximate
      A_phi = Acov_loc(3)

      ! Derivatives of vector potential - compute via numerical differentiation
      call compute_vector_potential_derivatives_gvec(field, s, theta, varphi, dA_theta_ds, dA_phi_ds)

    class default
      ! For other field types, compute or approximate these values (placeholder implementation)
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

  ! Interpolate complete VMEC data without field object
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

    call the existing VMEC routine
    call splint_vmec_data(s, theta, varphi, &
                          A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                          R, Z, alam, &
                          dR_ds, dR_dt, dR_dp, &
                          dZ_ds, dZ_dt, dZ_dp, &
                          dl_ds, dl_dt, dl_dp)
  end subroutine vmec_data_interpolate

end module vmec_field_adapter
