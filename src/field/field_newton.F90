module field_newton
  !> Module for field-agnostic Newton iterations
  !> Provides Newton solvers for coordinate transformations
  
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_base, only: magnetic_field_t
  
  implicit none
  private
  
  public :: newton_theta_from_canonical
  
contains
  
  !> Newton iteration to find field-specific theta from canonical theta
  !> Solves: theta + Lambda(s, theta, phi) = vartheta_canonical
  !> where Lambda is the stream function that depends on the field type
  subroutine newton_theta_from_canonical(field, s, vartheta_canonical, varphi, &
                                          theta, converged, iter_count)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: s              ! Normalized flux coordinate
    real(dp), intent(in) :: vartheta_canonical  ! Canonical poloidal angle
    real(dp), intent(in) :: varphi         ! Toroidal angle
    real(dp), intent(inout) :: theta       ! Field-specific theta (initial guess on input)
    logical, intent(out) :: converged      ! Convergence flag
    integer, intent(out), optional :: iter_count  ! Number of iterations
    
    ! Local variables
    real(dp), parameter :: epserr = 1.0e-12_dp  ! Convergence tolerance
    integer, parameter :: max_iter = 100
    real(dp) :: alam, dl_dt, deltheta
    integer :: iter
    
    ! Initialize
    converged = .false.
    
    ! Newton iteration loop
    do iter = 1, max_iter
      ! Get stream function and its theta derivative for current field
      call get_stream_function(field, s, theta, varphi, alam, dl_dt)
      
      ! Newton update: solve vartheta = theta + Lambda
      ! f(theta) = theta + Lambda(theta) - vartheta = 0
      ! f'(theta) = 1 + dLambda/dtheta
      deltheta = (vartheta_canonical - theta - alam) / (1.0_dp + dl_dt)
      theta = theta + deltheta
      
      ! Check convergence
      if (abs(deltheta) < epserr) then
        converged = .true.
        exit
      end if
    end do
    
    if (present(iter_count)) iter_count = iter
    
    ! Warn if not converged
    if (.not. converged) then
      print *, 'WARNING: Newton iteration for theta did not converge'
      print *, '  s =', s, ', vartheta =', vartheta_canonical, ', varphi =', varphi
      print *, '  Final theta =', theta, ', deltheta =', deltheta
    end if
    
  end subroutine newton_theta_from_canonical
  
  
  !> Get stream function Lambda and its derivative for a given field
  !> This abstracts the field-specific stream function evaluation
  subroutine get_stream_function(mag_field, s, theta, varphi, alam, dl_dt)
    use field, only: vmec_field_t
#ifdef GVEC_AVAILABLE
    use field, only: gvec_field_t
    use vmec_field_adapter, only: vmec_lambda_interpolate_with_field
#else
    use vmec_field_eval, only: vmec_lambda_interpolate_with_field
#endif
    
    class(magnetic_field_t), intent(in) :: mag_field
    real(dp), intent(in) :: s, theta, varphi
    real(dp), intent(out) :: alam   ! Stream function Lambda
    real(dp), intent(out) :: dl_dt  ! dLambda/dtheta
    
    ! Dispatch based on field type
    select type (mag_field)
    type is (vmec_field_t)
      ! For VMEC, use the existing lambda interpolation
      call vmec_lambda_interpolate_with_field(mag_field, s, theta, varphi, alam, dl_dt)
    
#ifdef GVEC_AVAILABLE
    type is (gvec_field_t)
      ! For GVEC, Lambda is available as LA
      call vmec_lambda_interpolate_with_field(mag_field, s, theta, varphi, alam, dl_dt)
#endif
    
    class default
      ! For other fields, stream function may not be defined
      alam = 0.0_dp
      dl_dt = 0.0_dp
    end select
    
  end subroutine get_stream_function
  
end module field_newton