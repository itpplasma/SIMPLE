module spline_3d_interpolation
    !> Shared 3D spline interpolation module for canonical and Boozer coordinates
    !> Provides a common interface for 3D tensor product spline evaluation
    !> 
    !> This is a simplified version that only handles function evaluation,
    !> not derivatives. Derivatives require the derf1/derf2 arrays from
    !> canonical_coordinates_mod which creates circular dependencies.
    
    use, intrinsic :: iso_fortran_env, only: dp => real64
    
    implicit none
    private
    
    public :: spline_3d_evaluate, spline_3d_evaluate_vector
    
contains

    !> Evaluate 3D tensor product spline for a single quantity
    !> Uses nested Horner's method for efficient evaluation
    subroutine spline_3d_evaluate(coeffs, ns_s, ns_tp, &
                                  ds, dtheta, dphi, &
                                  is, i_theta, i_phi, &
                                  result)
        real(dp), intent(in) :: coeffs(:,:,:,:)  ! (ns_s+1, ns_tp+1, ns_tp+1, grid)
        integer, intent(in) :: ns_s, ns_tp
        real(dp), intent(in) :: ds, dtheta, dphi
        integer, intent(in) :: is, i_theta, i_phi
        real(dp), intent(out) :: result
        
        ! Local variables for nested interpolation
        real(dp) :: stp_temp(ns_tp+1, ns_tp+1)
        real(dp) :: sp_temp(ns_tp+1)
        integer :: k, nstp
        
        nstp = ns_tp + 1
        
        ! Interpolation over s (radial direction)
        stp_temp(1:nstp, 1:nstp) = coeffs(ns_s+1, :, :, is)
        
        do k = ns_s, 1, -1
            stp_temp(1:nstp, 1:nstp) = coeffs(k, :, :, is) + ds * stp_temp(1:nstp, 1:nstp)
        end do
        
        ! Interpolation over theta (poloidal direction)
        sp_temp(1:nstp) = stp_temp(nstp, 1:nstp)
        
        do k = ns_tp, 1, -1
            sp_temp(1:nstp) = stp_temp(k, 1:nstp) + dtheta * sp_temp(1:nstp)
        end do
        
        ! Interpolation over phi (toroidal direction)
        result = sp_temp(nstp)
        
        do k = ns_tp, 1, -1
            result = sp_temp(k) + dphi * result
        end do
        
    end subroutine spline_3d_evaluate

    !> Evaluate 3D tensor product spline for multiple quantities (vector)
    !> This is optimized for cases where multiple physical quantities 
    !> (e.g., sqg, B_theta, B_phi) are stored in the same tensor
    subroutine spline_3d_evaluate_vector(coeffs, nquant, ns_s, ns_tp, &
                                        ds, dtheta, dphi, &
                                        is, i_theta, i_phi, &
                                        results)
        real(dp), intent(in) :: coeffs(:,:,:,:,:)  ! (nquant, ns_s+1, ns_tp+1, ns_tp+1, grid)
        integer, intent(in) :: nquant, ns_s, ns_tp
        real(dp), intent(in) :: ds, dtheta, dphi
        integer, intent(in) :: is, i_theta, i_phi
        real(dp), intent(out) :: results(:)  ! (nquant)
        
        ! Local variables for nested interpolation
        real(dp) :: stp_temp(nquant, ns_tp+1, ns_tp+1)
        real(dp) :: sp_temp(nquant, ns_tp+1)
        integer :: k, nstp
        
        nstp = ns_tp + 1
        
        ! Interpolation over s (radial direction)
        stp_temp(:, 1:nstp, 1:nstp) = coeffs(:, ns_s+1, :, :, is)
        
        do k = ns_s, 1, -1
            stp_temp(:, 1:nstp, 1:nstp) = coeffs(:, k, :, :, is) + ds * stp_temp(:, 1:nstp, 1:nstp)
        end do
        
        ! Interpolation over theta (poloidal direction)
        sp_temp(:, 1:nstp) = stp_temp(:, nstp, 1:nstp)
        
        do k = ns_tp, 1, -1
            sp_temp(:, 1:nstp) = stp_temp(:, k, 1:nstp) + dtheta * sp_temp(:, 1:nstp)
        end do
        
        ! Interpolation over phi (toroidal direction)
        results(:) = sp_temp(:, nstp)
        
        do k = ns_tp, 1, -1
            results(:) = sp_temp(:, k) + dphi * results(:)
        end do
        
    end subroutine spline_3d_evaluate_vector

end module spline_3d_interpolation