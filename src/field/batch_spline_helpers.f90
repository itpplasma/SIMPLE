!> Module providing shared helper functions for batch spline operations
module batch_spline_helpers
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, evaluate_batch_splines_3d
    implicit none
    
    public :: evaluate_batch_splines_3d_single
    
contains
    
    !> Evaluate single quantity from batch splines
    subroutine evaluate_batch_splines_3d_single(spl, x, iq, y)
        type(BatchSplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: iq
        real(dp), intent(out) :: y
        
        real(dp) :: y_batch(spl%num_quantities)
        
        call evaluate_batch_splines_3d(spl, x, y_batch)
        y = y_batch(iq)
    end subroutine evaluate_batch_splines_3d_single
    
end module batch_spline_helpers