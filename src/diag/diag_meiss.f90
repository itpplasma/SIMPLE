module diag_meiss
!> Diagnostic routines for field_can_meiss.f90
!> Provides visualization and analysis tools for canonical coordinate transformations

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_meiss, only: rh_can, grid_indices_t, n_r, xmin, xmax, h_r
use fortplot, only: figure, plot, savefig, xlabel, ylabel, title

implicit none
private

public :: plot_rh_can_vs_rc

contains

subroutine plot_rh_can_vs_rc(i_th_in, i_phi_in, filename)
    !> Plot rh_can over r_c for fixed i_th and i_phi indices
    !> 
    !> @param i_th_in  Theta grid index (default: 1)
    !> @param i_phi_in Phi grid index (default: 1)  
    !> @param filename Output filename (default: "diag_meiss.pdf")
    
    integer, intent(in), optional :: i_th_in, i_phi_in
    character(len=*), intent(in), optional :: filename
    
    ! Local variables
    integer :: i_th, i_phi
    integer :: i, n_points
    real(dp), dimension(:), allocatable :: r_c_array, dz1_vals, dz2_vals
    real(dp), dimension(2) :: z, dz
    character(len=100) :: output_file, output_file1, output_file2
    character(len=50) :: plot_title
    
    ! Set defaults
    i_th = 1
    i_phi = 1
    if (present(i_th_in)) i_th = i_th_in
    if (present(i_phi_in)) i_phi = i_phi_in
    
    output_file = "diag_meiss.pdf"
    if (present(filename)) output_file = trim(filename)
    
    ! Create separate filenames for the two plots
    output_file1 = trim(output_file(1:len_trim(output_file)-4)) // "_dz1.pdf"
    output_file2 = trim(output_file(1:len_trim(output_file)-4)) // "_dz2.pdf"
    
    ! Create r_c array - use same grid as field_can_meiss
    n_points = n_r
    allocate(r_c_array(n_points))
    allocate(dz1_vals(n_points))
    allocate(dz2_vals(n_points))
    
    ! Fill r_c array
    do i = 1, n_points
        r_c_array(i) = xmin(1) + (xmax(1) - xmin(1)) * real(i-1, dp) / real(n_points-1, dp)
    end do
    
    ! Initialize z (lam_phi=0, chi_gauge=0 for diagnostic purposes)
    z = [0.0_dp, 0.0_dp]
    
    ! Compute both components of rh_can
    do i = 1, n_points
        call rh_can(r_c_array(i), z, dz, i_th, i_phi)
        dz1_vals(i) = dz(1)  ! dz(1) = -hr/hp
        dz2_vals(i) = dz(2)  ! dz(2) = Ar + Ap*dz(1)
    end do
    
    ! Create two separate plots for better visualization
    
    ! First plot: dz(1) = -hr/hp
    call figure()
    call plot(r_c_array, dz1_vals, label="dz(1) = -hr/hp", linestyle="b-")
    call xlabel("r_c")
    call ylabel("dz(1) = -hr/hp")
    write(plot_title, '(A,I0,A,I0,A)') "dz(1) vs r_c (i_th=", i_th, ", i_phi=", i_phi, ")"
    call title(trim(plot_title))
    
    ! Save first plot
    call savefig(trim(output_file1))
    
    ! Second plot: dz(2) = Ar + Ap*dz(1)  
    call figure()
    call plot(r_c_array, dz2_vals, label="dz(2) = Ar + Ap*dz(1)", linestyle="r-")
    call xlabel("r_c")
    call ylabel("dz(2) = Ar + Ap*dz(1)")
    write(plot_title, '(A,I0,A,I0,A)') "dz(2) vs r_c (i_th=", i_th, ", i_phi=", i_phi, ")"
    call title(trim(plot_title))
    
    ! Save second plot
    call savefig(trim(output_file2))
    
    ! Cleanup
    deallocate(r_c_array, dz1_vals, dz2_vals)
    
    print *, "Diagnostic plots saved to: ", trim(output_file1), " and ", trim(output_file2)
    
end subroutine plot_rh_can_vs_rc

end module diag_meiss