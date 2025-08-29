module diag_albert
!> Diagnostic routines for Albert canonical coordinate system
!> Provides contour plots of vector potential and magnetic field strength

use, intrinsic :: iso_fortran_env, only: dp => real64
use pyplot_module, only: pyplot
use field_can_albert, only: Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc, &
    n_r, n_th, n_phi, xmin, xmax

implicit none
private

public :: plot_albert_contours

contains

subroutine plot_albert_contours()
    !> Generate contour plots of Aph_of_xc, hth_of_xc, hph_of_xc, and Bmod_of_xc
    !> over theta and phi for three radial slices (inner, middle, outer)

    type(pyplot) :: plt
    integer :: i_r_inner, i_r_middle, i_r_outer
    real(dp), dimension(:), allocatable :: th_array, ph_array
    real(dp), dimension(:,:), allocatable :: contour_data
    character(len=100) :: filename
    character(len=80) :: plot_title_str
    real(dp) :: s_inner, s_middle, s_outer
    integer :: i_th, i_ph

    ! Define radial slice indices - inner (25%), middle (50%), outer (75%)
    i_r_inner = max(1, n_r / 4)
    i_r_middle = n_r / 2
    i_r_outer = max(1, 3 * n_r / 4)
    
    ! Calculate corresponding s values
    s_inner = xmin(1) + (xmax(1) - xmin(1)) * real(i_r_inner - 1, dp) / &
        real(n_r - 1, dp)
    s_middle = xmin(1) + (xmax(1) - xmin(1)) * real(i_r_middle - 1, dp) / &
        real(n_r - 1, dp)
    s_outer = xmin(1) + (xmax(1) - xmin(1)) * real(i_r_outer - 1, dp) / &
        real(n_r - 1, dp)

    ! Allocate coordinate arrays
    allocate(th_array(n_th))
    allocate(ph_array(n_phi))
    allocate(contour_data(n_th, n_phi))

    ! Create coordinate arrays
    do i_th = 1, n_th
        th_array(i_th) = xmin(2) + (xmax(2) - xmin(2)) * &
            real(i_th - 1, dp) / real(n_th - 1, dp)
    end do

    do i_ph = 1, n_phi
        ph_array(i_ph) = xmin(3) + (xmax(3) - xmin(3)) * &
            real(i_ph - 1, dp) / real(n_phi - 1, dp)
    end do

    ! Generate contour plots for each radial slice and field component
    ! Order: Aph, hth, hph, Bmod for each radial slice
    
    ! Inner slice - Aph_of_xc contour
    contour_data = Aph_of_xc(i_r_inner, :, :)
    filename = "albert_Aph_inner_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc contour at s=", s_inner, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_aph_inner.py')

    ! Inner slice - hth_of_xc contour
    contour_data = hth_of_xc(i_r_inner, :, :)
    filename = "albert_hth_inner_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hth_of_xc contour at s=", s_inner, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_hth_inner.py')

    ! Inner slice - hph_of_xc contour
    contour_data = hph_of_xc(i_r_inner, :, :)
    filename = "albert_hph_inner_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hph_of_xc contour at s=", s_inner, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_hph_inner.py')

    ! Inner slice - Bmod_of_xc contour
    contour_data = Bmod_of_xc(i_r_inner, :, :)
    filename = "albert_Bmod_inner_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc contour at s=", s_inner, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_bmod_inner.py')

    ! Middle slice - Aph_of_xc contour
    contour_data = Aph_of_xc(i_r_middle, :, :)
    filename = "albert_Aph_middle_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc contour at s=", s_middle, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_aph_middle.py')

    ! Middle slice - hth_of_xc contour
    contour_data = hth_of_xc(i_r_middle, :, :)
    filename = "albert_hth_middle_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hth_of_xc contour at s=", s_middle, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_hth_middle.py')

    ! Middle slice - hph_of_xc contour
    contour_data = hph_of_xc(i_r_middle, :, :)
    filename = "albert_hph_middle_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hph_of_xc contour at s=", s_middle, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_hph_middle.py')

    ! Middle slice - Bmod_of_xc contour
    contour_data = Bmod_of_xc(i_r_middle, :, :)
    filename = "albert_Bmod_middle_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc contour at s=", s_middle, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_bmod_middle.py')

    ! Outer slice - Aph_of_xc contour
    contour_data = Aph_of_xc(i_r_outer, :, :)
    filename = "albert_Aph_outer_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc contour at s=", s_outer, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_aph_outer.py')

    ! Outer slice - hth_of_xc contour
    contour_data = hth_of_xc(i_r_outer, :, :)
    filename = "albert_hth_outer_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hth_of_xc contour at s=", s_outer, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_hth_outer.py')

    ! Outer slice - hph_of_xc contour
    contour_data = hph_of_xc(i_r_outer, :, :)
    filename = "albert_hph_outer_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hph_of_xc contour at s=", s_outer, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_hph_outer.py')

    ! Outer slice - Bmod_of_xc contour
    contour_data = Bmod_of_xc(i_r_outer, :, :)
    filename = "albert_Bmod_outer_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc contour at s=", s_outer, " (Albert)"
    
    call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
        title=trim(plot_title_str), figsize=[10,8])
    call plt%add_contour(th_array, ph_array, contour_data, &
        linestyle='-', colorbar=.true.)
    call plt%savefig(trim(filename), pyfile='albert_bmod_outer.py')

    ! Cleanup
    deallocate(th_array, ph_array, contour_data)

    print *, "Albert coordinate diagnostic plots generated successfully!"
    print *, "Files created:"
    print *, "  Inner slice:"
    print *, "    albert_Aph_inner_contour.png, albert_hth_inner_contour.png"
    print *, "    albert_hph_inner_contour.png, albert_Bmod_inner_contour.png"
    print *, "  Middle slice:"
    print *, "    albert_Aph_middle_contour.png, albert_hth_middle_contour.png"
    print *, "    albert_hph_middle_contour.png, albert_Bmod_middle_contour.png"
    print *, "  Outer slice:"
    print *, "    albert_Aph_outer_contour.png, albert_hth_outer_contour.png"
    print *, "    albert_hph_outer_contour.png, albert_Bmod_outer_contour.png"
    print *, "Also generated Python files for reproducibility."

end subroutine plot_albert_contours

end module diag_albert