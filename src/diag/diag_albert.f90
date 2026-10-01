#ifdef USE_FORTPLOT
module diag_albert
!> Diagnostic routines for Albert canonical coordinate system
!> Provides contour plots of vector potential and magnetic field strength

use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortplot, only: figure_t
use field_can_albert, only: Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc, &
    n_r, n_th, n_phi, xmin, xmax

implicit none
private

    public :: plot_albert_contours, make_albert_contour

contains

subroutine plot_albert_contours()
    !> Generate contour plots of Aph_of_xc, hth_of_xc, hph_of_xc, and Bmod_of_xc
    !> over theta and phi for three radial slices (inner, middle, outer)

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
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Inner slice - hth_of_xc contour
    contour_data = hth_of_xc(i_r_inner, :, :)
    filename = "albert_hth_inner_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hth_of_xc contour at s=", s_inner, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Inner slice - hph_of_xc contour
    contour_data = hph_of_xc(i_r_inner, :, :)
    filename = "albert_hph_inner_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hph_of_xc contour at s=", s_inner, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Inner slice - Bmod_of_xc contour
    contour_data = Bmod_of_xc(i_r_inner, :, :)
    filename = "albert_Bmod_inner_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc contour at s=", s_inner, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Middle slice - Aph_of_xc contour
    contour_data = Aph_of_xc(i_r_middle, :, :)
    filename = "albert_Aph_middle_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc contour at s=", s_middle, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Middle slice - hth_of_xc contour
    contour_data = hth_of_xc(i_r_middle, :, :)
    filename = "albert_hth_middle_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hth_of_xc contour at s=", s_middle, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Middle slice - hph_of_xc contour
    contour_data = hph_of_xc(i_r_middle, :, :)
    filename = "albert_hph_middle_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hph_of_xc contour at s=", s_middle, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Middle slice - Bmod_of_xc contour
    contour_data = Bmod_of_xc(i_r_middle, :, :)
    filename = "albert_Bmod_middle_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc contour at s=", s_middle, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Outer slice - Aph_of_xc contour
    contour_data = Aph_of_xc(i_r_outer, :, :)
    filename = "albert_Aph_outer_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc contour at s=", s_outer, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Outer slice - hth_of_xc contour
    contour_data = hth_of_xc(i_r_outer, :, :)
    filename = "albert_hth_outer_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hth_of_xc contour at s=", s_outer, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Outer slice - hph_of_xc contour
    contour_data = hph_of_xc(i_r_outer, :, :)
    filename = "albert_hph_outer_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "hph_of_xc contour at s=", s_outer, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

    ! Outer slice - Bmod_of_xc contour
    contour_data = Bmod_of_xc(i_r_outer, :, :)
    filename = "albert_Bmod_outer_contour.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc contour at s=", s_outer, " (Albert)"
    
    call save_albert_contour(th_array, ph_array, contour_data, &
        trim(filename), trim(plot_title_str))

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

end subroutine plot_albert_contours


    subroutine make_albert_contour(theta, phi, data, plot_title, fig)
        ! data follows the Albert field convention: data(theta_index, phi_index).
        real(dp), contiguous, intent(in) :: theta(:), phi(:)
        real(dp), intent(in) :: data(:, :)
        character(len=*), intent(in) :: plot_title
        type(figure_t), intent(inout) :: fig
        real(dp), allocatable :: image_data(:, :)

        call fig%initialize(width=1000, height=800)
        call fig%set_xlabel('theta')
        call fig%set_ylabel('phi')
        call fig%set_title(plot_title)
        call fig%grid(enabled=.true.)
        ! Fortplot stores grid values in image order: z(phi_index, theta_index).
        image_data = transpose(data)
        call fig%add_contour(theta, phi, image_data)
        call fig%colorbar()
    end subroutine make_albert_contour

    subroutine save_albert_contour(theta, phi, data, filename, plot_title)
        real(dp), contiguous, intent(in) :: theta(:), phi(:)
        real(dp), intent(in) :: data(:, :)
        character(len=*), intent(in) :: filename, plot_title
        type(figure_t) :: fig

        call make_albert_contour(theta, phi, data, plot_title, fig)
        call fig%savefig(filename)
    end subroutine save_albert_contour

end module diag_albert
#else
module diag_albert
    implicit none
    private
    public :: plot_albert_contours
contains
    subroutine plot_albert_contours()
        print *, "Warning: plot_albert_contours requires fortplot "// &
            "(disabled for this build)"
    end subroutine plot_albert_contours
end module diag_albert
#endif
