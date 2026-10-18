program test_diag_albert
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use diag_albert, only: make_albert_contour, plot_albert_contours
    use fortplot, only: figure_t
    use fortplot_plot_data, only: plot_data_t, PLOT_TYPE_CONTOUR
    use field_can_albert, only: Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc, &
                               n_r, n_th, n_phi, xmin, xmax
    implicit none

    call check_grid([-1.0_dp, 0.5_dp, 2.0_dp], &
                    [1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp], .false., &
                    'albert_rectangular.png', 'Rectangular analytical field')
    call check_grid([-2.0_dp, 0.5_dp, 3.0_dp], &
                    [-2.0_dp, 0.5_dp, 3.0_dp], .true., &
                    'albert_square.png', 'Square antisymmetric field')
    call check_all_diagnostics()
    print *, 'Albert numerical figure checks passed.'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write(error_unit, '(A)') 'Albert diagnostic failure: '//message
            error stop 1
        end if
    end subroutine require

    pure real(dp) function grid_field(theta, phi, square) result(value)
        real(dp), intent(in) :: theta, phi
        logical, intent(in) :: square

        if (square) then
            value = theta - phi + 0.125_dp*(theta**2 - phi**2)
        else
            value = 2.0_dp*theta - 3.0_dp*phi + &
                    0.5_dp*theta*phi + 0.25_dp*theta**2
        end if
    end function grid_field

    subroutine check_grid(theta, phi, square, filename, plot_title)
        real(dp), contiguous, intent(in) :: theta(:), phi(:)
        logical, intent(in) :: square
        character(len=*), intent(in) :: filename, plot_title
        type(figure_t), target :: fig
        type(plot_data_t), pointer :: plots(:)
        real(dp) :: data(size(theta), size(phi)), expected
        integer :: i, j

        do j = 1, size(phi)
            do i = 1, size(theta)
                data(i, j) = grid_field(theta(i), phi(j), square)
            end do
        end do
        call make_albert_contour(theta, phi, data, plot_title, fig)

        call require(fig%get_plot_count() == 1, 'expected one contour field')
        plots => fig%get_plots()
        call require(associated(plots), 'figure has no stored plot data')
        call require(size(plots) >= 1, 'stored plot data is empty')
        call require(plots(1)%plot_type == PLOT_TYPE_CONTOUR, &
                     'diagnostic is not a contour plot')
        call require(.not. plots(1)%fill_contours, 'contours unexpectedly filled')
        call require(plots(1)%use_color_levels, 'contour lines have no level colors')
        call require(allocated(plots(1)%x_grid), 'theta coordinates missing')
        call require(allocated(plots(1)%y_grid), 'phi coordinates missing')
        call require(allocated(plots(1)%z_grid), 'field values missing')
        call require(size(plots(1)%x_grid) == size(theta), 'theta size changed')
        call require(size(plots(1)%y_grid) == size(phi), 'phi size changed')
        call require(size(plots(1)%z_grid, 1) == size(phi), &
                     'stored field first dimension must follow phi')
        call require(size(plots(1)%z_grid, 2) == size(theta), &
                     'stored field second dimension must follow theta')
        call require(maxval(abs(plots(1)%x_grid - theta)) < 1.0e-12_dp, &
                     'theta coordinates changed')
        call require(maxval(abs(plots(1)%y_grid - phi)) < 1.0e-12_dp, &
                     'phi coordinates changed')
        do i = 1, size(theta)
            do j = 1, size(phi)
                expected = grid_field(theta(i), phi(j), square)
                call require(abs(plots(1)%z_grid(j, i) - expected) < 1.0e-12_dp, &
                             'analytical field is at the wrong coordinates')
            end do
        end do
        if (square) then
            call require(abs(plots(1)%z_grid(3, 1) + 5.625_dp) < 1.0e-12_dp, &
                         'square upper-left corner was transposed')
            call require(abs(plots(1)%z_grid(1, 3) - 5.625_dp) < 1.0e-12_dp, &
                         'square lower-right corner was transposed')
        end if

        call require(fig%get_width() == 1000, 'figure width changed')
        call require(fig%get_height() == 800, 'figure height changed')
        call require(allocated(fig%xlabel), 'theta label missing')
        call require(allocated(fig%ylabel), 'phi label missing')
        call require(allocated(fig%title), 'figure title missing')
        call require(fig%xlabel == 'theta', 'theta label changed')
        call require(fig%ylabel == 'phi', 'phi label changed')
        call require(fig%title == plot_title, 'figure title changed')
        call require(fig%state%grid_enabled, 'diagnostic grid is disabled')
        call require(fig%state%colorbar_enabled, 'diagnostic colorbar is disabled')
        call require(fig%state%colorbar_plot_index == 1, &
                     'colorbar is not attached to the contour field')
        call require(plots(1)%show_colorbar, 'contour excludes its colorbar')
        call fig%savefig(filename)
        nullify (plots)
    end subroutine check_grid

    pure real(dp) function component_field(component, s, theta, phi) result(value)
        integer, intent(in) :: component
        real(dp), intent(in) :: s, theta, phi
        real(dp) :: c

        c = real(component, dp)
        value = 100.0_dp*c + 5.0_dp*s + (c + s)*theta + &
                (0.5_dp*c + 2.0_dp*s)*phi + &
                (0.2_dp*c - s)*theta*phi + 0.1_dp*c*theta**2
    end function component_field

    subroutine check_all_diagnostics()
        real(dp), parameter :: theta(3) = [-1.0_dp, 0.5_dp, 2.0_dp]
        real(dp), parameter :: phi(5) = [1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp]
        real(dp), parameter :: radial(8) = &
            [0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp, 0.7_dp, 0.8_dp, 0.9_dp]
        real(dp), parameter :: expected_slices(3) = [0.3_dp, 0.5_dp, 0.7_dp]
        character(len=4), parameter :: names(4) = &
            [character(len=4) :: 'Aph', 'hth', 'hph', 'Bmod']
        character(len=11), parameter :: fields(4) = &
            [character(len=11) :: 'Aph_of_xc', 'hth_of_xc', &
                                 'hph_of_xc', 'Bmod_of_xc']
        character(len=6), parameter :: slices(3) = &
            [character(len=6) :: 'inner', 'middle', 'outer']
        type(figure_t) :: expected_figure
        real(dp) :: expected_data(3, 5)
        character(len=100) :: filename, plot_title
        integer :: i, j, k, component, slice

        n_r = 8
        n_th = 3
        n_phi = 5
        xmin = [0.2_dp, -1.0_dp, 1.0_dp]
        xmax = [0.9_dp, 2.0_dp, 3.0_dp]
        allocate (Aph_of_xc(8, 3, 5), hth_of_xc(8, 3, 5), &
                  hph_of_xc(8, 3, 5), Bmod_of_xc(8, 3, 5))
        do j = 1, 5
            do i = 1, 3
                do k = 1, 8
                    Aph_of_xc(k, i, j) = component_field(1, radial(k), theta(i), phi(j))
                    hth_of_xc(k, i, j) = component_field(2, radial(k), theta(i), phi(j))
                    hph_of_xc(k, i, j) = component_field(3, radial(k), theta(i), phi(j))
                    Bmod_of_xc(k, i, j) = &
                        component_field(4, radial(k), theta(i), phi(j))
                end do
            end do
        end do

        call plot_albert_contours()

        ! Expected slices use radial indices 2, 4, 6, independent of the index formula.
        do slice = 1, 3
            do component = 1, 4
                do j = 1, 5
                    do i = 1, 3
                        expected_data(i, j) = component_field(component, &
                            expected_slices(slice), theta(i), phi(j))
                    end do
                end do
                write(plot_title, '(A,F5.3,A)') &
                    trim(fields(component))//' contour at s=', &
                    expected_slices(slice), ' (Albert)'
                filename = 'expected_albert_'//trim(names(component))//'_'// &
                           trim(slices(slice))//'_contour.png'
                call make_albert_contour(theta, phi, expected_data, &
                                         trim(plot_title), expected_figure)
                call expected_figure%savefig(trim(filename))
            end do
        end do
        deallocate (Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc)
    end subroutine check_all_diagnostics

end program test_diag_albert
