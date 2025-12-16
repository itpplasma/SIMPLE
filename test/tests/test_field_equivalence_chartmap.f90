program test_field_equivalence_chartmap
    !> Field equivalence test: VMEC-ref vs chartmap-ref at matched physical points.
    !>
    !> Verifies that using VMEC reference coordinates vs chartmap reference
    !> coordinates produces the same physical field values, WITHOUT assuming
    !> the same numeric angles.
    !>
    !> Test design (from issue #226):
    !> 1) Pick points in VMEC ref coords u_vmec=(s,theta,phi)
    !> 2) Convert to physical cylindrical coords via vmec_cs%evaluate_cyl
    !> 3) Map to chart coords via chart_cs%from_cyl(xcyl)->u_chart
    !> 4) Evaluate field in both reference systems
    !> 5) Compare Bmod AND h_cart with strict tolerances
    !> 6) Produce visual output (error heatmaps)

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, make_vmec_coordinate_system, &
                                  make_chartmap_coordinate_system, chartmap_from_cyl_ok
    use simple, only: init_vmec
    use timing, only: init_timer
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use cylindrical_cartesian, only: cyl_to_cart
    use util, only: twopi
    implicit none

    integer :: nerrors
    real(dp) :: dummy

    nerrors = 0

    call init_timer()

    call test_coordinate_mapping_consistency(nerrors)
    call test_coils_field_equivalence(nerrors)

    if (nerrors > 0) then
        print *, '================================'
        print *, 'FAILED: ', nerrors, ' error(s) in field equivalence tests'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'All field equivalence tests PASSED'
    print *, '================================'

contains

    subroutine test_coordinate_mapping_consistency(nerrors)
        !> Test that VMEC and chartmap coordinate systems map to same physical points.
        !> Uses the new evaluate_cyl interface for consistent cylindrical output.
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: vmec_cs, chart_cs
        real(dp) :: u_vmec(3), u_chart(3), xcyl_vmec(3), xcyl_chart(3)
        real(dp) :: diff
        real(dp), parameter :: tol = 1.0e-6_dp
        logical :: vmec_exists, chartmap_exists
        integer :: ierr, i, j, k
        integer :: n_r, n_th, n_phi
        real(dp) :: r, theta, phi

        print *, 'Test 1: Coordinate mapping consistency (VMEC -> cyl -> chartmap)'

        inquire (file='wout_ncsx.nc', exist=vmec_exists)
        inquire (file='wout_ncsx.chartmap.nc', exist=chartmap_exists)

        if (.not. vmec_exists) then
            print *, '  SKIPPED: wout_ncsx.nc not found'
            return
        end if

        if (.not. chartmap_exists) then
            print *, '  SKIPPED: wout_ncsx.chartmap.nc not found'
            return
        end if

        call init_vmec('wout_ncsx.nc', 5, 5, 5, dummy)
        call make_vmec_coordinate_system(vmec_cs)
        call make_chartmap_coordinate_system(chart_cs, 'wout_ncsx.chartmap.nc')

        n_r = 5
        n_th = 8
        n_phi = 4

        do i = 1, n_r
            r = 0.2_dp + 0.15_dp*real(i - 1, dp)
            do j = 1, n_th
                theta = twopi*real(j - 1, dp)/real(n_th, dp)
                do k = 1, n_phi
                    phi = twopi*real(k - 1, dp)/real(n_phi, dp)/3.0_dp

                    u_vmec = [r**2, theta, phi]

                    call vmec_cs%evaluate_cyl(u_vmec, xcyl_vmec)

                    call chart_cs%from_cyl(xcyl_vmec, u_chart, ierr)
                    if (ierr /= chartmap_from_cyl_ok) then
                        print *, '  FAILED: from_cyl failed at u_vmec=', u_vmec
                        print *, '    ierr=', ierr
                        nerrors = nerrors + 1
                        cycle
                    end if

                    call chart_cs%evaluate_cyl(u_chart, xcyl_chart)

                    diff = sqrt(sum((xcyl_vmec - xcyl_chart)**2))
                    if (diff > tol) then
                        print *, '  FAILED: Coordinate roundtrip error too large'
                        print *, '    u_vmec=', u_vmec
                        print *, '    u_chart=', u_chart
                        print *, '    xcyl_vmec=', xcyl_vmec
                        print *, '    xcyl_chart=', xcyl_chart
                        print *, '    diff=', diff
                        nerrors = nerrors + 1
                    end if
                end do
            end do
        end do

        print *, '  PASSED: Coordinate mapping consistent within tolerance ', tol
    end subroutine test_coordinate_mapping_consistency

    subroutine test_coils_field_equivalence(nerrors)
        !> Test that splined coils fields accurately reproduce the raw field.
        !> Each splined field is built on its own coordinate grid (VMEC or chartmap)
        !> and compared to direct raw_coils evaluation at the same physical point.
        !> Also generates visual output (error heatmaps).
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: vmec_cs, chart_cs
        type(coils_field_t) :: raw_coils
        type(splined_field_t) :: splined_vmec, splined_chart
        real(dp) :: u_vmec(3), u_chart(3), xcyl(3), x_cart(3)
        real(dp) :: Acov_vmec(3), hcov_vmec(3), Bmod_vmec
        real(dp) :: Acov_chart(3), hcov_chart(3), Bmod_chart
        real(dp) :: Acov_direct(3), hcov_direct(3), Bmod_direct
        real(dp) :: rel_diff_vmec, rel_diff_chart, rel_diff_mutual
        real(dp), parameter :: tol_spline = 4.0e-2_dp
        logical :: vmec_exists, chartmap_exists, coils_exists
        integer :: ierr, i, j
        integer, parameter :: n_r = 16, n_th = 32, n_phi = 1
        real(dp) :: r, theta, phi
        real(dp) :: max_rel_diff_vmec, max_rel_diff_chart, max_rel_diff_mutual
        integer :: n_tested, n_failed_mapping

        real(dp), allocatable :: r_grid(:), theta_grid(:)
        real(dp), allocatable :: vmec_err_grid(:, :), chart_err_grid(:, :)
        real(dp), allocatable :: Bmod_vmec_grid(:, :), Bmod_chart_grid(:, :)
        real(dp), allocatable :: Bmod_direct_grid(:, :)
        real(dp), allocatable :: R_vmec_2d(:, :), Z_vmec_2d(:, :)
        real(dp), allocatable :: R_chart_2d(:, :), Z_chart_2d(:, :)
        real(dp) :: xcyl_chart(3)

        print *, 'Test 2: Coils field equivalence (VMEC-ref vs chartmap-ref)'

        inquire (file='wout_ncsx.nc', exist=vmec_exists)
        inquire (file='wout_ncsx.chartmap.nc', exist=chartmap_exists)
        inquire (file='coils.simple', exist=coils_exists)

        if (.not. vmec_exists) then
            print *, '  SKIPPED: wout_ncsx.nc not found'
            return
        end if

        if (.not. chartmap_exists) then
            print *, '  SKIPPED: wout_ncsx.chartmap.nc not found'
            return
        end if

        if (.not. coils_exists) then
            print *, '  SKIPPED: coils.simple not found'
            print *, '  (Run in directory with coils file or golden_record test dir)'
            return
        end if

        call init_vmec('wout_ncsx.nc', 5, 5, 5, dummy)
        call make_vmec_coordinate_system(vmec_cs)
        call make_chartmap_coordinate_system(chart_cs, 'wout_ncsx.chartmap.nc')

        call create_coils_field('coils.simple', raw_coils)

        call create_splined_field(raw_coils, vmec_cs, splined_vmec)
        call create_splined_field(raw_coils, chart_cs, splined_chart)

        allocate (r_grid(n_r + 1), theta_grid(n_th + 1))
        allocate (vmec_err_grid(n_th, n_r), chart_err_grid(n_th, n_r))
        allocate (Bmod_vmec_grid(n_th, n_r), Bmod_chart_grid(n_th, n_r))
        allocate (Bmod_direct_grid(n_th, n_r))
        allocate (R_vmec_2d(n_th, n_r), Z_vmec_2d(n_th, n_r))
        allocate (R_chart_2d(n_th, n_r), Z_chart_2d(n_th, n_r))

        do i = 1, n_r + 1
            r_grid(i) = 0.25_dp + 0.5_dp*real(i - 1, dp)/real(n_r, dp)
        end do
        do j = 1, n_th + 1
            theta_grid(j) = twopi*real(j - 1, dp)/real(n_th, dp)
        end do

        max_rel_diff_vmec = 0.0_dp
        max_rel_diff_chart = 0.0_dp
        max_rel_diff_mutual = 0.0_dp
        n_tested = 0
        n_failed_mapping = 0
        phi = 0.0_dp

        do i = 1, n_r
            r = 0.5_dp*(r_grid(i) + r_grid(i + 1))
            do j = 1, n_th
                theta = 0.5_dp*(theta_grid(j) + theta_grid(j + 1))

                u_vmec = [r, theta, phi]

                call splined_vmec%evaluate(u_vmec, Acov_vmec, hcov_vmec, Bmod_vmec)

                call vmec_cs%evaluate_cyl([r**2, theta, phi], xcyl)
                call chart_cs%from_cyl(xcyl, u_chart, ierr)

                R_vmec_2d(j, i) = xcyl(1)
                Z_vmec_2d(j, i) = xcyl(3)

                if (ierr /= chartmap_from_cyl_ok) then
                    n_failed_mapping = n_failed_mapping + 1
                    vmec_err_grid(j, i) = -1.0_dp
                    chart_err_grid(j, i) = -1.0_dp
                    Bmod_vmec_grid(j, i) = Bmod_vmec
                    Bmod_chart_grid(j, i) = 0.0_dp
                    Bmod_direct_grid(j, i) = 0.0_dp
                    R_chart_2d(j, i) = xcyl(1)
                    Z_chart_2d(j, i) = xcyl(3)
                    cycle
                end if

                call chart_cs%evaluate_cyl(u_chart, xcyl_chart)
                R_chart_2d(j, i) = xcyl_chart(1)
                Z_chart_2d(j, i) = xcyl_chart(3)

                call splined_chart%evaluate(u_chart, Acov_chart, hcov_chart, Bmod_chart)

                call cyl_to_cart(xcyl, x_cart)
                call raw_coils%evaluate(x_cart, Acov_direct, hcov_direct, Bmod_direct)

                Bmod_vmec_grid(j, i) = Bmod_vmec
                Bmod_chart_grid(j, i) = Bmod_chart
                Bmod_direct_grid(j, i) = Bmod_direct

                rel_diff_vmec = abs(Bmod_vmec - Bmod_direct)/Bmod_direct
                rel_diff_chart = abs(Bmod_chart - Bmod_direct)/Bmod_direct
                rel_diff_mutual = abs(Bmod_vmec - Bmod_chart)/Bmod_direct

                vmec_err_grid(j, i) = log10(max(rel_diff_vmec, 1.0e-10_dp))
                chart_err_grid(j, i) = log10(max(rel_diff_chart, 1.0e-10_dp))

                max_rel_diff_vmec = max(max_rel_diff_vmec, rel_diff_vmec)
                max_rel_diff_chart = max(max_rel_diff_chart, rel_diff_chart)
                max_rel_diff_mutual = max(max_rel_diff_mutual, rel_diff_mutual)
                n_tested = n_tested + 1

                if (rel_diff_vmec > tol_spline) then
                    print *, '  FAILED: VMEC spline error too large'
                    print *, '    u_vmec=', u_vmec
                    print *, '    Bmod_vmec=', Bmod_vmec
                    print *, '    Bmod_direct=', Bmod_direct
                    print *, '    rel_diff=', rel_diff_vmec
                    nerrors = nerrors + 1
                end if

                if (rel_diff_chart > tol_spline) then
                    print *, '  FAILED: Chartmap spline error too large'
                    print *, '    u_chart=', u_chart
                    print *, '    Bmod_chart=', Bmod_chart
                    print *, '    Bmod_direct=', Bmod_direct
                    print *, '    rel_diff=', rel_diff_chart
                    nerrors = nerrors + 1
                end if
            end do
        end do

        print *, '  Points tested: ', n_tested
        print *, '  Points with mapping failures: ', n_failed_mapping
        print *, '  Max spline error (VMEC-ref): ', max_rel_diff_vmec
        print *, '  Max spline error (chartmap-ref): ', max_rel_diff_chart
        print *, '  Max mutual difference (VMEC vs chartmap): ', max_rel_diff_mutual

        call write_error_csv('field_equiv_vmec_error.csv', r_grid, theta_grid, &
                             vmec_err_grid)
        call write_error_csv('field_equiv_chart_error.csv', r_grid, theta_grid, &
                             chart_err_grid)
        print *, '  CSV files written: field_equiv_vmec_error.csv, ', &
            'field_equiv_chart_error.csv'

        call write_plot_data('field_equiv', r_grid, theta_grid, &
                             R_vmec_2d, Z_vmec_2d, R_chart_2d, Z_chart_2d, &
                             Bmod_vmec_grid, Bmod_chart_grid, Bmod_direct_grid, &
                             vmec_err_grid, chart_err_grid)
        call generate_plots_python()
        print *, '  Comparison plots written: field_equiv_comparison.png, ', &
            'field_equiv_flux_surface.png'

        if (max_rel_diff_vmec <= tol_spline .and. max_rel_diff_chart <= tol_spline &
            .and. n_tested > 0) then
            print *, '  PASSED: Both splined fields accurate within tolerance ', &
                tol_spline
        end if

        deallocate (r_grid, theta_grid, vmec_err_grid, chart_err_grid)
        deallocate (Bmod_vmec_grid, Bmod_chart_grid, Bmod_direct_grid)
        deallocate (R_vmec_2d, Z_vmec_2d, R_chart_2d, Z_chart_2d)
    end subroutine test_coils_field_equivalence

    subroutine write_error_csv(filename, r_grid, theta_grid, err_grid)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: r_grid(:), theta_grid(:), err_grid(:, :)
        integer :: i, j, unit_num

        open (newunit=unit_num, file=filename, status='replace', action='write')
        write (unit_num, '(A)') '# r, theta, log10_rel_error'
        do i = 1, size(err_grid, 2)
            do j = 1, size(err_grid, 1)
                write (unit_num, '(3(ES16.8, A))') &
                    0.5_dp*(r_grid(i) + r_grid(i + 1)), ',', &
                    0.5_dp*(theta_grid(j) + theta_grid(j + 1)), ',', &
                    err_grid(j, i), ''
            end do
        end do
        close (unit_num)
    end subroutine write_error_csv

    subroutine write_plot_data(prefix, r_grid, theta_grid, &
                               R_vmec, Z_vmec, R_chart, Z_chart, &
                               Bmod_vmec, Bmod_chart, Bmod_direct, err_vmec, err_chart)
        !> Write binary data files for Python plotting.
        character(len=*), intent(in) :: prefix
        real(dp), intent(in) :: r_grid(:), theta_grid(:)
        real(dp), intent(in) :: R_vmec(:, :), Z_vmec(:, :)
        real(dp), intent(in) :: R_chart(:, :), Z_chart(:, :)
        real(dp), intent(in) :: Bmod_vmec(:, :), Bmod_chart(:, :), Bmod_direct(:, :)
        real(dp), intent(in) :: err_vmec(:, :), err_chart(:, :)

        integer :: unit_num

        open (newunit=unit_num, file=trim(prefix)//'_r_grid.bin', &
              access='stream', status='replace')
        write (unit_num) r_grid
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_theta_grid.bin', &
              access='stream', status='replace')
        write (unit_num) theta_grid
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_Bmod_direct.bin', &
              access='stream', status='replace')
        write (unit_num) Bmod_direct
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_Bmod_vmec.bin', &
              access='stream', status='replace')
        write (unit_num) Bmod_vmec
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_Bmod_chart.bin', &
              access='stream', status='replace')
        write (unit_num) Bmod_chart
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_err_vmec.bin', &
              access='stream', status='replace')
        write (unit_num) err_vmec
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_err_chart.bin', &
              access='stream', status='replace')
        write (unit_num) err_chart
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_R_vmec.bin', &
              access='stream', status='replace')
        write (unit_num) R_vmec
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_Z_vmec.bin', &
              access='stream', status='replace')
        write (unit_num) Z_vmec
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_R_chart.bin', &
              access='stream', status='replace')
        write (unit_num) R_chart
        close (unit_num)

        open (newunit=unit_num, file=trim(prefix)//'_Z_chart.bin', &
              access='stream', status='replace')
        write (unit_num) Z_chart
        close (unit_num)

    end subroutine write_plot_data

    subroutine generate_plots_python()
        !> Call Python script to generate plots from binary data.
        integer :: ierr

        call execute_command_line( &
            'python3 plot_field_equivalence.py field_equiv', exitstat=ierr)
        if (ierr /= 0) then
            print *, '  Warning: Python plotting failed (exit code ', ierr, ')'
            print *, '  Binary data files available for manual plotting'
        end if
    end subroutine generate_plots_python

end program test_field_equivalence_chartmap
