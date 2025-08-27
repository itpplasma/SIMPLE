module batch_spline_migration
    ! Module to facilitate migration from individual splines to batch splines
    ! Provides utility functions and compatibility wrappers
    
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: SplineData3D, BatchSplineData3D, &
        construct_splines_3d, construct_batch_splines_3d, &
        evaluate_splines_3d, evaluate_batch_splines_3d, &
        evaluate_splines_3d_der, evaluate_batch_splines_3d_der, &
        evaluate_splines_3d_der2, evaluate_batch_splines_3d_der2, &
        destroy_splines_3d, destroy_batch_splines_3d
    
    implicit none
    
    ! Configuration flags for gradual migration
    logical :: USE_BATCH_MEISS = .true.
    logical :: USE_BATCH_ALBERT = .true.
    logical :: USE_BATCH_COILS = .true.
    
    ! Performance monitoring
    integer :: n_individual_calls = 0
    integer :: n_batch_calls = 0
    real(dp) :: time_individual = 0.0d0
    real(dp) :: time_batch = 0.0d0
    
contains

    subroutine migrate_meiss_to_batch(spl_Ath, spl_Aph, spl_hth, spl_hph, spl_Bmod, &
                                      spl_batch)
        ! Migrate existing Meiss individual splines to batch format
        type(SplineData3D), intent(in) :: spl_Ath, spl_Aph, spl_hth, spl_hph, spl_Bmod
        type(BatchSplineData3D), intent(out) :: spl_batch
        
        real(dp), allocatable :: field_batch(:,:,:,:)
        real(dp) :: x_max(3)
        integer :: n1, n2, n3
        
        ! Get dimensions from first spline
        n1 = size(spl_Ath%coeff, 4)
        n2 = size(spl_Ath%coeff, 5)
        n3 = size(spl_Ath%coeff, 6)
        
        ! Allocate batch array
        allocate(field_batch(n1, n2, n3, 5))
        
        ! Copy coefficient data (using highest order coefficients)
        field_batch(:,:,:,1) = spl_Ath%coeff(spl_Ath%order(1), 0, 0, :, :, :)
        field_batch(:,:,:,2) = spl_Aph%coeff(spl_Aph%order(1), 0, 0, :, :, :)
        field_batch(:,:,:,3) = spl_hth%coeff(spl_hth%order(1), 0, 0, :, :, :)
        field_batch(:,:,:,4) = spl_hph%coeff(spl_hph%order(1), 0, 0, :, :, :)
        field_batch(:,:,:,5) = spl_Bmod%coeff(spl_Bmod%order(1), 0, 0, :, :, :)
        
        ! Calculate x_max from x_min, h_step, and num_points
        x_max = spl_Ath%x_min + spl_Ath%h_step * (spl_Ath%num_points - 1)
        
        ! Construct batch spline
        call construct_batch_splines_3d(spl_Ath%x_min, x_max, field_batch, &
            spl_Ath%order, spl_Ath%periodic, spl_batch)
        
        deallocate(field_batch)
        
        print *, 'Migrated Meiss field splines to batch format (5 components)'
    end subroutine migrate_meiss_to_batch

    
    subroutine evaluate_meiss_compatible(use_batch, &
        spl_Ath, spl_Aph, spl_hth, spl_hph, spl_Bmod, spl_batch, &
        x, f_Ath, f_Aph, f_hth, f_hph, f_Bmod, &
        df_Ath, df_Aph, df_hth, df_hph, df_Bmod)
        ! Compatible evaluation that can use either individual or batch splines
        logical, intent(in) :: use_batch
        type(SplineData3D), intent(in), optional :: spl_Ath, spl_Aph, spl_hth, spl_hph, spl_Bmod
        type(BatchSplineData3D), intent(in), optional :: spl_batch
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: f_Ath, f_Aph, f_hth, f_hph, f_Bmod
        real(dp), intent(out), optional :: df_Ath(3), df_Aph(3), df_hth(3), df_hph(3), df_Bmod(3)
        
        real(dp) :: y_batch(5), dy_batch(3, 5)
        real(dp) :: t_start, t_end
        
        if (use_batch) then
            if (.not. present(spl_batch)) then
                error stop 'Batch spline required when use_batch=.true.'
            end if
            
            call cpu_time(t_start)
            
            if (present(df_Ath)) then
                call evaluate_batch_splines_3d_der(spl_batch, x, y_batch, dy_batch)
                f_Ath = y_batch(1); df_Ath = dy_batch(:, 1)
                f_Aph = y_batch(2); df_Aph = dy_batch(:, 2)
                f_hth = y_batch(3); df_hth = dy_batch(:, 3)
                f_hph = y_batch(4); df_hph = dy_batch(:, 4)
                f_Bmod = y_batch(5); df_Bmod = dy_batch(:, 5)
            else
                call evaluate_batch_splines_3d(spl_batch, x, y_batch)
                f_Ath = y_batch(1)
                f_Aph = y_batch(2)
                f_hth = y_batch(3)
                f_hph = y_batch(4)
                f_Bmod = y_batch(5)
            end if
            
            call cpu_time(t_end)
            time_batch = time_batch + (t_end - t_start)
            n_batch_calls = n_batch_calls + 1
            
        else
            if (.not. present(spl_Ath)) then
                error stop 'Individual splines required when use_batch=.false.'
            end if
            
            call cpu_time(t_start)
            
            if (present(df_Ath)) then
                call evaluate_splines_3d_der(spl_Ath, x, f_Ath, df_Ath)
                call evaluate_splines_3d_der(spl_Aph, x, f_Aph, df_Aph)
                call evaluate_splines_3d_der(spl_hth, x, f_hth, df_hth)
                call evaluate_splines_3d_der(spl_hph, x, f_hph, df_hph)
                call evaluate_splines_3d_der(spl_Bmod, x, f_Bmod, df_Bmod)
            else
                call evaluate_splines_3d(spl_Ath, x, f_Ath)
                call evaluate_splines_3d(spl_Aph, x, f_Aph)
                call evaluate_splines_3d(spl_hth, x, f_hth)
                call evaluate_splines_3d(spl_hph, x, f_hph)
                call evaluate_splines_3d(spl_Bmod, x, f_Bmod)
            end if
            
            call cpu_time(t_end)
            time_individual = time_individual + (t_end - t_start)
            n_individual_calls = n_individual_calls + 1
        end if
    end subroutine evaluate_meiss_compatible

    
    subroutine report_performance_stats()
        ! Report performance statistics
        real(dp) :: avg_time_individual = 0.0_dp, avg_time_batch = 0.0_dp, speedup
        
        print *, '========================================'
        print *, 'Batch Spline Performance Report'
        print *, '========================================'
        
        if (n_individual_calls > 0) then
            avg_time_individual = time_individual / n_individual_calls
            print *, 'Individual splines:'
            print *, '  Total calls:', n_individual_calls
            print *, '  Total time:', time_individual, 's'
            print *, '  Avg time per call:', avg_time_individual*1000.0d0, 'ms'
        end if
        
        if (n_batch_calls > 0) then
            avg_time_batch = time_batch / n_batch_calls
            print *, 'Batch splines:'
            print *, '  Total calls:', n_batch_calls
            print *, '  Total time:', time_batch, 's'
            print *, '  Avg time per call:', avg_time_batch*1000.0d0, 'ms'
        end if
        
        if (n_individual_calls > 0 .and. n_batch_calls > 0) then
            speedup = avg_time_individual / avg_time_batch
            print *, '========================================'
            print *, 'Speedup: ', speedup, 'x'
            print *, 'Time saved:', (time_individual - time_batch), 's'
            print *, 'Efficiency gain:', (1.0d0 - avg_time_batch/avg_time_individual)*100.0d0, '%'
        end if
        
        print *, '========================================'
    end subroutine report_performance_stats

    
    subroutine verify_batch_equivalence(spl_individual, spl_batch, n_test_points)
        ! Verify that batch and individual splines produce identical results
        type(SplineData3D), dimension(:), intent(in) :: spl_individual
        type(BatchSplineData3D), intent(in) :: spl_batch
        integer, intent(in) :: n_test_points
        
        real(dp) :: x_test(3), y_individual(size(spl_individual)), y_batch(size(spl_individual))
        real(dp) :: max_error, avg_error, error
        integer :: i, j, n_comp
        logical :: all_match
        
        n_comp = size(spl_individual)
        max_error = 0.0d0
        avg_error = 0.0d0
        all_match = .true.
        
        do i = 1, n_test_points
            ! Generate random test point
            call random_number(x_test)
            x_test(1) = spl_individual(1)%x_min(1) + x_test(1) * &
                       (spl_individual(1)%h_step(1) * (spl_individual(1)%num_points(1) - 1))
            x_test(2) = spl_individual(1)%x_min(2) + x_test(2) * &
                       (spl_individual(1)%h_step(2) * (spl_individual(1)%num_points(2) - 1))
            x_test(3) = spl_individual(1)%x_min(3) + x_test(3) * &
                       (spl_individual(1)%h_step(3) * (spl_individual(1)%num_points(3) - 1))
            
            ! Evaluate individual splines
            do j = 1, n_comp
                call evaluate_splines_3d(spl_individual(j), x_test, y_individual(j))
            end do
            
            ! Evaluate batch spline
            call evaluate_batch_splines_3d(spl_batch, x_test, y_batch)
            
            ! Check equivalence
            do j = 1, n_comp
                error = abs(y_individual(j) - y_batch(j))
                max_error = max(max_error, error)
                avg_error = avg_error + error
                
                if (error > 1.0d-12) then
                    all_match = .false.
                    print *, 'Mismatch at point', i, 'component', j, 'error:', error
                end if
            end do
        end do
        
        avg_error = avg_error / (n_test_points * n_comp)
        
        print *, '========================================'
        print *, 'Batch Spline Equivalence Verification'
        print *, '========================================'
        print *, 'Test points:', n_test_points
        print *, 'Components:', n_comp
        print *, 'Max error:', max_error
        print *, 'Avg error:', avg_error
        if (all_match) then
            print *, 'PASSED: Batch and individual splines are equivalent'
        else
            print *, 'FAILED: Discrepancies found'
        end if
        print *, '========================================'
    end subroutine verify_batch_equivalence

    
    subroutine set_batch_mode(meiss, albert, coils)
        ! Configure which modules use batch splines
        logical, intent(in), optional :: meiss, albert, coils
        
        if (present(meiss)) USE_BATCH_MEISS = meiss
        if (present(albert)) USE_BATCH_ALBERT = albert
        if (present(coils)) USE_BATCH_COILS = coils
        
        print *, 'Batch spline configuration:'
        print *, '  Meiss:', USE_BATCH_MEISS
        print *, '  Albert:', USE_BATCH_ALBERT
        print *, '  Coils:', USE_BATCH_COILS
    end subroutine set_batch_mode

end module batch_spline_migration