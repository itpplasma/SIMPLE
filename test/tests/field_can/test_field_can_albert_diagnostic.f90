program test_field_can_albert_diagnostic
    !> Test Albert coordinate initialization produces valid field values.
    !> Verifies physical constraints on field components and coordinate transforms.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use magfie_sub, only: ALBERT
    use velo_mod, only: isw_field_type
    use field, only: vmec_field_t, create_vmec_field
    use field_can_albert, only: init_albert, psi_inner, psi_outer, &
        psi_of_x, Ath_norm, dpsi_dr_positive
    use field_can_meiss, only: spl_field_batch, xmin, xmax, n_r, n_th, n_phi
    use interpolate, only: evaluate_batch_splines_3d
    use params, only: coord_input

    implicit none

    type(tracer_t) :: norb
    type(vmec_field_t) :: magfie
    integer :: i_r, i_th, i_phi, n_failed
    real(dp) :: x(3), x_spl(3), y_batch(5)
    real(dp) :: Bmod_min, Bmod_max
    logical :: file_exists

    n_failed = 0
    isw_field_type = ALBERT

    print *, 'Test: Albert coordinate field initialization'

    inquire(file='wout.nc', exist=file_exists)
    if (.not. file_exists) then
        print *, 'FAILED: Required VMEC file (wout.nc) not found'
        error stop 1
    end if

    coord_input = 'wout.nc'
    call create_vmec_field(magfie)
    call init_field(norb, 'wout.nc', 5, 5, 3, 0)

    ! Test 1: Grid dimensions must be positive
    print *, 'Test 1: Grid dimensions'
    if (n_r <= 0 .or. n_th <= 0 .or. n_phi <= 0) then
        print *, '  FAILED: Grid dimensions must be positive'
        print *, '    n_r=', n_r, ' n_th=', n_th, ' n_phi=', n_phi
        n_failed = n_failed + 1
    end if

    ! Test 2: psi boundaries must be ordered correctly
    print *, 'Test 2: psi boundaries'
    if (psi_inner >= psi_outer) then
        print *, '  FAILED: psi_inner must be less than psi_outer'
        print *, '    psi_inner=', psi_inner, ' psi_outer=', psi_outer
        n_failed = n_failed + 1
    end if

    ! Test 3: Ath_norm must be nonzero
    print *, 'Test 3: Ath normalization'
    if (abs(Ath_norm) < 1.0e-15_dp) then
        print *, '  FAILED: Ath_norm must be nonzero, got ', Ath_norm
        n_failed = n_failed + 1
    end if

    ! Test 4: psi_of_x must be monotonic in r (either increasing or decreasing)
    print *, 'Test 4: psi_of_x monotonicity'
    if (dpsi_dr_positive) then
        if (psi_of_x(n_r, n_th/2, n_phi/2) <= psi_of_x(1, n_th/2, n_phi/2)) then
            print *, '  FAILED: psi_of_x should increase with r when dpsi_dr_positive=.true.'
            n_failed = n_failed + 1
        end if
    else
        if (psi_of_x(n_r, n_th/2, n_phi/2) >= psi_of_x(1, n_th/2, n_phi/2)) then
            print *, '  FAILED: psi_of_x should decrease with r when dpsi_dr_positive=.false.'
            n_failed = n_failed + 1
        end if
    end if

    ! Test 5: Bmod from splines must be positive everywhere
    print *, 'Test 5: Bmod positivity across grid'
    Bmod_min = huge(1.0_dp)
    Bmod_max = 0.0_dp
    do i_phi = 1, n_phi, max(1, n_phi/4)
        do i_th = 1, n_th, max(1, n_th/4)
            do i_r = 1, n_r, max(1, n_r/4)
                x(1) = xmin(1) + (i_r-1)*(xmax(1)-xmin(1))/(n_r-1)
                x(2) = xmin(2) + (i_th-1)*(xmax(2)-xmin(2))/(n_th-1)
                x(3) = xmin(3) + (i_phi-1)*(xmax(3)-xmin(3))/(n_phi-1)
                ! Swap coordinates: physics [r, th, phi] -> spline [phi, th, r]
                x_spl(1) = x(3)
                x_spl(2) = x(2)
                x_spl(3) = x(1)
                call evaluate_batch_splines_3d(spl_field_batch, x_spl, y_batch)
                ! y_batch(5) is Bmod
                Bmod_min = min(Bmod_min, y_batch(5))
                Bmod_max = max(Bmod_max, y_batch(5))
            end do
        end do
    end do

    if (Bmod_min <= 0.0_dp) then
        print *, '  FAILED: Bmod must be positive everywhere, min=', Bmod_min
        n_failed = n_failed + 1
    end if

    if (Bmod_max > 200000.0_dp) then
        print *, '  FAILED: Bmod unreasonably large (CGS), max=', Bmod_max
        n_failed = n_failed + 1
    end if

    ! Test 6: xmin < xmax for all coordinates
    print *, 'Test 6: Coordinate bounds ordering'
    if (any(xmin >= xmax)) then
        print *, '  FAILED: xmin must be less than xmax'
        print *, '    xmin=', xmin
        print *, '    xmax=', xmax
        n_failed = n_failed + 1
    end if

    if (n_failed == 0) then
        print *, '================================'
        print *, 'All Albert diagnostic tests PASSED'
        print *, '================================'
        stop 0
    else
        print *, '================================'
        print *, n_failed, ' tests FAILED'
        print *, '================================'
        error stop 1
    end if

end program test_field_can_albert_diagnostic
