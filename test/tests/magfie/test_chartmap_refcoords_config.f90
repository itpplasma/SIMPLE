program test_chartmap_refcoords_config
    !> Regression test for issue #206:
    !> Ensure chartmap can be used as reference coordinates (coord_input)
    !> while still using a VMEC equilibrium for the field (field_input).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: field_input, coord_input
    use velo_mod, only: isw_field_type
    use magfie_sub, only: REFCOORDS, init_magfie, magfie

    implicit none

    type(tracer_t) :: tracer
    logical :: has_wout, has_chartmap
    real(dp) :: u(3)
    real(dp) :: bmod, sqrtg
    real(dp) :: bder(3), hcov(3), hctr(3), hcurl(3)

    inquire (file='wout.nc', exist=has_wout)
    inquire (file='wout.chartmap.nc', exist=has_chartmap)
    if (.not. has_wout) then
        print *, 'FAIL: wout.nc not found in test directory'
        error stop 1
    end if
    if (.not. has_chartmap) then
        print *, 'FAIL: wout.chartmap.nc not found in test directory'
        error stop 1
    end if

    field_input = 'wout.nc'
    coord_input = 'wout.chartmap.nc'
    isw_field_type = REFCOORDS

    call init_field(tracer, coord_input, 5, 5, 5, 0)
    call init_magfie(REFCOORDS)

    u = [0.5_dp, 0.0_dp, 0.0_dp]
    call magfie(u, bmod, sqrtg, bder, hcov, hctr, hcurl)

    if (bmod <= 0.0_dp) then
        print *, 'FAIL: bmod must be positive, got ', bmod
        error stop 1
    end if
    if (sqrtg /= sqrtg) then
        print *, 'FAIL: sqrtg is NaN'
        error stop 1
    end if
    if (sqrtg == 0.0_dp) then
        print *, 'FAIL: sqrtg must be non-zero, got ', sqrtg
        error stop 1
    end if

    print *, 'PASS: chartmap coord_input + VMEC field_input initializes and evaluates'
end program test_chartmap_refcoords_config
