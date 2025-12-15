program test_magfie_coils
    !> Test coils field evaluation comparing VMEC field, raw Biot-Savart,
    !> and splined field with real assertions.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use magfie_sub, only: VMEC
    use velo_mod, only: isw_field_type
    use field_base, only: magnetic_field_t
    use field_vmec, only: vmec_field_t
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use magfie_sub, only: magfie_vmec
    use util, only: twopi
    use cylindrical_cartesian, only: cyl_to_cart

    implicit none

    type(vmec_field_t) :: vmec_field
    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: splined_coils
    real(dp) :: dummy, x(3), Acov(3), hcov(3), Bmod
    real(dp) :: Acov_raw(3), hcov_raw(3), Bmod_raw
    real(dp) :: Acov_spl(3), hcov_spl(3), Bmod_spl
    integer :: n_failed

    logical :: wout_exists, coils_exists

    n_failed = 0
    isw_field_type = VMEC

    inquire(file='wout.nc', exist=wout_exists)
    inquire(file='coils.simple', exist=coils_exists)
    if (.not. wout_exists) then
        print *, 'FAILED: Required VMEC file (wout.nc) not found'
        error stop 1
    end if
    if (.not. coils_exists) then
        print *, 'FAILED: Required coils file (coils.simple) not found'
        error stop 1
    end if

    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates('wout.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, splined_coils)

    x = [0.5d0, 0.3d0, 0.2d0]

    ! Test 1: VMEC field returns physically valid values
    print *, 'Test 1: VMEC field evaluation'
    call vmec_field%evaluate(x, Acov, hcov, Bmod)

    if (Bmod <= 0.0_dp) then
        print *, '  FAILED: vmec_field Bmod must be positive, got ', Bmod
        n_failed = n_failed + 1
    end if

    if (Bmod > 200000.0_dp) then
        print *, '  FAILED: vmec_field Bmod unreasonably large (CGS) ', Bmod
        n_failed = n_failed + 1
    end if

    ! Test 2: Splined coils field returns positive Bmod
    print *, 'Test 2: Splined coils field evaluation'
    call splined_coils%evaluate(x, Acov_spl, hcov_spl, Bmod_spl)

    if (Bmod_spl <= 0.0_dp) then
        print *, '  FAILED: splined_coils Bmod must be positive, got ', Bmod_spl
        n_failed = n_failed + 1
    end if

    ! Test 3: Raw Biot-Savart at reference point
    print *, 'Test 3: Raw Biot-Savart evaluation'
    call evaluate_raw_at_ref(raw_coils, x, Acov_raw, hcov_raw, Bmod_raw)

    if (Bmod_raw <= 0.0_dp) then
        print *, '  FAILED: raw_coils Bmod must be positive, got ', Bmod_raw
        n_failed = n_failed + 1
    end if

    ! Test 4: Splined and raw should agree within 1%
    print *, 'Test 4: Splined vs raw Biot-Savart agreement'
    if (abs(Bmod_spl - Bmod_raw) / Bmod_raw > 0.01_dp) then
        print *, '  FAILED: Bmod_spl and Bmod_raw differ by more than 1%'
        print *, '    Bmod_spl = ', Bmod_spl, ' Bmod_raw = ', Bmod_raw
        n_failed = n_failed + 1
    end if

    ! Test 5: hcov should be normalized (|h| ~ 1 since h = B/|B|)
    print *, 'Test 5: hcov normalization'
    if (abs(sqrt(sum(hcov**2)) - 1.0_dp) > 0.1_dp) then
        print *, '  FAILED: vmec_field hcov not normalized, |h| = ', sqrt(sum(hcov**2))
        n_failed = n_failed + 1
    end if

    ! Test 6: magfie_vmec consistency
    print *, 'Test 6: magfie_vmec evaluation'
    call test_magfie(n_failed)

    if (n_failed == 0) then
        print *, '================================'
        print *, 'All magfie_coils tests PASSED'
        print *, '================================'
        stop 0
    else
        print *, '================================'
        print *, n_failed, ' tests FAILED'
        print *, '================================'
        error stop 1
    end if

contains

    subroutine evaluate_raw_at_ref(field, x_spline, Acov, hcov, Bmod)
        !> Helper to evaluate raw coils at a spline coordinate point.
        !> x_spline is in (r, theta, phi) where r = sqrt(s).
        !> ref_coords expects VMEC coords (s, theta, phi).
        type(coils_field_t), intent(in) :: field
        real(dp), intent(in) :: x_spline(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod

        real(dp) :: x_vmec(3), x_cyl(3), x_cart(3)

        ! Convert spline coords (r, theta, phi) to VMEC coords (s, theta, phi)
        x_vmec = [x_spline(1)**2, x_spline(2), x_spline(3)]
        call ref_coords%evaluate_cyl(x_vmec, x_cyl)
        call cyl_to_cart(x_cyl, x_cart)
        call field%evaluate(x_cart, Acov, hcov, Bmod)
    end subroutine evaluate_raw_at_ref


    subroutine test_magfie(n_failed)
        integer, intent(inout) :: n_failed
        real(dp) :: bmod, sqrtg
        real(dp), dimension(3) :: bder, hcovar, hctrvr, hcurl

        call magfie_vmec(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

        ! Bmod must be positive
        if (bmod <= 0.0_dp) then
            print *, '  FAILED: magfie_vmec bmod must be positive, got ', bmod
            n_failed = n_failed + 1
        end if

        ! sqrtg (Jacobian) must be positive for proper orientation
        if (sqrtg <= 0.0_dp) then
            print *, '  FAILED: magfie_vmec sqrtg must be positive, got ', sqrtg
            n_failed = n_failed + 1
        end if
    end subroutine test_magfie

end program test_magfie_coils
