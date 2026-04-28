program test_chartmap_scaling
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_sub, only: load_boozer_from_chartmap, reset_boozer_batch_splines, splint_boozer_coord
    use field_boozer_chartmap, only: boozer_chartmap_field_t, create_boozer_chartmap_field
    use new_vmec_stuff_mod, only: vmec_B_scale, vmec_RZ_scale
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use scaled_chartmap_coordinates, only: scaled_chartmap_coordinate_system_t

    implicit none

    type(boozer_chartmap_field_t), allocatable :: base_field, scaled_field
    real(dp), parameter :: rz_scale = 3.0_dp
    real(dp), parameter :: b_scale = 5.0_dp
    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: u(3), x_cart_base(3), x_cart_scaled(3)
    real(dp) :: ref_cart_base(3), ref_cart_scaled(3)
    real(dp) :: Acov_base(3), hcov_base(3), Bmod_base
    real(dp) :: Acov_scaled(3), hcov_scaled(3), Bmod_scaled
    real(dp) :: A_theta_can, A_phi_can, dA_theta_can, dA_phi_can, d2A_phi_can, d3A_phi_can
    real(dp) :: B_theta_can, dB_theta_can, d2B_theta_can
    real(dp) :: B_phi_can, dB_phi_can, d2B_phi_can
    real(dp) :: Bmod_can, B_r_can, dBmod_can(3), dB_r_can(3), d2Bmod_can(6), d2B_r_can(6)
    real(dp) :: u_back(3), x_scaled_input(3)
    real(dp) :: flux_scale
    character(len=*), parameter :: filename = 'test_boozer_chartmap.nc'
    integer :: ierr, nfail

    nfail = 0
    flux_scale = b_scale*rz_scale*rz_scale
    u = [0.5_dp, 0.7_dp, 0.3_dp]

    vmec_B_scale = 1.0_dp
    vmec_RZ_scale = 1.0_dp
    call create_boozer_chartmap_field(filename, base_field)
    call base_field%coords%evaluate_cart(u, x_cart_base)
    call base_field%evaluate(u, Acov_base, hcov_base, Bmod_base)
    call init_reference_coordinates(filename)
    call ref_coords%evaluate_cart(u, ref_cart_base)

    vmec_B_scale = b_scale
    vmec_RZ_scale = rz_scale
    call create_boozer_chartmap_field(filename, scaled_field)
    call scaled_field%coords%evaluate_cart(u, x_cart_scaled)
    call scaled_field%evaluate(u, Acov_scaled, hcov_scaled, Bmod_scaled)
    call init_reference_coordinates(filename)
    call ref_coords%evaluate_cart(u, ref_cart_scaled)

    if (maxval(abs(x_cart_scaled - rz_scale*x_cart_base)) > tol*max(1.0_dp, maxval(abs(x_cart_scaled)))) then
        print *, 'FAIL: field coords cart scaling mismatch'
        nfail = nfail + 1
    end if

    if (maxval(abs(ref_cart_scaled - rz_scale*ref_cart_base)) > tol*max(1.0_dp, maxval(abs(ref_cart_scaled)))) then
        print *, 'FAIL: reference coords cart scaling mismatch'
        nfail = nfail + 1
    end if

    if (abs(Bmod_scaled - b_scale*Bmod_base) > tol*max(1.0_dp, abs(Bmod_scaled))) then
        print *, 'FAIL: Bmod scaling mismatch'
        nfail = nfail + 1
    end if

    if (abs(scaled_field%torflux - flux_scale*base_field%torflux) > tol*max(1.0_dp, abs(scaled_field%torflux))) then
        print *, 'FAIL: torflux scaling mismatch'
        nfail = nfail + 1
    end if

    if (maxval(abs(Acov_scaled(2:3) - flux_scale*Acov_base(2:3))) > tol*max(1.0_dp, maxval(abs(Acov_scaled(2:3))))) then
        print *, 'FAIL: Acov angular scaling mismatch'
        nfail = nfail + 1
    end if

    if (maxval(abs(hcov_scaled(2:3) - rz_scale*hcov_base(2:3))) > tol*max(1.0_dp, maxval(abs(hcov_scaled(2:3))))) then
        print *, 'FAIL: hcov scaling mismatch'
        nfail = nfail + 1
    end if

    call load_boozer_from_chartmap(filename)
    call splint_boozer_coord(u(1)**2, u(2), u(3), 0, &
                             A_theta_can, A_phi_can, dA_theta_can, dA_phi_can, d2A_phi_can, d3A_phi_can, &
                             B_theta_can, dB_theta_can, d2B_theta_can, &
                             B_phi_can, dB_phi_can, d2B_phi_can, &
                             Bmod_can, dBmod_can, d2Bmod_can, B_r_can, dB_r_can, d2B_r_can)

    if (abs(A_theta_can - Acov_scaled(2)) > tol*max(1.0_dp, abs(A_theta_can))) then
        print *, 'FAIL: canonical A_theta scaling mismatch'
        nfail = nfail + 1
    end if

    if (abs(A_phi_can - Acov_scaled(3)) > tol*max(1.0_dp, abs(A_phi_can))) then
        print *, 'FAIL: canonical A_phi scaling mismatch'
        nfail = nfail + 1
    end if

    if (abs(B_theta_can - hcov_scaled(2)*Bmod_scaled) > tol*max(1.0_dp, abs(B_theta_can))) then
        print *, 'FAIL: canonical B_theta scaling mismatch'
        nfail = nfail + 1
    end if

    if (abs(B_phi_can - hcov_scaled(3)*Bmod_scaled) > tol*max(1.0_dp, abs(B_phi_can))) then
        print *, 'FAIL: canonical B_phi scaling mismatch'
        nfail = nfail + 1
    end if

    if (abs(Bmod_can - Bmod_scaled) > 5.0e-3_dp*max(1.0_dp, abs(Bmod_can))) then
        print *, 'FAIL: canonical Bmod scaling mismatch'
        print *, '  canonical Bmod =', Bmod_can
        print *, '  direct Bmod    =', Bmod_scaled
        nfail = nfail + 1
    end if

    select type (cs => ref_coords)
    type is (scaled_chartmap_coordinate_system_t)
        x_scaled_input = ref_cart_scaled
        call cs%from_cart(x_scaled_input, u_back, ierr)
        if (ierr /= 0) then
            print *, 'FAIL: scaled chartmap from_cart ierr=', ierr
            nfail = nfail + 1
        else if (maxval(abs(u_back - u)) > 1.0e-8_dp) then
            print *, 'FAIL: scaled chartmap from_cart roundtrip mismatch'
            nfail = nfail + 1
        end if
    class default
        print *, 'FAIL: reference coordinates were not wrapped in scaled chartmap type'
        nfail = nfail + 1
    end select

    vmec_B_scale = 1.0_dp
    vmec_RZ_scale = 1.0_dp
    call reset_boozer_batch_splines

    if (nfail /= 0) then
        print *, nfail, 'chartmap scaling tests failed'
        error stop
    end if

    print *, 'PASS: chartmap scaling applies to field and coordinates'
end program test_chartmap_scaling
