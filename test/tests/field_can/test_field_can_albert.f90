program test_field_can_albert

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use util, only: twopi

    use simple, only: tracer_t
    use simple_main, only: init_field
    use magfie_sub, only: ALBERT
    use velo_mod, only: isw_field_type
    use field, only: vmec_field_t, create_vmec_field
    use field_can_meiss, only: n_r, n_th, n_phi, transformation_relerr, xmax, &
        get_meiss_coordinates, ref_to_integ_meiss
    use field_can_mod, only: eval_field => evaluate, field_can_t, field_can_from_name
    use field_can_albert, only: ref_to_integ_albert, integ_to_ref_albert, &
        r_of_xc, psi_of_x
    use params, only: canonical_grid_nr, canonical_grid_ntheta, &
        canonical_grid_nphi, canonical_ode_relerr, field_input, coord_input

    implicit none

    type(tracer_t) :: norb
    type(vmec_field_t) :: magfie

    isw_field_type = ALBERT
    call create_vmec_field(magfie)
    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    canonical_grid_nr = 8
    canonical_grid_ntheta = 9
    canonical_grid_nphi = 10
    canonical_ode_relerr = 1.0e-8_dp

    print *, 'init_field'
    call init_field(norb, 'wout.nc', 5, 5, 3, 0)
    if (any([n_r, n_th, n_phi] /= [8, 9, 10]) .or. &
        transformation_relerr /= 1.0e-8_dp) then
        error stop 'configured Albert map controls were not applied'
    end if
    call assert_albert_map

    call assert_albert_to_meiss_switch

    canonical_grid_nr = 7
    canonical_grid_ntheta = 8
    canonical_grid_nphi = 9
    canonical_ode_relerr = 5.0e-9_dp
    call init_field(norb, 'wout.nc', 5, 5, 3, 0)
    if (any([n_r, n_th, n_phi] /= [7, 8, 9]) .or. &
        transformation_relerr /= 5.0e-9_dp) then
        error stop 'reinitialized Albert map controls were not applied'
    end if
    call assert_albert_map

contains

    subroutine assert_albert_map
        real(dp), parameter :: tolerance = 5.0e-2_dp
        real(dp) :: xref(3), xref_back(3), xinteg(3), xinteg_outer(3)
        real(dp) :: error(3)
        type(field_can_t) :: field_value, periodic_value

        xref = [0.25_dp, 1.1_dp, 0.3_dp]
        call ref_to_integ_albert(xref, xinteg)
        call integ_to_ref_albert(xinteg, xref_back)
        error(1) = abs(xref_back(1) - xref(1))
        error(2) = abs(modulo(xref_back(2) - xref(2) + 0.5_dp*twopi, twopi) - &
            0.5_dp*twopi)
        error(3) = abs(modulo(xref_back(3) - xref(3) + 0.5_dp*twopi, twopi) - &
            0.5_dp*twopi)
        if (maxval(error) > tolerance) error stop 'Albert coordinate roundtrip failed'

        call ref_to_integ_albert([0.36_dp, xref(2), xref(3)], xinteg_outer)
        if (xinteg_outer(1) <= xinteg(1)) error stop 'Albert radial orientation failed'

        call eval_field(field_value, xinteg(1), xinteg(2), xinteg(3), 1)
        call eval_field(periodic_value, xinteg(1), xinteg(2), &
            xinteg(3) + xmax(3), 1)
        if (.not. ieee_is_finite(field_value%Bmod) .or. field_value%Bmod <= 0.0_dp .or. &
            .not. all(ieee_is_finite(field_value%dAth)) .or. &
            .not. all(ieee_is_finite(field_value%dAph)) .or. &
            .not. ieee_is_finite(field_value%hth) .or. &
            .not. ieee_is_finite(field_value%hph)) then
            error stop 'Albert field evaluation is not finite'
        end if
        if (abs(periodic_value%Bmod - field_value%Bmod) > tolerance) then
            error stop 'Albert field-period seam failed'
        end if
    end subroutine assert_albert_map

    subroutine assert_albert_to_meiss_switch
        real(dp) :: xinteg(3)
        type(field_can_t) :: field_value

        call field_can_from_name('meiss', magfie, 7, 8, 9, 5.0e-9_dp)
        if (allocated(r_of_xc) .or. allocated(psi_of_x)) then
            error stop 'Albert resources survived a Meiss mode switch'
        end if
        call get_meiss_coordinates
        call ref_to_integ_meiss([0.25_dp, 1.1_dp, 0.3_dp], xinteg)
        call eval_field(field_value, xinteg(1), xinteg(2), xinteg(3), 1)
        if (.not. ieee_is_finite(field_value%Bmod) .or. field_value%Bmod <= 0.0_dp) then
            error stop 'Meiss evaluation failed after Albert mode switch'
        end if
    end subroutine assert_albert_to_meiss_switch

end program test_field_can_albert
