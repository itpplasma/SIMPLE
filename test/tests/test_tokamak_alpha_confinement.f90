program test_tokamak_alpha_confinement
    !> Verify analytical Grad-Shafranov field integration via geoflux and
    !> canonical (Meiss) coordinates remain consistent for symplectic use.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use params, only: read_config, isw_field_type
    use field, only: field_from_file, MagneticField
    use field_geoflux, only: GeofluxField, geoflux_ready
    use field_can_mod, only: init_field_can, can_to_ref, ref_to_can, evaluate, &
        FieldCan_init
    use field_can_base, only: FieldCan
    use magfie_sub, only: MEISS
    implicit none

    character(len=256) :: config_file
    class(MagneticField), allocatable :: field_obj
    real(dp), dimension(3) :: x_ref, x_can, x_round
    real(dp), dimension(3) :: x_eval, acov, hcov, sqg
    real(dp) :: bmod_geo, rel_err, phi_err, twopi
    type(FieldCan) :: f_can
    real(dp), parameter :: tol = 1.0e-8_dp

    twopi = 2.0_dp*acos(-1.0_dp)

    config_file = '../../../examples/tokamak_alpha_confinement/simple.in'
    call read_config(config_file)

    if (isw_field_type /= MEISS) then
        error stop 'Tokamak config must enable Meiss canonical field'
    end if

    call field_from_file('analytical', field_obj, .true.)

    if (.not. geoflux_ready) then
        error stop 'Analytical geoflux field not initialized'
    end if

    select type(field_obj)
    type is (GeofluxField)
        ! Expected path: geoflux coordinates populated from analytical GS
    class default
        error stop 'field_from_file did not return GeofluxField'
    end select

    call init_field_can(MEISS, field_obj)
    call FieldCan_init(f_can)

    x_ref = [0.3_dp, 0.1_dp, 0.25_dp]

    call ref_to_can(x_ref, x_can)
    call can_to_ref(x_can, x_round)

    if (abs(x_round(1) - x_ref(1)) > tol) then
        error stop 'ref/can radial mismatch exceeds tolerance'
    end if

    if (abs(x_round(2) - x_ref(2)) > 1.0e-7_dp) then
        error stop 'ref/can poloidal mismatch exceeds tolerance'
    end if

    phi_err = modulo(abs(x_round(3) - x_ref(3)), twopi)
    phi_err = min(phi_err, twopi - phi_err)
    if (phi_err > 1.0e-6_dp) then
        error stop 'ref/can toroidal mismatch exceeds tolerance'
    end if

    x_eval = [sqrt(x_ref(1)), x_ref(2), x_ref(3)]
    call field_obj%evaluate(x_eval, acov, hcov, bmod_geo, sqg)

    if (bmod_geo <= 0.0_dp) then
        error stop 'Geoflux evaluation returned non-positive |B|'
    end if

    if (sqg(1) <= 0.0_dp) then
        error stop 'Geoflux Jacobian must be positive'
    end if

    call evaluate(f_can, x_can(1), x_can(2), x_can(3), 0)

    if (f_can%Bmod <= 0.0_dp) then
        error stop 'Canonical evaluation returned non-positive |B|'
    end if

    rel_err = abs(f_can%Bmod - bmod_geo) / bmod_geo

    if (rel_err > 5.0e-3_dp) then
        error stop 'Canonical and reference fields differ beyond tolerance'
    end if

    print *, 'Tokamak analytical field and canonical transforms consistent.'

end program test_tokamak_alpha_confinement
