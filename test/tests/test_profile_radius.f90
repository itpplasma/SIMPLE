program test_profile_radius
    !> Regression for #571: the collision operator must see the plasma
    !> profiles at the normalised toroidal flux s of the marker, whatever
    !> radius integ_coords integrates in. Each case places markers at known s,
    !> maps them to the tracing state and compares the collision coefficients
    !> at profile_s(z) with those at s.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use velo_mod, only: isw_field_type
    use magfie_sub, only: CANFLUX, VMEC, BOOZER, MEISS, ALBERT, REFCOORDS, &
                          SPECTRE, set_magfie_spectre_field
    use field_can_mod, only: ref_to_integ, init_field_can
    use params, only: field_input, coord_input, canonical_grid_nr, &
                      canonical_grid_ntheta, canonical_grid_nphi, canonical_ode_relerr
    use profile_radius, only: init_profile_radius, set_profile_radius, &
                              profile_s, RADIUS_RHO_TOR
    use simple_profiles, only: profile_type, active_profile, TWO_POWER
    use simple_profiles, only: Te_scale, Ti1_scale, Ti2_scale, ni1_scale, ni2_scale
    use simple_profiles, only: Te_p1, Te_p2, Ti1_p1, Ti1_p2, Ti2_p1, Ti2_p2
    use simple_profiles, only: ni1_p1, ni1_p2, ni2_p1, ni2_p2
    use collis_alp, only: init_collision_profiles, get_local_coeffs, nsorts

    implicit none

    real(dp), parameter :: s_known(3) = [0.09_dp, 0.25_dp, 0.64_dp]
    real(dp), parameter :: theta = 0.7_dp, phi = 0.3_dp
    character(len=64) :: mode
    character(len=1024) :: file
    type(tracer_t) :: norb

    call get_command_argument(1, mode)
    call get_command_argument(2, file)
    call steep_collision_profiles

    canonical_grid_nr = 16
    canonical_grid_ntheta = 17
    canonical_grid_nphi = 18
    canonical_ode_relerr = 1.0e-10_dp

    select case (trim(mode))
    case ('rho_lcfs')
        ! Extended VMEC chartmaps scale rho = rho_lcfs*sqrt(s) inside the LCFS.
        call set_profile_radius(RADIUS_RHO_TOR, .false., 0.75_dp)
        call check_states('rho_lcfs', 0.75_dp*sqrt(s_known), .false., 1.0e-14_dp)
    case ('vmec')
        call init_vmec_case(VMEC, 'wout.nc')
        call check_states('vmec', s_known, .false., 1.0e-14_dp)
    case ('canflux')
        call init_vmec_case(CANFLUX, 'wout.nc')
        call check_states('canflux', s_known, .true., 1.0e-14_dp)
    case ('boozer')
        call init_vmec_case(BOOZER, 'wout.nc')
        call check_states('boozer', s_known, .true., 1.0e-14_dp)
    case ('meiss')
        call init_vmec_case(MEISS, 'wout.nc')
        call check_states('meiss', s_known, .true., 1.0e-10_dp, legacy_differs=.true.)
    case ('albert')
        call init_vmec_case(ALBERT, 'wout.nc')
        ! psi/Ath_norm is not s (#571 item 4); the Albert inverse on this
        ! coarse grid round-trips s to ~1e-3 (test_field_can_albert: 5e-2).
        call check_states('albert', s_known, .true., 2.0e-3_dp, legacy_differs=.true.)
    case ('refcoords')
        call init_vmec_case(REFCOORDS, trim(file))
        call check_states('refcoords', sqrt(s_known), .false., 1.0e-14_dp, &
                          legacy_differs=.true.)
    case ('meiss_chartmap')
        call init_meiss_chartmap(trim(file))
        call check_states('meiss_chartmap', sqrt(s_known), .true., 1.0e-10_dp, &
                          legacy_differs=.true.)
    case ('spectre')
        call check_spectre(trim(file))
    case default
        error stop 'usage: test_profile_radius.x <case> [file]'
    end select
    print *, 'PASS ', trim(mode)

contains

    subroutine steep_collision_profiles
        real(dp) :: v0

        profile_type = 'two_power'
        active_profile = TWO_POWER
        Te_scale = 2.0e4_dp; Te_p1 = 1.0_dp; Te_p2 = 2.0_dp
        Ti1_scale = 2.0e4_dp; Ti1_p1 = 1.0_dp; Ti1_p2 = 2.0_dp
        Ti2_scale = 2.0e4_dp; Ti2_p1 = 1.0_dp; Ti2_p2 = 2.0_dp
        ni1_scale = 1.0e20_dp; ni1_p1 = 1.0_dp; ni1_p2 = 3.0_dp
        ni2_scale = 1.0e20_dp; ni2_p1 = 1.0_dp; ni2_p2 = 3.0_dp
        v0 = 1.0_dp
        call init_collision_profiles(2.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 3.5e6_dp, v0)
    end subroutine steep_collision_profiles

    subroutine init_vmec_case(field_type, coord_file)
        integer, intent(in) :: field_type
        character(*), intent(in) :: coord_file

        isw_field_type = field_type
        field_input = 'wout.nc'
        coord_input = coord_file
        call init_field(norb, 'wout.nc', 5, 5, 3, 0)
        call init_profile_radius
    end subroutine init_vmec_case

    subroutine init_meiss_chartmap(chartmap_file)
        !> MEISS on a field splined onto a rho_tor chartmap: the coils with
        !> REFCOORDS chartmap path of the WP-STEL runs.
        use field_vmec, only: vmec_field_t, create_vmec_field
        use field_splined, only: splined_field_t, create_splined_field
        use reference_coordinates, only: ref_coords
        character(*), intent(in) :: chartmap_file
        type(vmec_field_t) :: vmec_field
        type(splined_field_t) :: splined

        call init_vmec_case(REFCOORDS, chartmap_file)
        call create_vmec_field(vmec_field)
        call create_splined_field(vmec_field, ref_coords, splined)
        isw_field_type = MEISS
        call init_field_can(MEISS, splined, canonical_grid_nr, &
                            canonical_grid_ntheta, canonical_grid_nphi, &
                            canonical_ode_relerr)
        call init_profile_radius
    end subroutine init_meiss_chartmap

    subroutine check_states(name, r_ref, to_integ, s_tol, legacy_differs)
        character(*), intent(in) :: name
        real(dp), intent(in) :: r_ref(:)
        logical, intent(in) :: to_integ
        real(dp), intent(in) :: s_tol
        logical, intent(in), optional :: legacy_differs
        real(dp) :: xref(3), z(3)
        integer :: i

        do i = 1, size(s_known)
            xref = [r_ref(i), theta, phi]
            if (to_integ) then
                call ref_to_integ(xref, z)
            else
                z = xref
            end if
            call compare(name, s_known(i), z, s_tol)
            if (present(legacy_differs)) then
                if (legacy_differs .and. abs(z(1) - s_known(i)) < 1.0e-3_dp) then
                    print *, name, ': z(1) = ', z(1), ' s = ', s_known(i)
                    error stop 'case does not discriminate z(1) from s'
                end if
            end if
        end do
    end subroutine check_states

    subroutine compare(name, s, z, s_tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: s, z(3), s_tol
        real(dp) :: s_seen
        real(dp), dimension(nsorts) :: ef, vr, er, ef_seen, vr_seen, er_seen

        s_seen = profile_s(z)
        call get_local_coeffs(s, ef, vr, er)
        call get_local_coeffs(s_seen, ef_seen, vr_seen, er_seen)
        print '(a,1x,a,3es24.16)', 'case', name, s, z(1), s_seen
        if (abs(s_seen - s) > s_tol) then
            print *, name, ': profile radius ', s_seen, ' expected ', s
            error stop 'collision profile radius differs from s'
        end if
        if (any(abs(ef_seen - ef) > 50.0_dp*s_tol*abs(ef)) .or. &
            any(abs(vr_seen - vr) > 50.0_dp*s_tol*abs(vr)) .or. &
            any(abs(er_seen - er) > 50.0_dp*s_tol*abs(er))) then
            print *, name, ': coefficients at s = ', s, ef, vr, er
            print *, name, ': coefficients seen   ', ef_seen, vr_seen, er_seen
            error stop 'collision coefficients differ from the profile at s'
        end if
    end subroutine compare

    subroutine check_spectre(h5file)
        !> Stacked chart rho_g: interfaces carry the cumulative tflux exactly;
        !> inside the axis volume libneo's toroidal flux label is the oracle.
        use field_spectre, only: spectre_field_t, create_spectre_field
        character(*), intent(in) :: h5file
        type(spectre_field_t) :: field
        real(dp) :: label, tflux(3)
        integer :: ierr, k

        call create_spectre_field(field, h5file, ierr)
        if (ierr /= 0) error stop 'cannot load SPECTRE fixture'
        if (field%data%Mvol /= 3) error stop 'SPECTRE fixture must have 3 volumes'
        call set_magfie_spectre_field(field)
        isw_field_type = SPECTRE
        call init_profile_radius

        tflux = field%data%tflux/field%data%tflux(3)
        call compare('spectre', 0.0_dp, [0.0_dp, theta, phi], 1.0e-12_dp)
        do k = 1, 3
            call compare('spectre', tflux(k), [real(k, dp), theta, phi], 1.0e-12_dp)
        end do
        call field%axis_toroidal_flux_label(0.5_dp, label, ierr)
        if (ierr /= 0) error stop 'SPECTRE axis flux label failed'
        call compare('spectre', label, [0.5_dp, theta, phi], 1.0e-12_dp)
    end subroutine check_spectre

end program test_profile_radius
