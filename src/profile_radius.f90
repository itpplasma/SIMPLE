module profile_radius
    !> Map the tracing state to the radius of the plasma profiles.
    !>
    !> &profiles are functions of the normalised toroidal flux s. The radial
    !> state component z(1) is s only for integ_coords TEST, CANFLUX, VMEC and
    !> BOOZER (VMEC input or Boozer chartmap). MEISS and ALBERT integrate in
    !> their own radius and map back to the reference coordinates of their
    !> source field; REFCOORDS integrates in the reference coordinates; SPECTRE
    !> integrates in the stacked volume chart rho_g in [0, Mvol].
    !>
    !> With &config profile_label_file, REFCOORDS takes the profile radius from
    !> a tabulated label at the physical position of the marker instead
    !> (module profile_label), e.g. for geometric chartmaps without flux
    !> meaning.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: init_profile_radius, set_profile_radius, profile_s
    public :: RADIUS_S, RADIUS_RHO_TOR, RADIUS_SPECTRE, RADIUS_NONE, RADIUS_LABEL

    !> Reference radius is s.
    integer, parameter :: RADIUS_S = 1
    !> Reference radius is rho = rho_scale*sqrt(s).
    integer, parameter :: RADIUS_RHO_TOR = 2
    !> SPECTRE stacked chart rho_g: enclosed toroidal flux / tflux(Mvol).
    integer, parameter :: RADIUS_SPECTRE = 3
    !> No toroidal-flux meaning (geometric or poloidal radius).
    integer, parameter :: RADIUS_NONE = 4
    !> Tabulated profile label at the physical position (profile_label_file).
    integer, parameter :: RADIUS_LABEL = 5

    abstract interface
        function radius_i(z) result(r)
            import :: dp
            real(dp), intent(in) :: z(3)
            real(dp) :: r
        end function radius_i
    end interface

    !> Set by set_profile_radius; called from profile_s only when needed.
    procedure(radius_i), pointer, save :: reference_radius => null()
    procedure(radius_i), pointer, save :: spectre_s_of_rho => null()
    procedure(radius_i), pointer, save :: label_s => null()

    integer, save :: radius_kind = RADIUS_S
    logical, save :: map_to_reference = .false.
    real(dp), save :: rho_scale = 1.0_dp

contains

    subroutine set_profile_radius(kind, use_integ_to_ref, rho_lcfs)
        use field_can_meiss, only: integ_to_ref_radius_meiss
        integer, intent(in) :: kind
        logical, intent(in) :: use_integ_to_ref
        real(dp), intent(in), optional :: rho_lcfs

        if (kind < RADIUS_S .or. kind > RADIUS_LABEL) &
            error stop 'set_profile_radius: unknown radius kind'
        radius_kind = kind
        map_to_reference = use_integ_to_ref
        reference_radius => null()
        spectre_s_of_rho => null()
        label_s => null()
        if (use_integ_to_ref) reference_radius => integ_to_ref_radius_meiss
        if (kind == RADIUS_SPECTRE) spectre_s_of_rho => spectre_s
        if (kind == RADIUS_LABEL) then
            if (use_integ_to_ref) &
                error stop 'set_profile_radius: a profile label needs '// &
                'reference coordinates'
            label_s => label_s_refcoords
        end if
        rho_scale = 1.0_dp
        if (present(rho_lcfs)) then
            if (.not. (rho_lcfs > 0.0_dp .and. rho_lcfs <= 1.0_dp)) &
                error stop 'set_profile_radius: rho_lcfs must be in (0, 1]'
            rho_scale = rho_lcfs
        end if
    end subroutine set_profile_radius

    !> Select the mapping for the active integ_coords and field. Call after
    !> the field and canonical coordinates are initialised.
    subroutine init_profile_radius
        use velo_mod, only: isw_field_type
        use field_can_albert, only: integ_to_ref_radius_albert
        use magfie_sub, only: TEST, CANFLUX, VMEC, BOOZER, MEISS, ALBERT, &
                              REFCOORDS, SPECTRE, spectre_field
        use field_can_meiss, only: field_noncan
        use field_vmec, only: vmec_field_t
        use field_splined, only: splined_field_t
        use reference_coordinates, only: ref_coords
        use params, only: coord_input, profile_label_file
        use profile_label, only: load_profile_label

        if (len_trim(profile_label_file) > 0) then
            ! The label replaces the chart radius; it needs the physical position,
            ! which REFCOORDS states give through the reference coordinates.
            if (isw_field_type /= REFCOORDS) &
                error stop 'init_profile_radius: profile_label_file needs '// &
                'integ_coords = REFCOORDS'
            if (.not. allocated(ref_coords)) &
                error stop 'init_profile_radius: no reference coordinates'
            call load_profile_label(profile_label_file)
            call set_profile_radius(RADIUS_LABEL, .false.)
            return
        end if

        select case (isw_field_type)
        case (TEST, CANFLUX, VMEC, BOOZER)
            ! Boozer chartmaps integrate in s = rho^2 (ref_to_integ_boozer_chartmap).
            call set_profile_radius(RADIUS_S, .false.)
        case (REFCOORDS)
            if (.not. allocated(ref_coords)) &
                error stop 'init_profile_radius: no reference coordinates'
            call from_coordinates(ref_coords, .false.)
        case (MEISS, ALBERT)
            if (.not. allocated(field_noncan)) &
                error stop 'init_profile_radius: canonical source field not set'
            select type (f => field_noncan)
            type is (vmec_field_t)
                call set_profile_radius(RADIUS_S, .true.)
            type is (splined_field_t)
                call from_coordinates(f%coords, .true.)
            class default
                call set_profile_radius(RADIUS_NONE, .true.)
            end select
        case (SPECTRE)
            if (.not. allocated(spectre_field)) &
                error stop 'init_profile_radius: SPECTRE field not set'
            call set_profile_radius(RADIUS_SPECTRE, .false.)
        case default
            error stop 'init_profile_radius: unknown integ_coords'
        end select
        if (map_to_reference .and. isw_field_type == ALBERT) &
            reference_radius => integ_to_ref_radius_albert

        if (radius_kind == RADIUS_NONE .and. .not. profiles_uniform()) then
            print *, 'init_profile_radius: the radial coordinate of integ_coords =', &
                isw_field_type, ' with coord_input = ', trim(coord_input), &
                ' has no toroidal-flux meaning.'
            error stop 'radially varying &profiles need a toroidal-flux radius; '// &
                'use a rho_tor/psi_tor_norm chartmap, VMEC, a profile_label_file '// &
                'or radially uniform profiles'
        end if

    contains

        subroutine from_coordinates(cs, use_integ_to_ref)
            use libneo_coordinates, only: coordinate_system_t, RHO_TOR, &
                                          PSI_TOR_NORM, chartmap_coordinate_system_t
            use chartmap_metadata, only: chartmap_metadata_t, read_chartmap_metadata
            class(coordinate_system_t), intent(in) :: cs
            logical, intent(in) :: use_integ_to_ref
            type(chartmap_metadata_t) :: meta

            select type (cs)
            class is (chartmap_coordinate_system_t)
                select case (cs%rho_convention)
                case (RHO_TOR)
                    call read_chartmap_metadata(coord_input, meta)
                    ! Extended VMEC chartmaps put the LCFS at rho_lcfs < 1 and
                    ! scale rho = rho_lcfs*sqrt(s) inside it (libneo chartmap.py).
                    select case (trim(meta%rho_convention))
                    case ('vmec_extended', 'vmec_to_wall')
                        call set_profile_radius(RADIUS_RHO_TOR, use_integ_to_ref, &
                                                meta%rho_lcfs)
                    case default
                        call set_profile_radius(RADIUS_RHO_TOR, use_integ_to_ref)
                    end select
                case (PSI_TOR_NORM)
                    call set_profile_radius(RADIUS_S, use_integ_to_ref)
                case default
                    call set_profile_radius(RADIUS_NONE, use_integ_to_ref)
                end select
            class default
                ! VMEC reference coordinates: the radius is s.
                call set_profile_radius(RADIUS_S, use_integ_to_ref)
            end select
        end subroutine from_coordinates

    end subroutine init_profile_radius

    logical function profiles_uniform()
        use simple_profiles, only: get_plasma_params
        integer, parameter :: n = 101
        real(dp) :: p0(5), p(5)
        integer :: i

        call get_plasma_params(0.0_dp, p0(1), p0(2), p0(3), p0(4), p0(5))
        profiles_uniform = .true.
        do i = 2, n
            call get_plasma_params(real(i - 1, dp)/real(n - 1, dp), &
                                   p(1), p(2), p(3), p(4), p(5))
            if (any(p /= p0)) then
                profiles_uniform = .false.
                return
            end if
        end do
    end function profiles_uniform

    !> Normalised toroidal flux s at the tracing state z(1:3).
    !>
    !> This runs at every collision. It must not use or call procedures of
    !> modules that pull in IEEE arithmetic: gfortran then wraps the function in
    !> _gfortran_ieee_procedure_entry/exit (about 0.15 us per call, +11% of a
    !> collisional run). The expensive branches are therefore reached through
    !> the procedure pointers set in set_profile_radius.
    function profile_s(z) result(s)
        real(dp), intent(in) :: z(3)
        real(dp) :: s
        real(dp) :: r

        if (map_to_reference) then
            r = reference_radius(z)
        else
            r = z(1)
        end if
        select case (radius_kind)
        case (RADIUS_RHO_TOR)
            s = (r/rho_scale)**2
        case (RADIUS_SPECTRE)
            ! SPECTRE integrates in rho_g itself, so r = z(1) here.
            s = spectre_s_of_rho(z)
        case (RADIUS_LABEL)
            s = label_s(z)
        case default
            ! RADIUS_NONE is accepted only for radially uniform profiles.
            s = r
        end select
    end function profile_s

    !> Profile label at the physical position of the REFCOORDS state z.
    function label_s_refcoords(z) result(s)
        use reference_coordinates, only: ref_coords
        use profile_label, only: profile_label_at
        real(dp), intent(in) :: z(3)
        real(dp) :: s
        real(dp) :: x(3)

        call ref_coords%evaluate_cyl(z, x)
        s = profile_label_at(x)
    end function label_s_refcoords

    function spectre_s(z) result(s)
        use magfie_sub, only: spectre_field
        use spectre_basis, only: eval_spectre_toroidal_flux_fraction
        real(dp), intent(in) :: z(3)
        real(dp) :: s
        real(dp) :: rho_g, rho, s_local, fraction, derivative, inner
        integer :: lvol, mvol, ierr

        rho_g = z(1)
        mvol = spectre_field%data%Mvol
        rho = min(max(rho_g, 0.0_dp), real(mvol, dp))
        lvol = min(int(rho) + 1, mvol)
        s_local = 2.0_dp*(rho - real(lvol - 1, dp)) - 1.0_dp
        call eval_spectre_toroidal_flux_fraction(spectre_field%data, lvol, s_local, &
                                                 fraction, derivative, ierr)
        if (ierr /= 0) error stop 'profile_s: SPECTRE toroidal flux unavailable'
        ! tflux holds the cumulative toroidal flux enclosed by each interface.
        inner = 0.0_dp
        if (lvol > 1) inner = spectre_field%data%tflux(lvol - 1)
        s = (inner + fraction*(spectre_field%data%tflux(lvol) - inner)) &
            /spectre_field%data%tflux(mvol)
    end function spectre_s

end module profile_radius
