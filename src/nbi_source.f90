module nbi_source
    !> Neutral beam deposition along pencil beamlets and guiding-centre births.
    !> Attenuation uses the Suzuki et al. (1998) stopping cross section. A
    !> beamlet of power P with power fractions f_j at energies E/j injects
    !> R_j = f_j P/(E/j) neutrals per second, of which
    !> dR = R_j [exp(-tau(l)) - exp(-tau(l+dl))] ionize on [l, l+dl], with
    !> d tau/dl = n_e sigma (T. H. Stix, Plasma Phys. 14 (1972) 367, Sec. II).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use finite_check, only: is_finite
    use beam_stopping, only: suzuki_stopping, suzuki_impurity_supported
    implicit none
    private
    public :: nbi_beamlet_t, nbi_plasma_t, nbi_segment_t, read_nbi_beamlets, &
              beamlet_deposition, guiding_centre_shift, nbi_plasma_lookup, &
              ev_joule, amu_kg

    real(dp), parameter :: ev_joule = 1.602176634e-19_dp
    real(dp), parameter :: amu_kg = 1.66053906660e-27_dp

    type :: nbi_beamlet_t
        character(len=32) :: name = ''
        integer :: mass = 2                 ! beam isotope mass number
        real(dp) :: energy_ev = 0.0_dp      ! full-energy component [eV]
        real(dp) :: power_w = 0.0_dp        ! neutral power [W]
        real(dp) :: fraction(3) = [1.0_dp, 0.0_dp, 0.0_dp] ! power at E, E/2, E/3
        real(dp) :: origin(3) = 0.0_dp      ! Cartesian start point [m]
        real(dp) :: direction(3) = 0.0_dp   ! unit vector
    end type nbi_beamlet_t

    !> Local plasma on the beam line: inside the plasma, n_e [m^-3], T_e [eV],
    !> Z_eff and the impurity charge used by the stopping cross section.
    type :: nbi_plasma_t
        logical :: inside = .false.
        real(dp) :: ne = 0.0_dp, te = 0.0_dp, zeff = 1.0_dp
        integer :: z_impurity = 0
    end type nbi_plasma_t

    !> Ionization on one path segment of one energy component.
    type :: nbi_segment_t
        integer :: beamlet = 0, component = 0
        real(dp) :: l0 = 0.0_dp, l1 = 0.0_dp ! path length [m]
        real(dp) :: rate = 0.0_dp            ! ionizations per second
    end type nbi_segment_t

    abstract interface
        function nbi_plasma_lookup(x) result(plasma)
            import :: dp, nbi_plasma_t
            real(dp), intent(in) :: x(3)     ! Cartesian position [m]
            type(nbi_plasma_t) :: plasma
        end function nbi_plasma_lookup
    end interface

contains

    subroutine read_nbi_beamlets(filename, beamlets)
        !> One beamlet per line: name mass energy_eV power_W f_full f_half
        !> f_third x y z dx dy dz, Cartesian in metres; '#' starts a comment.
        character(len=*), intent(in) :: filename
        type(nbi_beamlet_t), allocatable, intent(out) :: beamlets(:)
        type(nbi_beamlet_t) :: b
        type(nbi_beamlet_t), allocatable :: tmp(:)
        character(len=1024) :: line
        integer :: unit, ios, n
        real(dp) :: norm

        allocate (beamlets(0))
        open (newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'nbi source: cannot open beamlet file'
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            read (line, *, iostat=ios) b%name, b%mass, b%energy_ev, b%power_w, &
                b%fraction, b%origin, b%direction
            if (ios /= 0) error stop 'nbi source: malformed beamlet line'
            if (.not. all(is_finite([b%energy_ev, b%power_w, b%fraction, b%origin, &
                                     b%direction]))) &
                error stop 'nbi source: nonfinite beamlet entry'
    if (b%mass < 1 .or. b%mass > 3) error stop 'nbi source: beam mass must be 1, 2 or 3'
      if (b%energy_ev <= 0.0_dp .or. b%power_w < 0.0_dp .or. any(b%fraction < 0.0_dp)) &
                error stop 'nbi source: nonpositive energy or negative power'
            if (abs(sum(b%fraction) - 1.0_dp) > 1.0e-6_dp) &
                error stop 'nbi source: power fractions must sum to one'
            norm = norm2(b%direction)
            if (norm <= 0.0_dp) error stop 'nbi source: zero beam direction'
            b%direction = b%direction/norm
            n = size(beamlets)
            allocate (tmp(n + 1))
            tmp(1:n) = beamlets
            tmp(n + 1) = b
            call move_alloc(tmp, beamlets)
        end do
        close (unit)
        if (size(beamlets) == 0) error stop 'nbi source: no beamlets in file'
    end subroutine read_nbi_beamlets

    subroutine beamlet_deposition(beamlet, component, plasma_at, dl, lmax, segments, &
                                  injected, shine)
        !> Ionization rates of one energy component along the beamlet path,
        !> integrating the optical depth with the midpoint rule on steps dl [m].
        type(nbi_beamlet_t), intent(in) :: beamlet
        integer, intent(in) :: component
        procedure(nbi_plasma_lookup) :: plasma_at
        real(dp), intent(in) :: dl, lmax
        type(nbi_segment_t), allocatable, intent(out) :: segments(:)
        real(dp), intent(out) :: injected, shine
        type(nbi_plasma_t) :: p
        type(nbi_segment_t), allocatable :: tmp(:)
        real(dp) :: energy, tau, dtau, survive0, survive1, l
        integer :: i, nstep, n

    if (component < 1 .or. component > 3) error stop 'nbi source: component must be 1-3'
        if (dl <= 0.0_dp .or. lmax <= dl) error stop 'nbi source: invalid path step'
        energy = beamlet%energy_ev/component
        injected = beamlet%fraction(component)*beamlet%power_w/(energy*ev_joule)
        allocate (segments(0))
        shine = injected
        if (injected <= 0.0_dp) return
        nstep = ceiling(lmax/dl)
        tau = 0.0_dp
        do i = 1, nstep
            l = (i - 0.5_dp)*dl
            p = plasma_at(beamlet%origin + l*beamlet%direction)
            if (.not. p%inside) cycle
            dtau = p%ne*suzuki_stopping(energy*1.0e-3_dp/beamlet%mass, p%ne, p%te, &
                                        beamlet%mass, p%zeff, p%z_impurity)*dl
            survive0 = exp(-tau)
            tau = tau + dtau
            survive1 = exp(-tau)
            n = size(segments)
            allocate (tmp(n + 1))
            tmp(1:n) = segments
            tmp(n + 1) = nbi_segment_t(0, component, l - 0.5_dp*dl, l + 0.5_dp*dl, &
                                       injected*(survive0 - survive1))
            call move_alloc(tmp, segments)
        end do
        shine = injected*exp(-tau)
    end subroutine beamlet_deposition

    pure function guiding_centre_shift(x, v, bhat, omega) result(xgc)
        !> Guiding centre of a gyrating ion at x with velocity v in a field
        !> along bhat, with signed gyrofrequency omega = qB/m:
        !> X = x + v x bhat/omega.
        real(dp), intent(in) :: x(3), v(3), bhat(3), omega
        real(dp) :: xgc(3)

        xgc = x + [v(2)*bhat(3) - v(3)*bhat(2), v(3)*bhat(1) - v(1)*bhat(3), &
                   v(1)*bhat(2) - v(2)*bhat(1)]/omega
    end function guiding_centre_shift

end module nbi_source
