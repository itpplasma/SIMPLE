program test_nbi_birth_equilibrium
    !> Neutral beam births in a VMEC equilibrium. Radial injection through the
    !> axis gives near-perpendicular pitch, tangential injection near-parallel
    !> pitch, and every guiding centre lies one gyroradius from its birth point.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_api, only: simple_reset_config, simple_init, simple_sample
    use params, only: field_input, coord_input, startmode, isw_field_type, integmode, &
                      nper, npoiper, ntimstep, deterministic, n_d, n_e, &
                      particle_energy_eV, nbi_file, ntestpart, contr_pp
    use reference_coordinates, only: ref_coords
    use field_can_mod, only: ref_to_integ
    use magfie_sub, only: magfie
    use nbi_birth, only: nbi_marker_position
    use nbi_source, only: ev_joule
    use util, only: p_mass
    use timing, only: init_timer
    use spline_vmec_sub, only: vmec_field
    implicit none
    real(dp), parameter :: e_si = 1.602176634e-19_dp, energy = 1.0e6_dp
    real(dp) :: xgc(3), xi(3), bmod, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)
    real(dp), allocatable :: z(:, :)
    real(dp) :: at, ap, dat, dap, iota, sqg, alam, dl_ds, dl_dt, dl_dp, bctr_t, bctr_p
    real(dp) :: bco_s, bco_t, bco_p, e_cov(3, 3), bcart(3), direction(3), worst_pitch
    real(dp) :: rho, speed, lam, worst, injected, ionized, shine, value
    character(1024) :: filename, mode, key, line
    integer :: i, unit, ios

    call get_command_argument(1, filename)
    call get_command_argument(2, mode)
    call init_timer
    call simple_reset_config
    field_input = trim(filename)
    coord_input = trim(filename)
    isw_field_type = 1
    integmode = 0
    startmode = 8
    deterministic = .true.
    nper = 1
    npoiper = 4
    ntimstep = 2
    n_d = 2
    n_e = 1
    particle_energy_eV = energy
    ntestpart = 400
    contr_pp = -1.0e10_dp
    nbi_file = 'beamlets.dat'
    open (newunit=unit, file='beamlets.dat', status='replace')
    if (trim(mode) == 'radial') then
        write (unit, '(a)') 'R 2 1.0e6 1.0e6 1.0 0.0 0.0 14.0 0.0 0.0 -1.0 0.0 0.0'
    else
 write (unit, '(a)') 'T 2 1.0e6 1.0e6 1.0 0.0 0.0 10.2265706505403 -8.0 0.0 0.0 1.0 0.0'
    end if
    close (unit)
    call simple_init
    allocate (z(5, ntestpart))
    call simple_sample(ntestpart, z)

    injected = -1.0_dp
    ionized = -1.0_dp
    shine = -1.0_dp
    open (newunit=unit, file='nbi_source.dat', status='old')
    do
        read (unit, '(a)', iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        read (line, *, iostat=ios) key, value
        if (ios /= 0) cycle
        if (key == 'injected_rate') injected = value
        if (key == 'ionized_rate') ionized = value
        if (key == 'shine_rate') shine = value
    end do
    close (unit)
if (min(injected, ionized, shine) < 0.0_dp) error stop 'rates missing in nbi_source.dat'
    if (abs((ionized + shine)/injected - 1.0_dp) > 1.0e-12_dp) &
        error stop 'ionized plus shine-through differs from injected'

    speed = sqrt(2.0_dp*energy*ev_joule/(2.0e-3_dp*p_mass))
    worst = 0.0_dp
    worst_pitch = 0.0_dp
    do i = 1, size(z, 2)
        lam = z(5, i)
        if (trim(mode) == 'radial' .and. abs(lam) > 0.2_dp) &
            error stop 'radial injection gives a large pitch'
        ! Independent field direction: VMEC contravariant components, moved
        ! from the straight-field-line angle to the VMEC angle, on libneo's
        ! covariant basis.
        call vmec_field(z(1, i), z(2, i), z(3, i), at, ap, dat, dap, iota, sqg, alam, &
                        dl_ds, dl_dt, dl_dp, bctr_t, bctr_p, bco_s, bco_t, bco_p)
        call ref_coords%covariant_basis(z(1:3, i), e_cov)
       bcart = (bctr_t - bctr_p*dl_dp)/(1.0_dp + dl_dt)*e_cov(:, 2) + bctr_p*e_cov(:, 3)
        direction = [0.0_dp, 1.0_dp, 0.0_dp]
        if (trim(mode) == 'radial') direction = [-1.0_dp, 0.0_dp, 0.0_dp]
   worst_pitch = max(worst_pitch, abs(lam - dot_product(direction, bcart)/norm2(bcart)))
        call ref_coords%evaluate_cart(z(1:3, i), xgc)
        call ref_to_integ(z(1:3, i), xi)
        call magfie(xi, bmod, sqrtg, bder, hcov, hctr, hcurl)
        ! Gyroradius in m from B in Gauss.
   rho = 2.0e-3_dp*p_mass*speed*sqrt(max(0.0_dp, 1.0_dp - lam**2))/(e_si*bmod*1.0e-4_dp)
        if (rho < 1.0e-3_dp) cycle
  worst = max(worst, abs(norm2(1.0e-2_dp*xgc - nbi_marker_position(:, i))/rho - 1.0_dp))
    end do
    print '(a,a,a,es10.3)', trim(mode), ': worst relative gyroradius error ', '', worst
    ! The shift uses the field at the particle, the check the field at the
    ! guiding centre; they differ at order rho/L.
    if (worst > 0.05_dp) error stop 'guiding centre is not one gyroradius from birth'
    print '(a,es10.3)', 'worst pitch difference to the independent field ', worst_pitch
    ! The pitch uses the field at the guiding centre in both routes.
    if (worst_pitch > 1.0e-3_dp) error stop 'pitch differs from the independent field'
    print '(a)', 'PASS'
end program test_nbi_birth_equilibrium
