module nbi_birth
    !> Guiding-centre markers born from neutral beam ionization in a VMEC
    !> equilibrium. Each marker carries an equal share of the total ionization
    !> rate and the toroidal angular momentum its neutral brought in.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use nbi_source, only: nbi_beamlet_t, nbi_plasma_t, nbi_segment_t, &
                          read_nbi_beamlets, beamlet_deposition, guiding_centre_shift, &
                          ev_joule
    use util, only: p_mass
    implicit none
    private
    public :: sample_nbi_birth, vmec_beam_plasma, nbi_marker_angmom, nbi_marker_rate, &
              nbi_marker_position, field_direction

    !> Toroidal angular momentum of each marker's neutral [kg m^2/s], about the
    !> machine axis along +Z, and the ionizations per second it represents.
    real(dp), allocatable :: nbi_marker_angmom(:), nbi_marker_position(:, :)
    real(dp) :: nbi_marker_rate = 0.0_dp

    real(dp), parameter :: e_cgs = 4.80320471e-10_dp, c_cgs = 2.99792458e10_dp
    ! Beam ions have mass A m_p, as in SIMPLE's orbit equations.
    real(dp), parameter :: proton_kg = 1.0e-3_dp*p_mass

contains

    function vmec_beam_plasma(x) result(p)
        !> Plasma at Cartesian x [m] from the flux label s of the profiles and
        !> SIMPLE profiles. Species 2 with charge above one acts as the impurity.
        use reference_coordinates, only: ref_coords
        use simple_profiles, only: get_plasma_params
        use params, only: Z1, Z2
        real(dp), intent(in) :: x(3)
        type(nbi_plasma_t) :: p
        real(dp) :: u(3), s, te, ti1, ti2, n1, n2
        integer :: ierr

        call ref_coords%from_cyl([100.0_dp*hypot(x(1), x(2)), atan2(x(2), x(1)), &
                                  100.0_dp*x(3)], u, ierr)
        p%inside = .false.
        if (ierr /= 0) return
        s = flux_label(u)
        p%inside = s >= 0.0_dp .and. s < 1.0_dp
        if (.not. p%inside) return
        call get_plasma_params(s, te, ti1, ti2, n1, n2)
        p%te = te
        p%ne = Z1*n1 + Z2*n2
        p%zeff = (Z1**2*n1 + Z2**2*n2)/p%ne
        p%z_impurity = 0
        if (Z2 > 1.5_dp .and. n2 > 0.0_dp) p%z_impurity = nint(Z2)
    end function vmec_beam_plasma

    function flux_label(u) result(s)
        !> Normalised toroidal flux of the reference point u. The chart radius
        !> of a chartmap is rho, and beyond the LCFS of an extended chart it
        !> exceeds rho_lcfs; only VMEC reference coordinates give s directly.
        use field_can_mod, only: ref_to_integ
        use profile_radius, only: profile_s
        real(dp), intent(in) :: u(3)
        real(dp) :: s
        real(dp) :: xi(3)

        call ref_to_integ(u, xi)
        s = profile_s(xi)
    end function flux_label

    subroutine field_direction(u, bhat, bmod)
        !> Unit field vector in Cartesian components at reference point u, from
        !> the contravariant unit vector and central differences of the map.
        use reference_coordinates, only: ref_coords
        use field_can_mod, only: integ_to_ref, ref_to_integ
        use magfie_sub, only: magfie
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: bhat(3), bmod
        real(dp), parameter :: h(3) = [1.0e-5_dp, 1.0e-5_dp, 1.0e-5_dp]
        real(dp) :: xi(3), xp(3), xm(3), up(3), um(3), cp(3), cm(3), e(3, 3)
        real(dp) :: sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)
        integer :: i

        call ref_to_integ(u, xi)
        call magfie(xi, bmod, sqrtg, bder, hcov, hctr, hcurl)
        do i = 1, 3
            xp = xi
            xm = xi
            xp(i) = xp(i) + h(i)
            xm(i) = xm(i) - h(i)
            ! Within h of the axis the backward point would lie at negative s.
            if (i == 1) xm(1) = max(xm(1), 0.0_dp)
            call integ_to_ref(xp, up)
            call integ_to_ref(xm, um)
            call ref_coords%evaluate_cart(up, cp)
            call ref_coords%evaluate_cart(um, cm)
            e(:, i) = (cp - cm)/(xp(i) - xm(i))
        end do
        bhat = matmul(e, hctr)
        bhat = bhat/norm2(bhat)
    end subroutine field_direction

    subroutine sample_nbi_birth(zstart, integration_coords)
        use reference_coordinates, only: ref_coords
        use field_can_mod, only: ref_to_integ
        use params, only: nbi_file, nbi_path_step, nbi_path_length, n_d, n_e, &
                          particle_energy_eV_effective, marker_weights_input, sbeg, &
                          num_surf
        use simple_mpi, only: mpi_is_root, mpi_ranks, mpi_broadcast_r64
        real(dp), intent(out) :: zstart(:, :)
        logical, intent(in) :: integration_coords
        type(nbi_beamlet_t), allocatable :: beamlets(:)
        type(nbi_segment_t), allocatable :: seg(:), all_seg(:)
        real(dp), allocatable :: cdf(:)
        integer, allocatable :: trial_ordinal(:)
        real(dp) :: injected, shine, total_injected, total_shine, total_ionized
     real(dp) :: injected_power, ionized_power, injected_torque, rnd(2), l, energy, mass
        real(dp) :: x(3), v(3), u(3), ugc(3), xgc(3), bhat(3), bmod, omega, speed
        real(dp), allocatable :: rate_table(:, :)
        real(dp) :: prompt_angmom, prompt_rate, prompt_torque
        integer :: ib, j, k, n, ipart, lo, hi, mid, ierr, unit, ntrial

        n = size(zstart, 2)
        if (n < 1) error stop 'nbi births require at least one marker'
        call read_nbi_beamlets(trim(nbi_file), beamlets)
        if (any(beamlets%mass /= nint(n_d)) .or. n_e /= 1) &
            error stop 'nbi births require n_d equal to the beam mass and n_e = 1'
  if (maxval(beamlets%energy_ev) > particle_energy_eV_effective*(1.0_dp + 1.0e-12_dp)) &
           error stop 'nbi births require particle energy at least the full beam energy'
        allocate (all_seg(0), rate_table(4, size(beamlets)*3))
        total_injected = 0.0_dp
        total_shine = 0.0_dp
        injected_power = 0.0_dp
        ionized_power = 0.0_dp
        injected_torque = 0.0_dp
        do ib = 1, size(beamlets)
            do j = 1, 3
             call beamlet_deposition(beamlets(ib), j, vmec_beam_plasma, nbi_path_step, &
                                        nbi_path_length, seg, injected, shine)
                seg%beamlet = ib
                all_seg = [all_seg, seg]
                energy = beamlets(ib)%energy_ev/j
                mass = beamlets(ib)%mass*proton_kg
                speed = sqrt(2.0_dp*energy*ev_joule/mass)
                total_injected = total_injected + injected
                total_shine = total_shine + shine
                injected_power = injected_power + injected*energy*ev_joule
                ionized_power = ionized_power + (injected - shine)*energy*ev_joule
                ! Angular momentum flux of the incident neutrals about +Z.
                injected_torque = injected_torque + injected*mass*speed* &
                                  (beamlets(ib)%origin(1)*beamlets(ib)%direction(2) - &
                                   beamlets(ib)%origin(2)*beamlets(ib)%direction(1))
            rate_table(:, 3*(ib - 1) + j) = [real(ib, dp), real(j, dp), injected, shine]
            end do
        end do
        total_ionized = sum(all_seg%rate)
   if (total_ionized <= 0.0_dp) error stop 'nbi births: no ionization inside the plasma'
        allocate (cdf(0:size(all_seg)))
        cdf(0) = 0.0_dp
        do k = 1, size(all_seg)
            cdf(k) = cdf(k - 1) + all_seg(k)%rate
        end do
        cdf = cdf/cdf(size(all_seg))

        if (allocated(nbi_marker_angmom)) deallocate (nbi_marker_angmom)
        if (allocated(nbi_marker_position)) deallocate (nbi_marker_position)
        allocate (nbi_marker_angmom(n), nbi_marker_position(3, n))
        allocate (trial_ordinal(n))
        ! Births whose guiding centre lies outside the last closed surface are
        ! prompt losses. They are counted and replaced by a new draw.
        ntrial = 0
        prompt_angmom = 0.0_dp
        ipart = 0
        do while (ipart < n)
            ntrial = ntrial + 1
        if (ntrial > 100*n) error stop 'nbi births: almost all births are prompt losses'
            call random_number(rnd)
            lo = 0
            hi = size(all_seg)
            do while (hi - lo > 1)
                mid = (lo + hi)/2
                if (cdf(mid) <= rnd(1)) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            k = hi
            ib = all_seg(k)%beamlet
            j = all_seg(k)%component
            l = all_seg(k)%l0 + rnd(2)*(all_seg(k)%l1 - all_seg(k)%l0)
            x = beamlets(ib)%origin + l*beamlets(ib)%direction
            energy = beamlets(ib)%energy_ev/j
            mass = beamlets(ib)%mass*proton_kg
            speed = sqrt(2.0_dp*energy*ev_joule/mass)
            v = speed*beamlets(ib)%direction
            ! Guiding centre in cm and cm/s, Gaussian units.
            call ref_coords%from_cyl([100.0_dp*hypot(x(1), x(2)), atan2(x(2), x(1)), &
                                      100.0_dp*x(3)], u, ierr)
            if (ierr /= 0) then
                prompt_angmom = prompt_angmom + mass*(x(1)*v(2) - x(2)*v(1))
                cycle
            end if
            if (flux_label(u) >= 1.0_dp) then
                ! The segment midpoint was inside, this point is not.
                prompt_angmom = prompt_angmom + mass*(x(1)*v(2) - x(2)*v(1))
                cycle
            end if
            call field_direction(u, bhat, bmod)
            omega = e_cgs*bmod/(1.0e3_dp*mass*c_cgs)
            xgc = guiding_centre_shift(100.0_dp*x, 100.0_dp*v, bhat, omega)
            if (.not. all(ieee_is_finite(xgc))) &
                error stop 'NBI birth: non-finite guiding centre'
      call ref_coords%from_cyl([hypot(xgc(1), xgc(2)), atan2(xgc(2), xgc(1)), xgc(3)], &
                                     ugc, ierr)
            if (ierr /= 0) then
                prompt_angmom = prompt_angmom + mass*(x(1)*v(2) - x(2)*v(1))
                cycle
            end if
            if (flux_label(ugc) >= 1.0_dp) then
                prompt_angmom = prompt_angmom + mass*(x(1)*v(2) - x(2)*v(1))
                cycle
            end if
            ipart = ipart + 1
            trial_ordinal(ipart) = ntrial
            nbi_marker_angmom(ipart) = mass*(x(1)*v(2) - x(2)*v(1))
            nbi_marker_position(:, ipart) = x
            call field_direction(ugc, bhat, bmod)
            if (integration_coords) then
                call ref_to_integ(ugc, zstart(1:3, ipart))
            else
                zstart(1:3, ipart) = ugc
            end if
            zstart(4, ipart) = sqrt(energy/particle_energy_eV_effective)
            zstart(5, ipart) = dot_product(v, bhat)/speed
        end do
        nbi_marker_rate = total_ionized/ntrial
        if (mpi_ranks > 1) then
            ! The tracing driver broadcasts rank 0's start states.
            block
                real(dp) :: rate(1)
                rate(1) = nbi_marker_rate
                call mpi_broadcast_r64(rate, 1)
                nbi_marker_rate = rate(1)
            end block
            call mpi_broadcast_r64(nbi_marker_angmom, n)
            call mpi_broadcast_r64(nbi_marker_position, size(nbi_marker_position))
        end if
        prompt_rate = total_ionized*(ntrial - n)/ntrial
        prompt_torque = prompt_angmom*total_ionized/ntrial
        if (allocated(marker_weights_input)) deallocate (marker_weights_input)
        allocate (marker_weights_input(n), source=nbi_marker_rate)
        num_surf = 2
        sbeg(1:2) = [0.0_dp, 1.0_dp]

        if (.not. mpi_is_root()) return
        open (newunit=unit, file='nbi_source.dat', status='replace')
        write (unit, '(a)') '# Neutral beam deposition, rates in 1/s, powers in W,'
        write (unit, '(a)') '# torque about +Z in N m. Each marker carries marker_rate.'
        write (unit, '(a,es24.16)') 'injected_rate ', total_injected
        write (unit, '(a,es24.16)') 'ionized_rate ', total_ionized
        write (unit, '(a,es24.16)') 'shine_rate ', total_shine
        write (unit, '(a,es24.16)') 'injected_power ', injected_power
        write (unit, '(a,es24.16)') 'ionized_power ', ionized_power
        write (unit, '(a,es24.16)') 'injected_torque ', injected_torque
        write (unit, '(a,es24.16)') 'prompt_loss_rate ', prompt_rate
        write (unit, '(a,es24.16)') 'prompt_loss_torque ', prompt_torque
        write (unit, '(a,es24.16)') 'marker_rate ', nbi_marker_rate
        write (unit, '(a,i0)') 'marker_count ', n
        write (unit, '(a)') '# beamlet component injected_rate shine_rate'
        do k = 1, size(rate_table, 2)
            write (unit, '(2i6,2es24.16)') nint(rate_table(1:2, k)), rate_table(3:4, k)
        end do
        close (unit)
        open (newunit=unit, file='nbi_births.dat', status='replace')
        write (unit, '(a)') '# marker neutral_angular_momentum_kg_m2_s birth_x_m birth_y_m birth_z_m'
        do ipart = 1, n
            write (unit, '(i10,4es24.16)') ipart, nbi_marker_angmom(ipart), &
                nbi_marker_position(:, ipart)
        end do
        close (unit)
        open (newunit=unit, file='nbi_birth_trials.dat', status='replace')
        write (unit, '(a)') '# marker trial_ordinal'
        do ipart = 1, n
            write (unit, '(i0,1x,i0)') ipart, trial_ordinal(ipart)
        end do
        close (unit)
        print *, 'NBI ionized power [W] = ', ionized_power, ' (nbi_source.dat)'
    end subroutine sample_nbi_birth

end module nbi_birth
