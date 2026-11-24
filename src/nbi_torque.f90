module nbi_torque
    !> Toroidal torque of neutral beam ions on the bulk plasma, from the
    !> guiding-centre angular momentum L = R m v_par B_phi / abs(B) about +Z
    !> and the reference poloidal flux -A_phi(s) in the VMEC flux gauge.
    !> The bulk receives L_neutral - L at birth (the radial current of the
    !> ionization displacement), -dL at each collision, (q/c) dA_phi along the
    !> orbit (j x B of the radial fast-ion current), and the remaining L at
    !> thermalization. In an axisymmetric field L + (q/c) A_phi is conserved
    !> between collisions. In 3D the field ledger records the residual
    !> -dL - (q/c) dA_phi relative to the reference flux torque.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: nbi_tally_enabled, nbi_tally_init, nbi_tally_begin, nbi_tally_drift, &
              nbi_tally_collision, nbi_tally_end, nbi_tally_write, nbi_gc_momentum, &
              nbi_tally_total

    logical :: nbi_tally_enabled = .false.
    integer :: nbins = 0
    ! Per-bin angular momentum transfer [kg m^2/s per marker].
    real(dp), allocatable :: dep_birth(:), dep_coll(:), dep_jxb(:), dep_therm(:), &
                             dep_field(:), dep_neutron(:)
    real(dp), allocatable :: last_l(:), last_a(:), last_s(:), neutral_l(:)
    real(dp) :: lost_l = 0.0_dp, alive_l = 0.0_dp, marker_rate = 0.0_dp
    integer :: n_lost = 0, n_therm = 0, n_alive = 0, n_begin = 0

    ! g cm^2/s to kg m^2/s
    real(dp), parameter :: to_si = 1.0e-7_dp

contains

    subroutine nbi_tally_init(n_bins, angmom, rate)
        use fusion_source, only: init_hermite
        integer, intent(in) :: n_bins
        real(dp), intent(in) :: angmom(:), rate
        integer :: n

        n = size(angmom)
        nbins = n_bins
        if (allocated(dep_birth)) deallocate (dep_birth, dep_coll, dep_jxb, dep_therm, &
                                              dep_field, dep_neutron, last_l, last_a, &
                                              last_s, neutral_l)
        allocate (dep_birth(nbins), dep_coll(nbins), dep_jxb(nbins), dep_therm(nbins), &
                  dep_field(nbins), dep_neutron(nbins), source=0.0_dp)
        allocate (last_l(n), last_a(n), last_s(n), source=0.0_dp)
        allocate (neutral_l(n), source=angmom)
        lost_l = 0.0_dp
        alive_l = 0.0_dp
        n_lost = 0
        n_therm = 0
        n_alive = 0
        marker_rate = rate
        nbi_tally_enabled = .true.
        ! Quadrature nodes before the parallel orbit loop.
        call init_hermite
    end subroutine nbi_tally_init

    subroutine nbi_gc_momentum(z, l, a, s)
        !> Axial momentum and reference flux term [kg m^2/s] at state z.
        use field_can_mod, only: integ_to_ref
        use spline_vmec_sub, only: splint_vmec_data, compute_field_components
        use params, only: v0, n_d, n_e
        use util, only: c, e_charge, p_mass
        real(dp), intent(in) :: z(5)
        real(dp), intent(out) :: l, a, s
        real(dp) :: u(3), at, ap, dat, dap, iota, sqg, alam, dl_ds, dl_dt, dl_dp
        real(dp) :: r, height, dr_ds, dr_dt, dr_dp, dz_ds, dz_dt, dz_dp
        real(dp) :: bctr_t, bctr_p, bco_s, bco_t, bco_p, bmod

        call integ_to_ref(z(1:3), u)
        call splint_vmec_data(u(1), u(2), u(3), ap, at, dap, dat, iota, &
                              r, height, alam, dr_ds, dr_dt, dr_dp, dz_ds, dz_dt, &
                              dz_dp, dl_ds, dl_dt, dl_dp)
        call compute_field_components(r, dr_ds, dr_dt, dr_dp, dz_ds, dz_dt, dz_dp, &
                                      dat, dap, dl_ds, dl_dt, dl_dp, sqg, bctr_t, &
                                      bctr_p, bco_s, bco_t, bco_p)
        bmod = sqrt(bctr_t*bco_t + bctr_p*bco_p)
        ! The covariant VMEC component includes R and Z changes at fixed s, theta.
        l = n_d*p_mass*v0*z(4)*z(5)*r**2*bctr_p/bmod*to_si
        a = n_e*e_charge/c*ap*to_si
        s = u(1)
    end subroutine nbi_gc_momentum

    pure integer function bin_of(s)
        real(dp), intent(in) :: s
        bin_of = min(nbins, max(1, int(s*nbins) + 1))
    end function bin_of

    subroutine add(array, s, value)
        real(dp), intent(inout) :: array(:)
        real(dp), intent(in) :: s, value
        integer :: k
        k = bin_of(s)
        !$omp atomic update
        array(k) = array(k) + value
    end subroutine add

    subroutine nbi_tally_begin(ipart, z)
        integer, intent(in) :: ipart
        real(dp), intent(in) :: z(5)
        call nbi_gc_momentum(z, last_l(ipart), last_a(ipart), last_s(ipart))
        call add(dep_birth, last_s(ipart), neutral_l(ipart) - last_l(ipart))
        !$omp atomic update
        n_begin = n_begin + 1
    end subroutine nbi_tally_begin

    subroutine nbi_tally_drift(ipart, z, dt)
        integer, intent(in) :: ipart
        real(dp), intent(in) :: z(5), dt
        real(dp) :: l, a, s, smid

        call nbi_gc_momentum(z, l, a, s)
        smid = 0.5_dp*(s + last_s(ipart))
        call add(dep_jxb, smid, a - last_a(ipart))
        call add(dep_field, smid, -(l - last_l(ipart)) - (a - last_a(ipart)))
        last_l(ipart) = l
        last_a(ipart) = a
        last_s(ipart) = s
        call add(dep_neutron, s, beam_target_neutrons(z, s, dt))
    end subroutine nbi_tally_drift

    function beam_target_neutrons(z, s, dt) result(yield)
        !> Expected D(d,n)3He reactions on Maxwellian deuterium during the
        !> accepted drift duration dt [s]: n_D <sigma u> dt.
        use params, only: n_d, v0, am1, am2
        use simple_profiles, only: get_plasma_params
        use fusion_source, only: dd_beam_maxwellian_rate
        use util, only: p_mass
        real(dp), intent(in) :: z(5), s, dt
        real(dp) :: yield, te, ti1, ti2, n1, n2, nd, ti

        yield = 0.0_dp
        if (nint(n_d) /= 2) return
        call get_plasma_params(s, te, ti1, ti2, n1, n2)
        nd = 0.0_dp
        ti = 0.0_dp
        if (nint(am2) == 2) then
            nd = nd + n2
            ti = ti2
        end if
        if (nint(am1) == 2) then
            nd = nd + n1
            ti = ti1
        end if
        ! Speed in m/s; SIMPLE's ions have mass n_d m_p, p_mass in g.
        yield = nd*dd_beam_maxwellian_rate(1.0e-2_dp*v0*z(4), ti, 2.0e-3_dp*p_mass)*dt
    end function beam_target_neutrons

    subroutine nbi_tally_collision(ipart, z)
        !> Collisions change only the velocity at fixed position.
        integer, intent(in) :: ipart
        real(dp), intent(in) :: z(5)
        real(dp) :: l, a, s

        call nbi_gc_momentum(z, l, a, s)
        call add(dep_coll, s, -(l - last_l(ipart)))
        last_l(ipart) = l
    end subroutine nbi_tally_collision

    subroutine nbi_tally_end(ipart, thermalized, lost)
        integer, intent(in) :: ipart
        logical, intent(in) :: thermalized, lost
        if (thermalized) then
            call add(dep_therm, last_s(ipart), last_l(ipart))
            !$omp atomic update
            n_therm = n_therm + 1
        else if (lost) then
            !$omp atomic update
            lost_l = lost_l + last_l(ipart)
            !$omp atomic update
            n_lost = n_lost + 1
        else
            !$omp atomic update
            alive_l = alive_l + last_l(ipart)
            !$omp atomic update
            n_alive = n_alive + 1
        end if
    end subroutine nbi_tally_end

    pure function nbi_tally_total() result(total)
        !> Bulk torque, field residual, lost and still-fast angular momentum
        !> rates [N m], the injected rate of the markers, and the beam-target
        !> neutron rate [1/s].
        real(dp) :: total(6)
        total(1) = marker_rate*(sum(dep_birth) + sum(dep_coll) + sum(dep_jxb) + &
                                sum(dep_therm))
        total(2) = marker_rate*sum(dep_field)
        total(3) = marker_rate*lost_l
        total(4) = marker_rate*alive_l
        total(5) = marker_rate*sum(neutral_l)
        total(6) = marker_rate*sum(dep_neutron)
    end function nbi_tally_total

    subroutine nbi_tally_write(filename)
        !> Torque per radial bin in N m and the bin volume from the VMEC metric.
        use reference_coordinates, only: ref_coords
        character(len=*), intent(in) :: filename
        integer, parameter :: ns = 8, nt = 64, np = 128
        real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)
        real(dp) :: volume(nbins), g(3, 3), ginv(3, 3), sqrtg, u(3), total(6)
        integer :: k, i, j, m, unit

        volume = 0.0_dp
        do k = 1, nbins
            do i = 1, ns
                u(1) = (k - 1 + (i - 0.5_dp)/ns)/nbins
                do j = 1, nt
                    u(2) = twopi*(j - 0.5_dp)/nt
                    do m = 1, np
                        u(3) = twopi*(m - 0.5_dp)/np
                        call ref_coords%metric_tensor(u, g, ginv, sqrtg)
                        volume(k) = volume(k) + abs(sqrtg)
                    end do
                end do
            end do
            volume(k) = volume(k)*twopi**2/(nbins*ns*nt*np)*1.0e-6_dp
        end do
        total = nbi_tally_total()
        open (newunit=unit, file=filename, status='replace')
        write (unit, '(a)') '# Beam torque on the bulk plasma about +Z in N m per bin of VMEC s.'
        write (unit, '(a,es24.16)') '# bulk_total ', total(1)
        write (unit, '(a,es24.16)') '# field_total ', total(2)
        write (unit, '(a,es24.16)') '# lost_angular_momentum_rate ', total(3)
        write (unit, '(a,es24.16)') '# still_fast_angular_momentum_rate ', total(4)
       write (unit, '(a,es24.16)') '# injected_ionized_angular_momentum_rate ', total(5)
        write (unit, '(a,es24.16)') '# beam_target_neutron_rate_per_s ', total(6)
      write (unit, '(a,4i8)') '# markers begun thermalized lost still_fast ', n_begin, &
            n_therm, n_lost, n_alive
     write (unit, '(a)') '# s_lo s_hi volume_m3 birth collisional jxb thermalized field neutrons_per_s'
        do k = 1, nbins
           write (unit, '(2f10.5,7es20.11)') real(k - 1, dp)/nbins, real(k, dp)/nbins, &
                volume(k), marker_rate*[dep_birth(k), dep_coll(k), dep_jxb(k), &
                                        dep_therm(k), dep_field(k), dep_neutron(k)]
        end do
        close (unit)
    end subroutine nbi_tally_write

end module nbi_torque
