module thermal_ntv
    !> Linear delta-f tally of the thermal NTV torque density.
    !>
    !> The background is the canonical Maxwellian F_c(P, H), the local Maxwellian
    !> with s replaced by the flux label of the canonical toroidal momentum P. It
    !> solves the unperturbed collisionless equation, so w = delta f / F_c = 0 at
    !> birth is the steady state without the perturbation. A local Maxwellian
    !> F_M(s) is not: with dw = -A ds the weights carry the bounce displacement
    !> s(0) - <s> of their birth, which only collisions remove. Along each
    !> collisionless step the weight follows
    !>   dw = -(A_1 + A_2 E/T) dP / (dP/ds),
    !> with A_1 = n'/n + e Phi'/T - 3/2 T'/T and A_2 = T'/T, ' = d/ds, the
    !> derivative of ln F_c at fixed total energy. Collisions move markers at
    !> fixed weight; the source C[F_c] of the collisional axisymmetric response
    !> is not included. The torque on a particle from the field is
    !> the change of canonical toroidal momentum P = (q/c) psi* along a
    !> collisionless step, zero in an axisymmetric field. The torque density on
    !> the species is n <w dP/dt>, estimated per marker by
    !>   c_i = int w dP,  d_i = int dt
    !> over the steps after window_start with the marker in the bin [s_lo, s_hi],
    !> and n sum c_i / sum d_i.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: ntv_tally_enabled, ntv_tally_init, ntv_tally_begin, ntv_tally_drift, &
              ntv_tally_collision, ntv_tally_write, ntv_tally_weight, ntv_tally_marker

    logical :: ntv_tally_enabled = .false.
    real(dp) :: a1 = 0, a2 = 0, energy_ratio = 0, s_lo = 0, s_hi = 1, t_start = 0
    real(dp), allocatable :: w(:), last_p(:), last_s(:), contrib(:), dwell(:)

    ! g cm^2/s to kg m^2/s
    real(dp), parameter :: to_si = 1.0e-7_dp

contains

    subroutine ntv_tally_init(n_markers, drive, temperature_ev, bin, window_start)
        use params, only: particle_energy_eV_effective, integmode, isw_field_type
        use magfie_sub, only: BOOZER
        integer, intent(in) :: n_markers
        real(dp), intent(in) :: drive(2), temperature_ev, bin(2), window_start
        if (temperature_ev <= 0) error stop 'thermal_ntv: temperature must be > 0'
        if (integmode <= 0) error stop 'thermal_ntv: needs a symplectic integmode'
        if (isw_field_type /= BOOZER) error stop 'thermal_ntv: needs Boozer coordinates'
        energy_ratio = particle_energy_eV_effective/temperature_ev
        a1 = drive(1)
        a2 = drive(2)
        s_lo = bin(1)
        s_hi = bin(2)
        t_start = window_start
        if (allocated(w)) deallocate (w, last_p, last_s, contrib, dwell)
        allocate (w(n_markers), last_p(n_markers), last_s(n_markers), &
                  contrib(n_markers), dwell(n_markers), source=0.0_dp)
        ntv_tally_enabled = .true.
    end subroutine ntv_tally_init

    subroutine momentum_and_s(z, p, s)
        !> Canonical toroidal momentum [kg m^2/s] and flux label of state z,
        !> with psi* as in encode_symplectic_state.
        use field_can_mod, only: eval_field => evaluate, field_can_t
        use parmot_mod, only: ro0
        real(dp), intent(in) :: z(5)
        real(dp), intent(out) :: p, s
        type(field_can_t) :: f

        call eval_field(f, z(1), z(2), z(3), 0)
        p = momentum_unit()*(f%Aph + ro0*z(4)*z(5)*f%hph)
        s = flux_label(z)
    end subroutine momentum_and_s

    real(dp) function momentum_unit()
        !> (q/c) in kg m^2/s per G cm.
        use params, only: n_e
        use util, only: c, e_charge
        momentum_unit = n_e*e_charge/c*to_si
    end function momentum_unit

    real(dp) function momentum_per_s(z)
        !> dP/ds = (q/c) dA_phi/ds at the marker [kg m^2/s].
        use field_can_mod, only: eval_field => evaluate, field_can_t
        real(dp), intent(in) :: z(5)
        type(field_can_t) :: f
        call eval_field(f, z(1), z(2), z(3), 0)
        momentum_per_s = momentum_unit()*f%dAph(1)
    end function momentum_per_s

    pure real(dp) function flux_label(z)
        !> s, the first Boozer integration coordinate (ntv_tally_init).
        real(dp), intent(in) :: z(5)
        flux_label = z(1)
    end function flux_label

    subroutine ntv_tally_begin(ipart, z)
        integer, intent(in) :: ipart
        real(dp), intent(in) :: z(5)
        w(ipart) = 0
        contrib(ipart) = 0
        dwell(ipart) = 0
        call momentum_and_s(z, last_p(ipart), last_s(ipart))
    end subroutine ntv_tally_begin

    subroutine ntv_tally_drift(ipart, z, time, pphi)
        !> One collisionless step ending at time [s]. pphi is the integrator's
        !> normalized canonical momentum, conserved exactly in axisymmetry;
        !> psi* rebuilt from z is not, since the integrator's v_par belongs to
        !> its intermediate point.
        use params, only: dtaumin, v0
        use parmot_mod, only: ro0
        integer, intent(in) :: ipart
        real(dp), intent(in) :: z(5), time, pphi
        real(dp) :: p, s, dw, w_mid

        p = momentum_unit()*ro0/sqrt(2.0_dp)*pphi
        s = flux_label(z)
        dw = -(a1 + a2*energy_ratio*z(4)**2)*(p - last_p(ipart))/momentum_per_s(z)
        w_mid = w(ipart) + 0.5_dp*dw
        w(ipart) = w(ipart) + dw
        if (time >= t_start .and. 0.5_dp*(s + last_s(ipart)) >= s_lo .and. &
            0.5_dp*(s + last_s(ipart)) < s_hi) then
            contrib(ipart) = contrib(ipart) + w_mid*(p - last_p(ipart))
            dwell(ipart) = dwell(ipart) + dtaumin/v0
        end if
        last_p(ipart) = p
        last_s(ipart) = s
    end subroutine ntv_tally_drift

    subroutine ntv_tally_collision(ipart, z)
        !> Collisions change P without a torque from the field.
        integer, intent(in) :: ipart
        real(dp), intent(in) :: z(5)
        real(dp) :: s
        call momentum_and_s(z, last_p(ipart), s)
    end subroutine ntv_tally_collision

    real(dp) function ntv_tally_weight(ipart)
        integer, intent(in) :: ipart
        ntv_tally_weight = w(ipart)
    end function ntv_tally_weight

    subroutine ntv_tally_marker(ipart, c_i, d_i)
        integer, intent(in) :: ipart
        real(dp), intent(out) :: c_i, d_i
        c_i = contrib(ipart)
        d_i = dwell(ipart)
    end subroutine ntv_tally_marker

    subroutine ntv_tally_write(filename)
        !> One row per marker: c_i [kg m^2/s], d_i [s], final weight.
        !> Torque density per unit density [N m] is sum c / sum d.
        character(*), intent(in) :: filename
        integer :: unit, i
        open (newunit=unit, file=filename, status='replace')
        write (unit, '(a,4es16.8)') '# A1 A2 s_lo s_hi ', a1, a2, s_lo, s_hi
        write (unit, '(a,es16.8)') '# sum c / sum d [N m per particle] ', &
            sum(contrib)/max(sum(dwell), tiny(1.0_dp))
        do i = 1, size(w)
            write (unit, '(3es24.15)') contrib(i), dwell(i), w(i)
        end do
        close (unit)
    end subroutine ntv_tally_write

end module thermal_ntv
