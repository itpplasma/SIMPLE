module field_can_boozer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: apply_radial_electric_potential, &
                              count_field_evaluation, field_can_t, twopi

    implicit none

    public :: integ_to_ref_boozer_chartmap, ref_to_integ_boozer_chartmap
    public :: evaluate_boozer_regular

contains

    subroutine evaluate_boozer_regular(x, y, phi, mode, rf)
        !> Regular near-axis Boozer field at X = sqrt(s) cos(theta_B),
        !> Y = sqrt(s) sin(theta_B), phi_B (libneo boozer_regular). mode 0, 1, 2
        !> returns values, first and second derivatives in (X, Y, phi); every
        !> component is C^2 through the magnetic axis. Requires the field to be
        !> initialised with field_can_mod regular_axis_field set.
        use boozer_regular, only: regular_field_t, eval_boozer_regular
        real(dp), intent(in) :: x, y, phi
        integer, intent(in) :: mode
        type(regular_field_t), intent(inout) :: rf

        call eval_boozer_regular(x, y, phi, mode, rf)

        ! Only mode 2 carries second derivatives, the d2 class of the counter.
        if (mode >= 2) then
            call count_field_evaluation(1)
        else
            call count_field_evaluation(0)
        end if
    end subroutine evaluate_boozer_regular

    subroutine integ_to_ref_boozer_chartmap(xinteg, xref)
        !> Convert integrator coords (s, theta_B, phi_B) to
        !> chartmap reference coords (rho, theta_B, phi_B).
        real(dp), intent(in) :: xinteg(3)
        real(dp), intent(out) :: xref(3)

        xref(1) = sqrt(max(xinteg(1), 0.0_dp))
        xref(2) = xinteg(2)
        xref(3) = xinteg(3)
    end subroutine integ_to_ref_boozer_chartmap

    subroutine ref_to_integ_boozer_chartmap(xref, xinteg)
        !> Convert chartmap reference coords (rho, theta_B, phi_B) to
        !> integrator coords (s, theta_B, phi_B).
        real(dp), intent(in) :: xref(3)
        real(dp), intent(out) :: xinteg(3)

        xinteg(1) = xref(1)**2
        xinteg(2) = xref(2)
        xinteg(3) = xref(3)
    end subroutine ref_to_integ_boozer_chartmap

    subroutine evaluate_boozer(f, r, th_c, ph_c, mode_secders)
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: r, th_c, ph_c
        integer, intent(in) :: mode_secders

        call eval_field_booz(f, r, th_c, ph_c, mode_secders)
        call apply_radial_electric_potential(f, abs(r))
        if (r < 0.0_dp) then
            ! Axis continuation (#370): the potential is a flux function, so
            ! its value lives at |s|, while its radial derivative carries the
            ! sign of the signed label.
            f%dPhie(1) = -f%dPhie(1)
        end if

        call count_field_evaluation(mode_secders)
    end subroutine evaluate_boozer

    subroutine integ_to_ref_boozer(xinteg, xref)
        use boozer_sub, only: boozer_to_vmec
        real(dp), intent(in) :: xinteg(3)
        real(dp), intent(out) :: xref(3)

        xref(1) = xinteg(1)
        call boozer_to_vmec( &
            xinteg(1), mod(xinteg(2), twopi), mod(xinteg(3), twopi), xref(2), xref(3))
        xref(2) = mod(xref(2), twopi)
        xref(3) = mod(xref(3), twopi)
    end subroutine integ_to_ref_boozer

    subroutine ref_to_integ_boozer(xref, xinteg)
        use boozer_sub, only: vmec_to_boozer
        use new_vmec_stuff_mod, only: nper
        real(dp), intent(in) :: xref(3)
        real(dp), intent(out) :: xinteg(3)
        real(dp) :: phi_ref, fper

        xinteg(1) = xref(1)
        phi_ref = mod(xref(3), twopi)
        call vmec_to_boozer( &
            xref(1), mod(xref(2), twopi), phi_ref, xinteg(2), xinteg(3))
        ! vmec_to_boozer returns phi_B in the first field period; restore the
        ! whole periods of the reference angle so the start keeps its sector.
        fper = twopi/real(nper, dp)
        xinteg(3) = xinteg(3) + fper*real(nint((phi_ref - xinteg(3))/fper), dp)
        xinteg(2) = mod(xinteg(2), twopi)
        xinteg(3) = mod(xinteg(3), twopi)
    end subroutine ref_to_integ_boozer

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine eval_field_booz(f, r, th_c, ph_c, mode_secders)
        use boozer_sub, only: splint_boozer_coord, boozer_state
        !
        ! Evaluates magnetic field in Boozer canonical coordinates (r, th_c, ph_c)
        ! and stores results in variable f
        ! Works for A_th linear in r (toroidal flux as radial variable)
        !
        ! Sign convention: A_theta = torflux*r is the toroidal flux psi_tor (on the
        ! poloidal covariant component), with torflux = -phi_vmec/(2*pi) > 0 since
        ! VMEC is left-handed (signgs = -1). A_phi = -chi is the toroidal covariant
        ! component carrying the poloidal flux. See libneo spline_vmec_data /
        ! vmecinm_m and gvec_simple/CONVENTIONS.md.
        !
        ! mode_secders = 0: no second derivatives
        ! mode_secders = 1: second derivatives only in d/dr^2
        ! mode_secders = 2: all second derivatives, including mixed
        !

        type(field_can_t), intent(inout) :: f
        double precision, intent(in) :: r, th_c, ph_c
        integer, intent(in) :: mode_secders

        double precision :: d3Aphdr3, dummy, dummy3(3), dummy6(6), &
            Bth, Bph, dBth, dBph, d2Bth, d2Bph
        double precision :: r_cont, th_cont
        logical :: crossed

        ! Count evaluations for profiling

        ! initialize to zero - no angular derivatives will be set due to straight field
        ! line Ath(r) Aph(r)
        f%dAth = 0d0
        f%dAph = 0d0

        ! initialize all 2nd derivatives to zero, as the mode decides, which one to use
        f%d2Ath = 0d0
        f%d2Aph = 0d0
        f%d2hth = 0d0
        f%d2hph = 0d0
        f%d2Bmod = 0d0

        ! Axis continuation (#370): a negative radial label s denotes the
        ! physical point (|s|, theta+pi) on the opposite ray. Evaluate the
        ! field there and fold the sign of the signed label into the odd
        ! radial derivatives, so the extended Hamiltonian is consistent
        ! through the axis and Newton iterates may cross it like any other
        ! point.
        crossed = r < 0.0_dp
        if (crossed) then
            r_cont = -r
            th_cont = th_c + 0.5_dp*twopi
        else
            r_cont = r
            th_cont = th_c
        end if

        call splint_boozer_coord(r_cont, th_cont, ph_c, mode_secders, &
                                 f%Ath, f%Aph, f%dAth(1), f%dAph(1), f%d2Aph(1), &
                                 d3Aphdr3, &
                                 Bth, dBth, d2Bth, &
                                 Bph, dBph, d2Bph, &
                                 f%Bmod, f%dBmod, f%d2Bmod, dummy, dummy3, dummy6)

        if (crossed) then
            f%dAth(1) = -f%dAth(1)
            f%dAph(1) = -f%dAph(1)
            dBth = -dBth
            dBph = -dBph
            f%dBmod(1) = -f%dBmod(1)
            f%d2Bmod(2:3) = -f%d2Bmod(2:3)
        end if

        if (mode_secders == 2) then
            call boozer_h_all(f, Bth, dBth, d2Bth, Bph, dBph, d2Bph)
        else if (mode_secders > 0) then
            call boozer_h_rr(f, Bth, dBth, d2Bth, Bph, dBph, d2Bph)
        else
            call boozer_h_first(f, Bth, dBth, Bph, dBph)
        end if
    end subroutine eval_field_booz

    subroutine eval_field_booz_device(f, r, th_c, ph_c, mode_secders)
        ! Device-only counterpart of eval_field_booz. Initialization checks and
        ! host diagnostic counters stay in splint_boozer_coord's host wrapper.
        ! The radial-mixed Newton evaluation is gpu_eval_field_rmix in simple_gpu.
        use boozer_sub, only: splint_boozer_coord_device

        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: r, th_c, ph_c
        integer, intent(in) :: mode_secders
        !$acc routine seq

        real(dp) :: d3Aphdr3, dummy, dummy3(3), dummy6(6), &
                    Bth, Bph, dBth, dBph, d2Bth, d2Bph
        real(dp) :: r_cont, th_cont
        logical :: crossed

        f%dAth = 0.0_dp
        f%dAph = 0.0_dp
        f%d2Ath = 0.0_dp
        f%d2Aph = 0.0_dp
        f%d2hth = 0.0_dp
        f%d2hph = 0.0_dp

        ! Axis continuation (#370), mirroring eval_field_booz: a negative
        ! radial label denotes the physical point (|s|, theta+pi).
        crossed = r < 0.0_dp
        if (crossed) then
            r_cont = -r
            th_cont = th_c + 0.5_dp*twopi
        else
            r_cont = r
            th_cont = th_c
        end if

        call splint_boozer_coord_device(r_cont, th_cont, ph_c, mode_secders, &
                                        f%Ath, f%Aph, f%dAth(1), f%dAph(1), &
                                        f%d2Aph(1), d3Aphdr3, &
                                        Bth, dBth, d2Bth, Bph, dBph, d2Bph, &
                                        f%Bmod, f%dBmod, f%d2Bmod, &
                                        dummy, dummy3, dummy6)

        if (crossed) then
            f%dAth(1) = -f%dAth(1)
            f%dAph(1) = -f%dAph(1)
            dBth = -dBth
            dBph = -dBph
            f%dBmod(1) = -f%dBmod(1)
            f%d2Bmod(2:3) = -f%d2Bmod(2:3)
        end if

        if (mode_secders == 2) then
            call boozer_h_all(f, Bth, dBth, d2Bth, Bph, dBph, d2Bph)
        else if (mode_secders > 0) then
            call boozer_h_rr(f, Bth, dBth, d2Bth, Bph, dBph, d2Bph)
        else
            call boozer_h_first(f, Bth, dBth, Bph, dBph)
        end if
    end subroutine eval_field_booz_device

    include '../generated/boozer_h_derivatives.inc'

end module field_can_boozer
