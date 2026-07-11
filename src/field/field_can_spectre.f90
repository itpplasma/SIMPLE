module field_can_spectre
    !> Per-volume Meiss canonical coordinates for SPECTRE equilibria (#439).
    !>
    !> SPECTRE's A_s = 0 gauge leaves the covariant h_s /= 0 (metric off-diagonals),
    !> so the phase-space 1-form keeps a rho_par*B_s ds term. Each SPECTRE volume is
    !> smooth in its own stacked-rho slab r in [lvol-1, lvol]; a global chart would
    !> integrate the Meiss construction ODE through the interface current sheets.
    !> This module runs the existing Meiss machinery once per volume on the identity
    !> radial scaling (the SPECTRE chart is its own radial coordinate) and harvests
    !> the resulting batch splines into a per-volume slot.
    !>
    !> Construction stays in the field's SI covariant units. The Gaussian-CGS
    !> conversion SIMPLE integrates in is applied only here at evaluation, uniformly
    !> to values and derivatives because r, theta, zeta are all dimensionless:
    !> Bmod Tesla->Gauss, covariant h meter->cm, covariant A (flux-like, [B]*L^2)
    !> Tesla*m^2 -> Gauss*cm^2.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    use field_spectre, only: spectre_field_t
    use field_can_base, only: field_can_t, n_field_evaluations, twopi
    use field_can_meiss, only: meiss_volume_t, init_meiss, get_meiss_coordinates, &
                               harvest_meiss_volume, cleanup_meiss
    use coordinate_scaling, only: identity_scaling_t
    use interpolate, only: evaluate_batch_splines_3d, &
                           evaluate_batch_splines_3d_der, &
                           evaluate_batch_splines_3d_der2
    use magfie_sub, only: TESLA_TO_GAUSS, M_TO_CM

    implicit none
    private

    public :: construct_spectre_coordinates, evaluate_spectre_can, &
              integ_to_ref_spectre, ref_to_integ_spectre, cleanup_spectre, &
              spectre_volumes, spectre_mvol

    integer, parameter :: N_R_VOL = 48, N_TH_VOL = 48, N_PHI_VOL = 32
    real(dp), parameter :: AXIS_OFFSET = 1.0e-3_dp
    real(dp), parameter :: A_SCALE = TESLA_TO_GAUSS*M_TO_CM**2
    real(dp), parameter :: H_SCALE = M_TO_CM
    real(dp), parameter :: B_SCALE = TESLA_TO_GAUSS

    type(meiss_volume_t), allocatable :: spectre_volumes(:)
    integer :: spectre_mvol = 0

contains

    subroutine construct_spectre_coordinates(field_noncan)
        !> Build one Meiss chart per SPECTRE volume on r in [lvol-1, lvol]. The axis
        !> volume starts at AXIS_OFFSET so the construction never touches the r = 0
        !> coordinate singularity.
        class(magnetic_field_t), intent(in) :: field_noncan

        integer :: lvol
        real(dp) :: rmin, rmax

        call cleanup_spectre

        select type (fld => field_noncan)
        type is (spectre_field_t)
            spectre_mvol = fld%data%Mvol
        class default
            error stop 'construct_spectre_coordinates: field is not spectre_field_t'
        end select

        allocate (spectre_volumes(spectre_mvol))

        do lvol = 1, spectre_mvol
            rmin = real(lvol - 1, dp)
            if (lvol == 1) rmin = rmin + AXIS_OFFSET
            rmax = real(lvol, dp)

            call init_meiss(field_noncan, N_R_VOL, N_TH_VOL, N_PHI_VOL, rmin, rmax, &
                            0.0_dp, twopi, identity_scaling_t())
            call get_meiss_coordinates
            call harvest_meiss_volume(spectre_volumes(lvol))
        end do

        call cleanup_meiss
    end subroutine construct_spectre_coordinates


    integer function volume_index(r) result(lvol)
        real(dp), intent(in) :: r

        lvol = min(int(r) + 1, spectre_mvol)
        lvol = max(1, lvol)
    end function volume_index


    pure function clamp_to_volume(r, vol) result(rc)
        real(dp), intent(in) :: r
        type(meiss_volume_t), intent(in) :: vol
        real(dp) :: rc

        ! Newton iterates inside a symplectic step may probe just past the volume
        ! edge. The splines carry no valid polynomial extension there, so clamp the
        ! evaluation to the construction domain edge. This is NOT the direct-basis
        ! path where clamping was rejected: the accepted state leaving the volume
        ! still triggers a boundary-stop in the caller.
        rc = max(vol%rmin, min(vol%rmax, r))
    end function clamp_to_volume


    subroutine evaluate_spectre_can(f, r, th_c, ph_c, mode_secders)
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: r, th_c, ph_c
        integer, intent(in) :: mode_secders

        integer :: lvol
        real(dp) :: x(3)
        real(dp) :: y(5), dy(3, 5), d2y(6, 5)

        n_field_evaluations = n_field_evaluations + 1

        lvol = volume_index(r)
        x = [clamp_to_volume(r, spectre_volumes(lvol)), th_c, ph_c]

        if (mode_secders > 0) then
            call evaluate_batch_splines_3d_der2(spectre_volumes(lvol)%spl_field, &
                                                x, y, dy, d2y)
            f%d2Ath = d2y(:, 1)*A_SCALE
            f%d2Aph = d2y(:, 2)*A_SCALE
            f%d2hth = d2y(:, 3)*H_SCALE
            f%d2hph = d2y(:, 4)*H_SCALE
            f%d2Bmod = d2y(:, 5)*B_SCALE
        else
            call evaluate_batch_splines_3d_der(spectre_volumes(lvol)%spl_field, &
                                               x, y, dy)
        end if

        f%Ath = y(1)*A_SCALE
        f%Aph = y(2)*A_SCALE
        f%hth = y(3)*H_SCALE
        f%hph = y(4)*H_SCALE
        f%Bmod = y(5)*B_SCALE

        f%dAth = dy(:, 1)*A_SCALE
        f%dAph = dy(:, 2)*A_SCALE
        f%dhth = dy(:, 3)*H_SCALE
        f%dhph = dy(:, 4)*H_SCALE
        f%dBmod = dy(:, 5)*B_SCALE
    end subroutine evaluate_spectre_can


    subroutine integ_to_ref_spectre(xinteg, xref)
        real(dp), intent(in) :: xinteg(3)
        real(dp), intent(out) :: xref(3)

        integer :: lvol
        real(dp) :: x(3), y(2)

        lvol = volume_index(xinteg(1))
        x = [clamp_to_volume(xinteg(1), spectre_volumes(lvol)), xinteg(2), xinteg(3)]
        call evaluate_batch_splines_3d(spectre_volumes(lvol)%spl_transform, x, y)

        xref(1) = xinteg(1)
        xref(2) = modulo(xinteg(2), twopi)
        xref(3) = modulo(xinteg(3) + y(1), twopi)
    end subroutine integ_to_ref_spectre


    subroutine ref_to_integ_spectre(xref, xinteg)
        real(dp), intent(in) :: xref(3)
        real(dp), intent(out) :: xinteg(3)

        real(dp), parameter :: TOL = 1d-12
        integer, parameter :: MAX_ITER = 16

        integer :: lvol, i
        real(dp) :: x(3), y(2), dy(3, 2), phi_prev

        lvol = volume_index(xref(1))
        xinteg(1) = xref(1)
        xinteg(2) = modulo(xref(2), twopi)
        xinteg(3) = modulo(xref(3), twopi)

        do i = 1, MAX_ITER
            x = [clamp_to_volume(xinteg(1), spectre_volumes(lvol)), xinteg(2), &
                 xinteg(3)]
            call evaluate_batch_splines_3d_der(spectre_volumes(lvol)%spl_transform, &
                                               x, y, dy)
            phi_prev = xinteg(3)
            xinteg(3) = phi_prev - (phi_prev + y(1) - xref(3))/(1d0 + dy(3, 1))
            if (abs(xinteg(3) - phi_prev) < TOL) return
        end do
        print *, 'WARNING: ref_to_integ_spectre did not converge after', &
            MAX_ITER, 'iterations'
    end subroutine ref_to_integ_spectre


    subroutine cleanup_spectre
        if (allocated(spectre_volumes)) deallocate (spectre_volumes)
        spectre_mvol = 0
    end subroutine cleanup_spectre

end module field_can_spectre
