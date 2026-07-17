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

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
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
    use diag_counters, only: count_event, EVT_SPECTRE_REF_INVERSE_MAXIT, &
        EVT_SPECTRE_INVALID_STATE

    implicit none
    private

    public :: construct_spectre_coordinates, evaluate_spectre_can, &
              integ_to_ref_spectre, ref_to_integ_spectre, cleanup_spectre, &
              spectre_volumes, spectre_mvol, set_spectre_volume_lock, &
              set_spectre_construction_grid

    !> Per-volume Meiss construction grid. phi = -1 selects the automatic policy
    !> below. set_spectre_construction_grid
    !> threads the namelist knob (params: spectre_ncon_r/th/phi) in before the build,
    !> kept as a setter because field_can_* cannot use params (dependency cycle).
    integer :: ncon_r = 48, ncon_th = 48, ncon_phi = -1, ncon_order = 5
    integer :: ncon_ode_step_limit = 1000000
    real(dp) :: ncon_ode_relerr = 1.0e-2_dp

    !> Automatic phi resolution. Axisymmetric SPECTRE
    !> fields (all n = 0) are phi-invariant, so the construction only needs a minimal
    !> periodic-quintic phi grid; AXISYM_NPHI is the dominant memory saver for tokamak
    !> cases (arrays and splines scale linearly in n_phi). It stays above the quintic
    !> order+1 = 6 so spl_five_per is well conditioned, and preserves the symplectic
    !> energy invariants; a positive spectre_ncon_phi overrides it (raise to NPHI_FULL
    !> for high-order convergence, lower for RK45-only or memory-bound runs).
    integer, parameter :: AXISYM_NPHI = 8, NPHI_FULL = 32
    real(dp), parameter :: AXIS_OFFSET = 1.0e-3_dp
    !> Sheet-adjacent regularization band, the canonical analogue of the RK45
    !> path's drift_band (spectre_orbit): the guiding-center 1-form degenerates
    !> (dpth/dr = 0, i.e. Bstar_par -> 0) in a thin unphysical layer next to
    !> each interface, where every state on the degeneracy manifold is a
    !> spurious root of the implicit steppers as h -> 0. Inside the band A and
    !> Bmod are extended linearly in r from the band edge while the radial slope
    !> of h_theta, h_phi is blended quadratically to zero over one band width
    !> (C1 everywhere, so the implicit-step residuals stay smooth); at and
    !> beyond the interface dpth/dr = sqrtg*Bmod/(ro0*h_phi) > 0 independent of
    !> pitch, so the canonical structure is non-degenerate by construction and
    !> all conservation properties survive. The interface jump itself still
    !> uses the true both-side |B| via magfie. The width is chosen so the
    !> band-edge h gradients are subcritical (ro0*vpar*curl h < B) for 3.5 MeV
    !> alphas on the test equilibria: a monotone blend of a subcritical slope
    !> stays subcritical, so no degeneracy re-enters the blend zone.
    real(dp), parameter :: EDGE_BAND = 5.0e-2_dp
    real(dp), parameter :: A_SCALE = TESLA_TO_GAUSS*M_TO_CM**2
    real(dp), parameter :: H_SCALE = M_TO_CM
    real(dp), parameter :: B_SCALE = TESLA_TO_GAUSS

    type(meiss_volume_t), allocatable :: spectre_volumes(:)
    integer :: spectre_mvol = 0

    !> Home-volume lock for multi-volume symplectic stepping (#441). Newton
    !> iterates of an implicit step may probe past an interior interface, where
    !> the int(r) dispatch would evaluate the neighbour volume's discontinuous
    !> field; with the lock set, evaluation stays on the home volume's splines,
    !> clamped at its edge (the canonical analogue of velo_can_clamped on the
    !> RK45 path). Each traced marker owns one thread.
    integer :: lock_lvol = 0
    !$omp threadprivate(lock_lvol)

contains

    elemental logical function ieee_is_finite(value)
        real(dp), intent(in) :: value

        integer(int64), parameter :: exponent_mask = &
            int(z'7FF0000000000000', int64)
        integer(int64) :: bits

        bits = transfer(value, bits)
        ieee_is_finite = iand(bits, exponent_mask) /= exponent_mask
    end function ieee_is_finite

    subroutine set_spectre_volume_lock(lvol)
        !> lvol > 0 pins field evaluation to that volume; 0 restores dispatch on
        !> int(r).
        integer, intent(in) :: lvol

        lock_lvol = lvol
    end subroutine set_spectre_volume_lock

    integer function active_volume(r) result(lvol)
        real(dp), intent(in) :: r

        if (lock_lvol > 0) then
            lvol = min(lock_lvol, spectre_mvol)
        else
            lvol = volume_index(r)
        end if
    end function active_volume

    subroutine set_spectre_construction_grid(n_r, n_th, n_phi, spline_order, &
                                             ode_step_limit, ode_relerr)
        !> Set the per-volume Meiss construction grid before the build (n_phi <= 0
        !> selects the automatic phi policy). Called from init_spectre_field.
        integer, intent(in) :: n_r, n_th, n_phi
        integer, intent(in), optional :: spline_order
        integer, intent(in), optional :: ode_step_limit
        real(dp), intent(in), optional :: ode_relerr

        ncon_r = n_r
        ncon_th = n_th
        ncon_phi = n_phi
        ncon_order = 5
        if (present(spline_order)) ncon_order = spline_order
        ncon_ode_step_limit = 1000000
        if (present(ode_step_limit)) ncon_ode_step_limit = ode_step_limit
        ncon_ode_relerr = 1.0e-2_dp
        if (present(ode_relerr)) ncon_ode_relerr = ode_relerr
    end subroutine set_spectre_construction_grid

    subroutine construct_spectre_coordinates(field_noncan, ierr)
        !> Build one Meiss chart per SPECTRE volume on r in [lvol-1, lvol]. The axis
        !> volume starts at AXIS_OFFSET so the construction never touches the r = 0
        !> coordinate singularity.
        class(magnetic_field_t), intent(in) :: field_noncan
        integer, intent(out), optional :: ierr

        integer :: lvol, nphi_eff, meiss_ierr
        real(dp) :: rmin, rmax
        logical :: axisym

        call cleanup_spectre

        select type (fld => field_noncan)
        type is (spectre_field_t)
            spectre_mvol = fld%data%Mvol
            axisym = all(fld%data%in == 0)
        class default
            error stop 'construct_spectre_coordinates: field is not spectre_field_t'
        end select

        if (ncon_phi > 0) then
            nphi_eff = ncon_phi
        else if (axisym) then
            nphi_eff = AXISYM_NPHI
        else
            nphi_eff = NPHI_FULL
        end if

        print '(A,I0,A,I0,A,L1,A,I0,A,L1,A,I0,A,I0,A,ES9.2)', &
            ' spectre_construction_grid: n_r=', ncon_r, ' n_th=', ncon_th, &
            ' axisym=', axisym, ' n_phi=', nphi_eff, &
            ' phi_auto=', ncon_phi <= 0, ' order=', ncon_order, &
            ' ode_steps=', ncon_ode_step_limit, ' ode_relerr=', ncon_ode_relerr

        allocate (spectre_volumes(spectre_mvol))

        do lvol = 1, spectre_mvol
            rmin = real(lvol - 1, dp)
            if (lvol == 1) rmin = rmin + AXIS_OFFSET
            rmax = real(lvol, dp)

            call init_meiss(field_noncan, ncon_r, ncon_th, nphi_eff, rmin, rmax, &
                            0.0_dp, twopi, identity_scaling_t(), &
                            transformation_relerr_=ncon_ode_relerr, &
                            spline_order_=ncon_order, &
                            ode_step_limit_=ncon_ode_step_limit)
            call get_meiss_coordinates(meiss_ierr)
            if (meiss_ierr /= 0) then
                call cleanup_spectre
                if (present(ierr)) then
                    ierr = lvol
                    return
                end if
                error stop 'SPECTRE canonical construction ODE failed'
            end if
            call harvest_meiss_volume(spectre_volumes(lvol))
        end do

        call cleanup_meiss
        if (present(ierr)) ierr = 0
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

        ! Angle-transform evaluation: the batch splines carry no valid
        ! polynomial extension past the construction domain, so pin the radius
        ! to the edge. Field evaluation instead goes through band_clamp and
        ! extend_band (see EDGE_BAND).
        rc = max(vol%rmin, min(vol%rmax, r))
    end function clamp_to_volume


    pure function band_clamp(r, vol) result(rc)
        !> Evaluation radius for the canonical field: pinned EDGE_BAND inside the
        !> volume; evaluate_spectre_can extends the values over the remaining
        !> distance (see EDGE_BAND and extend_band).
        real(dp), intent(in) :: r
        type(meiss_volume_t), intent(in) :: vol
        real(dp) :: rc

        rc = max(vol%rmin + EDGE_BAND, min(vol%rmax - EDGE_BAND, r))
    end function band_clamp


    pure subroutine extend_band(y, dy, d2y, dr_out, secders)
        !> Radial extension of the harvested field over the signed distance
        !> dr_out from the band edge (see EDGE_BAND). Per quantity the value
        !> follows y + y_r*G(dr_out): G = dr_out (linear, C-infinity) for Ath,
        !> Aph, Bmod; for hth, hph the slope G' blends quadratically from 1 to 0
        !> over one band width and stays 0 beyond, so the h-components go over
        !> smoothly (C1) into radially frozen values at the interface. All first
        !> derivatives are updated consistently with the extended values; only
        !> curvatures unavailable without third derivatives stay frozen, which
        !> degrades the Newton Jacobian alone.
        real(dp), intent(inout) :: y(5), dy(3, 5), d2y(6, 5)
        real(dp), intent(in) :: dr_out
        logical, intent(in) :: secders

        integer :: k
        real(dp) :: a, s, gval, gp, gpp, y1

        a = abs(dr_out)
        s = sign(1.0_dp, dr_out)

        do k = 1, 5
            if (k == 3 .or. k == 4) then
                if (a < EDGE_BAND) then
                    gval = dr_out - s*0.5_dp*a*a/EDGE_BAND
                    gp = 1.0_dp - a/EDGE_BAND
                    gpp = -s/EDGE_BAND
                else
                    gval = s*0.5_dp*EDGE_BAND
                    gp = 0.0_dp
                    gpp = 0.0_dp
                end if
            else
                gval = dr_out
                gp = 1.0_dp
                gpp = 0.0_dp
            end if

            y1 = dy(1, k)
            y(k) = y(k) + y1*gval
            if (secders) then
                dy(2, k) = dy(2, k) + d2y(2, k)*gval
                dy(3, k) = dy(3, k) + d2y(3, k)*gval
                d2y(1, k) = y1*gpp
                d2y(2, k) = d2y(2, k)*gp
                d2y(3, k) = d2y(3, k)*gp
            end if
            dy(1, k) = y1*gp
        end do
    end subroutine extend_band


    subroutine evaluate_spectre_can(f, r, th_c, ph_c, mode_secders)
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: r, th_c, ph_c
        integer, intent(in) :: mode_secders

        integer :: lvol, k
        real(dp) :: x(3), dr_out
        real(dp) :: y(5), dy(3, 5), d2y(6, 5)

        n_field_evaluations = n_field_evaluations + 1

        lvol = active_volume(r)
        x = [band_clamp(r, spectre_volumes(lvol)), th_c, ph_c]
        dr_out = r - x(1)

        if (mode_secders > 0) then
            call evaluate_batch_splines_3d_der2(spectre_volumes(lvol)%spl_field, &
                                                x, y, dy, d2y)
            if (dr_out /= 0.0_dp) then
                call extend_band(y, dy, d2y, dr_out, .true.)
            end if
            f%d2Ath = d2y(:, 1)*A_SCALE
            f%d2Aph = d2y(:, 2)*A_SCALE
            f%d2hth = d2y(:, 3)*H_SCALE
            f%d2hph = d2y(:, 4)*H_SCALE
            f%d2Bmod = d2y(:, 5)*B_SCALE
        else
            call evaluate_batch_splines_3d_der(spectre_volumes(lvol)%spl_field, &
                                               x, y, dy)
            if (dr_out /= 0.0_dp) then
                call extend_band(y, dy, d2y, dr_out, .false.)
            end if
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

        lvol = active_volume(xinteg(1))
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
        real(dp) :: x(3), y(2), dy(3, 2), phi, phi_next, target
        real(dp) :: residual, correction, best_phi, best_residual, denom

        if (.not. all(ieee_is_finite(xref))) then
            xinteg = xref
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if

        lvol = active_volume(xref(1))
        xinteg(1) = xref(1)
        xinteg(2) = modulo(xref(2), twopi)
        target = modulo(xref(3), twopi)
        phi = target
        best_phi = phi
        best_residual = huge(1.0_dp)

        do i = 1, MAX_ITER
            x = [clamp_to_volume(xinteg(1), spectre_volumes(lvol)), xinteg(2), &
                 modulo(phi, twopi)]
            call evaluate_batch_splines_3d_der(spectre_volumes(lvol)%spl_transform, &
                                               x, y, dy)
            if (.not. ieee_is_finite(y(1)) .or. &
                .not. ieee_is_finite(dy(3, 1))) exit
            residual = modulo(phi + y(1) - target + 0.5_dp*twopi, twopi) - &
                0.5_dp*twopi
            if (abs(residual) < best_residual) then
                best_residual = abs(residual)
                best_phi = phi
            end if
            denom = 1.0_dp + dy(3, 1)
            if (abs(denom) <= 100.0_dp*epsilon(1.0_dp)) exit
            correction = max(-0.5_dp*twopi, min(0.5_dp*twopi, residual/denom))
            phi_next = phi - correction
            if (.not. ieee_is_finite(phi_next)) exit
            phi = phi_next
            if (abs(correction) < TOL) then
                xinteg(3) = phi + twopi*anint((xref(3) - phi)/twopi)
                return
            end if
        end do
        xinteg(3) = best_phi + twopi*anint((xref(3) - best_phi)/twopi)
        call count_event(EVT_SPECTRE_REF_INVERSE_MAXIT)
    end subroutine ref_to_integ_spectre


    subroutine cleanup_spectre
        if (allocated(spectre_volumes)) deallocate (spectre_volumes)
        spectre_mvol = 0
    end subroutine cleanup_spectre

end module field_can_spectre
