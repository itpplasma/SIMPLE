program test_albert_transform_diagnostic
    !> Diagnostic test to isolate error sources in Albert coordinate transforms.
    !>
    !> The Albert transform chain:
    !>   ref_to_integ: s -> r=sqrt(s) -> Ath(r,th,ph) -> psi = Ath/Ath_norm
    !>   integ_to_ref: psi -> r(psi,th,ph) via spline -> s = r^2
    !>
    !> Error sources:
    !>   1. Meiss transform (shown to be ~1e-16, negligible)
    !>   2. Ath spline interpolation (Meiss spl_field_batch)
    !>   3. Psi grid range restriction (safe side approach)
    !>   4. r(psi) inversion via Lagrange interpolation
    !>   5. r_of_xc spline interpolation

    use, intrinsic :: iso_fortran_env, only: dp => real64

    use field, only: vmec_field_t, create_vmec_field
    use simple, only: init_vmec
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, cleanup_meiss, &
        spl_field_batch, xmin, xmax, n_r, n_th, n_phi, twopi
    use field_can_albert, only: get_albert_coordinates, psi_inner, psi_outer, &
        psi_of_x, Ath_norm, r_of_xc, spl_r_batch, &
        integ_to_ref_albert, ref_to_integ_albert
    use interpolate, only: evaluate_batch_splines_3d

    implicit none

    real(dp) :: fper
    type(vmec_field_t) :: magfie

    print *, 'Initializing VMEC...'
    call init_vmec('wout.nc', 5, 5, 5, fper)
    call create_vmec_field(magfie)

    print *, ''
    print *, 'Initializing Meiss and Albert coordinates (32^3 grid)...'
    call init_meiss(magfie, 32, 32, 32, 0.01d0, 0.99d0, 0.0d0, twopi)
    call get_albert_coordinates()

    print *, ''
    call diagnose_psi_range()
    call diagnose_transform_steps()
    call diagnose_grid_resolution_effect()

    call cleanup_meiss()

contains

    subroutine diagnose_psi_range()
        !> Check the psi grid range restriction.
        real(dp) :: psi_min_inner, psi_max_inner, psi_min_outer, psi_max_outer
        real(dp) :: psi_range_full, psi_range_safe, coverage

        print *, '=== Psi Grid Range Analysis ==='

        ! Analyze psi variation at inner and outer radii
        psi_min_inner = minval(psi_of_x(1, :, :))
        psi_max_inner = maxval(psi_of_x(1, :, :))
        psi_min_outer = minval(psi_of_x(n_r, :, :))
        psi_max_outer = maxval(psi_of_x(n_r, :, :))

        print *, '  At r_min (inner radius):'
        print *, '    psi range:', psi_min_inner, 'to', psi_max_inner
        print *, '    variation:', psi_max_inner - psi_min_inner

        print *, '  At r_max (outer radius):'
        print *, '    psi range:', psi_min_outer, 'to', psi_max_outer
        print *, '    variation:', psi_max_outer - psi_min_outer

        psi_range_full = psi_max_outer - psi_min_inner
        psi_range_safe = psi_outer - psi_inner
        coverage = psi_range_safe / psi_range_full * 100d0

        print *, '  Safe psi range: [', psi_inner, ',', psi_outer, ']'
        print *, '  Full psi range: [', psi_min_inner, ',', psi_max_outer, ']'
        print *, '  Coverage:', coverage, '%'
        print *, ''
    end subroutine diagnose_psi_range


    subroutine diagnose_transform_steps()
        !> Trace through transform steps to identify error accumulation.
        use field_can_meiss, only: ref_to_integ_meiss, integ_to_ref_meiss

        real(dp) :: x_ref(3), x_meiss(3), x_meiss_back(3), x_ref_back(3)
        real(dp) :: x_albert(3)
        real(dp) :: psi_forward, r_inverse
        real(dp) :: y_ath(5), y_r(1)
        real(dp) :: err_meiss, err_ath_spline, err_r_spline, err_total

        print *, '=== Transform Step Analysis (at mid-radius) ==='

        ! Test point: s=0.5, theta=pi, phi=0.5
        x_ref = [0.5d0, 3.14159d0, 0.5d0]

        print *, '  Starting point (ref coords): s=', x_ref(1), &
                 ' th=', x_ref(2), ' ph=', x_ref(3)

        ! Step 1: ref -> meiss (should be exact)
        call ref_to_integ_meiss(x_ref, x_meiss)
        print *, '  Step 1 - ref_to_integ_meiss:'
        print *, '    r =', x_meiss(1), ' (expected', sqrt(x_ref(1)), ')'

        ! Step 2: Evaluate Ath spline at meiss coords
        call evaluate_batch_splines_3d(spl_field_batch, x_meiss, y_ath)
        psi_forward = y_ath(1) / Ath_norm
        print *, '  Step 2 - Ath spline evaluation:'
        print *, '    Ath =', y_ath(1), ' psi = Ath/Ath_norm =', psi_forward

        ! Step 3: Evaluate r spline at (psi, th, ph) - this is the inverse
        x_albert = [psi_forward, x_meiss(2), x_meiss(3)]
        call evaluate_batch_splines_3d(spl_r_batch, x_albert, y_r)
        r_inverse = y_r(1)
        print *, '  Step 3 - r_of_xc spline evaluation:'
        print *, '    r_inverse =', r_inverse, ' (expected', x_meiss(1), ')'
        err_r_spline = abs(r_inverse - x_meiss(1))
        print *, '    Error in r:', err_r_spline

        ! Step 4: meiss -> ref
        x_meiss_back = [r_inverse, x_meiss(2), x_meiss(3)]
        call integ_to_ref_meiss(x_meiss_back, x_ref_back)
        print *, '  Step 4 - integ_to_ref_meiss:'
        print *, '    s_back =', x_ref_back(1), ' (expected', x_ref(1), ')'

        ! Full roundtrip via Albert functions
        print *, ''
        print *, '  Full roundtrip via Albert functions:'
        call ref_to_integ_albert(x_ref, x_albert)
        call integ_to_ref_albert(x_albert, x_ref_back)
        err_total = abs(x_ref(1) - x_ref_back(1))
        print *, '    x_ref =', x_ref(1)
        print *, '    x_albert =', x_albert(1)
        print *, '    x_ref_back =', x_ref_back(1)
        print *, '    Total roundtrip error in s:', err_total
        print *, ''
    end subroutine diagnose_transform_steps


    subroutine diagnose_grid_resolution_effect()
        !> Test how grid resolution affects transform accuracy.
        real(dp) :: x_ref(3), x_albert(3), x_ref_back(3)
        real(dp) :: errors(5), s_values(5)
        integer :: i

        print *, '=== Error vs Radial Position ==='

        s_values = [0.1d0, 0.25d0, 0.5d0, 0.75d0, 0.9d0]

        do i = 1, 5
            x_ref = [s_values(i), 1.5d0, 0.5d0]
            call ref_to_integ_albert(x_ref, x_albert)
            call integ_to_ref_albert(x_albert, x_ref_back)
            errors(i) = abs(x_ref(1) - x_ref_back(1))
            print *, '  s =', s_values(i), ' -> error =', errors(i)
        end do

        print *, ''
        print *, '  Observation: Error tends to be larger near boundaries'
        print *, '  This is due to psi variation with (theta, phi) at fixed r'
        print *, ''
    end subroutine diagnose_grid_resolution_effect

end program test_albert_transform_diagnostic
