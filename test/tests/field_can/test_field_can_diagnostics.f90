program test_field_can_diagnostics
    !> Comprehensive diagnostic test for ALL coordinate system modes.
    !> Dumps intermediate values at every stage for golden record comparison.
    !>
    !> Supported coordinate systems:
    !>   - CANFLUX (1): Canonical flux coordinates
    !>   - BOOZER (2): Boozer coordinates
    !>   - MEISS (3): Meiss canonical coordinates
    !>   - ALBERT (4): Albert canonical coordinates
    !>
    !> Output files:
    !>   - field_can_diagnostic_<mode>.dat: Detailed intermediate values
    !>
    !> Usage: Run in directory with simple.in and wout.nc (and coils.simple for ALBERT)

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
        integ_coords, field_input
    use field_can_mod, only: field_can_t, evaluate, integ_to_ref, ref_to_integ, &
        init_field_can
    use magfie_sub, only: CANFLUX, BOOZER, MEISS, ALBERT

    implicit none

    type(tracer_t) :: norb
    integer :: mode, unit_id
    character(len=32) :: mode_name
    character(len=64) :: output_file
    logical :: success
    character(len=256) :: config_file

    config_file = 'simple.in'
    call read_config(config_file)

    mode = integ_coords
    call get_mode_name(mode, mode_name)

    write(output_file, '(A,I1,A)') 'field_can_diagnostic_', mode, '.dat'

    print *, '=============================================='
    print *, 'Field Can Diagnostics for mode: ', trim(mode_name)
    print *, '  integ_coords = ', mode
    print *, '  field_input = ', trim(field_input)
    print *, '  netcdffile = ', trim(netcdffile)
    print *, '=============================================='

    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)

    open(newunit=unit_id, file=trim(output_file), status='replace')

    write(unit_id, '(A)') '# Field Can Diagnostic Output'
    write(unit_id, '(A,A)') '# Mode: ', trim(mode_name)
    write(unit_id, '(A,I2)') '# integ_coords: ', mode
    write(unit_id, '(A,A)') '# field_input: ', trim(field_input)
    write(unit_id, '(A,A)') '# netcdffile: ', trim(netcdffile)
    write(unit_id, '(A)')

    select case (mode)
    case (CANFLUX)
        call dump_canflux_diagnostics(unit_id)
    case (BOOZER)
        call dump_boozer_diagnostics(unit_id)
    case (MEISS)
        call dump_meiss_diagnostics(unit_id)
    case (ALBERT)
        call dump_albert_diagnostics(unit_id)
    case default
        write(unit_id, '(A,I2)') '# ERROR: Unknown mode ', mode
    end select

    call dump_field_evaluation_diagnostics(unit_id, mode)
    call dump_coordinate_transform_diagnostics(unit_id, mode)

    close(unit_id)

    print *, 'Diagnostic output written to: ', trim(output_file)
    print *, 'SUCCESS: Diagnostics completed'

contains

    subroutine get_mode_name(mode, name)
        integer, intent(in) :: mode
        character(len=*), intent(out) :: name

        select case (mode)
        case (CANFLUX)
            name = 'CANFLUX'
        case (BOOZER)
            name = 'BOOZER'
        case (MEISS)
            name = 'MEISS'
        case (ALBERT)
            name = 'ALBERT'
        case default
            name = 'UNKNOWN'
        end select
    end subroutine get_mode_name


    subroutine dump_canflux_diagnostics(unit_id)
        use field_can_flux, only: twopi
        integer, intent(in) :: unit_id

        write(unit_id, '(A)') '=== CANFLUX Specific Diagnostics ==='
        write(unit_id, '(A,E25.17)') 'twopi = ', twopi
        write(unit_id, '(A)') '# CANFLUX uses VMEC coordinates directly'
        write(unit_id, '(A)') '# No intermediate grid transformation'
        write(unit_id, '(A)')
    end subroutine dump_canflux_diagnostics


    subroutine dump_boozer_diagnostics(unit_id)
        use field_can_boozer, only: twopi
        integer, intent(in) :: unit_id

        write(unit_id, '(A)') '=== BOOZER Specific Diagnostics ==='
        write(unit_id, '(A,E25.17)') 'twopi = ', twopi
        write(unit_id, '(A)') '# BOOZER transforms to Boozer coordinates'
        write(unit_id, '(A)')
    end subroutine dump_boozer_diagnostics


    subroutine dump_meiss_diagnostics(unit_id)
        use field_can_meiss, only: xmin, xmax, n_r, n_th, n_phi, order, &
            spl_field_batch
        use interpolate, only: evaluate_batch_splines_3d
        integer, intent(in) :: unit_id
        real(dp) :: x(3), x_spl(3), y_batch(5)
        integer :: i_r, i_th, i_phi

        write(unit_id, '(A)') '=== MEISS Specific Diagnostics ==='
        write(unit_id, '(A,I6)') 'n_r = ', n_r
        write(unit_id, '(A,I6)') 'n_th = ', n_th
        write(unit_id, '(A,I6)') 'n_phi = ', n_phi
        write(unit_id, '(A,I6)') 'order = ', order
        write(unit_id, '(A,3E25.17)') 'xmin = ', xmin
        write(unit_id, '(A,3E25.17)') 'xmax = ', xmax
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# Meiss spline values: i_r, i_th, i_phi, r, th, ph, Ath, Aph, hth, hph, Bmod'
        do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
            do i_th = 1, n_th, max(1, (n_th-1)/4)
                do i_r = 1, n_r, max(1, (n_r-1)/4)
                    x(1) = xmin(1) + (i_r - 1) * (xmax(1) - xmin(1)) / (n_r - 1)
                    x(2) = xmin(2) + (i_th - 1) * (xmax(2) - xmin(2)) / (n_th - 1)
                    x(3) = xmin(3) + (i_phi - 1) * (xmax(3) - xmin(3)) / (n_phi - 1)
                    ! Swap coordinates: physics [r, th, phi] -> spline [phi, th, r]
                    x_spl(1) = x(3)
                    x_spl(2) = x(2)
                    x_spl(3) = x(1)
                    call evaluate_batch_splines_3d(spl_field_batch, x_spl, y_batch)
                    write(unit_id, '(3I6,8E25.17)') i_r, i_th, i_phi, &
                        x(1), x(2), x(3), &
                        y_batch(1), y_batch(2), y_batch(3), y_batch(4), y_batch(5)
                end do
            end do
        end do
        write(unit_id, '(A)')
    end subroutine dump_meiss_diagnostics


    subroutine dump_albert_diagnostics(unit_id)
        use field_can_meiss, only: xmin, xmax, n_r, n_th, n_phi, order, &
            spl_field_batch
        use field_can_albert, only: psi_inner, psi_outer, psi_of_x, Ath_norm, &
            dpsi_dr_positive, r_of_xc, Aph_of_xc, hth_of_xc, hph_of_xc, &
            Bmod_of_xc, spl_albert_batch
        use interpolate, only: evaluate_batch_splines_3d
        integer, intent(in) :: unit_id
        real(dp) :: x(3), x_spl(3), y_batch_meiss(5), y_batch_albert(4), psi, h_psi
        integer :: i_r, i_th, i_phi

        write(unit_id, '(A)') '=== ALBERT Specific Diagnostics ==='

        write(unit_id, '(A)') '# Grid parameters (from Meiss base)'
        write(unit_id, '(A,I6)') 'n_r = ', n_r
        write(unit_id, '(A,I6)') 'n_th = ', n_th
        write(unit_id, '(A,I6)') 'n_phi = ', n_phi
        write(unit_id, '(A,I6)') 'order = ', order
        write(unit_id, '(A,3E25.17)') 'xmin = ', xmin
        write(unit_id, '(A,3E25.17)') 'xmax = ', xmax
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# Normalization and psi transformation'
        write(unit_id, '(A,E25.17)') 'Ath_norm = ', Ath_norm
        write(unit_id, '(A,E25.17)') 'psi_inner = ', psi_inner
        write(unit_id, '(A,E25.17)') 'psi_outer = ', psi_outer
        write(unit_id, '(A,L2)') 'dpsi_dr_positive = ', dpsi_dr_positive
        h_psi = (psi_outer - psi_inner) / dble(n_r - 1)
        write(unit_id, '(A,E25.17)') 'h_psi = ', h_psi
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# psi_of_x extremes'
        write(unit_id, '(A,E25.17)') 'maxval(psi_of_x(1,:,:)) = ', &
            maxval(psi_of_x(1,:,:))
        write(unit_id, '(A,E25.17)') 'minval(psi_of_x(1,:,:)) = ', &
            minval(psi_of_x(1,:,:))
        write(unit_id, '(A,E25.17)') 'maxval(psi_of_x(n_r,:,:)) = ', &
            maxval(psi_of_x(n_r,:,:))
        write(unit_id, '(A,E25.17)') 'minval(psi_of_x(n_r,:,:)) = ', &
            minval(psi_of_x(n_r,:,:))
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# Meiss spline samples: i_r, i_th, i_phi, r, th, ph, Ath, Aph, hth, hph, Bmod'
        do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
            do i_th = 1, n_th, max(1, (n_th-1)/4)
                do i_r = 1, n_r, max(1, (n_r-1)/4)
                    x(1) = xmin(1) + (i_r - 1) * (xmax(1) - xmin(1)) / (n_r - 1)
                    x(2) = xmin(2) + (i_th - 1) * (xmax(2) - xmin(2)) / (n_th - 1)
                    x(3) = xmin(3) + (i_phi - 1) * (xmax(3) - xmin(3)) / (n_phi - 1)
                    ! Swap coordinates: physics [r, th, phi] -> spline [phi, th, r]
                    x_spl(1) = x(3)
                    x_spl(2) = x(2)
                    x_spl(3) = x(1)
                    call evaluate_batch_splines_3d(spl_field_batch, x_spl, y_batch_meiss)
                    write(unit_id, '(3I6,8E25.17)') i_r, i_th, i_phi, &
                        x(1), x(2), x(3), y_batch_meiss(1), y_batch_meiss(2), &
                        y_batch_meiss(3), y_batch_meiss(4), y_batch_meiss(5)
                end do
            end do
        end do
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# psi_of_x grid samples: i_r, i_th, i_phi, psi'
        do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
            do i_th = 1, n_th, max(1, (n_th-1)/4)
                do i_r = 1, n_r, max(1, (n_r-1)/4)
                    write(unit_id, '(3I6,E25.17)') i_r, i_th, i_phi, &
                        psi_of_x(i_r, i_th, i_phi)
                end do
            end do
        end do
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# r_of_xc samples: i_psi, i_th, i_phi, r'
        do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
            do i_th = 1, n_th, max(1, (n_th-1)/4)
                do i_r = 1, n_r, max(1, (n_r-1)/4)
                    write(unit_id, '(3I6,E25.17)') i_r, i_th, i_phi, &
                        r_of_xc(i_r, i_th, i_phi)
                end do
            end do
        end do
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# Field components after psi transform: i_psi, i_th, i_phi, Aph, hth, hph, Bmod'
        do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
            do i_th = 1, n_th, max(1, (n_th-1)/4)
                do i_r = 1, n_r, max(1, (n_r-1)/4)
                    write(unit_id, '(3I6,4E25.17)') i_r, i_th, i_phi, &
                        Aph_of_xc(i_r, i_th, i_phi), hth_of_xc(i_r, i_th, i_phi), &
                        hph_of_xc(i_r, i_th, i_phi), Bmod_of_xc(i_r, i_th, i_phi)
                end do
            end do
        end do
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# Albert spline evaluation: psi, th, ph, Aph, hth, hph, Bmod'
        do i_phi = 1, 3
            do i_th = 1, 3
                do i_r = 1, 5
                    psi = psi_inner + (i_r - 1) * (psi_outer - psi_inner) / 4
                    x = [psi, &
                         xmin(2) + (i_th - 1) * (xmax(2) - xmin(2)) / 2, &
                         xmin(3) + (i_phi - 1) * (xmax(3) - xmin(3)) / 2]
                    ! Swap coordinates: physics [psi, th, phi] -> spline [phi, th, psi]
                    x_spl(1) = x(3)
                    x_spl(2) = x(2)
                    x_spl(3) = x(1)
                    call evaluate_batch_splines_3d(spl_albert_batch, x_spl, y_batch_albert)
                    write(unit_id, '(7E25.17)') x(1), x(2), x(3), &
                        y_batch_albert(1), y_batch_albert(2), &
                        y_batch_albert(3), y_batch_albert(4)
                end do
            end do
        end do
        write(unit_id, '(A)')

        write(unit_id, '(A)') '# Summary statistics'
        write(unit_id, '(A,E25.17)') 'sum(psi_of_x) = ', sum(psi_of_x)
        write(unit_id, '(A,E25.17)') 'sum(r_of_xc) = ', sum(r_of_xc)
        write(unit_id, '(A,E25.17)') 'sum(Aph_of_xc) = ', sum(Aph_of_xc)
        write(unit_id, '(A,E25.17)') 'sum(hth_of_xc) = ', sum(hth_of_xc)
        write(unit_id, '(A,E25.17)') 'sum(hph_of_xc) = ', sum(hph_of_xc)
        write(unit_id, '(A,E25.17)') 'sum(Bmod_of_xc) = ', sum(Bmod_of_xc)
        write(unit_id, '(A,E25.17)') 'minval(Bmod_of_xc) = ', minval(Bmod_of_xc)
        write(unit_id, '(A,E25.17)') 'maxval(Bmod_of_xc) = ', maxval(Bmod_of_xc)
        write(unit_id, '(A)')
    end subroutine dump_albert_diagnostics


    subroutine dump_field_evaluation_diagnostics(unit_id, mode)
        integer, intent(in) :: unit_id, mode
        type(field_can_t) :: f
        real(dp) :: test_points(3, 5)
        integer :: i

        write(unit_id, '(A)') '=== Field Evaluation at Test Points ==='
        write(unit_id, '(A)') '# x1, x2, x3, Ath, Aph, hth, hph, Bmod'

        test_points(:, 1) = [0.3_dp, 0.5_dp, 0.2_dp]
        test_points(:, 2) = [0.5_dp, 1.5_dp, 0.5_dp]
        test_points(:, 3) = [0.7_dp, 3.0_dp, 1.0_dp]
        test_points(:, 4) = [0.4_dp, 4.5_dp, 0.8_dp]
        test_points(:, 5) = [0.6_dp, 5.5_dp, 1.5_dp]

        do i = 1, 5
            call evaluate(f, test_points(1, i), test_points(2, i), &
                test_points(3, i), 0)
            write(unit_id, '(8E25.17)') test_points(1, i), test_points(2, i), &
                test_points(3, i), f%Ath, f%Aph, f%hth, f%hph, f%Bmod
        end do
        write(unit_id, '(A)')
    end subroutine dump_field_evaluation_diagnostics


    subroutine dump_coordinate_transform_diagnostics(unit_id, mode)
        integer, intent(in) :: unit_id, mode
        real(dp) :: xref(3), xinteg(3), xref_back(3)
        real(dp) :: test_ref(3, 5)
        integer :: i

        write(unit_id, '(A)') '=== Coordinate Transform Roundtrip ==='
        write(unit_id, '(A)') '# xref -> xinteg -> xref_back'
        write(unit_id, '(A)') '# ref1, ref2, ref3, integ1, integ2, integ3, back1, back2, back3'

        test_ref(:, 1) = [0.3_dp, 0.5_dp, 0.2_dp]
        test_ref(:, 2) = [0.5_dp, 1.5_dp, 0.5_dp]
        test_ref(:, 3) = [0.7_dp, 3.0_dp, 1.0_dp]
        test_ref(:, 4) = [0.4_dp, 4.5_dp, 0.8_dp]
        test_ref(:, 5) = [0.6_dp, 5.5_dp, 1.5_dp]

        do i = 1, 5
            xref = test_ref(:, i)
            call ref_to_integ(xref, xinteg)
            call integ_to_ref(xinteg, xref_back)
            write(unit_id, '(9E25.17)') xref(1), xref(2), xref(3), &
                xinteg(1), xinteg(2), xinteg(3), &
                xref_back(1), xref_back(2), xref_back(3)
        end do
        write(unit_id, '(A)')
    end subroutine dump_coordinate_transform_diagnostics

end program test_field_can_diagnostics
