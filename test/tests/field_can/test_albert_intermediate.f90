program test_albert_intermediate
    !> Comprehensive intermediate value check for the Albert coordinate pipeline.
    !> Reads configuration from simple.in (like golden record tests).
    !> Dumps values at EVERY stage to pinpoint any discrepancy sources.
    !>
    !> Pipeline stages:
    !>   1. Grid parameters (n_r, n_th, n_phi, xmin, xmax)
    !>   2. Meiss field initialization (spl_field_batch values)
    !>   3. psi_of_x grid computation (from Ath component)
    !>   4. Ath_norm normalization factor
    !>   5. psi_inner, psi_outer boundary computation
    !>   6. r_of_xc inverse mapping
    !>   7. Field components after psi transform (Aph, hth, hph, Bmod)_of_xc
    !>   8. Final Albert spline values (spl_albert_batch)

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, field_input
    use field_can_meiss, only: spl_field_batch, xmin, xmax, n_r, n_th, n_phi, order
    use field_can_albert, only: psi_inner, psi_outer, psi_of_x, Ath_norm, &
        dpsi_dr_positive, r_of_xc, Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc, &
        spl_albert_batch
    use interpolate, only: evaluate_batch_splines_3d

    implicit none

    type(tracer_t) :: norb
    integer :: i_r, i_th, i_phi, unit_id
    real(dp) :: x(3), y_batch_meiss(5), y_batch_albert(4)
    real(dp) :: r, th, ph, psi
    real(dp) :: h_psi
    character(len=256) :: config_file

    config_file = 'simple.in'
    call read_config(config_file)

    print *, 'Initializing field with Albert coordinates...'
    print *, '  field_input = ', trim(field_input)
    print *, '  netcdffile = ', trim(netcdffile)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)

    open(newunit=unit_id, file='albert_intermediate.dat', status='replace')

    ! STAGE 0: Grid parameters
    write(unit_id, '(A)') '=== STAGE 0: Grid parameters ==='
    write(unit_id, '(A,I6)') 'n_r = ', n_r
    write(unit_id, '(A,I6)') 'n_th = ', n_th
    write(unit_id, '(A,I6)') 'n_phi = ', n_phi
    write(unit_id, '(A,I6)') 'order = ', order
    write(unit_id, '(A,3E25.17)') 'xmin = ', xmin
    write(unit_id, '(A,3E25.17)') 'xmax = ', xmax
    write(unit_id, *)

    ! STAGE 1: Normalization constants computed in init_psi_grid
    write(unit_id, '(A)') '=== STAGE 1: Normalization constants ==='
    write(unit_id, '(A,E25.17)') 'Ath_norm = ', Ath_norm
    write(unit_id, '(A,E25.17)') 'psi_inner = ', psi_inner
    write(unit_id, '(A,E25.17)') 'psi_outer = ', psi_outer
    write(unit_id, '(A,L2)') 'dpsi_dr_positive = ', dpsi_dr_positive
    h_psi = (psi_outer - psi_inner) / dble(n_r - 1)
    write(unit_id, '(A,E25.17)') 'h_psi = ', h_psi
    write(unit_id, *)

    ! STAGE 2: Meiss spline values at grid points
    write(unit_id, '(A)') '=== STAGE 2: Meiss spl_field_batch at grid points ==='
    write(unit_id, '(A)') '# Format: i_r, i_th, i_phi, r, th, ph, Ath, Aph, hth,'&
        //' hph, Bmod'
    ! Sample every 1/4 of grid
    do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
        ph = xmin(3) + (i_phi - 1) * (xmax(3) - xmin(3)) / (n_phi - 1)
        do i_th = 1, n_th, max(1, (n_th-1)/4)
            th = xmin(2) + (i_th - 1) * (xmax(2) - xmin(2)) / (n_th - 1)
            do i_r = 1, n_r, max(1, (n_r-1)/4)
                r = xmin(1) + (i_r - 1) * (xmax(1) - xmin(1)) / (n_r - 1)
                x = [r, th, ph]
                call evaluate_batch_splines_3d(spl_field_batch, x, y_batch_meiss)
                write(unit_id, '(3I6,8E25.17)') i_r, i_th, i_phi, r, th, ph, &
                    y_batch_meiss(1), y_batch_meiss(2), y_batch_meiss(3), &
                    y_batch_meiss(4), y_batch_meiss(5)
            end do
        end do
    end do
    write(unit_id, *)

    ! STAGE 3: psi_of_x grid (raw, before normalization is applied for storage)
    ! Note: psi_of_x is already normalized by Ath_norm in init_psi_grid
    write(unit_id, '(A)') '=== STAGE 3: psi_of_x grid samples ==='
    write(unit_id, '(A)') '# Format: i_r, i_th, i_phi, psi'
    do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
        do i_th = 1, n_th, max(1, (n_th-1)/4)
            do i_r = 1, n_r, max(1, (n_r-1)/4)
                write(unit_id, '(3I6,E25.17)') i_r, i_th, i_phi, &
                    psi_of_x(i_r, i_th, i_phi)
            end do
        end do
    end do
    write(unit_id, *)

    ! STAGE 4: psi_of_x extremes used for boundary computation
    write(unit_id, '(A)') '=== STAGE 4: psi_of_x extremes for boundaries ==='
    write(unit_id, '(A,E25.17)') 'maxval(psi_of_x(1,:,:)) = ', &
        maxval(psi_of_x(1,:,:))
    write(unit_id, '(A,E25.17)') 'minval(psi_of_x(1,:,:)) = ', &
        minval(psi_of_x(1,:,:))
    write(unit_id, '(A,E25.17)') 'maxval(psi_of_x(n_r,:,:)) = ', &
        maxval(psi_of_x(n_r,:,:))
    write(unit_id, '(A,E25.17)') 'minval(psi_of_x(n_r,:,:)) = ', &
        minval(psi_of_x(n_r,:,:))
    write(unit_id, '(A,E25.17)') 'psi_of_x(n_r,n_th/2,n_phi/2) = ', &
        psi_of_x(n_r, n_th/2, n_phi/2)
    write(unit_id, '(A,E25.17)') 'psi_of_x(1,n_th/2,n_phi/2) = ', &
        psi_of_x(1, n_th/2, n_phi/2)
    write(unit_id, *)

    ! STAGE 5: r_of_xc transformation result
    write(unit_id, '(A)') '=== STAGE 5: r_of_xc samples ==='
    write(unit_id, '(A)') '# Format: i_psi, i_th, i_phi, r_of_xc'
    do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
        do i_th = 1, n_th, max(1, (n_th-1)/4)
            do i_r = 1, n_r, max(1, (n_r-1)/4)
                write(unit_id, '(3I6,E25.17)') i_r, i_th, i_phi, &
                    r_of_xc(i_r, i_th, i_phi)
            end do
        end do
    end do
    write(unit_id, *)

    ! STAGE 6: Field components after psi transform
    write(unit_id, '(A)') '=== STAGE 6: Field components after grid_r_to_psi ==='
    write(unit_id, '(A)') '# Format: i_psi, i_th, i_phi, Aph, hth, hph, Bmod'
    do i_phi = 1, n_phi, max(1, (n_phi-1)/4)
        do i_th = 1, n_th, max(1, (n_th-1)/4)
            do i_r = 1, n_r, max(1, (n_r-1)/4)
                write(unit_id, '(3I6,4E25.17)') i_r, i_th, i_phi, &
                    Aph_of_xc(i_r, i_th, i_phi), hth_of_xc(i_r, i_th, i_phi), &
                    hph_of_xc(i_r, i_th, i_phi), Bmod_of_xc(i_r, i_th, i_phi)
            end do
        end do
    end do
    write(unit_id, *)

    ! STAGE 7: Albert spline evaluation at canonical coordinate test points
    write(unit_id, '(A)') '=== STAGE 7: spl_albert_batch evaluation ==='
    write(unit_id, '(A)') '# Format: psi, th_c, ph_c, Aph, hth, hph, Bmod'
    do i_phi = 1, 3
        ph = xmin(3) + (i_phi - 1) * (xmax(3) - xmin(3)) / 2
        do i_th = 1, 3
            th = xmin(2) + (i_th - 1) * (xmax(2) - xmin(2)) / 2
            do i_r = 1, 5
                psi = psi_inner + (i_r - 1) * (psi_outer - psi_inner) / 4
                x = [psi, th, ph]
                call evaluate_batch_splines_3d(spl_albert_batch, x, y_batch_albert)
                write(unit_id, '(7E25.17)') psi, th, ph, &
                    y_batch_albert(1), y_batch_albert(2), &
                    y_batch_albert(3), y_batch_albert(4)
            end do
        end do
    end do
    write(unit_id, *)

    ! STAGE 8: Statistics
    write(unit_id, '(A)') '=== STAGE 8: Summary statistics ==='
    write(unit_id, '(A,E25.17)') 'sum(psi_of_x) = ', sum(psi_of_x)
    write(unit_id, '(A,E25.17)') 'sum(r_of_xc) = ', sum(r_of_xc)
    write(unit_id, '(A,E25.17)') 'sum(Aph_of_xc) = ', sum(Aph_of_xc)
    write(unit_id, '(A,E25.17)') 'sum(hth_of_xc) = ', sum(hth_of_xc)
    write(unit_id, '(A,E25.17)') 'sum(hph_of_xc) = ', sum(hph_of_xc)
    write(unit_id, '(A,E25.17)') 'sum(Bmod_of_xc) = ', sum(Bmod_of_xc)
    write(unit_id, '(A,E25.17)') 'minval(Bmod_of_xc) = ', minval(Bmod_of_xc)
    write(unit_id, '(A,E25.17)') 'maxval(Bmod_of_xc) = ', maxval(Bmod_of_xc)

    close(unit_id)
    print *, 'Intermediate values written to albert_intermediate.dat'
    print *, 'SUCCESS: All stages dumped for comparison'

end program test_albert_intermediate
