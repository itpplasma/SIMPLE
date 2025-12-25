program test_field_can_albert_coils_diagnostic
    !> Diagnostic test to dump intermediate Albert coordinate values for coils field.
    !> Reads configuration from simple.in (like golden record tests).
    !> Run with both main and current branch to compare intermediate values.
    !>
    !> Usage: Run in a directory containing simple.in with:
    !>   isw_field_type = 4       (Albert coordinates)
    !>   field_input = 'coils.simple'
    !>   netcdffile = 'wout.nc'

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
        field_input
    use field_can_albert, only: psi_inner, psi_outer, psi_of_x, Ath_norm, &
        dpsi_dr_positive, r_of_xc, Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc
    use field_can_meiss, only: spl_field_batch, xmin, xmax, n_r, n_th, n_phi
    use interpolate, only: evaluate_batch_splines_3d

    implicit none

    type(tracer_t) :: norb
    integer :: i_r, i_th, i_phi, unit_id
    real(dp) :: x(3), x_spl(3), y_batch(5)
    integer :: mid_r, mid_th, mid_phi
    character(len=256) :: config_file

    config_file = 'simple.in'
    call read_config(config_file)

    print *, 'Initializing field with Albert coordinates...'
    print *, '  field_input = ', trim(field_input)
    print *, '  netcdffile = ', trim(netcdffile)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)

    open(newunit=unit_id, file='albert_coils_diagnostic.dat', status='replace')

    write(unit_id, '(A)') '# Albert coordinate diagnostic output (coils version)'
    write(unit_id, '(A,A)') '# field_input: ', trim(field_input)
    write(unit_id, '(A,3E24.16)') '# xmin:', xmin
    write(unit_id, '(A,3E24.16)') '# xmax:', xmax
    write(unit_id, '(A,3I8)') '# n_r, n_th, n_phi:', n_r, n_th, n_phi
    write(unit_id, '(A,E24.16)') '# Ath_norm:', Ath_norm
    write(unit_id, '(A,E24.16)') '# psi_inner:', psi_inner
    write(unit_id, '(A,E24.16)') '# psi_outer:', psi_outer
    write(unit_id, '(A,L2)') '# dpsi_dr_positive:', dpsi_dr_positive
    write(unit_id, *)

    write(unit_id, '(A)') '# psi_of_x at corners and center:'
    write(unit_id, '(A,E24.16)') '# psi_of_x(1,1,1):', psi_of_x(1,1,1)
    write(unit_id, '(A,E24.16)') '# psi_of_x(n_r,1,1):', psi_of_x(n_r,1,1)
    mid_th = n_th/2
    mid_phi = n_phi/2
    write(unit_id, '(A,E24.16)') '# psi_of_x(1,n_th/2,n_phi/2):', &
        psi_of_x(1,mid_th,mid_phi)
    write(unit_id, '(A,E24.16)') '# psi_of_x(n_r,n_th/2,n_phi/2):', &
        psi_of_x(n_r,mid_th,mid_phi)
    write(unit_id, *)

    mid_r = n_r/2
    write(unit_id, '(A)') '# Albert transformed values at center (n_r/2,n_th/2,n_phi/2):'
    write(unit_id, '(A,E24.16)') '# r_of_xc:', r_of_xc(mid_r,mid_th,mid_phi)
    write(unit_id, '(A,E24.16)') '# Aph_of_xc:', Aph_of_xc(mid_r,mid_th,mid_phi)
    write(unit_id, '(A,E24.16)') '# hth_of_xc:', hth_of_xc(mid_r,mid_th,mid_phi)
    write(unit_id, '(A,E24.16)') '# hph_of_xc:', hph_of_xc(mid_r,mid_th,mid_phi)
    write(unit_id, '(A,E24.16)') '# Bmod_of_xc:', Bmod_of_xc(mid_r,mid_th,mid_phi)
    write(unit_id, *)

    write(unit_id, '(A)') '# Meiss spline values: i_r, i_th, i_phi, ' // &
        'r, th, ph, Ath, Aph, hth, hph, Bmod'
    do i_phi = 1, n_phi, max(1, n_phi/8)
        do i_th = 1, n_th, max(1, n_th/8)
            do i_r = 1, n_r, max(1, n_r/8)
                x(1) = xmin(1) + (i_r-1)*(xmax(1)-xmin(1))/(n_r-1)
                x(2) = xmin(2) + (i_th-1)*(xmax(2)-xmin(2))/(n_th-1)
                x(3) = xmin(3) + (i_phi-1)*(xmax(3)-xmin(3))/(n_phi-1)
                ! Swap coordinates: physics [r, th, phi] -> spline [phi, th, r]
                x_spl(1) = x(3)
                x_spl(2) = x(2)
                x_spl(3) = x(1)
                call evaluate_batch_splines_3d(spl_field_batch, x_spl, y_batch)
                write(unit_id, '(3I5,8E24.16)') i_r, i_th, i_phi, &
                    x(1), x(2), x(3), &
                    y_batch(1), y_batch(2), y_batch(3), y_batch(4), y_batch(5)
            end do
        end do
    end do

    close(unit_id)
    print *, 'Diagnostic output written to albert_coils_diagnostic.dat'

end program test_field_can_albert_coils_diagnostic
