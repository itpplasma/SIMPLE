program test_field_can_albert_diagnostic
    !> Diagnostic test to dump intermediate Albert coordinate values.
    !> Run with both main and current branch to compare.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: Tracer
    use simple_main, only: init_field
    use magfie_sub, only: ALBERT
    use velo_mod, only: isw_field_type
    use field, only: VmecField
    use field_can_albert, only: init_albert, psi_inner, psi_outer, &
        psi_of_x, Ath_norm, dpsi_dr_positive
    use field_can_meiss, only: spl_field_batch, xmin, xmax, n_r, n_th, n_phi
    use interpolate, only: evaluate_batch_splines_3d

    implicit none

    type(Tracer) :: norb
    class(VmecField), allocatable :: magfie
    integer :: i_r, i_th, i_phi, unit_id
    real(dp) :: x(3), y_batch(5)

    isw_field_type = ALBERT
    magfie = VmecField()

    print *, 'Initializing field with Albert coordinates...'
    call init_field(norb, 'wout.nc', 5, 5, 3, 0)

    ! Open diagnostic output file
    open(newunit=unit_id, file='albert_diagnostic.dat', status='replace')

    ! Write header with key parameters
    write(unit_id, '(A)') '# Albert coordinate diagnostic output'
    write(unit_id, '(A,3E20.12)') '# xmin:', xmin
    write(unit_id, '(A,3E20.12)') '# xmax:', xmax
    write(unit_id, '(A,3I8)') '# n_r, n_th, n_phi:', n_r, n_th, n_phi
    write(unit_id, '(A,E20.12)') '# Ath_norm:', Ath_norm
    write(unit_id, '(A,E20.12)') '# psi_inner:', psi_inner
    write(unit_id, '(A,E20.12)') '# psi_outer:', psi_outer
    write(unit_id, '(A,L2)') '# dpsi_dr_positive:', dpsi_dr_positive
    write(unit_id, *)

    ! Dump Meiss spline values at a subset of grid points
    write(unit_id, '(A)') '# Meiss spline values: i_r, i_th, i_phi, r, th, ph, Ath, Aph, hth, hph, Bmod'
    do i_phi = 1, n_phi, max(1, n_phi/8)
        do i_th = 1, n_th, max(1, n_th/8)
            do i_r = 1, n_r, max(1, n_r/8)
                x(1) = xmin(1) + (i_r-1)*(xmax(1)-xmin(1))/(n_r-1)
                x(2) = xmin(2) + (i_th-1)*(xmax(2)-xmin(2))/(n_th-1)
                x(3) = xmin(3) + (i_phi-1)*(xmax(3)-xmin(3))/(n_phi-1)
                call evaluate_batch_splines_3d(spl_field_batch, x, y_batch)
                write(unit_id, '(3I5,8E20.12)') i_r, i_th, i_phi, x(1), x(2), x(3), &
                    y_batch(1), y_batch(2), y_batch(3), y_batch(4), y_batch(5)
            end do
        end do
    end do

    ! Dump psi_of_x at corners
    write(unit_id, *)
    write(unit_id, '(A)') '# psi_of_x at corners:'
    write(unit_id, '(A,E20.12)') '# psi_of_x(1,1,1):', psi_of_x(1,1,1)
    write(unit_id, '(A,E20.12)') '# psi_of_x(n_r,1,1):', psi_of_x(n_r,1,1)
    write(unit_id, '(A,E20.12)') '# psi_of_x(1,n_th/2,n_phi/2):', psi_of_x(1,n_th/2,n_phi/2)
    write(unit_id, '(A,E20.12)') '# psi_of_x(n_r,n_th/2,n_phi/2):', psi_of_x(n_r,n_th/2,n_phi/2)

    close(unit_id)
    print *, 'Diagnostic output written to albert_diagnostic.dat'

end program test_field_can_albert_diagnostic
