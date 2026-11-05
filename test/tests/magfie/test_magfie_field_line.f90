program test_magfie_field_line
    !> magfie_field_line must return the field strength and contravariant
    !> direction of the active magfie, for VMEC (its own single evaluation)
    !> and for any other field (delegation).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use magfie_sub, only: magfie, magfie_field_line, init_magfie, VMEC, TEST

    implicit none

    real(dp), parameter :: tol_rel = 1.0e-13_dp
    real(dp) :: dummy
    integer :: n_failed

    n_failed = 0
    call init_vmec('wout.nc', 5, 5, 5, dummy)

    call init_magfie(VMEC)
    call compare_on_grid('VMEC', twopi/real(nper, dp), n_failed)
    call init_magfie(TEST)
    call compare_on_grid('TEST', twopi, n_failed)

    if (n_failed > 0) then
        print '(i0,a)', n_failed, ' magfie_field_line mismatches'
        error stop 1
    end if
    print *, 'test_magfie_field_line: passed'

contains

    subroutine compare_on_grid(label, phi_period, n_failed)
        character(*), intent(in) :: label
        real(dp), intent(in) :: phi_period
        integer, intent(inout) :: n_failed
        real(dp) :: x(3), bmod, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)
        real(dp) :: bmod_line, hctr_line(3), scale
        integer :: i_r, i_th, i_ph

        do i_r = 1, 4
            do i_th = 1, 5
                do i_ph = 1, 3
                    x = [0.1_dp + 0.25_dp*real(i_r - 1, dp), &
                         0.3_dp + 1.1_dp*real(i_th - 1, dp), &
                         phi_period*(0.1_dp + 0.3_dp*real(i_ph - 1, dp))]
                    if (label == 'TEST') x(1) = 0.4_dp*x(1)
                    call magfie(x, bmod, sqrtg, bder, hcov, hctr, hcurl)
                    call magfie_field_line(x, bmod_line, hctr_line)
                    scale = maxval(abs(hctr))
                    if (abs(bmod_line - bmod) > tol_rel*abs(bmod) .or. &
                        maxval(abs(hctr_line - hctr)) > tol_rel*scale) then
                        print '(a,a,3es14.6)', label, ' mismatch at ', x
                        n_failed = n_failed + 1
                    end if
                end do
            end do
        end do
    end subroutine compare_on_grid

end program test_magfie_field_line
