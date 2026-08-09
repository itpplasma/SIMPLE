program test_boozer_rk_symplectic_oracle
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap, only: load_boozer_from_chartmap
    use boozer_rk_tables, only: rk_tables_ready, &
        splint_boozer_symplectic_table_device
    use field_can_base, only: field_can_t
    use field_can_boozer, only: eval_field_booz

    implicit none

    integer, parameter :: npoint = 256
    real(dp), parameter :: raw_limit = 5.0e-2_dp
    real(dp), parameter :: hessian_limit = 2.0e-1_dp
    character(1024) :: filename
    type(field_can_t) :: f
    real(dp) :: s, theta, phi
    real(dp) :: aphi, daphi, d2aphi
    real(dp) :: btheta, dbtheta, d2btheta
    real(dp) :: bphi, dbphi, d2bphi, bmod
    real(dp) :: dbmod(3), d2bmod(6)
    real(dp) :: max_raw, max_hessian
    real(dp) :: raw_max(8), hessian_max(6)
    integer :: i, status, failures

    call get_command_argument(1, filename, status=status)
    if (status /= 0 .or. len_trim(filename) == 0) &
        error stop 'usage: test_boozer_rk_symplectic_oracle.x chartmap.nc'

    call load_boozer_from_chartmap(trim(filename))
    if (.not. rk_tables_ready) &
        error stop 'chartmap does not contain native RK tables'

    max_raw = 0.0_dp
    max_hessian = 0.0_dp
    raw_max = 0.0_dp
    hessian_max = 0.0_dp
    failures = 0
    do i = 1, npoint
        s = 0.015_dp + 0.96_dp*real(mod(37*i, 251), dp)/251.0_dp
        theta = 2.0_dp*acos(-1.0_dp)*real(mod(19*i, 257), dp)/257.0_dp
        phi = 2.0_dp*acos(-1.0_dp)*real(mod(43*i, 263), dp)/263.0_dp

        call eval_field_booz(f, s, theta, phi, 3)
        call splint_boozer_symplectic_table_device(s, theta, phi, 3, &
            aphi, daphi, d2aphi, btheta, dbtheta, d2btheta, &
            bphi, dbphi, d2bphi, bmod, dbmod, d2bmod)

        call compare_value(f%Aph, aphi, raw_max(1), failures)
        call compare_value(f%dAph(1), daphi, raw_max(2), failures)
        call compare_value(f%d2Aph(1), d2aphi, hessian_max(1), failures)
        call compare_value(f%hth, btheta/f%Bmod, raw_max(3), failures)
        call compare_value(f%hph, bphi/f%Bmod, raw_max(4), failures)
        call compare_value(f%Bmod, bmod, raw_max(5), failures)
        call compare_vector(f%dBmod, dbmod, raw_max(6:8), failures)
        call compare_vector(f%d2Bmod(1:3), d2bmod(1:3), hessian_max(2:4), failures)
        call compare_value(f%d2hth(1), &
            d2btheta/f%Bmod - 2.0_dp*dbtheta*dbmod(1)/f%Bmod**2 + &
            btheta/f%Bmod**2*(2.0_dp*dbmod(1)**2/f%Bmod - d2bmod(1)), &
            hessian_max(5), failures)
        call compare_value(f%d2hph(1), &
            d2bphi/f%Bmod - 2.0_dp*dbphi*dbmod(1)/f%Bmod**2 + &
            bphi/f%Bmod**2*(2.0_dp*dbmod(1)**2/f%Bmod - d2bmod(1)), &
            hessian_max(6), failures)
    end do

    max_raw = maxval(raw_max)
    max_hessian = maxval(hessian_max)

    write (*, '(a,8(es12.4,1x))') 'raw by quantity (Aphi,dAphi,hth,hph,Bmod,dB_s,dB_th,dB_ph) = ', raw_max
    write (*, '(a,6(es12.4,1x))') 'Hessian by quantity (d2Aphi,d2B_ss,d2B_st,d2B_sp,d2hth,d2hph) = ', hessian_max
    write (*, '(a,es12.4)') 'maximum raw relative difference = ', max_raw
    write (*, '(a,es12.4)') 'maximum Hessian relative difference = ', max_hessian
    if (failures > 0 .or. max_raw > raw_limit .or. max_hessian > hessian_limit) &
        error stop 'compact symplectic table differs from generic Boozer contract'
    print *, 'test_boozer_rk_symplectic_oracle: PASSED'

contains

    subroutine compare_value(left, right, maximum, failures)
        real(dp), intent(in) :: left, right
        real(dp), intent(inout) :: maximum
        integer, intent(inout) :: failures
        real(dp) :: scale, error

        if (.not. ieee_is_finite(left) .or. .not. ieee_is_finite(right)) then
            failures = failures + 1
            return
        end if
        scale = max(1.0_dp, abs(left), abs(right))
        error = abs(left - right)/scale
        maximum = max(maximum, error)
    end subroutine compare_value

    subroutine compare_vector(left, right, maximum, failures)
        real(dp), intent(in) :: left(:), right(:)
        real(dp), intent(inout) :: maximum(:)
        integer, intent(inout) :: failures
        integer :: j

        do j = 1, size(left)
            call compare_value(left(j), right(j), maximum(j), failures)
        end do
    end subroutine compare_vector

end program test_boozer_rk_symplectic_oracle
