program test_canflux_threads
    use, intrinsic :: iso_fortran_env, only: dp => real64
!$  use omp_lib, only: omp_set_dynamic, omp_set_num_threads, omp_get_num_threads
    use simple, only: init_vmec
    use field, only: vmec_field_t, create_vmec_field
    use get_can_sub, only: get_canonical_coordinates_with_field, can_to_vmec
    use field_can_flux, only: eval_field_can
    use field_can_base, only: field_can_t

    implicit none

    integer, parameter :: npoints = 256, nthreads = 4
    type(vmec_field_t) :: vmec_field
    real(dp) :: serial(7, npoints), threaded(7, npoints), rebuilt(7, npoints)
    real(dp) :: field_period
    integer :: i, observed_threads

!$  call omp_set_dynamic(.false.)
!$  call omp_set_num_threads(1)
    call init_vmec('wout.nc', 5, 5, 3, field_period)
    call create_vmec_field(vmec_field)
    call get_canonical_coordinates_with_field(vmec_field)
    do i = 1, npoints
        call sample(i, serial(:, i))
    end do

    observed_threads = 1
    !$omp parallel num_threads(nthreads)
    !$omp single
!$  observed_threads = omp_get_num_threads()
    !$omp end single
    !$omp do schedule(static, 1)
    do i = 1, npoints
        call sample(i, threaded(:, i))
    end do
    !$omp end do
    !$omp end parallel
    if (observed_threads /= nthreads) error stop 'expected four OpenMP threads'
    if (any(threaded /= serial)) then
        print *, 'conversion difference:', maxval(abs(threaded - serial))
        error stop 'threaded CANFLUX conversions differ from serial'
    end if

!$  call omp_set_num_threads(nthreads)
    call get_canonical_coordinates_with_field(vmec_field)
    do i = 1, npoints
        call sample(i, rebuilt(:, i))
    end do
    if (any(rebuilt /= serial)) then
        print *, 'construction difference:', maxval(abs(rebuilt - serial))
        error stop 'threaded CANFLUX construction differs from serial'
    end if
    print *, 'PASS'

contains

    subroutine sample(index, values)
        integer, intent(in) :: index
        real(dp), intent(out) :: values(7)
        type(field_can_t) :: f
        real(dp) :: s, theta, phi

        s = 0.05_dp + 0.9_dp*real(modulo(37*index, npoints), dp)/npoints
        theta = 6.0_dp*real(modulo(71*index, npoints), dp)/npoints
        phi = field_period*real(modulo(97*index, npoints), dp)/npoints
        call can_to_vmec(s, theta, phi, values(1), values(2))
        call eval_field_can(f, s, theta, phi, 0)
        values(3:) = [f%Bmod, f%Ath, f%Aph, f%hth, f%hph]
    end subroutine sample

end program test_canflux_threads
