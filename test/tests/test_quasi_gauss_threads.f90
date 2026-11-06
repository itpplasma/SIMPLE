program test_quasi_gauss_threads
    ! The quasi-Newton Gauss steppers must give every particle the same
    ! trajectory whether it is traced alone or next to others on concurrent
    ! threads. Each orbit is traced once on one thread as the reference and
    ! then again with all orbits shared among several threads; every step
    ! must agree bit for bit. s = 4 is the largest Gauss stage count.
!$  use omp_lib, only: omp_get_num_threads
    use orbit_symplectic, only: orbit_sympl_init
    use orbit_symplectic_base, only: symplectic_integrator_t, GAUSS2, GAUSS4
    use orbit_symplectic_quasi, only: f_quasi => f, si_quasi => si
    use orbit_quasi_implicit, only: timestep_rk_gauss_quasi
    use field_can_mod, only: field_can_t, eval_field => evaluate
    use simple, only: tracer_t, init_params
    use parmot_mod, only: ro0
    use simple_main, only: init_field
    use params, only: coord_input, field_input
    use new_vmec_stuff_mod, only: rmajor
    use util, only: twopi
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: NORB = 16, NT = 200, NTHREADS = 4
    type(tracer_t) :: norb_tracer
    real(dp) :: dt
    integer :: nfail, observed_threads

    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    call init_field(norb_tracer, 'wout.nc', 5, 5, 3, 2)
    call init_params(norb_tracer, 2, 4, 3.5d6, 64, 1, 1d-8)
    dt = twopi*rmajor*1.0d2/48

    nfail = 0
    observed_threads = 1
    call check_stages(2, GAUSS2)
    call check_stages(4, GAUSS4)
    if (observed_threads < 2) error stop 'parallel region ran on one thread'
    if (nfail > 0) error stop 'threaded quasi Gauss orbits differ from serial'
    print *, 'PASS'

contains

    subroutine check_stages(s, mode)
        integer, intent(in) :: s, mode
        type(symplectic_integrator_t) :: si0(NORB)
        type(field_can_t) :: f0(NORB)
        real(dp) :: serial(4, NT, NORB), threaded(4, NT, NORB)
        integer :: i, nsteps

        do i = 1, NORB
            call start_orbit(i, mode, si0(i), f0(i))
        end do

        do i = 1, NORB
            call trace(s, si0(i), f0(i), serial(:, :, i), nsteps)
            if (nsteps /= NT) then
                print *, 'Gauss s =', s, ' orbit', i, ' stopped after', nsteps
                nfail = nfail + 1
            end if
        end do

        !$omp parallel num_threads(NTHREADS) private(nsteps)
        !$omp single
!$      observed_threads = max(observed_threads, omp_get_num_threads())
        !$omp end single
        !$omp do schedule(static, 1)
        do i = 1, NORB
            call trace(s, si0(i), f0(i), threaded(:, :, i), nsteps)
        end do
        !$omp end do
        !$omp end parallel

        do i = 1, NORB
            if (any(threaded(:, :, i) /= serial(:, :, i))) then
                print *, 'Gauss s =', s, ' orbit', i, ' max |difference| =', &
                    maxval(abs(threaded(:, :, i) - serial(:, :, i)))
                nfail = nfail + 1
            end if
        end do
    end subroutine check_stages

    ! Passing orbits of different pitch and radius; set up serially because
    ! orbit_sympl_init writes module state.
    subroutine start_orbit(i, mode, si, f)
        integer, intent(in) :: i, mode
        type(symplectic_integrator_t), intent(out) :: si
        type(field_can_t), intent(out) :: f
        real(dp) :: z0(4), vpar0

        z0(1) = 0.3d0 + 0.02d0*i
        z0(2) = 0.7d0 + 0.3d0*i
        z0(3) = 0.1d0
        vpar0 = 0.6d0 + 0.02d0*i
        call eval_field(f, z0(1), z0(2), z0(3), 0)
        f%mu = 0.5d0**2*(1d0 - vpar0**2)/f%Bmod*2d0
        f%ro0 = ro0/sqrt(2d0)
        f%vpar = vpar0*sqrt(2d0)
        z0(4) = vpar0*f%hph + f%Aph/f%ro0
        call orbit_sympl_init(si, f, z0, dt, 1, 1d-12, mode)
    end subroutine start_orbit

    subroutine trace(s, si0, f0, out, nsteps)
        integer, intent(in) :: s
        type(symplectic_integrator_t), intent(in) :: si0
        type(field_can_t), intent(in) :: f0
        real(dp), intent(out) :: out(4, NT)
        integer, intent(out) :: nsteps
        integer :: ierr

        si_quasi = si0
        f_quasi = f0
        out = 0d0
        do nsteps = 0, NT - 1
            call timestep_rk_gauss_quasi(s, ierr)
            if (ierr /= 0) exit
            out(:, nsteps + 1) = si_quasi%z
        end do
    end subroutine trace

end program test_quasi_gauss_threads
