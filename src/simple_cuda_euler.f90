module simple_cuda_euler
    ! Bridge to the direct-CUDA symplectic Euler tracer (simple_cuda_euler.cu),
    ! the counterpart of simple_gpu trace_orbits_gpu_range on the Boozer
    ! quintic splines. The particle state crosses as structure-of-arrays.
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
    use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr, c_loc
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_mod, only: field_can_t
    use orbit_symplectic_base, only: symplectic_integrator_t, sympl_rmax
    use orbit_symplectic_base, only: sympl_predictor, sympl_stop_rule
    use orbit_symplectic_base, only: sympl_stop_kappa, symplectic_newton_warning_mode
    use boozer_sub, only: boozer_state, get_boozer_batch_splines
    use batch_interpolate_types, only: BatchSplineData1D, BatchSplineData3D

    implicit none
    private

    public :: trace_orbits_cuda_euler, cuda_euler_supported, cuda_euler_system_probe

    integer, parameter :: NSTATE = 25

    type, bind(c) :: field_c_t
        type(c_ptr) :: aphi_coeff
        integer(c_size_t) :: aphi_count
        integer(c_int) :: aphi_num_points, aphi_periodic
        real(c_double) :: aphi_x_min, aphi_h_step
        type(c_ptr) :: bcov_coeff
        integer(c_size_t) :: bcov_count
        integer(c_int) :: bcov_num_points, bcov_periodic
        real(c_double) :: bcov_x_min, bcov_h_step
        type(c_ptr) :: field_coeff
        integer(c_size_t) :: field_count
        integer(c_int) :: field_num_points(3), field_periodic(3)
        real(c_double) :: field_x_min(3), field_h_step(3), field_inv_h_step(3)
        real(c_double) :: field_period(3)
        real(c_double) :: torflux
        integer(c_int) :: nper
    end type field_c_t

    type, bind(c) :: controls_c_t
        real(c_double) :: rmax, stop_kappa
        integer(c_int) :: predictor, stop_rule, warning_mode
    end type controls_c_t

    interface
        function simple_cuda_euler_trace(n, field, controls, state, istate, nt, &
                                         ntau, lstep, ltime, nfev, profile_ms) &
            bind(c, name='simple_cuda_euler_trace') result(ierr)
            import :: c_int, c_double, c_int64_t, field_c_t, controls_c_t
            integer(c_int), value :: n, nt
            type(field_c_t), intent(in) :: field
            type(controls_c_t), intent(in) :: controls
            real(c_double), intent(inout) :: state(*)
            integer(c_int), intent(inout) :: istate(*)
            integer(c_int), intent(in) :: ntau(*)
            integer(c_int), intent(out) :: lstep(*)
            real(c_double), intent(out) :: ltime(*)
            integer(c_int64_t), intent(out) :: nfev(*)
            real(c_double), intent(out) :: profile_ms(5)
            integer(c_int) :: ierr
        end function simple_cuda_euler_trace

        function simple_cuda_euler_system_probe(n, field, input, output) &
            bind(c, name='simple_cuda_euler_system_probe') result(ierr)
            import :: c_int, c_double, field_c_t
            integer(c_int), value :: n
            type(field_c_t), intent(in) :: field
            real(c_double), intent(in) :: input(*)
            real(c_double), intent(out) :: output(*)
            integer(c_int) :: ierr
        end function simple_cuda_euler_system_probe
    end interface

contains

    !> True when the active Boozer field is the single-quantity quintic
    !> representation the generated CUDA kernels evaluate.
    logical function cuda_euler_supported()
        type(BatchSplineData1D), pointer :: aphi, bcov
        type(BatchSplineData3D), pointer :: field3d
        logical :: ready

        call get_boozer_batch_splines(aphi, bcov, field3d, ready)
        cuda_euler_supported = ready
        if (.not. ready) return
        cuda_euler_supported = aphi%order == 5 .and. aphi%num_quantities == 1
        cuda_euler_supported = cuda_euler_supported .and. bcov%order == 5 .and. &
                               bcov%num_quantities == 2
        cuda_euler_supported = cuda_euler_supported .and. all(field3d%order == 5) &
                               .and. field3d%num_quantities == 1
    end function cuda_euler_supported

    subroutine trace_orbits_cuda_euler(si, f, npart, ntimstep, ntau_macro, &
                                       loss_step, loss_time, zend, nfev, profile_ms)
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(inout) :: f(npart)
        integer, intent(in) :: npart, ntimstep
        integer, intent(in) :: ntau_macro(ntimstep)
        integer, intent(out) :: loss_step(npart), nfev(npart)
        real(dp), intent(out) :: loss_time(npart), zend(4, npart)
        real(dp), intent(out), optional :: profile_ms(5)

        type(field_c_t) :: field
        type(controls_c_t) :: controls
        real(c_double), allocatable :: state(:, :)
        integer(c_int), allocatable :: istate(:, :)
        integer(c_int64_t), allocatable :: nfev64(:)
        integer(c_int), allocatable :: ntau(:), lstep(:)
        real(c_double) :: profile(5)
        integer(c_int) :: n, nt
        integer :: i, ierr

        call fill_field(field)

        controls%rmax = sympl_rmax
        controls%stop_kappa = sympl_stop_kappa
        controls%predictor = sympl_predictor
        controls%stop_rule = sympl_stop_rule
        controls%warning_mode = merge(1, 0, symplectic_newton_warning_mode)

        allocate (state(npart, NSTATE), istate(npart, 2), nfev64(npart), &
                  lstep(npart))
        do i = 1, npart
            state(i, 1:4) = si(i)%z
            state(i, 5) = si(i)%pthold
            state(i, 6) = f(i)%pth
            state(i, 7) = si(i)%dt
            state(i, 8) = si(i)%atol
            state(i, 9) = si(i)%rtol
            state(i, 10) = si(i)%warning_factor
            state(i, 11) = f(i)%mu
            state(i, 12) = f(i)%ro0
            state(i, 13) = si(i)%hist_dt
            state(i, 14:17) = si(i)%hist_z
            state(i, 18:21) = si(i)%hist_dz(:, 1)
            state(i, 22:25) = si(i)%hist_dz(:, 2)
            istate(i, 1) = si(i)%ntau
            istate(i, 2) = si(i)%nhist
        end do
        ntau = ntau_macro

        n = int(npart, c_int)
        nt = int(ntimstep, c_int)
        ierr = simple_cuda_euler_trace(n, field, controls, state, istate, nt, ntau, &
                                       lstep, loss_time, nfev64, profile)
        if (ierr /= 0) then
            print *, 'simple_cuda_euler_trace: CUDA error ', ierr
            error stop 'trace_orbits_cuda_euler failed'
        end if

        do i = 1, npart
            si(i)%z = state(i, 1:4)
            si(i)%pthold = state(i, 5)
            f(i)%pth = state(i, 6)
            si(i)%hist_dt = state(i, 13)
            si(i)%hist_z = state(i, 14:17)
            si(i)%hist_dz(:, 1) = state(i, 18:21)
            si(i)%hist_dz(:, 2) = state(i, 22:25)
            si(i)%nhist = istate(i, 2)
            zend(:, i) = si(i)%z
            loss_step(i) = lstep(i)
            nfev(i) = int(nfev64(i))
        end do
        if (present(profile_ms)) profile_ms = profile
    end subroutine trace_orbits_cuda_euler

    !> Generated Euler system and its auxiliary field quantities at n points;
    !> input(9, n) and output(21, n) as documented in simple_cuda_euler.h.
    subroutine cuda_euler_system_probe(n, input, output)
        integer, intent(in) :: n
        real(dp), intent(in) :: input(9, n)
        real(dp), intent(out) :: output(21, n)

        type(field_c_t) :: field
        integer :: ierr

        call fill_field(field)
        ierr = simple_cuda_euler_system_probe(int(n, c_int), field, input, output)
        if (ierr /= 0) then
            print *, 'simple_cuda_euler_system_probe: CUDA error ', ierr
            error stop 'cuda_euler_system_probe failed'
        end if
    end subroutine cuda_euler_system_probe

    subroutine fill_field(field)
        type(field_c_t), intent(out) :: field
        type(BatchSplineData1D), pointer :: aphi, bcov
        type(BatchSplineData3D), pointer :: field3d
        logical :: ready

        call get_boozer_batch_splines(aphi, bcov, field3d, ready)
        if (.not. ready) error stop 'trace_orbits_cuda_euler: Boozer splines not built'
        if (.not. cuda_euler_supported()) error stop &
            'trace_orbits_cuda_euler: needs quintic Boozer splines with one |B|'

        field%aphi_coeff = c_loc(aphi%coeff)
        field%aphi_count = size(aphi%coeff, kind=c_size_t)
        field%aphi_num_points = aphi%num_points
        field%aphi_periodic = merge(1, 0, aphi%periodic)
        field%aphi_x_min = aphi%x_min
        field%aphi_h_step = aphi%h_step
        field%bcov_coeff = c_loc(bcov%coeff)
        field%bcov_count = size(bcov%coeff, kind=c_size_t)
        field%bcov_num_points = bcov%num_points
        field%bcov_periodic = merge(1, 0, bcov%periodic)
        field%bcov_x_min = bcov%x_min
        field%bcov_h_step = bcov%h_step
        field%field_coeff = c_loc(field3d%coeff)
        field%field_count = size(field3d%coeff, kind=c_size_t)
        field%field_num_points = field3d%num_points
        field%field_periodic = merge(1, 0, field3d%periodic)
        field%field_x_min = field3d%x_min
        field%field_h_step = field3d%h_step
        field%field_inv_h_step = field3d%inv_h_step
        field%field_period = field3d%period
        field%torflux = boozer_state%torflux
        field%nper = boozer_state%nper

    end subroutine fill_field

end module simple_cuda_euler
