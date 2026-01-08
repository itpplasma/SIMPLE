program test_chartmap_wall_losses
    !> Test chartmap wall boundary loss computations.
    !>
    !> Compares loss fractions between 5 coordinate systems:
    !>   Mode 1: VMEC coordinates with RK45 integration
    !>   Mode 2: Chartmap reference coordinates with RK45 integration
    !>   Mode 3: Boozer canonical coordinates
    !>   Mode 4: Meiss canonical from VMEC field
    !>   Mode 5: Meiss canonical from chartmap field
    !>
    !> Coordinate conventions:
    !>   - VMEC:    z(1) = s (normalized toroidal flux), boundary at s=1
    !>   - Chartmap: z(1) = rho = sqrt(s), boundary at rho=1
    !>   - Boozer:  z(1) = s (same as VMEC), only angles transformed
    !>   - Meiss:   z(1) = r with sqrt_s_scaling, so r = sqrt(s_ref)
    !>              For Meiss-VMEC: s_ref = s, boundary at r=1 means s=1
    !>              For Meiss-Chart: s_ref = rho, boundary at r=1 means rho=1
    !>
    !> Outputs NetCDF file for visualization of loss positions and times.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_max_threads
    use netcdf
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use magfie_sub, only: init_magfie, magfie, VMEC, BOOZER, MEISS, REFCOORDS, &
                          set_magfie_refcoords_field
    use alpha_lifetime_sub, only: orbit_timestep_axis
    use parmot_mod, only: ro0, rmu
    use orbit_symplectic_base, only: symplectic_integrator_t, MIDPOINT
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl
    use field_can_mod, only: field_can_t, eval_field => evaluate
    use field_vmec, only: vmec_field_t, create_vmec_field
    use field_splined, only: splined_field_t, create_splined_field
    use field_can_mod, only: init_field_can, integ_to_ref, ref_to_integ, &
                             field_can_from_id
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, cleanup_meiss
    use boozer_sub, only: vmec_to_boozer, boozer_to_vmec, get_boozer_coordinates
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use libneo_coordinates, only: coordinate_system_t, vmec_coordinate_system_t, &
                                  make_chartmap_coordinate_system

    implicit none

    character(len=*), parameter :: wout_file = 'wout.nc'
    character(len=1000) :: chartmap_file, chartmap_offset_file
    character(len=1000) :: output_nc

    type(vmec_field_t) :: vmec_field
    type(splined_field_t) :: splined_vmec, splined_chartmap
    class(coordinate_system_t), allocatable :: chartmap_coords

    integer, parameter :: n_particles = 64
    integer, parameter :: n_steps = 2000
    real(dp), parameter :: trace_time = 1.0e-3_dp
    real(dp), parameter :: relerr = 1.0e-10_dp
    real(dp), parameter :: s_start = 0.8_dp

    real(dp) :: dtau, tau, v0_alpha, fper, zeta_period
    real(dp) :: z0(5, n_particles)

    real(dp), allocatable :: loss_times_vmec(:), loss_times_chart(:), &
                             loss_times_boozer(:)
    real(dp), allocatable :: loss_times_meiss_vmec(:), loss_times_meiss_chart(:)
    real(dp), allocatable :: loss_pos_vmec(:, :), loss_pos_chart(:, :), &
                             loss_pos_boozer(:, :)
    real(dp), allocatable :: loss_pos_meiss_vmec(:, :), loss_pos_meiss_chart(:, :)
    integer, allocatable :: lost_step_vmec(:), lost_step_chart(:), lost_step_boozer(:)
    integer, allocatable :: lost_step_meiss_vmec(:), lost_step_meiss_chart(:)

    real(dp) :: confined_vmec, confined_chart, confined_boozer
    real(dp) :: confined_meiss_vmec, confined_meiss_chart
    real(dp) :: confined_boozer_mid, confined_meiss_vmec_mid, confined_meiss_chart_mid
    real(dp) :: t_vmec, t_chart, t_boozer, t_meiss_vmec, t_meiss_chart
    real(dp) :: t_boozer_mid, t_meiss_vmec_mid, t_meiss_chart_mid
    real(dp) :: t_start, t_end
    real(dp), allocatable :: loss_times_boozer_mid(:), loss_times_meiss_vmec_mid(:)
    real(dp), allocatable :: loss_times_meiss_chart_mid(:)
    real(dp), allocatable :: loss_pos_boozer_mid(:, :), loss_pos_meiss_vmec_mid(:, :)
    real(dp), allocatable :: loss_pos_meiss_chart_mid(:, :)
    integer, allocatable :: lost_step_boozer_mid(:), lost_step_meiss_vmec_mid(:)
    integer, allocatable :: lost_step_meiss_chart_mid(:)
    integer :: i, ierr, n_failed

    n_failed = 0

    print *, '========================================================'
    print *, 'Chartmap Wall Losses Test'
    print *, '========================================================'
    print *, 'Comparing loss fractions: VMEC vs Chartmap vs Meiss'
    print *, '========================================================'
    print '(A,I0)', ' OpenMP max threads: ', omp_get_max_threads()
    print *

    chartmap_file = 'wout.chartmap.nc'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, chartmap_file)
    end if

    output_nc = 'chartmap_wall_losses.nc'
    if (command_argument_count() >= 2) then
        call get_command_argument(2, output_nc)
    end if

    allocate (loss_times_vmec(n_particles), loss_times_chart(n_particles))
    allocate (loss_times_boozer(n_particles))
    allocate (loss_times_meiss_vmec(n_particles), loss_times_meiss_chart(n_particles))
    allocate (loss_pos_vmec(3, n_particles), loss_pos_chart(3, n_particles))
    allocate (loss_pos_boozer(3, n_particles))
    allocate (loss_pos_meiss_vmec(3, n_particles), loss_pos_meiss_chart(3, n_particles))
    allocate (lost_step_vmec(n_particles), lost_step_chart(n_particles))
    allocate (lost_step_boozer(n_particles))
    allocate (lost_step_meiss_vmec(n_particles), lost_step_meiss_chart(n_particles))
    allocate (loss_times_boozer_mid(n_particles))
    allocate (loss_times_meiss_vmec_mid(n_particles), &
              loss_times_meiss_chart_mid(n_particles))
    allocate (loss_pos_boozer_mid(3, n_particles))
    allocate (loss_pos_meiss_vmec_mid(3, n_particles))
    allocate (loss_pos_meiss_chart_mid(3, n_particles))
    allocate (lost_step_boozer_mid(n_particles))
    allocate (lost_step_meiss_vmec_mid(n_particles), &
              lost_step_meiss_chart_mid(n_particles))

    print *, 'Step 1: Initialize VMEC equilibrium'
    call init_vmec(wout_file, 5, 5, 5, fper)
    call init_reference_coordinates(wout_file)
    call set_physics_parameters(v0_alpha)
    zeta_period = twopi/real(nper, dp)

    tau = trace_time*v0_alpha
    dtau = tau/n_steps
    print '(A,ES12.4)', '  Trace time (s): ', trace_time
    print '(A,I0)', '  Number of particles: ', n_particles
    print '(A,I0)', '  Time steps: ', n_steps

    print *, 'Step 2: Load chartmap'
    print *, '  Using: ', trim(chartmap_file)

    print *, 'Step 3: Create field objects'
    call create_vmec_field(vmec_field)
    call create_splined_field(vmec_field, ref_coords, splined_vmec)
    call make_chartmap_coordinate_system(chartmap_coords, chartmap_file)
    call create_splined_field(vmec_field, chartmap_coords, splined_chartmap)

    print *, 'Step 4: Initialize particles'
    call init_particles(z0)

    print *, '========================================================'
    print *, 'Running loss computations...'
    print *, '========================================================'

    print *, 'Mode 1: VMEC coordinates (RK45)'
    call init_magfie(VMEC)
    call cpu_time(t_start)
    call trace_particles_vmec(z0, loss_times_vmec, loss_pos_vmec, lost_step_vmec)
    call cpu_time(t_end)
    t_vmec = t_end - t_start
    confined_vmec = count_confined(loss_times_vmec)
    print '(A,F6.2,A,F8.3,A)', '  Confined fraction: ', confined_vmec*100.0_dp, &
        '%  (', t_vmec, ' s)'

    print *, 'Mode 2: Chartmap reference coordinates (RK45)'
    call set_magfie_refcoords_field(splined_chartmap)
    call init_magfie(REFCOORDS)
    call cpu_time(t_start)
    call trace_particles_chartmap(z0, loss_times_chart, loss_pos_chart, lost_step_chart)
    call cpu_time(t_end)
    t_chart = t_end - t_start
    confined_chart = count_confined(loss_times_chart)
    print '(A,F6.2,A,F8.3,A)', '  Confined fraction: ', confined_chart*100.0_dp, &
        '%  (', t_chart, ' s)'

    print *, 'Mode 3: Boozer canonical coordinates (RK45)'
    call get_boozer_coordinates
    call init_magfie(BOOZER)
    call cpu_time(t_start)
    call trace_particles_boozer(z0, loss_times_boozer, loss_pos_boozer, &
                                lost_step_boozer)
    call cpu_time(t_end)
    t_boozer = t_end - t_start
    confined_boozer = count_confined(loss_times_boozer)
    print '(A,F6.2,A,F8.3,A)', '  Confined fraction: ', confined_boozer*100.0_dp, &
        '%  (', t_boozer, ' s)'

    print *, 'Mode 4: Meiss canonical from VMEC (RK45)'
    call init_meiss(vmec_field, n_r_=16, n_th_=17, n_phi_=18, &
                    rmin=0.1d0, rmax=0.99d0)
    call get_meiss_coordinates
    call field_can_from_id(MEISS)  ! Set up ref_to_integ/integ_to_ref pointers
    call init_magfie(MEISS)
    call cpu_time(t_start)
    call trace_particles_meiss_vmec(z0, loss_times_meiss_vmec, loss_pos_meiss_vmec, &
                                    lost_step_meiss_vmec)
    call cpu_time(t_end)
    t_meiss_vmec = t_end - t_start
    call cleanup_meiss
    confined_meiss_vmec = count_confined(loss_times_meiss_vmec)
    print '(A,F6.2,A,F8.3,A)', '  Confined fraction: ', confined_meiss_vmec*100.0_dp, &
        '%  (', t_meiss_vmec, ' s)'

    print *, 'Mode 5: Meiss canonical from chartmap (RK45)'
    call init_meiss(splined_chartmap, n_r_=16, n_th_=17, n_phi_=18, &
                    rmin=0.1d0, rmax=0.99d0)
    call get_meiss_coordinates
    call field_can_from_id(MEISS)  ! Set up ref_to_integ/integ_to_ref pointers
    call init_magfie(MEISS)
    call cpu_time(t_start)
    call trace_particles_meiss(z0, loss_times_meiss_chart, loss_pos_meiss_chart, &
                               lost_step_meiss_chart)
    call cpu_time(t_end)
    t_meiss_chart = t_end - t_start
    call cleanup_meiss
    confined_meiss_chart = count_confined(loss_times_meiss_chart)
    print '(A,F6.2,A,F8.3,A)', '  Confined fraction: ', confined_meiss_chart*100.0_dp, &
        '%  (', t_meiss_chart, ' s)'

    print *, 'Mode 6: Boozer canonical coordinates (MIDPOINT)'
    call get_boozer_coordinates
    call field_can_from_id(BOOZER)
    call cpu_time(t_start)
    call trace_particles_boozer_sympl(z0, loss_times_boozer_mid, loss_pos_boozer_mid, &
                                      lost_step_boozer_mid)
    call cpu_time(t_end)
    t_boozer_mid = t_end - t_start
    confined_boozer_mid = count_confined(loss_times_boozer_mid)
    print '(A,F6.2,A,F8.3,A)', '  Confined fraction: ', confined_boozer_mid*100.0_dp, &
        '%  (', t_boozer_mid, ' s)'

    print *, 'Mode 7: Meiss canonical from VMEC (MIDPOINT)'
    call init_meiss(vmec_field, n_r_=16, n_th_=17, n_phi_=18, &
                    rmin=0.1d0, rmax=0.99d0)
    call get_meiss_coordinates
    call field_can_from_id(MEISS)
    call cpu_time(t_start)
    call trace_particles_meiss_vmec_sympl(z0, loss_times_meiss_vmec_mid, &
                                          loss_pos_meiss_vmec_mid, &
                                          lost_step_meiss_vmec_mid)
    call cpu_time(t_end)
    t_meiss_vmec_mid = t_end - t_start
    call cleanup_meiss
    confined_meiss_vmec_mid = count_confined(loss_times_meiss_vmec_mid)
    print '(A,F6.2,A,F8.3,A)', '  Confined fraction: ', &
        confined_meiss_vmec_mid*100.0_dp, '%  (', t_meiss_vmec_mid, ' s)'

    print *, 'Mode 8: Meiss canonical from chartmap (MIDPOINT)'
    call init_meiss(splined_chartmap, n_r_=16, n_th_=17, n_phi_=18, &
                    rmin=0.1d0, rmax=0.99d0)
    call get_meiss_coordinates
    call field_can_from_id(MEISS)
    call cpu_time(t_start)
    call trace_particles_meiss_chart_sympl(z0, loss_times_meiss_chart_mid, &
                                           loss_pos_meiss_chart_mid, &
                                           lost_step_meiss_chart_mid)
    call cpu_time(t_end)
    t_meiss_chart_mid = t_end - t_start
    call cleanup_meiss
    confined_meiss_chart_mid = count_confined(loss_times_meiss_chart_mid)
    print '(A,F6.2,A,F8.3,A)', '  Confined fraction: ', &
        confined_meiss_chart_mid*100.0_dp, '%  (', t_meiss_chart_mid, ' s)'

    print *
    print *, '========================================================'
    print *, 'Results Summary'
    print *, '========================================================'
    print '(A)', '  RK45 integrator:'
    print '(A,F6.2,A)', '    VMEC confined:         ', confined_vmec*100.0_dp, '%'
    print '(A,F6.2,A)', '    Chartmap confined:     ', confined_chart*100.0_dp, '%'
    print '(A,F6.2,A)', '    Boozer confined:       ', confined_boozer*100.0_dp, '%'
    print '(A,F6.2,A)', '    Meiss(VMEC) confined:  ', confined_meiss_vmec*100.0_dp, '%'
    print '(A,F6.2,A)', '    Meiss(Chart) confined: ', &
        confined_meiss_chart*100.0_dp, '%'
    print '(A)', '  Midpoint symplectic:'
    print '(A,F6.2,A)', '    Boozer confined:       ', confined_boozer_mid*100.0_dp, '%'
    print '(A,F6.2,A)', '    Meiss(VMEC) confined:  ', &
        confined_meiss_vmec_mid*100.0_dp, '%'
    print '(A,F6.2,A)', '    Meiss(Chart) confined: ', &
        confined_meiss_chart_mid*100.0_dp, '%'

    print *, 'Step 5: Write output for visualization'
    call write_output_netcdf

    if (abs(confined_vmec - confined_chart) > 0.15_dp) then
        print *, 'WARNING: VMEC vs Chartmap confined fraction differs by > 15%'
        n_failed = n_failed + 1
    else if (abs(confined_vmec - confined_chart) > 0.05_dp) then
        print *, 'NOTE: VMEC vs Chartmap differs by 5-15% (boundary detection differs)'
    end if

    if (n_failed > 0) then
        print *, '========================================================'
        print *, 'FAILED: ', n_failed, ' comparison(s) exceeded tolerance'
        print *, '========================================================'
        error stop 1
    end if

    print *, '========================================================'
    print *, 'PASSED: Loss fractions consistent across coordinate systems'
    print *, '========================================================'

contains

    subroutine set_physics_parameters(v0_out)
        use util, only: ev, p_mass, c_light => c, e_charge

        real(dp), intent(out) :: v0_out
        real(dp) :: E_alpha, rlarm
        integer :: n_d, n_e

        n_d = 4
        n_e = 2
        E_alpha = 3.5d6

        v0_out = sqrt(2.0_dp*E_alpha*ev/(n_d*p_mass))
        rlarm = v0_out*n_d*p_mass*c_light/(n_e*e_charge)
        ro0 = rlarm
        rmu = 1.0d8

        print '(A,ES12.4)', '  Alpha energy (eV): ', E_alpha
        print '(A,ES12.4)', '  v0 (cm/s): ', v0_out
    end subroutine set_physics_parameters

    subroutine init_particles(z)
        real(dp), intent(out) :: z(5, n_particles)
        integer :: i
        real(dp) :: theta, phi, lambda

        do i = 1, n_particles
            z(1, i) = s_start
            theta = twopi*real(mod(i - 1, 8), dp)/8.0_dp
            phi = zeta_period*real((i - 1)/8, dp)/8.0_dp
            z(2, i) = theta
            z(3, i) = phi
            z(4, i) = 1.0_dp
            lambda = -1.0_dp + 2.0_dp*real(mod(i - 1, 4), dp)/3.0_dp
            z(5, i) = lambda
        end do
    end subroutine init_particles

    subroutine trace_particles_vmec(z0, loss_times, loss_pos, lost_step)
        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        real(dp) :: z(5)
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

!$omp parallel do default(shared) private(i,j,ierr,z) schedule(static)
        do i = 1, n_particles
            z = z0(:, i)
            do j = 1, n_steps
                call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
                if (ierr /= 0 .or. z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    loss_pos(:, i) = z(1:3)
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
!$omp end parallel do
    end subroutine trace_particles_vmec

    subroutine trace_particles_chartmap(z0, loss_times, loss_pos, lost_step)
        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        real(dp) :: z(5), z_chart(5)
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

!$omp parallel do default(shared) private(i,j,ierr,z,z_chart) schedule(static)
        do i = 1, n_particles
            call vmec_to_chartmap(z0(:, i), z_chart)
            z = z_chart
            do j = 1, n_steps
                call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
                if (ierr /= 0 .or. z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    loss_pos(:, i) = z(1:3)
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
!$omp end parallel do
    end subroutine trace_particles_chartmap

    subroutine trace_particles_boozer(z0, loss_times, loss_pos, lost_step)
        !> Trace particles using Boozer canonical coordinates.
        !> Boozer coords: z(1) = s (same as VMEC), angles transformed.
        !> Boundary at s >= 1.0
        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        real(dp) :: z(5), theta_b, phi_b, theta_v, phi_v
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

!$omp parallel do default(shared) private(i,j,ierr,z,theta_b,phi_b,theta_v,phi_v) schedule(static)
        do i = 1, n_particles
            z = z0(:, i)
            call vmec_to_boozer(z(1), z(2), z(3), theta_b, phi_b)
            z(2) = theta_b
            z(3) = phi_b
            do j = 1, n_steps
                call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
                if (ierr /= 0 .or. z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    call boozer_to_vmec(z(1), z(2), z(3), theta_v, phi_v)
                    loss_pos(1, i) = z(1)
                    loss_pos(2, i) = theta_v
                    loss_pos(3, i) = phi_v
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
!$omp end parallel do
    end subroutine trace_particles_boozer

    subroutine trace_particles_meiss(z0, loss_times, loss_pos, lost_step)
        !> Trace particles using Meiss coords built from chartmap field.
        !> Meiss coords: z(1) = r = sqrt(rho) where rho = sqrt(s).
        !> Boundary at r >= 1.0 (which means rho=1, i.e., s=1).
        use field_can_mod, only: ref_to_integ

        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        real(dp) :: z(5), z_chart(5), x_ref(3), x_integ(3)
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

!$omp parallel do default(shared) private(i,j,ierr,z,z_chart,x_ref,x_integ) schedule(static)
        do i = 1, n_particles
            call vmec_to_chartmap(z0(:, i), z_chart)
            x_ref = z_chart(1:3)
            call ref_to_integ(x_ref, x_integ)
            z(1:3) = x_integ
            z(4:5) = z_chart(4:5)

            do j = 1, n_steps
                call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
                ! Check integrator coord directly: r >= 1.0 means rho >= 1 (boundary)
                if (ierr /= 0 .or. z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    call integ_to_ref(z(1:3), x_ref)
                    loss_pos(:, i) = x_ref
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
!$omp end parallel do
    end subroutine trace_particles_meiss

    subroutine trace_particles_meiss_vmec(z0, loss_times, loss_pos, lost_step)
        !> Trace particles using Meiss coords built from VMEC field.
        !> Meiss coords: z(1) = r = sqrt(s) where s is VMEC flux coord.
        !> Boundary at r >= 1.0 (which means s >= 1).

        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        real(dp) :: z(5), x_ref(3), x_integ(3)
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

!$omp parallel do default(shared) private(i,j,ierr,z,x_ref,x_integ) schedule(static)
        do i = 1, n_particles
            x_ref = z0(1:3, i)
            call ref_to_integ(x_ref, x_integ)
            z(1:3) = x_integ
            z(4:5) = z0(4:5, i)

            do j = 1, n_steps
                call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
                ! Check integrator coord directly: r >= 1.0 means s >= 1 (boundary)
                if (ierr /= 0 .or. z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    call integ_to_ref(z(1:3), x_ref)
                    loss_pos(:, i) = x_ref
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
!$omp end parallel do
    end subroutine trace_particles_meiss_vmec

    subroutine trace_particles_boozer_sympl(z0, loss_times, loss_pos, lost_step)
        !> Trace particles using Boozer canonical coordinates with midpoint symplectic.
        use field_can_mod, only: field_can_t, field_can_init, eval_field => evaluate

        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        type(field_can_t) :: f
        type(symplectic_integrator_t) :: integ
        real(dp) :: z(4), z5(5), theta_b, phi_b, theta_v, phi_v
        real(dp) :: dtau_sympl
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

        dtau_sympl = dtau/sqrt(2.0_dp)

!$omp parallel do default(shared) private(i,j,ierr,f,integ,z,z5,theta_b,phi_b,theta_v,phi_v) schedule(static)
        do i = 1, n_particles
            z5 = z0(:, i)
            call vmec_to_boozer(z5(1), z5(2), z5(3), theta_b, phi_b)
            z5(2) = theta_b
            z5(3) = phi_b

            call eval_field(f, z5(1), z5(2), z5(3), 0)
            f%mu = 0.5_dp*z5(4)**2*(1.0_dp - z5(5)**2)/f%Bmod*2.0_dp
            f%ro0 = ro0/sqrt(2.0_dp)
            f%vpar = z5(4)*z5(5)*sqrt(2.0_dp)
            z(1:3) = z5(1:3)
            z(4) = f%vpar*f%hph + f%Aph/f%ro0
            call orbit_sympl_init(integ, f, z, dtau_sympl, 1, 1.0e-12_dp, MIDPOINT)

            do j = 1, n_steps
                ierr = 0
                call orbit_timestep_sympl(integ, f, ierr)
                if (ierr /= 0 .or. integ%z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    call boozer_to_vmec(integ%z(1), integ%z(2), integ%z(3), &
                                        theta_v, phi_v)
                    loss_pos(1, i) = integ%z(1)
                    loss_pos(2, i) = theta_v
                    loss_pos(3, i) = phi_v
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
!$omp end parallel do
    end subroutine trace_particles_boozer_sympl

    subroutine trace_particles_meiss_vmec_sympl(z0, loss_times, loss_pos, lost_step)
        !> Trace particles using Meiss coords from VMEC with midpoint symplectic.
        !> Uses npoiper2=256 for proper timestep based on major radius.
        use field_can_mod, only: field_can_t, field_can_init, eval_field => evaluate
        use new_vmec_stuff_mod, only: rmajor

        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        type(field_can_t) :: f
        type(symplectic_integrator_t) :: integ
        real(dp) :: z(4), z5(5), x_ref(3), x_integ(3)
        real(dp) :: rbig, dtaumin, dtau_sympl
        integer, parameter :: npoiper2 = 256
        integer :: i, j, ierr, ntau_substeps

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

        rbig = rmajor*1.0d2
        dtaumin = twopi*rbig/real(npoiper2, dp)
        ntau_substeps = max(1, ceiling(dtau/dtaumin))
        dtau_sympl = dtaumin/sqrt(2.0_dp)

!$omp parallel do default(shared) private(i,j,ierr,f,integ,z,z5,x_ref,x_integ) schedule(static)
        do i = 1, n_particles
            x_ref = z0(1:3, i)
            call ref_to_integ(x_ref, x_integ)
            z5(1:3) = x_integ
            z5(4:5) = z0(4:5, i)

            call eval_field(f, z5(1), z5(2), z5(3), 0)
            f%mu = 0.5_dp*z5(4)**2*(1.0_dp - z5(5)**2)/f%Bmod*2.0_dp
            f%ro0 = ro0/sqrt(2.0_dp)
            f%vpar = z5(4)*z5(5)*sqrt(2.0_dp)
            z(1:3) = z5(1:3)
            z(4) = f%vpar*f%hph + f%Aph/f%ro0
            call orbit_sympl_init(integ, f, z, dtau_sympl, ntau_substeps, &
                                  1.0e-12_dp, MIDPOINT)

            do j = 1, n_steps
                ierr = 0
                call orbit_timestep_sympl(integ, f, ierr)
                if (ierr /= 0 .or. integ%z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    call integ_to_ref(integ%z(1:3), x_ref)
                    loss_pos(:, i) = x_ref
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
!$omp end parallel do
    end subroutine trace_particles_meiss_vmec_sympl

    subroutine trace_particles_meiss_chart_sympl(z0, loss_times, loss_pos, lost_step)
        !> Trace particles using Meiss coords from chartmap with midpoint symplectic.
        !> Uses npoiper2=256 for proper timestep based on major radius.
        use field_can_mod, only: field_can_t, field_can_init, eval_field => evaluate
        use new_vmec_stuff_mod, only: rmajor

        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        type(field_can_t) :: f
        type(symplectic_integrator_t) :: integ
        real(dp) :: z(4), z5(5), z_chart(5), x_ref(3), x_integ(3)
        real(dp) :: rbig, dtaumin, dtau_sympl
        integer, parameter :: npoiper2 = 256
        integer :: i, j, ierr, ntau_substeps

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

        rbig = rmajor*1.0d2
        dtaumin = twopi*rbig/real(npoiper2, dp)
        ntau_substeps = max(1, ceiling(dtau/dtaumin))
        dtau_sympl = dtaumin/sqrt(2.0_dp)

!$omp parallel do default(shared) private(i,j,ierr,f,integ,z,z5,z_chart,x_ref,x_integ) schedule(static)
        do i = 1, n_particles
            call vmec_to_chartmap(z0(:, i), z_chart)
            x_ref = z_chart(1:3)
            call ref_to_integ(x_ref, x_integ)
            z5(1:3) = x_integ
            z5(4:5) = z_chart(4:5)

            call eval_field(f, z5(1), z5(2), z5(3), 0)
            f%mu = 0.5_dp*z5(4)**2*(1.0_dp - z5(5)**2)/f%Bmod*2.0_dp
            f%ro0 = ro0/sqrt(2.0_dp)
            f%vpar = z5(4)*z5(5)*sqrt(2.0_dp)
            z(1:3) = z5(1:3)
            z(4) = f%vpar*f%hph + f%Aph/f%ro0
            call orbit_sympl_init(integ, f, z, dtau_sympl, ntau_substeps, &
                                  1.0e-12_dp, MIDPOINT)

            do j = 1, n_steps
                ierr = 0
                call orbit_timestep_sympl(integ, f, ierr)
                if (ierr /= 0 .or. integ%z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    call integ_to_ref(integ%z(1:3), x_ref)
                    loss_pos(:, i) = x_ref
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
!$omp end parallel do
    end subroutine trace_particles_meiss_chart_sympl

    subroutine vmec_to_chartmap(z_vmec, z_chart)
        real(dp), intent(in) :: z_vmec(5)
        real(dp), intent(out) :: z_chart(5)

        real(dp) :: x_cyl(3)
        integer :: ierr

        call ref_coords%evaluate_cyl(z_vmec(1:3), x_cyl)
        call chartmap_coords%from_cyl(x_cyl, z_chart(1:3), ierr)
        if (ierr /= 0) then
            print *, 'vmec_to_chartmap: from_cyl failed ierr=', ierr
            error stop 1
        end if
        z_chart(4:5) = z_vmec(4:5)
    end subroutine vmec_to_chartmap

    function count_confined(loss_times) result(frac)
        real(dp), intent(in) :: loss_times(n_particles)
        real(dp) :: frac
        integer :: n_conf

        n_conf = count(loss_times >= trace_time*0.999_dp)
        frac = real(n_conf, dp)/real(n_particles, dp)
    end function count_confined

    subroutine write_output_netcdf
        integer :: ncid, dimid_part
        integer :: varid_lt_vmec, varid_lt_chart, varid_lt_boozer
        integer :: varid_lt_meiss_v, varid_lt_meiss_c
        integer :: varid_lp_vmec, varid_lp_chart, varid_lp_boozer
        integer :: varid_lp_meiss_v, varid_lp_meiss_c
        integer :: varid_ls_vmec, varid_ls_chart, varid_ls_boozer
        integer :: varid_ls_meiss_v, varid_ls_meiss_c
        integer :: varid_z0
        integer :: dimid_5, dimid_3
        integer :: status

        status = nf90_create(trim(output_nc), nf90_netcdf4, ncid)
        if (status /= nf90_noerr) error stop 'Failed to create NetCDF'

        status = nf90_def_dim(ncid, 'n_particles', n_particles, dimid_part)
        status = nf90_def_dim(ncid, 'dim5', 5, dimid_5)
        status = nf90_def_dim(ncid, 'dim3', 3, dimid_3)

        status = nf90_put_att(ncid, nf90_global, 'description', &
                              'Chartmap wall losses comparison')
        status = nf90_put_att(ncid, nf90_global, 'trace_time', trace_time)
        status = nf90_put_att(ncid, nf90_global, 'n_steps', n_steps)

        status = nf90_def_var(ncid, 'z0', nf90_double, [dimid_5, dimid_part], varid_z0)
        status = nf90_def_var(ncid, 'loss_time_vmec', nf90_double, [dimid_part], &
                              varid_lt_vmec)
        status = nf90_def_var(ncid, 'loss_time_chartmap', nf90_double, [dimid_part], &
                              varid_lt_chart)
        status = nf90_def_var(ncid, 'loss_time_boozer', nf90_double, [dimid_part], &
                              varid_lt_boozer)
        status = nf90_def_var(ncid, 'loss_time_meiss_vmec', nf90_double, [dimid_part], &
                              varid_lt_meiss_v)
        status = nf90_def_var(ncid, 'loss_time_meiss_chart', nf90_double, &
                              [dimid_part], &
                              varid_lt_meiss_c)
        status = nf90_def_var(ncid, 'loss_pos_vmec', nf90_double, &
                              [dimid_3, dimid_part], varid_lp_vmec)
        status = nf90_def_var(ncid, 'loss_pos_chartmap', nf90_double, &
                              [dimid_3, dimid_part], varid_lp_chart)
        status = nf90_def_var(ncid, 'loss_pos_boozer', nf90_double, &
                              [dimid_3, dimid_part], varid_lp_boozer)
        status = nf90_def_var(ncid, 'loss_pos_meiss_vmec', nf90_double, &
                              [dimid_3, dimid_part], varid_lp_meiss_v)
        status = nf90_def_var(ncid, 'loss_pos_meiss_chart', nf90_double, &
                              [dimid_3, dimid_part], varid_lp_meiss_c)
        status = nf90_def_var(ncid, 'lost_step_vmec', nf90_int, [dimid_part], &
                              varid_ls_vmec)
        status = nf90_def_var(ncid, 'lost_step_chartmap', nf90_int, [dimid_part], &
                              varid_ls_chart)
        status = nf90_def_var(ncid, 'lost_step_boozer', nf90_int, [dimid_part], &
                              varid_ls_boozer)
        status = nf90_def_var(ncid, 'lost_step_meiss_vmec', nf90_int, [dimid_part], &
                              varid_ls_meiss_v)
        status = nf90_def_var(ncid, 'lost_step_meiss_chart', nf90_int, [dimid_part], &
                              varid_ls_meiss_c)

        status = nf90_enddef(ncid)

        status = nf90_put_var(ncid, varid_z0, z0)
        status = nf90_put_var(ncid, varid_lt_vmec, loss_times_vmec)
        status = nf90_put_var(ncid, varid_lt_chart, loss_times_chart)
        status = nf90_put_var(ncid, varid_lt_boozer, loss_times_boozer)
        status = nf90_put_var(ncid, varid_lt_meiss_v, loss_times_meiss_vmec)
        status = nf90_put_var(ncid, varid_lt_meiss_c, loss_times_meiss_chart)
        status = nf90_put_var(ncid, varid_lp_vmec, loss_pos_vmec)
        status = nf90_put_var(ncid, varid_lp_chart, loss_pos_chart)
        status = nf90_put_var(ncid, varid_lp_boozer, loss_pos_boozer)
        status = nf90_put_var(ncid, varid_lp_meiss_v, loss_pos_meiss_vmec)
        status = nf90_put_var(ncid, varid_lp_meiss_c, loss_pos_meiss_chart)
        status = nf90_put_var(ncid, varid_ls_vmec, lost_step_vmec)
        status = nf90_put_var(ncid, varid_ls_chart, lost_step_chart)
        status = nf90_put_var(ncid, varid_ls_boozer, lost_step_boozer)
        status = nf90_put_var(ncid, varid_ls_meiss_v, lost_step_meiss_vmec)
        status = nf90_put_var(ncid, varid_ls_meiss_c, lost_step_meiss_chart)

        status = nf90_close(ncid)
        print *, '  Wrote: ', trim(output_nc)
    end subroutine write_output_netcdf

end program test_chartmap_wall_losses
