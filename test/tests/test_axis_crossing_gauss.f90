program test_axis_crossing_gauss
    !> Gauss and Lobatto steps through the magnetic axis must commit the chart
    !> switch (r, theta) -> (|r|, theta + pi) as the Euler and midpoint steps do
    !> (#370): an accepted state never has a negative radius, which the output
    !> cannot map to cylindrical coordinates, and p_theta stays the value of the
    !> stored endpoint. The orbit starts next to the QA Boozer axis.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi
    use simple, only: tracer_t, init_sympl, init_params
    use simple_main, only: init_field
    use params, only: field_input, coord_input
    use new_vmec_stuff_mod, only: rmajor
    use orbit_symplectic, only: orbit_timestep_sympl
    use orbit_symplectic_base, only: GAUSS1, GAUSS2, GAUSS3, GAUSS4, LOBATTO3
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_val
    use diag_counters, only: diag_counters_init, diag_counters_total, &
                             EVT_R_NEGATIVE
    use velo_mod, only: isw_field_type
    use magfie_sub, only: BOOZER

    implicit none

    integer, parameter :: max_steps = 2000
    real(dp), parameter :: pth_tol = 1.0e-10_dp

    type(tracer_t) :: norb
    type(field_can_t) :: fcheck
    real(dp) :: z5(5), dtau, defect
    integer :: ierr, kt, mode
    integer(8) :: ncross_before
    character(16) :: scheme

    call get_command_argument(1, scheme)
    select case (trim(scheme))
    case ('gauss1')
        mode = GAUSS1
    case ('gauss2')
        mode = GAUSS2
    case ('gauss3')
        mode = GAUSS3
    case ('gauss4')
        mode = GAUSS4
    case ('lobatto3')
        mode = LOBATTO3
    case default
        error stop 'expected gauss1..gauss4 or lobatto3'
    end select

    isw_field_type = BOOZER
    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    call init_field(norb, 'wout.nc', 5, 5, 5, 1)
    call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0e-12_dp)
    call diag_counters_init

    dtau = twopi*rmajor*1.0e2_dp/64.0_dp
    z5 = [1.0e-5_dp, 2.0_dp, 0.1_dp, 1.0_dp, 0.95_dp]
    call init_sympl(norb%si, norb%f, z5, dtau, dtau, 1.0e-12_dp, mode)

    do kt = 1, max_steps
        ncross_before = diag_counters_total(EVT_R_NEGATIVE)
        call orbit_timestep_sympl(norb%si, norb%f, ierr)
        if (ierr /= 0) then
            print *, 'FAIL: step', kt, 'returned', ierr
            error stop 1
        end if
        if (norb%si%z(1) < 0.0_dp) then
            print *, 'FAIL: accepted state at step', kt, 'has r =', norb%si%z(1)
            error stop 1
        end if
        if (diag_counters_total(EVT_R_NEGATIVE) /= ncross_before) exit
    end do
    if (kt > max_steps) then
        print *, 'FAIL: no axis crossing was committed'
        error stop 1
    end if

    fcheck = norb%f
    call eval_field(fcheck, norb%si%z(1), norb%si%z(2), norb%si%z(3), 0)
    call get_val(fcheck, norb%si%z(4))
    defect = abs(norb%f%pth - fcheck%pth)/max(abs(fcheck%pth), tiny(1.0_dp))
    print '(a,a,a,i0,a,es10.2,a,es10.2)', ' ', trim(scheme), ': crossing at step ', &
        kt, '  r = ', norb%si%z(1), '  relative p_theta defect: ', defect
    if (.not. (defect <= pth_tol)) then
        print *, 'FAIL: carried p_theta differs from the fresh endpoint value'
        error stop 1
    end if
    print *, 'test_axis_crossing_gauss PASSED'
end program test_axis_crossing_gauss
