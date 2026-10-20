program test_axis_crossing_midpoint
    !> A chart switch must carry p_theta evaluated at the stored endpoint.
    !> The QA Boozer orbit below reaches the existing warning-mode acceptance,
    !> not a converged midpoint root. This test checks endpoint consistency;
    !> it does not establish the accuracy of that acceptance or later steps.
    !> The same field witnesses both the single-volume and sympl_rmax > 1
    !> dispatch paths without requiring a SPECTRE equilibrium.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi
    use simple, only: tracer_t, init_sympl, init_params
    use simple_main, only: init_field
    use params, only: field_input, coord_input
    use new_vmec_stuff_mod, only: rmajor
    use orbit_symplectic, only: orbit_timestep_sympl
    use orbit_symplectic_base, only: MIDPOINT, sympl_rmax
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_val
    use diag_counters, only: diag_counters_init, diag_counters_total, &
                             EVT_R_NEGATIVE
    use velo_mod, only: isw_field_type
    use magfie_sub, only: BOOZER

    implicit none

    integer, parameter :: max_steps = 2000
    real(dp), parameter :: pth_tol = 1.0e-12_dp

    type(tracer_t) :: norb
    type(field_can_t) :: fcheck
    real(dp) :: z5(5), dtau, defect
    integer :: ierr, kt
    integer(8) :: ncross_before
    character(16) :: path

    isw_field_type = BOOZER
    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    call init_field(norb, 'wout.nc', 5, 5, 5, 1)
    call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0e-12_dp)
    call get_command_argument(1, path)
    select case (trim(path))
    case ('single')
        sympl_rmax = 1.0_dp
    case ('multivolume')
        sympl_rmax = 2.0_dp
    case default
        error stop 'expected single or multivolume dispatch path'
    end select

    call diag_counters_init

    ! A coarse step produces an accepted negative endpoint near the axis.
    dtau = twopi*rmajor*1.0e2_dp/64.0_dp
    z5 = [1.0e-5_dp, 2.0_dp, 0.1_dp, 1.0_dp, 0.95_dp]
    call init_sympl(norb%si, norb%f, z5, dtau, dtau, 1.0e-12_dp, MIDPOINT)

    do kt = 1, max_steps
        ncross_before = diag_counters_total(EVT_R_NEGATIVE)
        call orbit_timestep_sympl(norb%si, norb%f, ierr)
        if (ierr /= 0) then
            print *, 'FAIL: step', kt, 'returned', ierr, 'before a crossing'
            error stop 1
        end if
        if (diag_counters_total(EVT_R_NEGATIVE) /= ncross_before) exit
    end do
    if (kt > max_steps) then
        print *, 'FAIL: no accepted midpoint axis crossing was produced'
        error stop 1
    end if

    fcheck = norb%f
    call eval_field(fcheck, norb%si%z(1), norb%si%z(2), norb%si%z(3), 0)
    call get_val(fcheck, norb%si%z(4))
    defect = abs(norb%f%pth - fcheck%pth)/max(abs(fcheck%pth), tiny(1.0_dp))
    print '(a,i0,a,es10.2,a,es10.2)', ' crossing at step ', kt, '  r = ', &
        norb%si%z(1), '  relative p_theta defect: ', defect
    if (.not. (defect <= pth_tol)) then
        print *, 'FAIL: carried p_theta differs from the fresh endpoint value'
        error stop 1
    end if
    print *, 'test_axis_crossing_midpoint PASSED'
end program test_axis_crossing_midpoint
