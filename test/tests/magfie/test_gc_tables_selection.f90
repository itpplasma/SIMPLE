program test_gc_tables_selection
    ! RK guiding-centre runs in reference coordinates on a rho chartmap use the
    ! guiding-centre tables for a coil file and for a frozen-current file,
    ! recognised by content and without an input key. A frozen-current file
    ! holding only the same coils gives the same tables. A VMEC field_input
    ! keeps the splined field and drops tables left from the previous field.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: coord_input, field_input
    use velo_mod, only: isw_field_type
    use magfie_sub, only: REFCOORDS, gc_tables, has_magfie_refcoords_field
    use field_gc_tables, only: evaluate_gc_tables
    use field_coils, only: coils_field_t, create_coils_field
    use neo_frozen_current_field, only: frozen_current_field_t, frozen_current_patch_t
    use neo_frozen_current_io, only: save_frozen_current_field
    use fortfem_meridian_bspline_n, only: meridian_bspline_build
    use neo_guiding_centre_field, only: GC_COUNT
    implicit none

    real(dp), parameter :: u(3, 3) = reshape([0.2_dp, 0.3_dp, 0.1_dp, &
                                              0.5_dp, 2.0_dp, 0.4_dp, &
                                              0.8_dp, 4.5_dp, 0.9_dp], [3, 3])
    type(tracer_t) :: norb
    real(dp) :: from_coils(GC_COUNT, 3), from_frozen(GC_COUNT, 3)

    call write_frozen_coils('coils.simple', 'frozen_ncsx.field')
    coord_input = 'wout_ncsx.chartmap.nc'
    isw_field_type = REFCOORDS

    field_input = 'coils.simple'
    call init_field(norb, 'wout_ncsx.nc', 5, 5, 3, 0)
    if (.not. gc_tables%initialized) error stop 'coil file: no tables'
    call sample(from_coils)

    field_input = 'frozen_ncsx.field'
    call init_field(norb, 'wout_ncsx.nc', 5, 5, 3, 0)
    if (.not. gc_tables%initialized) error stop 'frozen-current file: no tables'
    call sample(from_frozen)
    if (any(abs(from_frozen - from_coils) > 1.0e-10_dp*abs(from_coils) + &
            1.0e-12_dp*maxval(abs(from_coils)))) then
        print *, maxval(abs(from_frozen - from_coils))
        error stop 'frozen-current and coil tables differ'
    end if

    field_input = 'wout_ncsx.nc'
    call init_field(norb, 'wout_ncsx.nc', 5, 5, 3, 0)
    if (gc_tables%initialized) error stop 'VMEC field_input: stale tables'
    if (.not. has_magfie_refcoords_field()) error stop 'VMEC field_input: no field'
    print *, 'PASSED: guiding-centre tables chosen from field_input'

contains

    subroutine sample(q)
        real(dp), intent(out) :: q(GC_COUNT, 3)
        integer :: k

        do k = 1, 3
            call evaluate_gc_tables(gc_tables, u(:, k), q(1, k), q(2, k), &
                                    q(3:5, k), q(6:8, k), q(9:11, k), q(12:14, k))
        end do
    end subroutine sample

    !> Frozen-current file whose plasma part vanishes on a patch around the
    !> chart, so its field is that of the coils alone.
    subroutine write_frozen_coils(coil_file, path)
        character(*), intent(in) :: coil_file, path
        type(coils_field_t) :: raw
        type(frozen_current_field_t) :: f
        type(frozen_current_patch_t) :: patch
        integer :: status

        call create_coils_field(coil_file, raw)
        f%coils = raw%coils
        patch%r0 = 30
        patch%r1 = 300
        patch%z0 = -150
        patch%z1 = 150
        patch%nr = 4
        patch%nz = 4
        patch%kmax = 0
        call meridian_bspline_build(patch%r0, patch%r1, patch%z0, patch%z1, patch%nr, &
                                    patch%nz, patch%degree, f%plasma%sp, status)
        if (status /= 0) error stop 'meridian spline build'
        allocate (f%plasma%u(f%plasma%sp%n1, 0:0))
        f%plasma%u = 0
        call save_frozen_current_field(path, patch, f)
    end subroutine write_frozen_coils

end program test_gc_tables_selection
