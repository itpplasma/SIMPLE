program test_refcoords_file_detection
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: detect_refcoords_file_type, &
        refcoords_file_chartmap, refcoords_file_vmec_wout, refcoords_file_unknown
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use simple, only: init_vmec
    use timing, only: init_timer
    use params, only: coord_input
    implicit none

    integer :: nerrors
    real(dp) :: dummy

    nerrors = 0

    call init_timer()

    call test_detect_vmec_file(nerrors)
    call test_detect_chartmap_file(nerrors)
    call test_init_reference_coords_vmec(nerrors)
    call test_init_reference_coords_chartmap(nerrors)

    if (nerrors > 0) then
        print *, '================================'
        print *, 'FAILED: ', nerrors, ' error(s) in refcoords file detection'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'All refcoords file detection tests PASSED'
    print *, '================================'

contains

    subroutine test_detect_vmec_file(nerrors)
        integer, intent(inout) :: nerrors

        integer :: file_type, ierr
        character(len=2048) :: message
        logical :: file_exists

        print *, 'Test 1: detect VMEC wout file type'

        inquire(file='wout.nc', exist=file_exists)
        if (.not. file_exists) then
            print *, '  FAIL: wout.nc not found in test directory'
            nerrors = nerrors + 1
            return
        end if

        call detect_refcoords_file_type('wout.nc', file_type, ierr, message)
        if (ierr /= 0) then
            print *, '  FAILED: detect_refcoords_file_type error: ', trim(message)
            nerrors = nerrors + 1
            return
        end if

        if (file_type /= refcoords_file_vmec_wout) then
            print *, '  FAILED: expected VMEC_WOUT file_type, got ', file_type
            nerrors = nerrors + 1
            return
        end if

        print *, '  PASSED: wout.nc detected as VMEC_WOUT'
    end subroutine test_detect_vmec_file

    subroutine test_detect_chartmap_file(nerrors)
        integer, intent(inout) :: nerrors

        integer :: file_type, ierr
        character(len=2048) :: message
        logical :: file_exists

        print *, 'Test 2: detect chartmap file type'

        inquire(file='wout.chartmap.nc', exist=file_exists)
        if (.not. file_exists) then
            print *, '  FAIL: wout.chartmap.nc not found in test directory'
            nerrors = nerrors + 1
            return
        end if

        call detect_refcoords_file_type('wout.chartmap.nc', file_type, ierr, message)
        if (ierr /= 0) then
            print *, '  FAILED: detect_refcoords_file_type error: ', trim(message)
            nerrors = nerrors + 1
            return
        end if

        if (file_type /= refcoords_file_chartmap) then
            print *, '  FAILED: expected CHARTMAP file_type, got ', file_type
            nerrors = nerrors + 1
            return
        end if

        print *, '  PASSED: wout.chartmap.nc detected as CHARTMAP'
    end subroutine test_detect_chartmap_file

    subroutine test_init_reference_coords_vmec(nerrors)
        integer, intent(inout) :: nerrors

        logical :: file_exists

        print *, 'Test 3: init_reference_coordinates with VMEC file'

        inquire(file='wout.nc', exist=file_exists)
        if (.not. file_exists) then
            print *, '  FAIL: wout.nc not found in test directory'
            nerrors = nerrors + 1
            return
        end if

        coord_input = 'wout.nc'
        call init_vmec('wout.nc', 5, 5, 5, dummy)
        call init_reference_coordinates(coord_input)

        if (.not. allocated(ref_coords)) then
            print *, '  FAILED: ref_coords not allocated after init'
            nerrors = nerrors + 1
            return
        end if

        print *, '  PASSED: ref_coords allocated for VMEC file'
    end subroutine test_init_reference_coords_vmec

    subroutine test_init_reference_coords_chartmap(nerrors)
        integer, intent(inout) :: nerrors

        logical :: file_exists
        real(dp) :: u(3), x(3)

        print *, 'Test 4: init_reference_coordinates with chartmap file'

        inquire(file='wout.chartmap.nc', exist=file_exists)
        if (.not. file_exists) then
            print *, '  FAIL: wout.chartmap.nc not found in test directory'
            nerrors = nerrors + 1
            return
        end if

        coord_input = 'wout.chartmap.nc'
        call init_reference_coordinates(coord_input)

        if (.not. allocated(ref_coords)) then
            print *, '  FAILED: ref_coords not allocated after init'
            nerrors = nerrors + 1
            return
        end if

        u = [0.5_dp, 0.0_dp, 0.0_dp]
        call ref_coords%evaluate_cyl(u, x)

        if (x(1) <= 0.0_dp) then
            print *, '  FAILED: chartmap evaluate_cyl returned invalid R=', x(1)
            nerrors = nerrors + 1
            return
        end if

        print *, '  PASSED: ref_coords allocated and functional for chartmap file'
        print *, '    evaluate_cyl([0.5, 0, 0]) -> x =', x
    end subroutine test_init_reference_coords_chartmap

end program test_refcoords_file_detection
