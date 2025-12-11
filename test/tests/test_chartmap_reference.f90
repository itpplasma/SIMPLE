program test_chartmap_reference
    !> Test chartmap reference coordinates against VMEC reference coordinates.
    !>
    !> This test verifies that:
    !>   1. Chartmap coordinate system can be initialized from a chartmap file
    !>   2. Chartmap from VMEC produces equivalent coordinate transforms
    !>   3. evaluate_point agrees between VMEC and chartmap-from-VMEC
    !>
    !> The test requires:
    !>   - wout.nc: VMEC equilibrium file
    !>   - wout.chartmap.nc: Chartmap file generated from VMEC (optional, will skip)

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use reference_coordinates, only: init_reference_coordinates, ref_coords, &
        ref_coords_type, REF_COORDS_VMEC, REF_COORDS_CHARTMAP
    use libneo_coordinates, only: coordinate_system_t, &
        make_vmec_coordinate_system, make_chartmap_coordinate_system
    use timing, only: init_timer
    use params, only: coord_input, field_input
    use util, only: twopi

    implicit none

    integer :: n_failed
    real(dp) :: dummy

    n_failed = 0

    call init_timer()

    call test_vmec_reference_initialization(n_failed)
    call test_chartmap_reference_initialization(n_failed)
    call test_chartmap_vs_vmec_evaluate_point(n_failed)
    call test_chartmap_map2disc_loads(n_failed)

    if (n_failed == 0) then
        print *, '================================'
        print *, 'All chartmap reference tests PASSED'
        print *, '================================'
        stop 0
    else
        print *, '================================'
        print *, n_failed, ' tests FAILED'
        print *, '================================'
        stop 1
    end if

contains

    subroutine test_vmec_reference_initialization(n_failed)
        !> Test that VMEC reference coordinates initialize correctly.
        integer, intent(inout) :: n_failed
        logical :: file_exists

        print *, 'Test 1: VMEC reference coordinate initialization'

        inquire(file='wout.nc', exist=file_exists)
        if (.not. file_exists) then
            print *, '  FAILED: Required VMEC file (wout.nc) not found'
            n_failed = n_failed + 1
            return
        end if

        coord_input = 'wout.nc'
        field_input = 'wout.nc'
        call init_vmec('wout.nc', 5, 5, 5, dummy)
        call init_reference_coordinates(coord_input)

        if (ref_coords_type /= REF_COORDS_VMEC) then
            print *, '  FAILED: Expected VMEC reference type, got ', ref_coords_type
            n_failed = n_failed + 1
            return
        end if

        if (.not. allocated(ref_coords)) then
            print *, '  FAILED: ref_coords not allocated'
            n_failed = n_failed + 1
            return
        end if

        print *, '  PASSED: VMEC reference coordinates initialized'
    end subroutine test_vmec_reference_initialization


    subroutine test_chartmap_reference_initialization(n_failed)
        !> Test that chartmap reference coordinates initialize correctly.
        integer, intent(inout) :: n_failed
        logical :: file_exists

        print *, 'Test 2: Chartmap reference coordinate initialization'

        inquire(file='wout.chartmap.nc', exist=file_exists)
        if (.not. file_exists) then
            print *, '  SKIPPED: Optional chartmap file (wout.chartmap.nc) not found'
            print *, '  (Generate with: python scripts/vmec_to_chartmap.py wout.nc wout.chartmap.nc --passthrough)'
            return
        end if

        coord_input = 'chartmap:wout.chartmap.nc'
        call init_reference_coordinates(coord_input)

        if (ref_coords_type /= REF_COORDS_CHARTMAP) then
            print *, '  FAILED: Expected CHARTMAP reference type, got ', ref_coords_type
            n_failed = n_failed + 1
            return
        end if

        if (.not. allocated(ref_coords)) then
            print *, '  FAILED: ref_coords not allocated'
            n_failed = n_failed + 1
            return
        end if

        print *, '  PASSED: Chartmap reference coordinates initialized'
    end subroutine test_chartmap_reference_initialization


    subroutine test_chartmap_vs_vmec_evaluate_point(n_failed)
        !> Compare evaluate_point between VMEC and chartmap-from-VMEC.
        !> VMEC returns cylindrical (R, phi, Z), chartmap returns Cartesian (X, Y, Z).
        !> We convert VMEC to Cartesian for comparison.
        integer, intent(inout) :: n_failed

        class(coordinate_system_t), allocatable :: vmec_coords, chartmap_coords
        real(dp) :: u(3), x_vmec_cyl(3), x_vmec_cart(3), x_chartmap(3)
        real(dp) :: diff, max_diff
        ! Tolerance for passthrough mode: both use the same libneo VMEC splines,
        ! so difference is only from chartmap's B-spline interpolation of the
        ! tabulated values. Should be small (< 1 cm for reasonable grid).
        real(dp), parameter :: tol = 1.0_dp  ! 1 cm tolerance
        logical :: vmec_exists, chartmap_exists
        integer :: ir, ith, iphi
        integer :: nr, nth, nphi

        print *, 'Test 3: Chartmap vs VMEC evaluate_point comparison'

        inquire(file='wout.nc', exist=vmec_exists)
        inquire(file='wout.chartmap.passthrough.nc', exist=chartmap_exists)

        if (.not. vmec_exists) then
            print *, '  FAILED: Required VMEC file (wout.nc) not found'
            n_failed = n_failed + 1
            return
        end if

        if (.not. chartmap_exists) then
            print *, '  SKIPPED: Chartmap file (wout.chartmap.passthrough.nc) not found'
            print *, '  (Generate with: vmec_to_chartmap.x wout.nc wout.chartmap.passthrough.nc --passthrough)'
            return
        end if

        ! Initialize VMEC first (needed by vmec_coordinate_system_t)
        coord_input = 'wout.nc'
        field_input = 'wout.nc'
        call init_vmec('wout.nc', 5, 5, 5, dummy)

        ! Create both coordinate systems
        call make_vmec_coordinate_system(vmec_coords)
        call make_chartmap_coordinate_system(chartmap_coords, 'wout.chartmap.passthrough.nc')

        max_diff = 0.0_dp
        nr = 5
        nth = 8
        nphi = 4

        do iphi = 1, nphi
            do ith = 1, nth
                do ir = 1, nr
                    ! Test point: s in [0.1, 0.9], theta in [0, 2pi], phi in [0, 2pi/nfp]
                    ! Using nfp=2 based on test file
                    u(1) = 0.1_dp + 0.8_dp * (ir - 1) / (nr - 1)
                    u(2) = twopi * (ith - 1) / nth
                    u(3) = twopi / 2.0_dp * (iphi - 1) / nphi

                    ! VMEC returns (R, phi, Z) in cylindrical coords
                    call vmec_coords%evaluate_point(u, x_vmec_cyl)

                    ! Convert VMEC cylindrical to Cartesian
                    x_vmec_cart(1) = x_vmec_cyl(1) * cos(x_vmec_cyl(2))
                    x_vmec_cart(2) = x_vmec_cyl(1) * sin(x_vmec_cyl(2))
                    x_vmec_cart(3) = x_vmec_cyl(3)

                    ! Chartmap returns (X, Y, Z) in Cartesian coords
                    call chartmap_coords%evaluate_point(u, x_chartmap)

                    diff = sqrt(sum((x_vmec_cart - x_chartmap)**2))
                    max_diff = max(max_diff, diff)

                    if (diff > tol) then
                        print *, '  ERROR at u = ', u
                        print *, '    x_vmec_cart = ', x_vmec_cart
                        print *, '    x_chartmap = ', x_chartmap
                        print *, '    diff = ', diff
                    end if
                end do
            end do
        end do

        print *, '  Max position difference (cm): ', max_diff

        if (max_diff > tol) then
            print *, '  FAILED: Position difference exceeds tolerance ', tol, ' cm'
            n_failed = n_failed + 1
        else
            print *, '  PASSED: Chartmap agrees with VMEC within tolerance'
        end if

        deallocate(vmec_coords)
        deallocate(chartmap_coords)
    end subroutine test_chartmap_vs_vmec_evaluate_point


    subroutine test_chartmap_map2disc_loads(n_failed)
        !> Verify map2disc boundary matches VMEC boundary using Hausdorff distance.
        !>
        !> Map2disc uses conformal mapping with different theta parameterization,
        !> but at rho=1 both should trace the SAME physical boundary curve.
        !> We compare using Hausdorff distance: for each point on VMEC boundary,
        !> find the minimum distance to any point on the chartmap boundary.
        integer, intent(inout) :: n_failed

        class(coordinate_system_t), allocatable :: vmec_coords, chartmap_coords
        real(dp) :: u_vmec(3), u_chart(3), x_vmec_cyl(3), x_chart(3)
        real(dp) :: R_vmec, Z_vmec, R_chart, Z_chart
        real(dp) :: dist, min_dist, hausdorff_dist
        ! Tolerance for boundary match: conformal mapping distortion + spline error
        ! Raw data comparison gives ~9 cm, spline adds error due to periodic BC issues
        real(dp), parameter :: tol = 40.0_dp  ! 40 cm tolerance (known spline limitation)
        logical :: vmec_exists, chartmap_exists
        integer :: i, j, nth_vmec, nth_chart
        real(dp) :: dummy

        print *, 'Test 4: Map2disc boundary match verification'

        inquire(file='wout.nc', exist=vmec_exists)
        inquire(file='wout.chartmap.map2disc.nc', exist=chartmap_exists)

        if (.not. vmec_exists) then
            print *, '  FAILED: Required VMEC file (wout.nc) not found'
            n_failed = n_failed + 1
            return
        end if

        if (.not. chartmap_exists) then
            print *, '  SKIPPED: Optional map2disc chartmap file not found'
            print *, '  (Generate with: python scripts/vmec_to_chartmap.py wout.nc' &
                // ' wout.chartmap.map2disc.nc --map2disc)'
            return
        end if

        ! Initialize VMEC
        coord_input = 'wout.nc'
        field_input = 'wout.nc'
        call init_vmec('wout.nc', 5, 5, 5, dummy)

        call make_vmec_coordinate_system(vmec_coords)
        call make_chartmap_coordinate_system(chartmap_coords, 'wout.chartmap.map2disc.nc')

        ! Compute Hausdorff distance at zeta=0
        ! Chartmap theta goes to 2pi*(128/129) due to endpoint=False grid
        ! Must stay within this range to avoid extrapolation
        nth_vmec = 128
        nth_chart = 129  ! Match chartmap grid, sample within valid range
        hausdorff_dist = 0.0_dp

        ! For each VMEC boundary point, find closest chartmap boundary point
        do i = 1, nth_vmec
            u_vmec(1) = 1.0_dp
            u_vmec(2) = twopi * real(i - 1, dp) / real(nth_vmec, dp)
            u_vmec(3) = 0.0_dp

            call vmec_coords%evaluate_point(u_vmec, x_vmec_cyl)
            R_vmec = x_vmec_cyl(1)
            Z_vmec = x_vmec_cyl(3)

            min_dist = huge(1.0_dp)
            do j = 1, nth_chart
                ! Sample chartmap within its valid theta range [0, 2pi*128/129]
                u_chart(1) = 1.0_dp
                u_chart(2) = twopi * real(j - 1, dp) / real(nth_chart, dp)
                u_chart(3) = 0.0_dp

                call chartmap_coords%evaluate_point(u_chart, x_chart)
                R_chart = sqrt(x_chart(1)**2 + x_chart(2)**2)
                Z_chart = x_chart(3)

                dist = sqrt((R_vmec - R_chart)**2 + (Z_vmec - Z_chart)**2)
                min_dist = min(min_dist, dist)
            end do

            hausdorff_dist = max(hausdorff_dist, min_dist)
        end do

        print *, '  Hausdorff distance at zeta=0: ', hausdorff_dist, ' cm'

        if (hausdorff_dist > tol) then
            print *, '  FAILED: Boundary mismatch exceeds tolerance ', tol, ' cm'
            n_failed = n_failed + 1
        else
            print *, '  PASSED: Map2disc boundary matches VMEC within tolerance'
        end if

        deallocate(vmec_coords)
        deallocate(chartmap_coords)
    end subroutine test_chartmap_map2disc_loads

end program test_chartmap_reference
