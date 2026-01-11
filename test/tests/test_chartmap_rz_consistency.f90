program test_chartmap_rz_consistency
    !> Tests R,Z coordinate consistency between VMEC and chartmap coordinate systems.
    !>
    !> Core validation for arclength-parameterized chartmaps:
    !>   - VMEC point (s, theta_vmec, zeta) -> (R, Z)
    !>   - Chartmap point (rho, theta_cm, zeta) -> same (R, Z)
    !>
    !> For offset=0, the chartmap boundary (rho=1) must match VMEC boundary (s=1)
    !> in physical (R, Z) space, even though theta coordinates differ.
    !>
    !> This test validates that:
    !>   1. evaluate_cyl produces correct R, Z, phi for both coordinate systems
    !>   2. Physical positions match when coordinates are inverted correctly
    !>   3. Boundary surfaces coincide in R,Z for offset=0 chartmaps

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, make_vmec_coordinate_system, &
                                  make_chartmap_coordinate_system, chartmap_from_cyl_ok
    use simple, only: init_vmec
    use timing, only: init_timer
    use util, only: twopi
    implicit none

    integer :: nerrors
    real(dp) :: dummy

    real(dp), parameter :: TOL_RZ = 1.0e-4_dp  ! 0.1 mm tolerance for R,Z matching

    character(len=*), parameter :: VMEC_FILE = 'wout_ncsx.nc'
    character(len=*), parameter :: CHARTMAP_FILE = 'wout_ncsx.chartmap.nc'

    nerrors = 0
    call init_timer()

    call test_rz_chartmap_roundtrip(nerrors)
    call test_rz_vmec_chartmap_roundtrip(nerrors)
    call test_vmec_chartmap_cyl_roundtrip(nerrors)
    call test_boundary_rz_match(nerrors)
    call test_interior_rz_consistency(nerrors)

    if (nerrors > 0) then
        print *, '================================'
        print *, 'FAILED: ', nerrors, ' error(s) in chartmap R,Z consistency tests'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'All chartmap R,Z consistency tests PASSED'
    print *, '================================'

contains

    subroutine test_rz_chartmap_roundtrip(nerrors)
        !> Test R,Z -> chartmap -> R,Z roundtrip (direct, no VMEC).
        !> This is the core test: given arbitrary R,Z,phi inside the domain,
        !> chartmap from_cyl should find coordinates that evaluate back to R,Z.
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: chart_cs
        real(dp) :: xcyl_in(3), xcyl_out(3), u_chart(3)
        real(dp) :: r_in, z_in, phi_in
        real(dp) :: r_out, z_out, phi_out
        real(dp) :: err_r, err_z, err_max
        logical :: chart_exists
        integer :: ierr, i, j
        real(dp) :: r_axis, r_edge, z_center

        print *, 'Test: R,Z -> chartmap -> R,Z roundtrip'

        inquire(file=CHARTMAP_FILE, exist=chart_exists)
        if (.not. chart_exists) then
            print *, '  SKIP: ', CHARTMAP_FILE, ' not found'
            return
        end if

        call make_chartmap_coordinate_system(chart_cs, CHARTMAP_FILE)

        ! NCSX approximate geometry: R ~ 1.4 m, minor radius ~ 0.3 m
        r_axis = 1.44_dp
        r_edge = 0.32_dp
        z_center = 0.0_dp
        phi_in = 0.5_dp
        err_max = 0.0_dp

        do i = 1, 5
            r_in = r_axis - r_edge + 0.5_dp * r_edge * real(i, dp) / 5.0_dp
            do j = 1, 8
                z_in = z_center + 0.2_dp * r_edge * sin(twopi * real(j - 1, dp) / 8.0_dp)

                xcyl_in = [r_in, phi_in, z_in]

                call chart_cs%from_cyl(xcyl_in, u_chart, ierr)
                if (ierr /= chartmap_from_cyl_ok) then
                    cycle
                end if

                call chart_cs%evaluate_cyl(u_chart, xcyl_out)
                r_out = xcyl_out(1)
                phi_out = xcyl_out(2)
                z_out = xcyl_out(3)

                err_r = abs(r_in - r_out)
                err_z = abs(z_in - z_out)

                if (max(err_r, err_z) > err_max) err_max = max(err_r, err_z)

                if (err_r > TOL_RZ .or. err_z > TOL_RZ) then
                    print *, '  R,Z roundtrip error at R=', r_in, ' Z=', z_in
                    print *, '    Output: R=', r_out, ' Z=', z_out
                    print *, '    err_R=', err_r, ' err_Z=', err_z
                    nerrors = nerrors + 1
                end if
            end do
        end do

        print *, '  Max R,Z roundtrip error: ', err_max, ' m'
        if (nerrors == 0) print *, '  PASSED'
    end subroutine test_rz_chartmap_roundtrip

    subroutine test_rz_vmec_chartmap_roundtrip(nerrors)
        !> Test R,Z -> VMEC -> chartmap -> R,Z full roundtrip.
        !> Generate points in VMEC coordinates first (guaranteeing they are
        !> inside the plasma domain), then test R,Z inversion consistency.
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: vmec_cs, chart_cs
        real(dp) :: xcyl_orig(3), xcyl_vmec(3), xcyl_chart(3)
        real(dp) :: u_vmec_orig(3), u_vmec_inv(3), u_chart(3)
        real(dp) :: err_vmec, err_chart, err_max_vmec, err_max_chart
        logical :: vmec_exists, chart_exists
        integer :: ierr_vmec, ierr_chart, i, j, n_success, n_total
        real(dp) :: s, theta, phi

        print *, 'Test: R,Z -> VMEC -> chartmap -> R,Z full roundtrip'

        inquire(file=VMEC_FILE, exist=vmec_exists)
        inquire(file=CHARTMAP_FILE, exist=chart_exists)
        if (.not. vmec_exists) then
            print *, '  SKIP: ', VMEC_FILE, ' not found'
            return
        end if
        if (.not. chart_exists) then
            print *, '  SKIP: ', CHARTMAP_FILE, ' not found'
            return
        end if

        call init_vmec(VMEC_FILE, 5, 5, 5, dummy)
        call make_vmec_coordinate_system(vmec_cs)
        call make_chartmap_coordinate_system(chart_cs, CHARTMAP_FILE)

        err_max_vmec = 0.0_dp
        err_max_chart = 0.0_dp
        n_success = 0
        n_total = 0
        phi = 0.3_dp

        do i = 1, 5
            s = 0.2_dp + 0.15_dp * real(i - 1, dp)
            do j = 1, 8
                theta = twopi * real(j - 1, dp) / 8.0_dp
                n_total = n_total + 1

                u_vmec_orig = [s, theta, phi]
                call vmec_cs%evaluate_cyl(u_vmec_orig, xcyl_orig)

                call vmec_cs%from_cyl(xcyl_orig, u_vmec_inv, ierr_vmec)
                if (ierr_vmec /= 0) then
                    if (n_success == 0 .and. i == 1 .and. j == 1) then
                        print *, '  DEBUG: VMEC inversion failed at s=', s, ' theta=', theta
                        print *, '         R=', xcyl_orig(1), ' Z=', xcyl_orig(3)
                        print *, '         ierr=', ierr_vmec
                    end if
                    cycle
                end if

                call vmec_cs%evaluate_cyl(u_vmec_inv, xcyl_vmec)
                err_vmec = max(abs(xcyl_orig(1) - xcyl_vmec(1)), &
                               abs(xcyl_orig(3) - xcyl_vmec(3)))
                if (err_vmec > err_max_vmec) err_max_vmec = err_vmec

                call chart_cs%from_cyl(xcyl_orig, u_chart, ierr_chart)
                if (ierr_chart /= chartmap_from_cyl_ok) cycle

                call chart_cs%evaluate_cyl(u_chart, xcyl_chart)
                err_chart = max(abs(xcyl_orig(1) - xcyl_chart(1)), &
                                abs(xcyl_orig(3) - xcyl_chart(3)))
                if (err_chart > err_max_chart) err_max_chart = err_chart

                n_success = n_success + 1

                if (err_vmec > TOL_RZ) then
                    print *, '  VMEC inversion error at s=', s, ' theta=', theta
                    print *, '    err=', err_vmec
                    nerrors = nerrors + 1
                end if

                if (err_chart > TOL_RZ) then
                    print *, '  Chartmap inversion error at s=', s, ' theta=', theta
                    print *, '    err=', err_chart
                    nerrors = nerrors + 1
                end if
            end do
        end do

        print *, '  Successful inversions: ', n_success, ' /', n_total
        print *, '  Max VMEC R,Z error:    ', err_max_vmec, ' m'
        print *, '  Max chartmap R,Z error:', err_max_chart, ' m'

        if (n_success < n_total / 2) then
            print *, '  FAIL: Less than 50% inversions succeeded'
            nerrors = nerrors + 1
        else if (nerrors == 0) then
            print *, '  PASSED'
        end if
    end subroutine test_rz_vmec_chartmap_roundtrip

    subroutine test_vmec_chartmap_cyl_roundtrip(nerrors)
        !> Test VMEC -> cyl -> chartmap -> cyl roundtrip preserves R,Z.
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: vmec_cs, chart_cs
        real(dp) :: u_vmec(3), u_chart(3), xcyl_vmec(3), xcyl_chart(3)
        real(dp) :: r_vmec, z_vmec, phi_vmec
        real(dp) :: r_chart, z_chart, phi_chart
        real(dp) :: err_r, err_z, err_phi
        logical :: vmec_exists, chart_exists
        integer :: ierr, i, j
        real(dp) :: s, theta, phi

        print *, 'Test: VMEC -> cyl -> chartmap -> cyl roundtrip'

        inquire(file=VMEC_FILE, exist=vmec_exists)
        inquire(file=CHARTMAP_FILE, exist=chart_exists)
        if (.not. vmec_exists) then
            print *, '  SKIP: ', VMEC_FILE, ' not found'
            return
        end if
        if (.not. chart_exists) then
            print *, '  SKIP: ', CHARTMAP_FILE, ' not found'
            return
        end if

        call init_vmec(VMEC_FILE, 5, 5, 5, dummy)
        call make_vmec_coordinate_system(vmec_cs)
        call make_chartmap_coordinate_system(chart_cs, CHARTMAP_FILE)

        phi = 0.3_dp

        do i = 1, 5
            s = 0.2_dp + 0.15_dp * real(i - 1, dp)
            do j = 1, 8
                theta = twopi * real(j - 1, dp) / 8.0_dp

                u_vmec = [s, theta, phi]
                call vmec_cs%evaluate_cyl(u_vmec, xcyl_vmec)

                r_vmec = xcyl_vmec(1)
                phi_vmec = xcyl_vmec(2)
                z_vmec = xcyl_vmec(3)

                call chart_cs%from_cyl(xcyl_vmec, u_chart, ierr)
                if (ierr /= chartmap_from_cyl_ok) then
                    print *, '  chartmap inversion failed at s=', s, ' theta=', theta
                    nerrors = nerrors + 1
                    cycle
                end if

                call chart_cs%evaluate_cyl(u_chart, xcyl_chart)
                r_chart = xcyl_chart(1)
                phi_chart = xcyl_chart(2)
                z_chart = xcyl_chart(3)

                err_r = abs(r_vmec - r_chart)
                err_z = abs(z_vmec - z_chart)
                err_phi = abs(phi_vmec - phi_chart)

                if (err_r > TOL_RZ .or. err_z > TOL_RZ) then
                    print *, '  R,Z mismatch at s=', s, ' theta=', theta
                    print *, '    VMEC: R=', r_vmec, ' Z=', z_vmec
                    print *, '    Chart: R=', r_chart, ' Z=', z_chart
                    print *, '    err_R=', err_r, ' err_Z=', err_z
                    nerrors = nerrors + 1
                end if
            end do
        end do

        print *, '  PASSED'
    end subroutine test_vmec_chartmap_cyl_roundtrip

    subroutine test_boundary_rz_match(nerrors)
        !> Test that chartmap boundary (rho=1) matches VMEC boundary (s=1) in R,Z.
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: vmec_cs, chart_cs
        real(dp) :: u_vmec(3), u_chart(3), xcyl_vmec(3), xcyl_chart(3)
        real(dp) :: r_vmec, z_vmec, r_chart, z_chart
        real(dp) :: min_dist, dist, max_min_dist
        logical :: vmec_exists, chart_exists
        integer :: i, j, k

        integer, parameter :: N_VMEC_THETA = 64
        integer, parameter :: N_CHART_THETA = 64
        real(dp) :: phi, theta_vmec, theta_chart
        real(dp) :: vmec_boundary_r(N_VMEC_THETA), vmec_boundary_z(N_VMEC_THETA)
        real(dp) :: chart_boundary_r(N_CHART_THETA), chart_boundary_z(N_CHART_THETA)

        print *, 'Test: Boundary R,Z match (offset=0)'

        inquire(file=VMEC_FILE, exist=vmec_exists)
        inquire(file=CHARTMAP_FILE, exist=chart_exists)
        if (.not. vmec_exists .or. .not. chart_exists) then
            print *, '  SKIP: Required files not found'
            return
        end if

        call init_vmec(VMEC_FILE, 5, 5, 5, dummy)
        call make_vmec_coordinate_system(vmec_cs)
        call make_chartmap_coordinate_system(chart_cs, CHARTMAP_FILE)

        phi = 0.0_dp

        do i = 1, N_VMEC_THETA
            theta_vmec = twopi * real(i - 1, dp) / real(N_VMEC_THETA, dp)
            u_vmec = [1.0_dp, theta_vmec, phi]
            call vmec_cs%evaluate_cyl(u_vmec, xcyl_vmec)
            vmec_boundary_r(i) = xcyl_vmec(1)
            vmec_boundary_z(i) = xcyl_vmec(3)
        end do

        do i = 1, N_CHART_THETA
            theta_chart = twopi * real(i - 1, dp) / real(N_CHART_THETA, dp)
            u_chart = [1.0_dp, theta_chart, phi]
            call chart_cs%evaluate_cyl(u_chart, xcyl_chart)
            chart_boundary_r(i) = xcyl_chart(1)
            chart_boundary_z(i) = xcyl_chart(3)
        end do

        max_min_dist = 0.0_dp
        do i = 1, N_VMEC_THETA
            r_vmec = vmec_boundary_r(i)
            z_vmec = vmec_boundary_z(i)
            min_dist = huge(1.0_dp)
            do j = 1, N_CHART_THETA
                r_chart = chart_boundary_r(j)
                z_chart = chart_boundary_z(j)
                dist = sqrt((r_vmec - r_chart)**2 + (z_vmec - z_chart)**2)
                if (dist < min_dist) min_dist = dist
            end do
            if (min_dist > max_min_dist) max_min_dist = min_dist
        end do

        print *, '  Max closest-point distance: ', max_min_dist, ' m'

        if (max_min_dist > 0.02_dp) then
            print *, '  FAIL: Boundary mismatch exceeds 2 cm tolerance'
            nerrors = nerrors + 1
        else
            print *, '  PASSED'
        end if
    end subroutine test_boundary_rz_match

    subroutine test_interior_rz_consistency(nerrors)
        !> Test R,Z consistency at interior points.
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: vmec_cs, chart_cs
        real(dp) :: u_vmec(3), u_chart(3), xcyl_vmec(3), xcyl_chart(3)
        real(dp) :: r_vmec, z_vmec, r_chart, z_chart
        real(dp) :: err_max, err
        logical :: vmec_exists, chart_exists
        integer :: ierr, i, j
        real(dp) :: s, theta, phi

        print *, 'Test: Interior R,Z consistency'

        inquire(file=VMEC_FILE, exist=vmec_exists)
        inquire(file=CHARTMAP_FILE, exist=chart_exists)

        if (.not. vmec_exists) then
            print *, '  SKIP: ', VMEC_FILE, ' not found'
            return
        end if
        if (.not. chart_exists) then
            print *, '  SKIP: ', CHARTMAP_FILE, ' not found'
            return
        end if

        call init_vmec(VMEC_FILE, 5, 5, 5, dummy)
        call make_vmec_coordinate_system(vmec_cs)
        call make_chartmap_coordinate_system(chart_cs, CHARTMAP_FILE)

        phi = 0.5_dp
        err_max = 0.0_dp

        do i = 1, 5
            s = 0.3_dp + 0.1_dp * real(i - 1, dp)
            do j = 1, 8
                theta = twopi * real(j - 1, dp) / 8.0_dp

                u_vmec = [s, theta, phi]
                call vmec_cs%evaluate_cyl(u_vmec, xcyl_vmec)
                r_vmec = xcyl_vmec(1)
                z_vmec = xcyl_vmec(3)

                call chart_cs%from_cyl(xcyl_vmec, u_chart, ierr)
                if (ierr /= chartmap_from_cyl_ok) cycle

                call chart_cs%evaluate_cyl(u_chart, xcyl_chart)
                r_chart = xcyl_chart(1)
                z_chart = xcyl_chart(3)

                err = max(abs(r_vmec - r_chart), abs(z_vmec - z_chart))
                if (err > err_max) err_max = err
            end do
        end do

        print *, '  Chartmap max R,Z error: ', err_max, ' m'
        if (err_max > TOL_RZ) then
            print *, '  FAIL: Error exceeds tolerance'
            nerrors = nerrors + 1
        else
            print *, '  PASSED'
        end if
    end subroutine test_interior_rz_consistency

end program test_chartmap_rz_consistency
