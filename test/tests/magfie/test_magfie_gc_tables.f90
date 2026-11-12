program test_magfie_gc_tables
    ! Guiding-centre tables of the NCSX coil field on its chartmap against the
    ! exact quantities at random points between the nodes, on the chart grid
    ! and on a twice subdivided one, and against the splined reference route.
    ! The chart is a cubic spline, so its tangent basis has kinks at the chart
    ! knots and the tables converge at least at second order there.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi
    use field_coils, only: coils_field_t, create_coils_field
    use field_gc_tables, only: gc_tables_t, build_gc_tables_coils, evaluate_gc_tables
    use field_splined, only: splined_field_t, create_splined_field
    use magfie_sub, only: magfie_splined_refcoords, gc_tables, init_magfie, &
                          REFCOORDS, axis_radial_coordinate
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system, RHO_TOR, UNKNOWN
    use neo_biotsavart, only: compute_field_jacobian
    use neo_guiding_centre_field, only: guiding_centre_quantities, GC_COUNT
    implicit none

    integer, parameter :: npts = 200
    type(coils_field_t) :: raw
    class(coordinate_system_t), allocatable :: chart
    type(gc_tables_t) :: coarse, fine
    type(splined_field_t) :: splined
    real(dp) :: err_route(GC_COUNT)
    logical :: radial_is_rho
    integer :: status
    real(dp) :: u(3, npts), exact(GC_COUNT, npts), err_coarse(GC_COUNT), &
                err_fine(GC_COUNT), period
    integer :: k
    call make_chartmap_coordinate_system(chart, 'wout_ncsx.chartmap.nc')
    select type (chart)
    class is (chartmap_coordinate_system_t)
        chart%rho_convention = RHO_TOR
        period = twopi/real(chart%num_field_periods, dp)
    end select
    call create_coils_field('coils.simple', raw)
    call random_seed(put=[(20261108 + k, k=1, 64)])
    call random_number(u)
    u(1, :) = 0.02_dp + 0.9_dp*u(1, :)
    u(2, :) = twopi*u(2, :)
    u(3, :) = period*u(3, :)
    do k = 1, npts
        call exact_quantities(u(:, k), exact(:, k))
    end do
    ! With the tables installed, the RK axis map must use the chart's rho.
    call init_magfie(REFCOORDS)
    call build_gc_tables_coils(raw%coils, chart, 3, gc_tables, n_r=17, n_th=33, &
                               n_phi=65)
    call axis_radial_coordinate(radial_is_rho, status)
    if (status /= 0 .or. .not. radial_is_rho) error stop &
        'tables must use the rho axis map'
    call build_gc_tables_coils(raw%coils, chart, 3, coarse)
    call build_gc_tables_coils(raw%coils, chart, 3, fine, n_r=33, n_th=65, n_phi=129)
    err_coarse = table_error(coarse)
    err_fine = table_error(fine)
    call create_splined_field(raw, chart, splined)
    err_route = route_error()
    print '(a)', ' quantity  tables chart grid  tables 2x   splined route'
    do k = 1, GC_COUNT
        print '(i9, 3es12.3)', k, err_coarse(k), err_fine(k), err_route(k)
    end do
    ! A chart with a geometric radius (convention UNKNOWN, as map2disc wall
    ! charts) is tabulated like rho_tor.
    select type (chart)
    class is (chartmap_coordinate_system_t)
        chart%rho_convention = UNKNOWN
    end select
    call build_gc_tables_coils(raw%coils, chart, 3, coarse)
    if (any(abs(table_error(coarse) - err_coarse) > 0.0_dp)) error stop &
        'geometric-radius chart tables differ'
    if (minval(err_coarse/err_fine) < 2.5_dp) error stop 'tables do not converge'
    if (maxval(err_fine) > 5.0e-3_dp) error stop 'tables too far from the exact field'
    print *, 'PASSED: guiding-centre tables converge to the exact coil field'

contains

    subroutine exact_quantities(u, q)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: q(GC_COUNT)
        real(dp) :: x(3), e(3, 3), a(3), da(3, 3), b(3), db(3, 3)

        call chart%evaluate_cart(u, x)
        call chart%covariant_basis(u, e)
        call compute_field_jacobian(raw%coils, x, a, da, b, db)
        call guiding_centre_quantities(e, b, db, q)
    end subroutine exact_quantities

    function table_error(tab) result(err)
        type(gc_tables_t), intent(in) :: tab
        real(dp) :: err(GC_COUNT), q(GC_COUNT)
        integer :: i

        err = 0
        do i = 1, npts
            call evaluate_gc_tables(tab, u(:, i), q(1), q(2), q(3:5), q(6:8), &
                                    q(9:11), q(12:14))
            ! compare the regular products, as tabulated
            q(9:14) = q(2)*q(9:14)
            err = max(err, abs(q - regular(exact(:, i)))/regular_scale())
        end do
    end function table_error

    function route_error() result(err)
        real(dp) :: err(GC_COUNT), q(GC_COUNT)
        integer :: i

        err = 0
        do i = 1, npts
            call magfie_splined_refcoords(splined, u(:, i), q(1), q(2), q(3:5), &
                                          q(6:8), q(9:11), q(12:14))
            q(9:14) = q(2)*q(9:14)
            err = max(err, abs(q - regular(exact(:, i)))/regular_scale())
        end do
    end function route_error

    pure function regular(q) result(r)
        real(dp), intent(in) :: q(GC_COUNT)
        real(dp) :: r(GC_COUNT)

        r = q
        r(9:14) = q(2)*q(9:14)
    end function regular

    function regular_scale() result(s)
        real(dp) :: s(GC_COUNT)
        integer :: i

        s = 0
        do i = 1, npts
            s = max(s, abs(regular(exact(:, i))))
        end do
    end function regular_scale

end program test_magfie_gc_tables
