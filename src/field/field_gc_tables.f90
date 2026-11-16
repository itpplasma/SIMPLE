module field_gc_tables
    !> Guiding-centre field tables on a rho chartmap. The chart nodes hold |B|,
    !> sqrt(g), d_j log|B|, h_j, sqrt(g) h^j and sqrt(g) (curl h)^j, evaluated
    !> exactly from the Cartesian field and its Jacobian (libneo
    !> guiding_centre_quantities); a field call is one batch-spline evaluation.
    !> The contravariant parts are stored times sqrt(g), which keeps them
    !> regular at the axis, and divided by the interpolated sqrt(g) on use.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                           evaluate_batch_splines_3d, destroy_batch_splines_3d
    use libneo_coordinates, only: coordinate_system_t, &
                                  chartmap_coordinate_system_t, PSI_TOR_NORM, &
                                  PSI_POL_NORM
    use neo_biotsavart, only: coils_t, compute_field_jacobian
    use neo_guiding_centre_field, only: guiding_centre_quantities, GC_COUNT, &
                                  GC_BMOD, GC_SQRTG, GC_BDER, GC_HCOV, GC_HCTR, GC_HCURL
    use field_splined, only: configure_spline_grid

    implicit none
    private

    public :: gc_tables_t, gc_field_source_t, build_gc_tables, &
              build_gc_tables_coils, evaluate_gc_tables, clear_gc_tables

    ! Cubic, the order the tables were validated and timed with.
    integer, parameter :: gc_tables_order = 3
    ! sqrt(g) -> 0 at the axis, and the contravariant parts are divided by it.
    ! Interpolated across the first radial cells it erred by O(h**2) relative
    ! to a vanishing value, up to the wrong sign (#799). Within axis_cells
    ! cells, a table axis_refine times finer in rho and the chart's sqrt(g)
    ! are used instead.
    integer, parameter :: axis_cells = 2, axis_refine = 4

    type :: gc_tables_t
        type(BatchSplineData3D) :: spl
        !> The same quantities, axis_refine times finer in rho, below axis_rho.
        type(BatchSplineData3D) :: axis_spl
        real(dp) :: axis_rho = -huge(1.0_dp)
        !> The chart, for sqrt(g) below axis_rho.
        class(coordinate_system_t), allocatable :: coords
        logical :: initialized = .false.
    contains
        final :: clear_gc_tables
    end type gc_tables_t

    !> A Cartesian field B(x) with its Jacobian dB_i/dx_j (cm, G). Sources
    !> outside this module (field_gc_tables_input) extend it, so that the
    !> tracing modules using the tables never see their types.
    type, abstract :: gc_field_source_t
    contains
        procedure(field_jacobian_i), deferred :: field_jacobian
    end type gc_field_source_t

    abstract interface
        subroutine field_jacobian_i(self, x, b, db)
            import :: dp, gc_field_source_t
            class(gc_field_source_t), intent(in) :: self
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: b(3), db(3, 3)
        end subroutine field_jacobian_i
    end interface

    type, extends(gc_field_source_t) :: coil_source_t
        type(coils_t) :: coils
    contains
        procedure :: field_jacobian => coil_field_jacobian
    end type coil_source_t

contains

    !> Tables of the coil field (coordinates in cm, currents in statA, as
    !> create_coils_field) on the grid of the chartmap coords.
    subroutine build_gc_tables_coils(coils, coords, tab, n_r, n_th, n_phi)
        type(coils_t), intent(in) :: coils
        class(coordinate_system_t), intent(in) :: coords
        type(gc_tables_t), intent(inout) :: tab
        integer, intent(in), optional :: n_r, n_th, n_phi

        call build_gc_tables(coil_source_t(coils), coords, tab, n_r, n_th, n_phi)
    end subroutine build_gc_tables_coils

    subroutine coil_field_jacobian(self, x, b, db)
        class(coil_source_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3), db(3, 3)
        real(dp) :: a(3), da(3, 3)

        call compute_field_jacobian(self%coils, x, a, da, b, db)
    end subroutine coil_field_jacobian

    subroutine build_gc_tables(source, coords, tab, n_r, n_th, n_phi)
        class(gc_field_source_t), intent(in) :: source
        class(coordinate_system_t), intent(in) :: coords
        type(gc_tables_t), intent(inout) :: tab
        integer, intent(in), optional :: n_r, n_th, n_phi

        real(dp) :: xmin(3), xmax(3), step(3), axis_max(3)
        integer :: nr, nth, nph

        ! The tables are splined in the chart's own radius, which the RK axis
        ! map treats as rho (magfie axis_radial_coordinate); a normalised-flux
        ! radius would need the sqrt(s) scaling.
        select type (coords)
        class is (chartmap_coordinate_system_t)
            if (coords%rho_convention == PSI_TOR_NORM .or. &
                coords%rho_convention == PSI_POL_NORM) error stop &
                'build_gc_tables: the chartmap radius must not be normalised flux'
        class default
            error stop 'build_gc_tables: needs a chartmap reference chart'
        end select
        call configure_spline_grid(coords, n_r_in=n_r, n_th_in=n_th, &
                                   n_phi_in=n_phi, n_r=nr, n_th=nth, n_phi=nph, &
                                   xmin=xmin, xmax=xmax)
        step = (xmax - xmin)/real([nr, nth, nph] - 1, dp)
        call clear_gc_tables(tab)
        allocate (tab%coords, source=coords)
        call tabulate(source, coords, xmin, xmax, [nr, nth, nph], tab%spl)
        tab%axis_rho = xmin(1) + axis_cells*step(1)
        axis_max = xmax
        axis_max(1) = xmin(1) + (axis_cells + 1)*step(1)
        call tabulate(source, coords, xmin, axis_max, &
                      [axis_refine*(axis_cells + 1) + 1, nth, nph], tab%axis_spl)
        tab%initialized = .true.
    end subroutine build_gc_tables

    subroutine tabulate(source, coords, xmin, xmax, n, spl)
        class(gc_field_source_t), intent(in) :: source
        class(coordinate_system_t), intent(in) :: coords
        real(dp), intent(in) :: xmin(3), xmax(3)
        integer, intent(in) :: n(3)
        type(BatchSplineData3D), intent(out) :: spl

        real(dp), allocatable :: y(:, :, :, :)
        real(dp) :: step(3), u(3), x(3), e(3, 3), q(GC_COUNT), b(3), db(3, 3)
        integer :: i, j, k

        step = (xmax - xmin)/real(n - 1, dp)
        allocate (y(n(1), n(2), n(3), GC_COUNT))
        !$omp parallel do collapse(2) default(shared) &
        !$omp private(i, j, k, u, x, e, b, db, q)
        do k = 1, n(3)
            do j = 1, n(2)
                do i = 1, n(1)
                    u = xmin + step*real([i, j, k] - 1, dp)
                    call coords%evaluate_cart(u, x)
                    call coords%covariant_basis(u, e)
                    call source%field_jacobian(x, b, db)
                    call guiding_centre_quantities(e, b, db, q)
                    q(GC_HCTR:GC_HCTR + 2) = q(GC_SQRTG)*q(GC_HCTR:GC_HCTR + 2)
                    q(GC_HCURL:GC_HCURL + 2) = q(GC_SQRTG)*q(GC_HCURL:GC_HCURL + 2)
                    y(i, j, k, :) = q
                end do
            end do
        end do
        !$omp end parallel do
        call construct_batch_splines_3d(xmin, xmax, y, spread(gc_tables_order, 1, 3), &
                                        [.false., .true., .true.], spl)
    end subroutine tabulate

    !> magfie outputs at x = (rho, theta, zeta) of the chart.
    subroutine evaluate_gc_tables(tab, x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        type(gc_tables_t), intent(in) :: tab
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: y(GC_COUNT), e(3, 3)

        if (x(1) < tab%axis_rho) then
            call evaluate_batch_splines_3d(tab%axis_spl, x, y)
            call tab%coords%covariant_basis(x, e)
            y(GC_SQRTG) = dot_product(e(:, 1), [e(2, 2)*e(3, 3) - e(3, 2)*e(2, 3), &
                                                e(3, 2)*e(1, 3) - e(1, 2)*e(3, 3), &
                                                e(1, 2)*e(2, 3) - e(2, 2)*e(1, 3)])
        else
            call evaluate_batch_splines_3d(tab%spl, x, y)
        end if
        bmod = y(GC_BMOD)
        sqrtg = y(GC_SQRTG)
        bder = y(GC_BDER:GC_BDER + 2)
        hcovar = y(GC_HCOV:GC_HCOV + 2)
        hctrvr = y(GC_HCTR:GC_HCTR + 2)/sqrtg
        hcurl = y(GC_HCURL:GC_HCURL + 2)/sqrtg
    end subroutine evaluate_gc_tables

    subroutine clear_gc_tables(tab)
        type(gc_tables_t), intent(inout) :: tab

        if (tab%initialized) then
            call destroy_batch_splines_3d(tab%spl)
            call destroy_batch_splines_3d(tab%axis_spl)
        end if
        if (allocated(tab%coords)) deallocate (tab%coords)
        tab%initialized = .false.
    end subroutine clear_gc_tables

end module field_gc_tables
