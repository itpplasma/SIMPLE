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
              build_gc_tables_coils, evaluate_gc_tables

    type :: gc_tables_t
        type(BatchSplineData3D) :: spl
        logical :: initialized = .false.
    contains
        final :: gc_tables_cleanup
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
    subroutine build_gc_tables_coils(coils, coords, order, tab, n_r, n_th, n_phi)
        type(coils_t), intent(in) :: coils
        class(coordinate_system_t), intent(in) :: coords
        integer, intent(in) :: order
        type(gc_tables_t), intent(inout) :: tab
        integer, intent(in), optional :: n_r, n_th, n_phi

        call build_gc_tables(coil_source_t(coils), coords, order, tab, n_r, n_th, n_phi)
    end subroutine build_gc_tables_coils

    subroutine coil_field_jacobian(self, x, b, db)
        class(coil_source_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3), db(3, 3)
        real(dp) :: a(3), da(3, 3)

        call compute_field_jacobian(self%coils, x, a, da, b, db)
    end subroutine coil_field_jacobian

    subroutine build_gc_tables(source, coords, order, tab, n_r, n_th, n_phi)
        class(gc_field_source_t), intent(in) :: source
        class(coordinate_system_t), intent(in) :: coords
        integer, intent(in) :: order
        type(gc_tables_t), intent(inout) :: tab
        integer, intent(in), optional :: n_r, n_th, n_phi

        real(dp), allocatable :: y(:, :, :, :)
        real(dp) :: xmin(3), xmax(3), step(3), u(3), x(3), e(3, 3), q(GC_COUNT)
        real(dp) :: b(3), db(3, 3)
        integer :: nr, nth, nph, i, j, k

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
        allocate (y(nr, nth, nph, GC_COUNT))
        !$omp parallel do collapse(2) default(shared) &
        !$omp private(i, j, k, u, x, e, b, db, q)
        do k = 1, nph
            do j = 1, nth
                do i = 1, nr
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
        if (tab%initialized) call destroy_batch_splines_3d(tab%spl)
        call construct_batch_splines_3d(xmin, xmax, y, [order, order, order], &
                                        [.false., .true., .true.], tab%spl)
        tab%initialized = .true.
    end subroutine build_gc_tables

    !> magfie outputs at x = (rho, theta, zeta) of the chart.
    subroutine evaluate_gc_tables(tab, x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        type(gc_tables_t), intent(in) :: tab
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: y(GC_COUNT)

        call evaluate_batch_splines_3d(tab%spl, x, y)
        bmod = y(GC_BMOD)
        sqrtg = y(GC_SQRTG)
        bder = y(GC_BDER:GC_BDER + 2)
        hcovar = y(GC_HCOV:GC_HCOV + 2)
        hctrvr = y(GC_HCTR:GC_HCTR + 2)/sqrtg
        hcurl = y(GC_HCURL:GC_HCURL + 2)/sqrtg
    end subroutine evaluate_gc_tables

    subroutine gc_tables_cleanup(tab)
        type(gc_tables_t), intent(inout) :: tab

        if (tab%initialized) call destroy_batch_splines_3d(tab%spl)
        tab%initialized = .false.
    end subroutine gc_tables_cleanup

end module field_gc_tables
