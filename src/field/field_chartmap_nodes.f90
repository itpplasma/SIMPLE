module field_chartmap_nodes
    !> A field-on-chartmap file (libneo_chartmap_field: physical cylindrical A
    !> and B at the nodes of a chartmap) as a splined field on that chart. The
    !> covariant components are formed at the file's own nodes and splined
    !> there, without resampling, so coord_input must be the same chart.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_splined, only: splined_field_t, create_splined_field_from_nodes
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  PSI_TOR_NORM, PSI_POL_NORM
    use libneo_chartmap_field, only: chartmap_field_t, read_chartmap_field

    implicit none
    private

    public :: create_splined_field_from_chartmap

    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

contains

    subroutine create_splined_field_from_chartmap(path, ref_coords, field)
        character(*), intent(in) :: path
        class(coordinate_system_t), intent(in) :: ref_coords
        type(splined_field_t), intent(out) :: field

        type(chartmap_field_t) :: f
        real(dp), allocatable :: values(:, :, :, :)
        integer :: nr, nt, nz

        call read_chartmap_field(path, f)
        nr = size(f%rho)
        nt = size(f%theta)
        nz = size(f%zeta)
        call check_same_chart(path, f, ref_coords)
        allocate (values(nr, nt + 1, nz + 1, 7))
        call nodal_components(f, ref_coords, values(:, :nt, :nz, :))
        values(:, nt + 1, :nz, :) = values(:, 1, :nz, :)
        values(:, :, nz + 1, :) = values(:, :, 1, :)
        call create_splined_field_from_nodes(ref_coords, &
                                             [f%rho(1), f%theta(1), f%zeta(1)], &
                                             [f%rho(nr), f%theta(1) + twopi, &
                                              f%zeta(1) + twopi/f%num_field_periods], &
                                             values, field)
    end subroutine create_splined_field_from_chartmap

    !> A_i = A . e_i; the direction is stored as (B_R, B^phi, B_Z) like the
    !> sampled Cartesian sources (field_splined, regular_cartesian_direction).
    subroutine nodal_components(f, chart, values)
        type(chartmap_field_t), intent(in) :: f
        class(coordinate_system_t), intent(in) :: chart
        real(dp), intent(out) :: values(:, :, :, :)

        real(dp) :: u(3), x(3), e(3, 3), a(3), r, c, s
        integer :: i, j, k

        !$omp parallel do collapse(2) private(i, u, x, e, a, r, c, s)
        do k = 1, size(f%zeta)
            do j = 1, size(f%theta)
                do i = 1, size(f%rho)
                    u = [f%rho(i), f%theta(j), f%zeta(k)]
                    call chart%evaluate_cart(u, x)
                    call chart%covariant_basis(u, e)
                    r = hypot(x(1), x(2))
                    c = x(1)/r
                    s = x(2)/r
                    a = [c*f%a(1, i, j, k) - s*f%a(2, i, j, k), &
                         s*f%a(1, i, j, k) + c*f%a(2, i, j, k), f%a(3, i, j, k)]
                    values(i, j, k, 1:3) = matmul(a, e)
                    values(i, j, k, 4) = f%b(1, i, j, k)
                    values(i, j, k, 5) = f%b(2, i, j, k)/r
                    values(i, j, k, 6) = f%b(3, i, j, k)
                    values(i, j, k, 7) = norm2(f%b(:, i, j, k))
                end do
            end do
        end do
        !$omp end parallel do
        do i = 1, 3
            values(:, :, :, i) = values(:, :, :, i) - values(1, 1, 1, i)
        end do
    end subroutine nodal_components

    subroutine check_same_chart(path, f, chart)
        character(*), intent(in) :: path
        type(chartmap_field_t), intent(in) :: f
        class(coordinate_system_t), intent(in) :: chart

        select type (cs => chart)
        class is (chartmap_coordinate_system_t)
            if (cs%nrho /= size(f%rho) .or. cs%ntheta /= size(f%theta) .or. &
                cs%nzeta /= size(f%zeta) .or. &
                cs%num_field_periods /= f%num_field_periods) then
                print *, 'field_from_file: the grid of ', trim(path), &
                    ' differs from the coord_input chart; use the same chart'
                error stop
            end if
            if (cs%rho_convention == PSI_TOR_NORM .or. &
                cs%rho_convention == PSI_POL_NORM) then
                print *, 'field_from_file: ', trim(path), ' has a flux radius; ', &
                    'field-on-chartmap input needs a rho or geometric chart'
                error stop
            end if
        class default
            print *, 'field_from_file: ', trim(path), ' is a field on a chartmap; ', &
                'set coord_input to the same chartmap'
            error stop
        end select
    end subroutine check_same_chart

end module field_chartmap_nodes
