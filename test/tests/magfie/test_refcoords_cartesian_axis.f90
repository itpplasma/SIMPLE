module analytic_cartesian_axis_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    implicit none
    type, extends(magnetic_field_t) :: analytic_field_t
        real(dp) :: current = 0.0_dp
    contains
        procedure :: evaluate => analytic_evaluate
    end type analytic_field_t
contains
    subroutine analytic_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        class(analytic_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)
        real(dp) :: b(3)

        ! curl(A) = (-c*y, c*x, 1); c=0 is a uniform vacuum field.
        if (.not. allocated(self%coords)) error stop 'missing source coordinates'
        Acov = [-0.5_dp*x(2), 0.5_dp*x(1), &
                -0.5_dp*self%current*(x(1)**2 + x(2)**2)]
        b = [-self%current*x(2), self%current*x(1), 1.0_dp]
        Bmod = sqrt(sum(b*b))
        hcov = b/Bmod
        if (present(sqgBctr)) sqgBctr = b
    end subroutine analytic_evaluate
end module analytic_cartesian_axis_field

program test_refcoords_cartesian_axis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use analytic_cartesian_axis_field, only: analytic_field_t
    use cartesian_coordinates, only: cartesian_coordinate_system_t
    use simple, only: init_vmec
    use new_vmec_stuff_mod, only: nper
    use libneo_coordinates, only: chartmap_coordinate_system_t
    use reference_coordinates, only: ref_coords, init_reference_coordinates
    use field_splined, only: splined_field_t, create_splined_field
    use magfie_sub, only: init_magfie, magfie, REFCOORDS, set_magfie_refcoords_field
    implicit none
    real(dp), parameter :: radii(5) = [0.05_dp, 0.01_dp, 0.002_dp, &
                                       0.001_dp, 0.0002_dp]
    type(analytic_field_t) :: source
    type(splined_field_t) :: field
    real(dp) :: dummy, x(3), u(3), xyz(3), basis(3, 3), hcart(3)
    real(dp) :: bmod, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)
    real(dp) :: a(3), direct(3), bsource, direction_error, curl_error, expected
    integer :: geometry, kind, i, failures

    ! Public NCSX geometry only. Coarse field sampling must not turn a smooth
    ! analytic Cartesian direction into a singular vector near the axis.
    call init_vmec('wout_ncsx.nc', 5, 5, 1, dummy)
    allocate (cartesian_coordinate_system_t :: source%coords)
    failures = 0
    do geometry = 1, 2
        if (geometry == 1) then
            call init_reference_coordinates('wout_ncsx.nc')
        else
            call init_reference_coordinates('wout_ncsx.chartmap.nc')
        end if
        do kind = 1, 2
            source%current = real(kind - 1, dp)*1.0e-3_dp
         call create_splined_field(source, ref_coords, field, n_r=16, n_th=33, n_phi=32)
            call set_magfie_refcoords_field(field)
            call init_magfie(REFCOORDS)
            do i = 1, size(radii)
                x = [radii(i), 0.23_dp + 0.37_dp*i, &
                     0.12_dp + (0.07_dp*i + i - 3)*(2.0_dp*acos(-1.0_dp)/nper)]
                u = x
                select type (coords => ref_coords)
                class is (chartmap_coordinate_system_t)
                    call coords%evaluate_cart(u, xyz)
                    call coords%covariant_basis(u, basis)
                class default
                    u(1) = x(1)**2
                    call coords%evaluate_cart(u, xyz)
                    call coords%covariant_basis(u, basis)
                    basis(:, 1) = basis(:, 1)*2.0_dp*x(1)
                end select
                call source%evaluate(xyz, a, direct, bsource)
                call magfie(x, bmod, sqrtg, bder, hcov, hctr, hcurl)
                hcart = matmul(basis, hctr)
                direction_error = maxval(abs(hcart - direct))
                ! h.curl(h)=B.curl(B)/|B|^2; includes a non-vacuum case so
                ! replacing the parallel curl by zero cannot satisfy this oracle.
                expected = 2.0_dp*source%current/bsource**2
                curl_error = abs(dot_product(hcov, hcurl) - expected)
print *, 'geometry',geometry,'case',kind,'rho',x(1),'direction_error',direction_error, &
                 'parallel_curl_error', curl_error, 'unit_norm', dot_product(hcov, hctr)
                call check_derivative_contract(field, x, failures)
                if (direction_error > 1.0e-6_dp) failures = failures + 1
                if (curl_error > 2.0e-7_dp) failures = failures + 1
          if (abs(dot_product(hcov, hctr) - 1.0_dp) > 2.0e-6_dp) failures = failures + 1
            end do
        end do
    end do
    if (failures /= 0) error stop 'Cartesian direction lost near-axis regularity'
contains
    subroutine check_derivative_contract(f, point, failures)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: point(3)
        integer, intent(inout) :: failures
        real(dp) :: ac(3), hc(3), bm, da(3, 3), dh(3, 3), db(3)
        real(dp) :: ap(3), am(3), hp(3), hm(3), bp, bm2, shifted(3), step, fd(3)
        integer :: j
        call f%evaluate_with_der(point, ac, hc, bm, da, dh, db)
        do j = 1, 3
            step = 1.0e-5_dp
            if (j == 1) step = min(step, 0.01_dp*point(1))
            shifted = point; shifted(j) = point(j) + step
            call f%evaluate(shifted, ap, hp, bp)
            shifted(j) = point(j) - step
            call f%evaluate(shifted, am, hm, bm2)
            fd = (hp - hm)/(2*step)
            if (maxval(abs(dh(j, :) - fd)) > 1.0e-6_dp*max(1.0_dp, maxval(abs(fd)))) &
                failures = failures + 1
        end do
    end subroutine check_derivative_contract
end program test_refcoords_cartesian_axis
