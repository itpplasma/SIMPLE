program generate_test_boozer_chartmap
    !> Write an analytic circular-tokamak Boozer chartmap for the chartmap tests.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortio_netcdf_classic, only: NC_INT, NC_DOUBLE
    use fortio_netcdf_writer, only: classic_writer_t
    use fortio_status, only: fortio_status_t
    implicit none
    integer, parameter :: nfp = 2, n_rho = 17, n_theta = 33, n_zeta = 33
    real(dp), parameter :: b0 = 2.0e4_dp, r0 = 150.0_dp, a = 50.0_dp, iota = 0.5_dp
    real(dp), parameter :: rho_min = 1.0e-3_dp, epsilon = a/r0
    real(dp), parameter :: torflux = 0.5_dp*b0*a*a
    real(dp) :: twopi, rho(n_rho), s(n_rho), theta(n_theta), zeta(n_zeta)
    real(dp) :: a_phi(n_rho), b_theta(n_rho), b_phi(n_rho)
    real(dp), dimension(n_rho, n_theta, n_zeta) :: x, y, z, bmod
    type(classic_writer_t) :: nc
    type(fortio_status_t) :: st
    integer :: d_rho, d_s, d_theta, d_zeta, v(12), ir, it, iz
    character(len=1024) :: filename
    real(dp) :: radius

    filename = 'test_boozer_chartmap.nc'
    if (command_argument_count() > 0) call get_command_argument(1, filename)
    twopi = 2.0_dp*acos(-1.0_dp)
    do ir = 1, n_rho
        rho(ir) = rho_min + (1.0_dp - rho_min)*real(ir - 1, dp)/real(n_rho - 1, dp)
        s(ir) = rho_min**2 + (1.0_dp - rho_min**2)*real(ir - 1, dp)/real(n_rho - 1, dp)
    end do
    a_phi = -torflux*iota*s
    b_theta = b0*epsilon*iota
    b_phi = b0*r0
    theta = twopi*[(real(it, dp), it=0, n_theta - 1)]/n_theta
    zeta = twopi/nfp*[(real(iz, dp), iz=0, n_zeta - 1)]/n_zeta
    ! Geometry and Bmod on the endpoint-excluded angular grid.
    do iz = 1, n_zeta
        do it = 1, n_theta
            do ir = 1, n_rho
                radius = r0 + a*rho(ir)*cos(theta(it))
                x(ir, it, iz) = radius*cos(zeta(iz))
                y(ir, it, iz) = radius*sin(zeta(iz))
                z(ir, it, iz) = a*rho(ir)*sin(theta(it))
                bmod(ir, it, iz) = b0/(1.0_dp + rho(ir)*epsilon*cos(theta(it)))
            end do
        end do
    end do

    call nc%create(trim(filename), st); call check('create')
    call nc%define_dimension('rho', int(n_rho, int64), .false., d_rho, st); call check('dim rho')
    call nc%define_dimension('s', int(n_rho, int64), .false., d_s, st); call check('dim s')
    call nc%define_dimension('theta', int(n_theta, int64), .false., d_theta, st)
    call check('dim theta')
    call nc%define_dimension('zeta', int(n_zeta, int64), .false., d_zeta, st)
    call check('dim zeta')
    call nc%define_variable('rho', NC_DOUBLE, [d_rho], v(1), st); call check('rho')
    call nc%define_variable('s', NC_DOUBLE, [d_s], v(2), st); call check('s')
    call nc%define_variable('theta', NC_DOUBLE, [d_theta], v(3), st); call check('theta')
    call nc%define_variable('zeta', NC_DOUBLE, [d_zeta], v(4), st); call check('zeta')
    call nc%define_variable('x', NC_DOUBLE, [d_rho, d_theta, d_zeta], v(5), st); call check('x')
    call nc%define_variable('y', NC_DOUBLE, [d_rho, d_theta, d_zeta], v(6), st); call check('y')
    call nc%define_variable('z', NC_DOUBLE, [d_rho, d_theta, d_zeta], v(7), st); call check('z')
    call nc%define_variable('A_phi', NC_DOUBLE, [d_s], v(8), st); call check('A_phi')
    call nc%define_variable('B_theta', NC_DOUBLE, [d_rho], v(9), st); call check('B_theta')
    call nc%define_variable('B_phi', NC_DOUBLE, [d_rho], v(10), st); call check('B_phi')
    call nc%define_variable('Bmod', NC_DOUBLE, [d_rho, d_theta, d_zeta], v(11), st)
    call check('Bmod')
    call nc%define_variable('num_field_periods', NC_INT, [integer ::], v(12), st)
    call check('num_field_periods')
    do it = 5, 7
        call nc%put_attribute_text(v(it), 'units', 'cm', st); call check('units')
    end do
    call nc%put_attribute_text(v(8), 'radial_abscissa', 's', st); call check('radial_abscissa')
    call nc%put_attribute_text(-1, 'rho_convention', 'rho_tor', st); call check('rho_convention')
    call nc%put_attribute_text(-1, 'zeta_convention', 'boozer', st); call check('zeta_convention')
    call nc%put_attribute_r64(-1, 'rho_lcfs', [1.0_dp], st); call check('rho_lcfs')
    call nc%put_attribute_i32(-1, 'boozer_field', [1], st); call check('boozer_field')
    call nc%put_attribute_r64(-1, 'torflux', [torflux], st); call check('torflux')
    call nc%end_definition(st); call check('enddef')
    call nc%put_r64_1(v(1), rho, st); call check('put rho')
    call nc%put_r64_1(v(2), s, st); call check('put s')
    call nc%put_r64_1(v(3), theta, st); call check('put theta')
    call nc%put_r64_1(v(4), zeta, st); call check('put zeta')
    call nc%put_r64_3(v(5), x, st); call check('put x')
    call nc%put_r64_3(v(6), y, st); call check('put y')
    call nc%put_r64_3(v(7), z, st); call check('put z')
    call nc%put_r64_1(v(8), a_phi, st); call check('put A_phi')
    call nc%put_r64_1(v(9), b_theta, st); call check('put B_theta')
    call nc%put_r64_1(v(10), b_phi, st); call check('put B_phi')
    call nc%put_r64_3(v(11), bmod, st); call check('put Bmod')
    call nc%put_i32_scalar(v(12), nfp, st); call check('put num_field_periods')
    call nc%close_netcdf4(st); call check('close')
    print '(a,a)', 'Generated ', trim(filename)

contains

    subroutine check(context)
        character(len=*), intent(in) :: context
        if (.not. st%ok()) then
            print '(a,a,a,a)', 'NetCDF error in ', context, ': ', st%message
            error stop 1
        end if
    end subroutine check

end program generate_test_boozer_chartmap
