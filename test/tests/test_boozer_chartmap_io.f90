program test_boozer_chartmap_io
    !> Unit test for read_boozer_chartmap: verifies it parses Bmod on the
    !> endpoint-included field grid (theta_field/zeta_field), recovers the full
    !> periodic spans 2*pi (poloidal) and 2*pi/nfp (toroidal), and reads the
    !> optional rmajor attribute (with a clean fallback when absent).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use boozer_chartmap_io, only: boozer_chartmap_data_t, read_boozer_chartmap

    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: twopi = 2.0_dp * pi
    real(dp), parameter :: tol = 1.0e-12_dp
    integer, parameter :: n_rho = 3, n_theta_geom = 4, n_phi_geom = 4
    integer, parameter :: n_theta_field = 5, n_phi_field = 5, nfp = 2
    real(dp), parameter :: rmajor_val = 7.5_dp, torflux_val = 0.7_dp

    character(len=*), parameter :: file_new = 'io_unit_new.nc'
    character(len=*), parameter :: file_old = 'io_unit_old.nc'
    type(boozer_chartmap_data_t) :: d
    integer :: nfail

    nfail = 0

    ! --- Real format: separate field grid + rmajor -----------------------
    call write_chartmap(file_new, with_field_dims=.true., with_rmajor=.true.)
    call read_boozer_chartmap(file_new, d)

    call expect_int(d%n_rho, n_rho, 'n_rho', nfail)
    call expect_int(d%n_theta, n_theta_field, 'n_theta (field grid)', nfail)
    call expect_int(d%n_phi, n_phi_field, 'n_phi (field grid)', nfail)
    call expect_int(d%nfp, nfp, 'nfp', nfail)
    call expect_real(d%torflux, torflux_val, 'torflux', nfail)
    if (.not. d%has_rmajor) then
        print *, 'FAIL: has_rmajor false for file with rmajor'
        nfail = nfail + 1
    end if
    call expect_real(d%rmajor, rmajor_val, 'rmajor', nfail)
    call expect_real(real(d%n_theta - 1, dp)*d%h_theta, twopi, &
        'poloidal period (n_theta-1)*h_theta', nfail)
    call expect_real(real(d%n_phi - 1, dp)*d%h_phi, twopi/real(nfp, dp), &
        'toroidal period (n_phi-1)*h_phi', nfail)
    ! Bmod must round-trip on the field grid.
    call expect_real(d%Bmod(2, 3, 4), bmod_model(2, 3, 4), 'Bmod(2,3,4)', nfail)
    call expect_real(d%Bmod(n_rho, n_theta_field, n_phi_field), &
        bmod_model(n_rho, n_theta_field, n_phi_field), 'Bmod(end)', nfail)
    call expect_real(d%A_phi(2), 2.0_dp, 'A_phi(2)', nfail)

    ! --- Legacy format: no field dims, no rmajor -------------------------
    call write_chartmap(file_old, with_field_dims=.false., with_rmajor=.false.)
    call read_boozer_chartmap(file_old, d)
    call expect_int(d%n_theta, n_theta_geom, 'legacy n_theta falls back to geometry', nfail)
    call expect_int(d%n_phi, n_phi_geom, 'legacy n_phi falls back to geometry', nfail)
    if (d%has_rmajor) then
        print *, 'FAIL: has_rmajor true for file without rmajor'
        nfail = nfail + 1
    end if

    if (nfail /= 0) then
        print *, nfail, ' boozer_chartmap_io tests failed'
        error stop 'test_boozer_chartmap_io failed'
    end if
    print *, 'PASS: read_boozer_chartmap parses field grid, period, and rmajor'

contains

    pure function bmod_model(ir, it, ip) result(b)
        integer, intent(in) :: ir, it, ip
        real(dp) :: b
        b = 1.0_dp + 0.1_dp*ir + 0.01_dp*it + 0.001_dp*ip
    end function bmod_model

    subroutine write_chartmap(filename, with_field_dims, with_rmajor)
        character(len=*), intent(in) :: filename
        logical, intent(in) :: with_field_dims, with_rmajor

        integer :: ncid, dim_rho, dim_theta, dim_zeta
        integer :: dim_theta_f, dim_zeta_f, bmod_dims(3)
        integer :: var_rho, var_theta, var_zeta, var_aphi, var_bth, var_bph
        integer :: var_bmod, var_nfp
        integer :: nt_field, np_field, ir, it, ip
        real(dp) :: h_theta, h_phi
        real(dp) :: rho(n_rho), theta(n_theta_geom), zeta(n_phi_geom)
        real(dp) :: aphi(n_rho), bth(n_rho), bph(n_rho)
        real(dp), allocatable :: bmod(:, :, :)

        ! Geometry grid steps chosen so the field grid spans the full period.
        h_theta = twopi / real(n_theta_field - 1, dp)
        h_phi = twopi / real(nfp, dp) / real(n_phi_field - 1, dp)
        do ir = 1, n_rho
            rho(ir) = 0.1_dp + 0.9_dp*real(ir - 1, dp)/real(n_rho - 1, dp)
            aphi(ir) = real(ir, dp)
            bth(ir) = 10.0_dp*real(ir, dp)
            bph(ir) = 100.0_dp*real(ir, dp)
        end do
        do it = 1, n_theta_geom
            theta(it) = real(it - 1, dp)*h_theta
        end do
        do ip = 1, n_phi_geom
            zeta(ip) = real(ip - 1, dp)*h_phi
        end do

        if (with_field_dims) then
            nt_field = n_theta_field
            np_field = n_phi_field
        else
            nt_field = n_theta_geom
            np_field = n_phi_geom
        end if
        allocate (bmod(n_rho, nt_field, np_field))
        do ip = 1, np_field
            do it = 1, nt_field
                do ir = 1, n_rho
                    bmod(ir, it, ip) = bmod_model(ir, it, ip)
                end do
            end do
        end do

        call nc(nf90_create(filename, nf90_clobber, ncid), 'create')
        call nc(nf90_def_dim(ncid, 'rho', n_rho, dim_rho), 'dim rho')
        call nc(nf90_def_dim(ncid, 'theta', n_theta_geom, dim_theta), 'dim theta')
        call nc(nf90_def_dim(ncid, 'zeta', n_phi_geom, dim_zeta), 'dim zeta')
        if (with_field_dims) then
            call nc(nf90_def_dim(ncid, 'theta_field', n_theta_field, dim_theta_f), &
                'dim theta_field')
            call nc(nf90_def_dim(ncid, 'zeta_field', n_phi_field, dim_zeta_f), &
                'dim zeta_field')
            bmod_dims = [dim_rho, dim_theta_f, dim_zeta_f]
        else
            bmod_dims = [dim_rho, dim_theta, dim_zeta]
        end if

        call nc(nf90_def_var(ncid, 'rho', nf90_double, [dim_rho], var_rho), 'var rho')
        call nc(nf90_def_var(ncid, 'theta', nf90_double, [dim_theta], var_theta), 'var theta')
        call nc(nf90_def_var(ncid, 'zeta', nf90_double, [dim_zeta], var_zeta), 'var zeta')
        call nc(nf90_def_var(ncid, 'A_phi', nf90_double, [dim_rho], var_aphi), 'var A_phi')
        call nc(nf90_def_var(ncid, 'B_theta', nf90_double, [dim_rho], var_bth), 'var B_theta')
        call nc(nf90_def_var(ncid, 'B_phi', nf90_double, [dim_rho], var_bph), 'var B_phi')
        call nc(nf90_def_var(ncid, 'Bmod', nf90_double, bmod_dims, var_bmod), 'var Bmod')
        call nc(nf90_def_var(ncid, 'num_field_periods', nf90_int, var_nfp), 'var nfp')
        call nc(nf90_put_att(ncid, nf90_global, 'torflux', torflux_val), 'att torflux')
        if (with_rmajor) &
            call nc(nf90_put_att(ncid, nf90_global, 'rmajor', rmajor_val), 'att rmajor')
        call nc(nf90_enddef(ncid), 'enddef')

        call nc(nf90_put_var(ncid, var_rho, rho), 'put rho')
        call nc(nf90_put_var(ncid, var_theta, theta), 'put theta')
        call nc(nf90_put_var(ncid, var_zeta, zeta), 'put zeta')
        call nc(nf90_put_var(ncid, var_aphi, aphi), 'put A_phi')
        call nc(nf90_put_var(ncid, var_bth, bth), 'put B_theta')
        call nc(nf90_put_var(ncid, var_bph, bph), 'put B_phi')
        call nc(nf90_put_var(ncid, var_bmod, bmod), 'put Bmod')
        call nc(nf90_put_var(ncid, var_nfp, nfp), 'put nfp')
        call nc(nf90_close(ncid), 'close')
        deallocate (bmod)
    end subroutine write_chartmap

    subroutine nc(status, loc)
        integer, intent(in) :: status
        character(len=*), intent(in) :: loc
        if (status /= nf90_noerr) then
            print *, 'NetCDF write error at ', trim(loc), ': ', trim(nf90_strerror(status))
            error stop
        end if
    end subroutine nc

    subroutine expect_int(got, want, label, nfail)
        integer, intent(in) :: got, want
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (got /= want) then
            print *, 'FAIL: ', label, ' got ', got, ' want ', want
            nfail = nfail + 1
        end if
    end subroutine expect_int

    subroutine expect_real(got, want, label, nfail)
        real(dp), intent(in) :: got, want
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (abs(got - want) > tol*max(1.0_dp, abs(want))) then
            print *, 'FAIL: ', label, ' got ', got, ' want ', want
            nfail = nfail + 1
        end if
    end subroutine expect_real

end program test_boozer_chartmap_io
