program test_field_can_meiss

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use util, only: twopi
    use simple, only: tracer_t
    use simple_main, only: init_field
    use velo_mod, only: isw_field_type
    use field, only: vmec_field_t, create_vmec_field
    use field_can_mod, only: eval_field => evaluate, field_can_t
    use magfie_sub, only: MEISS
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, &
                               spline_transformation, init_canonical_field_components, &
           xmin, xmax, h_r, h_phi, h_th, ah_cov_on_slice, n_r, n_phi, n_th, &
           spl_transform_batch, transformation_relerr, ref_to_integ_meiss, &
           integ_to_ref_meiss
    use interpolate, only: evaluate_batch_splines_3d
    use params, only: canonical_grid_nr, canonical_grid_ntheta, &
        canonical_grid_nphi, canonical_ode_relerr, field_input, coord_input, &
        read_config
    implicit none

    type(tracer_t) :: norb
    type(vmec_field_t) :: magfie
    integer :: config_unit

    isw_field_type = MEISS
    call create_vmec_field(magfie)

    open(newunit=config_unit, file='canonical-map-config.nml', status='replace')
    write(config_unit, '(a)') '&config'
    write(config_unit, '(a)') 'canonical_grid_nr=8'
    write(config_unit, '(a)') 'canonical_grid_ntheta=9'
    write(config_unit, '(a)') 'canonical_grid_nphi=10'
    write(config_unit, '(a)') 'canonical_ode_relerr=1.0e-8'
    write(config_unit, '(a)') '/'
    close(config_unit)
    call read_config('canonical-map-config.nml')
    open(newunit=config_unit, file='canonical-map-config.nml', status='old')
    close(config_unit, status='delete')
    field_input = 'wout.nc'
    coord_input = 'wout.nc'

    print *, 'init_field'
    call init_field(norb, 'wout.nc', 5, 5, 3, 0)
    if (n_r /= 8 .or. n_th /= 9 .or. n_phi /= 10) then
        error stop 'configured canonical map grid was not applied'
    end if
    if (transformation_relerr /= 1.0e-8_dp) then
        error stop 'configured canonical ODE tolerance was not applied'
    end if
    call assert_meiss_map
    print *, 'test_covar_components'
    call test_covar_components

    print *, 'field_can_meiss.write_transformation'
    call write_transformation('lam_chi.out')

    print *, 'test_evaluate_vmec'
    call test_evaluate_vmec

    print *, 'test_evaluate_meiss'
    call test_evaluate_meiss

    call test_repeated_initialization

contains

    subroutine test_repeated_initialization
        call init_meiss(magfie, 12, 13, 14, 0.01_dp, 0.9_dp, 0.2_dp, 6.0_dp, &
            transformation_relerr_=1.0e-8_dp)
        if (any([n_r, n_th, n_phi] /= [12, 13, 14])) then
            error stop 'custom canonical map grid was not applied'
        end if
        if (transformation_relerr /= 1.0e-8_dp) then
            error stop 'custom canonical ODE tolerance was not applied'
        end if
        call get_meiss_coordinates
        call assert_meiss_map

        call init_meiss(magfie)
        if (any([n_r, n_th, n_phi] /= [62, 63, 64])) then
            error stop 'canonical map grid defaults were not restored'
        end if
        if (transformation_relerr /= 1.0e-11_dp) then
            error stop 'canonical ODE tolerance default was not restored'
        end if
        if (any(xmin /= [1.0e-6_dp, 0.0_dp, 0.0_dp]) .or. &
            xmax(1) /= 1.0_dp .or. xmax(2) /= twopi) then
            error stop 'canonical map bounds defaults were not restored'
        end if
    end subroutine test_repeated_initialization

    subroutine assert_meiss_map
        real(dp), parameter :: tolerance = 1.0e-7_dp
        real(dp) :: xref(3), xref_back(3), xinteg(3), xinteg_outer(3)
        real(dp) :: error(3)
        type(field_can_t) :: field_value, periodic_value

        xref = [0.25_dp, 1.1_dp, 0.3_dp]
        call ref_to_integ_meiss(xref, xinteg)
        call integ_to_ref_meiss(xinteg, xref_back)
        error(1) = abs(xref_back(1) - xref(1))
        error(2) = abs(modulo(xref_back(2) - xref(2) + 0.5_dp*twopi, twopi) - &
            0.5_dp*twopi)
        error(3) = abs(modulo(xref_back(3) - xref(3) + 0.5_dp*twopi, twopi) - &
            0.5_dp*twopi)
        if (maxval(error) > tolerance) error stop 'Meiss coordinate roundtrip failed'

        call ref_to_integ_meiss([0.36_dp, xref(2), xref(3)], xinteg_outer)
        if (xinteg_outer(1) <= xinteg(1)) error stop 'Meiss radial orientation failed'

        call eval_field(field_value, xinteg(1), xinteg(2), xinteg(3), 1)
        call eval_field(periodic_value, xinteg(1), xinteg(2), &
            xinteg(3) + xmax(3), 1)
        if (.not. ieee_is_finite(field_value%Bmod) .or. field_value%Bmod <= 0.0_dp .or. &
            .not. all(ieee_is_finite(field_value%dAth)) .or. &
            .not. all(ieee_is_finite(field_value%dAph)) .or. &
            .not. ieee_is_finite(field_value%hth) .or. &
            .not. ieee_is_finite(field_value%hph)) then
            error stop 'Meiss field evaluation is not finite'
        end if
        if (abs(periodic_value%Bmod - field_value%Bmod) > tolerance) then
            error stop 'Meiss field-period seam failed'
        end if
    end subroutine assert_meiss_map

    subroutine test_covar_components
        real(dp) :: r, phi, th
        real(dp) :: Ar, Ap, hr, hp
        integer :: i_r, i_phi, i_th
        integer :: funit

        open (newunit=funit, file='covar_components.out')
        write (funit, *) '#', ' r', ' phi', ' th', ' Arcov', ' Apcov', ' Atcov', &
            ' hrcov', ' hpcov', ' htcov', ' Bmod'
        do i_phi = 1, n_phi
            phi = xmin(3) + h_phi*(i_phi - 1)
            do i_th = 1, n_th
                th = xmin(2) + h_th*(i_th - 1)
                do i_r = 1, n_r
                    r = xmin(1) + h_r*(i_r - 1)
                    call ah_cov_on_slice(r, phi, i_th, Ar, Ap, hr, hp)
                    write (funit, *) r, phi, th, Ar, Ap, 0d0, hr, hp, 0d0, 0d0
                end do
            end do
        end do
        close (funit)
    end subroutine test_covar_components

    subroutine write_transformation(filename)
        character(*), intent(in) :: filename

        integer :: funit
        integer :: i_r, i_th, i_phi
        real(dp) :: r, th, phi
        real(dp) :: transformation(2)

        open (newunit=funit, file=filename, status='unknown')
        write (funit, *) '#', ' r', ' phi', ' th', ' lam_phi', ' chi_gauge'

        do i_th = 1, n_th
            th = xmin(2) + h_th*(i_th - 1)
            do i_phi = 1, n_phi
                phi = xmin(3) + h_phi*(i_phi - 1)
                do i_r = 1, n_r
                    r = xmin(1) + h_r*(i_r - 1)
                    call evaluate_batch_splines_3d(spl_transform_batch, &
                                                   [r, th, phi], transformation)
                    write (funit, *) r, phi, th, transformation
                end do
            end do
        end do

        close (funit)
    end subroutine write_transformation

    subroutine test_evaluate_vmec
        real(dp) :: r, phi, th
        real(dp) :: Acov(3), hcov(3), Bmod
        integer :: i_r, i_phi, i_th
        integer :: funit

        open (newunit=funit, file='field_vmec.out')
        write (funit, *) '#', ' r', ' phi', ' th', ' Arcov', ' Apcov', ' Atcov', &
            ' hrcov', ' hpcov', ' htcov', ' Bmod'
        do i_th = 1, n_th
            th = xmin(2) + h_th*(i_th - 1)
            do i_phi = 1, n_phi
                phi = xmin(3) + h_phi*(i_phi - 1)
                do i_r = 1, n_r
                    r = xmin(1) + h_r*(i_r - 1)
                    call magfie%evaluate([r, th, phi], Acov, hcov, Bmod)
                    write (funit, *) r, phi, th, Acov(1), Acov(3), Acov(2), &
                        hcov(1), hcov(3), hcov(2), Bmod
                end do
            end do
        end do
        close (funit)
    end subroutine test_evaluate_vmec

    subroutine test_evaluate_meiss
        real(dp) :: r, phi, th
        type(field_can_t) :: f
        integer :: i_r, i_phi, i_th
        integer :: funit

        open (newunit=funit, file='field_can_meiss.out')
        write (funit, *) '#', ' r', ' phi', ' th', ' Arcov', ' Apcov', ' Atcov', &
            ' hrcov', ' hpcov', ' htcov', ' Bmod'
        do i_th = 1, n_th
            th = xmin(3) + h_th*(i_th - 1)
            do i_phi = 1, n_phi
                phi = xmin(2) + h_phi*(i_phi - 1)
                do i_r = 1, n_r
                    r = xmin(1) + h_r*(i_r - 1)
                    call eval_field(f, r, th, phi, 0)
               write (funit, *) r, phi, th, 0d0, f%Aph, f%Ath, 0d0, f%hph, f%hth, f%Bmod
                end do
            end do
        end do
        close (funit)
    end subroutine test_evaluate_meiss

end program test_field_can_meiss
