program test_beam_stopping
    !> Compare the Suzuki et al. (1998) fit with the independent C
    !> implementation in ASCOT5 (src/suzuki.c, itpplasma/ascot5 a6cc085a),
    !> evaluated at the same points, and check the impurity enhancement.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use beam_stopping, only: suzuki_stopping
    implicit none
    ! E [keV/amu], ne [m^-3], Te [eV], impurity Z, n_Z/n_e, ASCOT5 sigma [m^2]
    real(dp), parameter :: ref(6, 8) = reshape([ &
                     30.0_dp, 3e19_dp, 1000.0_dp, 0.0_dp, 0.0_dp, 5.4711179033e-20_dp, &
                     46.5_dp, 5e19_dp, 2500.0_dp, 0.0_dp, 0.0_dp, 3.9416465970e-20_dp, &
                      10.0_dp, 2e19_dp, 300.0_dp, 0.0_dp, 0.0_dp, 1.0517221776e-19_dp, &
                    30.0_dp, 3e19_dp, 1000.0_dp, 7.0_dp, 0.02_dp, 5.7926440205e-20_dp, &
                     46.5_dp, 6e19_dp, 800.0_dp, 6.0_dp, 0.01_dp, 4.4450763464e-20_dp, &
                   200.0_dp, 5e19_dp, 5000.0_dp, 6.0_dp, 0.01_dp, 1.4865516954e-20_dp, &
                   500.0_dp, 8e19_dp, 10000.0_dp, 0.0_dp, 0.0_dp, 6.6838449862e-21_dp, &
             15.5_dp, 4e19_dp, 1500.0_dp, 5.0_dp, 0.03_dp, 8.0549463431e-20_dp], [6, 8])
    real(dp) :: sigma, zeff, nz, nd, err
    integer :: i, z

    do i = 1, size(ref, 2)
        z = nint(ref(4, i))
        nz = ref(5, i)
        nd = 1.0_dp - z*nz
        zeff = (nd + z**2*nz)/(nd + z*nz)
        sigma = suzuki_stopping(ref(1, i), ref(2, i), ref(3, i), 2, zeff, z)
        err = abs(sigma/ref(6, i) - 1.0_dp)
    print '(a,i0,a,es16.8,a,es9.2)', 'case ', i, ' sigma = ', sigma, ' rel. error ', err
        ! The reference values carry ten significant digits.
       if (err > 1.0e-9_dp) error stop 'beam stopping differs from the ASCOT5 reference'
    end do
    ! Impurities increase stopping at these conditions (Suzuki Fig. 3).
    if (suzuki_stopping(30.0_dp, 3e19_dp, 1000.0_dp, 2, 2.0_dp, 6) <= &
        suzuki_stopping(30.0_dp, 3e19_dp, 1000.0_dp, 2, 1.0_dp, 6)) &
        error stop 'carbon does not increase beam stopping'
    print '(a)', 'PASS'
end program test_beam_stopping
