program test_boozer_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_boozer_chartmap, only: boozer_chartmap_field_t, &
                                     create_boozer_chartmap_field, &
                                     is_boozer_chartmap

    implicit none

    type(boozer_chartmap_field_t), allocatable :: field
    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    character(len=256) :: filename
    logical :: is_bc
    integer :: nfail

    nfail = 0
    filename = 'test_boozer_chartmap.nc'

    ! Test detection
    is_bc = is_boozer_chartmap(filename)
    if (.not. is_bc) then
        print *, 'FAIL: is_boozer_chartmap returned false'
        nfail = nfail + 1
    else
        print *, 'PASS: is_boozer_chartmap detected file correctly'
    end if

    ! Test loading
    call create_boozer_chartmap_field(filename, field)

    if (.not. field%initialized) then
        print *, 'FAIL: field not initialized'
        error stop
    end if

    if (field%nfp /= 2) then
        print *, 'FAIL: nfp =', field%nfp, ' expected 2'
        nfail = nfail + 1
    else
        print *, 'PASS: nfp = 2'
    end if

    if (field%torflux <= 0.0_dp) then
        print *, 'FAIL: torflux =', field%torflux, ' expected > 0'
        nfail = nfail + 1
    else
        print *, 'PASS: torflux =', field%torflux
    end if

    ! Test field evaluation at mid-radius, theta=0, phi=0
    ! x(1) = rho = sqrt(s), so rho=0.5 means s=0.25
    x = [0.5_dp, 0.0_dp, 0.0_dp]  ! rho=0.5, theta=0, phi=0
    call field%evaluate(x, Acov, hcov, Bmod)

    ! Bmod should be positive and reasonable (analytic tokamak ~ B0)
    if (Bmod <= 0.0_dp) then
        print *, 'FAIL: Bmod =', Bmod, ' expected > 0'
        nfail = nfail + 1
    else
        print *, 'PASS: Bmod =', Bmod
    end if

    ! A_theta = torflux * s = torflux * rho^2 = torflux * 0.25
    if (abs(Acov(2) - field%torflux * 0.25_dp) > 1.0e-8_dp * abs(Acov(2))) then
        print *, 'FAIL: Acov(2) =', Acov(2), ' expected', field%torflux * 0.25_dp
        nfail = nfail + 1
    else
        print *, 'PASS: Acov(2) = torflux * rho^2'
    end if

    ! hcov(1) should be 0 (no radial component in Boozer)
    if (abs(hcov(1)) > 1.0e-15_dp) then
        print *, 'FAIL: hcov(1) =', hcov(1), ' expected 0'
        nfail = nfail + 1
    else
        print *, 'PASS: hcov(1) = 0'
    end if

    ! hcov(2) and hcov(3) should be B_theta/Bmod and B_phi/Bmod
    if (abs(hcov(2)) < 1.0e-15_dp .or. abs(hcov(3)) < 1.0e-15_dp) then
        print *, 'FAIL: hcov(2) or hcov(3) is zero'
        nfail = nfail + 1
    else
        print *, 'PASS: hcov(2) =', hcov(2), ' hcov(3) =', hcov(3)
    end if

    ! Test at another point: rho=0.7, theta=1, phi=0.5
    x = [0.7_dp, 1.0_dp, 0.5_dp]
    call field%evaluate(x, Acov, hcov, Bmod)

    if (Bmod <= 0.0_dp) then
        print *, 'FAIL: Bmod at second point =', Bmod
        nfail = nfail + 1
    else
        print *, 'PASS: Bmod at second point =', Bmod
    end if

    ! Summary
    print *, ''
    if (nfail == 0) then
        print *, 'All tests passed.'
    else
        print *, nfail, ' tests failed.'
        error stop 'Test failures'
    end if

end program test_boozer_chartmap
