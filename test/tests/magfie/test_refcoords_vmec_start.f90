program test_refcoords_vmec_start
    !> REFCOORDS on a VMEC reference integrates in rho = sqrt(s), while start.dat,
    !> sbeg and all output hold VMEC s. A start at reference (s, theta, phi) must
    !> therefore sit where magfie_vmec sees the same |B|, map back to the same s,
    !> and give the collision profiles that s.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: coord_input, field_input
    use velo_mod, only: isw_field_type
    use magfie_sub, only: magfie, init_magfie, VMEC, REFCOORDS
    use field_can_mod, only: ref_to_integ, integ_to_ref
    use profile_radius, only: init_profile_radius, profile_s

    implicit none

    ! Spline error of |B| is far below this; s against sqrt(s) differs by percent.
    real(dp), parameter :: tol_bmod = 1.0e-4_dp
    real(dp), parameter :: tol_s = 1.0e-12_dp
    real(dp), parameter :: s_values(3) = [0.09_dp, 0.3_dp, 0.6_dp]
    real(dp), parameter :: angles(2, 2) = reshape([0.4_dp, 0.3_dp, 2.5_dp, 1.1_dp], &
                                                  [2, 2])

    type(tracer_t) :: norb
    real(dp) :: u(3), x(3), u_back(3)
    real(dp) :: bmod_vmec, bmod_ref, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)
    integer :: i, j, nfail

    nfail = 0
    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    isw_field_type = REFCOORDS
    call init_field(norb, 'wout.nc', 5, 5, 3, 0)
    call init_profile_radius

    do i = 1, size(s_values)
        do j = 1, size(angles, 2)
            u = [s_values(i), angles(1, j), angles(2, j)]

            call init_magfie(VMEC)
            call magfie(u, bmod_vmec, sqrtg, bder, hcov, hctr, hcurl)

            call ref_to_integ(u, x)
            call init_magfie(REFCOORDS)
            call magfie(x, bmod_ref, sqrtg, bder, hcov, hctr, hcurl)
            if (abs(bmod_ref - bmod_vmec) > tol_bmod*abs(bmod_vmec)) then
                write (error_unit, '(a,3es12.4,a,2es16.8)') '  |B| at u =', u, &
                    ': REFCOORDS, VMEC =', bmod_ref, bmod_vmec
                nfail = nfail + 1
            end if

            call integ_to_ref(x, u_back)
            if (any(abs(u_back - u) > tol_s)) then
                write (error_unit, '(a,3es12.4,a,3es12.4)') '  roundtrip of', u, &
                    ' gives', u_back
                nfail = nfail + 1
            end if

            if (abs(profile_s(x) - u(1)) > tol_s) then
                write (error_unit, '(a,es12.4,a,es12.4)') '  profile_s at s =', &
                    u(1), ' is', profile_s(x)
                nfail = nfail + 1
            end if
        end do
    end do

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, ' check(s) failed'
        error stop 1
    end if
    print *, 'test_refcoords_vmec_start: all checks passed'
end program test_refcoords_vmec_start
