program test_magfie_vmec_near_axis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use simple, only: init_vmec
    use magfie_sub, only: magfie_vmec
    implicit none
    real(dp) :: dummy, x(3), b, jac, db(3), cov(3), ctr(3), hc(3)
    real(dp) :: bp, bm, cp(3), cm(3), ignored(3), ignored_curl(3), jc, step, radial_b, radial_h(3)
    real(dp), parameter :: radii(3) = [1d-6, 1d-4, 1d-3]
    integer :: i, nfailed

    call init_vmec('wout_ncsx.nc', 5, 5, 1, dummy)
    nfailed = 0
    do i = 1, size(radii)
        x = [radii(i), 0.37d0, 0.13d0]
        call magfie_vmec(x, b, jac, db, cov, ctr, hc)
        if (.not. all(ieee_is_finite(db)) .or. .not. all(ieee_is_finite(hc))) then
            print *, 'FAILED: nonfinite near-axis derivatives', x(1)
            nfailed = nfailed + 1
            cycle
        end if
        ! A much narrower native-field stencil checks the physical derivative.
        step = x(1)*1d-4
        x(1) = radii(i) + step
        call magfie_vmec(x, bp, jc, ignored, cp, ctr, ignored_curl)
        x(1) = radii(i) - step
        call magfie_vmec(x, bm, jc, ignored, cm, ctr, ignored_curl)
        radial_b = (bp - bm)/(2d0*step*b)
        radial_h = (cp - cm)/(2d0*step)
        if (abs(db(1) - radial_b) > 0.02d0*max(abs(radial_b), 1d-8)) then
         print *, 'FAILED: near-axis radial field derivative', radii(i), db(1), radial_b
            nfailed = nfailed + 1
        end if
        ! Angular derivatives are unchanged; their difference eliminates them.
        ! Covariant curl terms with a radial derivative must stay finite.
        if (.not. all(ieee_is_finite(radial_h))) nfailed = nfailed + 1
    end do
    if (nfailed > 0) error stop 'near-axis VMEC derivative regression'
    print *, 'PASS: positive-flux near-axis derivatives'
end program
