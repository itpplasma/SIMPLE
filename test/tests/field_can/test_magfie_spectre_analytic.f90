program test_magfie_spectre_analytic
    !> Regression witness for the analytic SPECTRE drift derivatives (bder, hcurl)
    !> that replaced central finite differences in magfie_spectre. magfie_spectre
    !> now takes the derivatives from libneo spectre_evaluate_der (A Hessian +
    !> analytic metric derivative); magfie_spectre_fd below is the verbatim
    !> central-difference algorithm it replaced. At interior points of every
    !> volume the analytic result must equal the finite-difference reference to
    !> the FD truncation floor, so a wrong Hessian, chart chain factor, or metric
    !> derivative term fails this test.
    !>
    !> The reference is fed the analytic sqrtg so the comparison isolates the
    !> drift derivatives (bder, dh_i/dx); sqrtg, bmod, hcovar and hctrvr are
    !> computed by the same formulas on both paths and unchanged by the switch.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_sub, only: magfie, init_magfie, set_magfie_spectre_field, &
                          compute_hcurl, SPECTRE, M_TO_CM
    use field_spectre, only: spectre_field_t, create_spectre_field
    use util, only: twopi
    implicit none

    integer, parameter :: NPTS = 50
    real(dp), parameter :: TOL = 1.0d-6, PAD = 0.05_dp

    character(len=1024) :: h5file
    type(spectre_field_t) :: field
    integer :: ierr, mvol, lvol, ip
    real(dp) :: x(3), rho_lo, rho_hi, r3(3)
    real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: bder_fd(3), hcurl_fd(3)
    real(dp) :: errB, errH, maxErrB, maxErrH

    if (command_argument_count() < 1) then
        print *, 'usage: test_magfie_spectre_analytic <file.h5>'
        error stop 1
    end if
    call get_command_argument(1, h5file)

    call create_spectre_field(field, trim(h5file), ierr)
    if (ierr /= 0) error stop 'create_spectre_field failed'
    mvol = field%data%Mvol

    call set_magfie_spectre_field(field)
    call init_magfie(SPECTRE)

    call seed_rng()
    maxErrB = 0.0_dp
    maxErrH = 0.0_dp

    do lvol = 1, mvol
        rho_lo = real(lvol - 1, dp) + PAD
        rho_hi = real(lvol, dp) - PAD
        do ip = 1, NPTS
            call random_number(r3)
            x(1) = rho_lo + (rho_hi - rho_lo)*r3(1)
            x(2) = twopi*r3(2)
            x(3) = twopi*r3(3)

            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            call magfie_spectre_fd(field, x, sqrtg, bder_fd, hcurl_fd)

            errB = maxval(abs(bder - bder_fd))/(maxval(abs(bder_fd)) + tiny(1.0_dp))
            errH = maxval(abs(hcurl - hcurl_fd))/(maxval(abs(hcurl_fd)) + tiny(1.0_dp))
            maxErrB = max(maxErrB, errB)
            maxErrH = max(maxErrH, errH)
        end do
    end do

    print '(A,I0,A,I0,A)', 'magfie_spectre analytic vs FD: ', mvol, &
        ' volumes x ', NPTS, ' interior points'
    print '(A,ES12.4)', '  max rel err bder  = ', maxErrB
    print '(A,ES12.4)', '  max rel err hcurl = ', maxErrH
    if (maxErrB >= TOL .or. maxErrH >= TOL) then
        print '(A,ES10.2)', 'FAIL: analytic drift derivative mismatch above ', TOL
        error stop 1
    end if
    print *, 'analytic == FD to truncation PASS'

contains

    subroutine seed_rng()
        integer :: n
        integer, allocatable :: seed(:)

        call random_seed(size=n)
        allocate (seed(n))
        seed = 20240711
        call random_seed(put=seed)
    end subroutine seed_rng

    subroutine magfie_spectre_fd(fld, x, sqrtg_cm, bder, hcurl)
        !> Verbatim central-difference drift derivatives from magfie_spectre before
        !> the analytic switch: the radial stencil stays inside the current volume
        !> [floor(rho_g), floor(rho_g)+1] so the |B| interface jump never enters
        !> the difference. Uses the analytic sqrtg so only bder and dh_i/dx differ.
        type(spectre_field_t), intent(in) :: fld
        real(dp), intent(in) :: x(3), sqrtg_cm
        real(dp), intent(out) :: bder(3), hcurl(3)

        real(dp) :: Acov(3), hcov_p(3), hcov_m(3), Bmod_p, Bmod_m
        real(dp) :: dhcov(3, 3), xp(3), xm(3), denom, lo, hi
        real(dp), parameter :: h_fd = 1.0d-4
        integer :: j

        do j = 1, 3
            xp = x
            xm = x
            denom = 2.0d0*h_fd
            if (j == 1) then
                lo = real(floor(x(1)), dp)
                hi = lo + 1.0d0
                if (x(1) - h_fd < lo) then
                    xm(1) = x(1)
                    xp(1) = x(1) + denom
                else if (x(1) + h_fd > hi) then
                    xm(1) = x(1) - denom
                    xp(1) = x(1)
                else
                    xm(1) = x(1) - h_fd
                    xp(1) = x(1) + h_fd
                end if
            else
                xm(j) = x(j) - h_fd
                xp(j) = x(j) + h_fd
            end if
            call fld%evaluate(xp, Acov, hcov_p, Bmod_p)
            call fld%evaluate(xm, Acov, hcov_m, Bmod_m)
            bder(j) = (log(Bmod_p) - log(Bmod_m))/denom
            dhcov(:, j) = (hcov_p - hcov_m)*M_TO_CM/denom
        end do

        call compute_hcurl(sqrtg_cm, dhcov, hcurl)
    end subroutine magfie_spectre_fd

end program test_magfie_spectre_analytic
