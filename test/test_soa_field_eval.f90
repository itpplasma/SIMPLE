program test_soa_field_eval
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_boozer, only: eval_field_booz, eval_field_booz_many
    use field_can_base, only: field_can_t
    use boozer_sub, only: get_boozer_coordinates
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm
    use simple_main, only: init_field
    use simple, only: tracer_t

    implicit none

    real(dp), parameter :: RTOL = 5.0d-12
    type(tracer_t) :: norb
    type(field_can_t) :: f
    character(len=256) :: config_file
    integer :: i, j, npts, nerrors

    real(dp), allocatable :: r(:), th(:), ph(:)
    real(dp), allocatable :: Ath(:), Aph(:), dAth_dr(:), dAph_dr(:), d2Aph_dr2(:)
    real(dp), allocatable :: hth(:), hph(:), dhth(:,:), dhph(:,:)
    real(dp), allocatable :: d2hth_dr2(:), d2hph_dr2(:)
    real(dp), allocatable :: Bmod(:), dBmod(:,:), d2Bmod(:,:)

    print *, '=========================================='
    print *, 'Testing eval_field_booz_many'
    print *, '=========================================='

    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)
    call get_boozer_coordinates()

    npts = 100
    allocate(r(npts), th(npts), ph(npts))
    allocate(Ath(npts), Aph(npts), dAth_dr(npts), dAph_dr(npts), d2Aph_dr2(npts))
    allocate(hth(npts), hph(npts), dhth(3, npts), dhph(3, npts))
    allocate(d2hth_dr2(npts), d2hph_dr2(npts))
    allocate(Bmod(npts), dBmod(3, npts), d2Bmod(6, npts))

    do i = 1, npts
        r(i) = 0.1d0 + 0.8d0 * real(i-1, dp) / real(npts-1, dp)
        th(i) = 6.28318530718d0 * real(i, dp) / real(npts, dp)
        ph(i) = 3.14159265359d0 * real(i, dp) / real(npts, dp)
    end do

    call eval_field_booz_many(npts, r, th, ph, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod)

    nerrors = 0
    do i = 1, npts
        call eval_field_booz(f, r(i), th(i), ph(i), 1)

        if (reldiff(Ath(i), f%Ath) > RTOL) then
            print *, 'Mismatch Ath at i=', i, Ath(i), f%Ath
            nerrors = nerrors + 1
        end if
        if (reldiff(Aph(i), f%Aph) > RTOL) then
            print *, 'Mismatch Aph at i=', i, Aph(i), f%Aph
            nerrors = nerrors + 1
        end if
        if (reldiff(dAth_dr(i), f%dAth(1)) > RTOL) then
            print *, 'Mismatch dAth_dr at i=', i, dAth_dr(i), f%dAth(1)
            nerrors = nerrors + 1
        end if
        if (reldiff(dAph_dr(i), f%dAph(1)) > RTOL) then
            print *, 'Mismatch dAph_dr at i=', i, dAph_dr(i), f%dAph(1)
            nerrors = nerrors + 1
        end if
        if (reldiff(d2Aph_dr2(i), f%d2Aph(1)) > RTOL) then
            print *, 'Mismatch d2Aph_dr2 at i=', i, d2Aph_dr2(i), f%d2Aph(1)
            nerrors = nerrors + 1
        end if
        if (reldiff(hth(i), f%hth) > RTOL) then
            print *, 'Mismatch hth at i=', i, hth(i), f%hth
            nerrors = nerrors + 1
        end if
        if (reldiff(hph(i), f%hph) > RTOL) then
            print *, 'Mismatch hph at i=', i, hph(i), f%hph
            nerrors = nerrors + 1
        end if
        do j = 1, 3
            if (reldiff(dhth(j, i), f%dhth(j)) > RTOL) then
                print *, 'Mismatch dhth at i=', i, 'j=', j
                nerrors = nerrors + 1
            end if
            if (reldiff(dhph(j, i), f%dhph(j)) > RTOL) then
                print *, 'Mismatch dhph at i=', i, 'j=', j
                nerrors = nerrors + 1
            end if
        end do
        if (reldiff(d2hth_dr2(i), f%d2hth(1)) > RTOL) then
            print *, 'Mismatch d2hth_dr2 at i=', i, d2hth_dr2(i), f%d2hth(1)
            nerrors = nerrors + 1
        end if
        if (reldiff(d2hph_dr2(i), f%d2hph(1)) > RTOL) then
            print *, 'Mismatch d2hph_dr2 at i=', i, d2hph_dr2(i), f%d2hph(1)
            nerrors = nerrors + 1
        end if
        if (reldiff(Bmod(i), f%Bmod) > RTOL) then
            print *, 'Mismatch Bmod at i=', i, Bmod(i), f%Bmod
            nerrors = nerrors + 1
        end if
        do j = 1, 3
            if (reldiff(dBmod(j, i), f%dBmod(j)) > RTOL) then
                print *, 'Mismatch dBmod at i=', i, 'j=', j
                nerrors = nerrors + 1
            end if
        end do
        do j = 1, 6
            if (reldiff(d2Bmod(j, i), f%d2Bmod(j)) > RTOL) then
                print *, 'Mismatch d2Bmod at i=', i, 'j=', j
                nerrors = nerrors + 1
            end if
        end do
    end do

    if (nerrors > 0) then
        print *, 'FAILED: ', nerrors, ' errors found'
        stop 1
    end if

    print *, 'PASSED: eval_field_booz_many matches single-point version'

contains

    pure function reldiff(a, b) result(rd)
        real(dp), intent(in) :: a, b
        real(dp) :: rd, denom
        denom = max(abs(a), abs(b), 1.0d-30)
        rd = abs(a - b) / denom
    end function reldiff

end program test_soa_field_eval
