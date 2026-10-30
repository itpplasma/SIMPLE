program test_classifier_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use simple_api, only: simple_read_config, simple_init, simple_sample, &
                          simple_classify
    use classification, only: iaaa_bou
    use magfie_sub, only: magfie
    use field_can_mod, only: ref_to_integ
    use params, only: zstart, class_passing, trap_par, fper, cut_in_per
    use util, only: twopi
    implicit none

    character(1024) :: wout, arg, radius
    character(256) :: message
    integer :: unit, surfaces, i, j, row, iteration, status
    real(dp) :: starts(5, 8), x(3), b, bmax, bmin, theta, pitch, expected
    real(dp) :: best(2), candidate(2), step(2), value, sign_extremum
    real(dp) :: extrema(2), boundary(2, 1001)
    logical :: improved, passing

    call get_command_argument(1, wout)
    call get_command_argument(2, arg)
    read (arg, *) surfaces
    call get_command_argument(3, radius)
    if (len_trim(radius) == 0) radius = "0.2"
    open (newunit=unit, file='simple.in', status='replace')
    write (unit, '(a)') '&config', 'multharm=3', 'contr_pp=-1d10', &
        'trace_time=1d-5', 'ntestpart=8', 'nper=1', 'npoiper=4', &
        'npoiper2=32', 'ntimstep=2', 'sbeg='//trim(radius)//',0.4', &
        "netcdffile='"//trim(wout)//"'", 'isw_field_type=2', &
        'deterministic=.True.', 'class_plot=.True.', 'fast_class=.True.', &
        'notrace_passing=1', 'tcut=-1d0', 'cut_in_per=0d0'
    write (unit, '(a,i0)') 'num_surf=', surfaces
    write (unit, '(a)') '/'
    close (unit)
    call simple_read_config('simple.in', status, message)
    if (status /= 0) error stop trim(message)
    call simple_init()
    call simple_sample(8, starts)
    call simple_classify()
    flush (iaaa_bou)

    ! Independent derivative-free surface search: a mesh followed by pattern
    ! refinement, without the production Newton solver or radial extrema cache.
    x(1) = zstart(1, 1)
    do j = 1, 2
        sign_extremum = real(3 - 2*j, dp)
        value = -huge(1.0_dp)
        do i = 0, 127
            do row = 0, 127
                x(2:3) = [twopi*i/128, fper*row/128]
                b = field_strength(x)
                if (sign_extremum*b <= value) cycle
                value = sign_extremum*b
                best = x(2:3)
            end do
        end do
        step = [twopi, fper]/128
        do iteration = 1, 100
            improved = .false.
            do i = -1, 1
                do row = -1, 1
                    candidate = best + step*[i, row]
                    x(2:3) = candidate
                    b = field_strength(x)
                    if (sign_extremum*b <= value) cycle
                    value = sign_extremum*b
                    best = candidate
                    improved = .true.
                end do
            end do
            if (.not. improved) step = step/2
            if (maxval(step) < 1.0e-10_dp) exit
        end do
        extrema(j) = sign_extremum*value
    end do
    bmax = extrema(1)
    bmin = extrema(2)
    open (newunit=unit, file='fort.20000', status='old')
    do row = 1, 1001
        read (unit, *) boundary(:, row)
        if (.not. all(ieee_is_finite(boundary(:, row)))) then
            print *, 'Nonfinite boundary at row', row
            error stop 1
        end if
    end do
    close (unit)
    do row = 1, 1001
        theta = boundary(1, row)
        pitch = boundary(2, row)
        x = [zstart(1, 1), theta, cut_in_per*fper]
        b = field_strength(x)
        ! Conservation of energy and magnetic moment at the mirror threshold:
        ! (1 - pitch**2)/B(cut) = 1/Bmax(surface).
        expected = 1.0_dp - b/bmax
        if (abs(pitch**2 - expected) > 2.0e-7_dp) then
            print *, 'Mirror boundary mismatch at row', row, pitch**2, expected
            error stop 2
        end if
    end do
    do i = 1, 8
        if (abs(zstart(1, i) - zstart(1, 1)) > 1.0e-12_dp) cycle
        call ref_to_integ(zstart(1:3, i), x)
        x(3) = cut_in_per*fper
        b = field_strength(x)
        passing = zstart(5, i)**2 > 1.0_dp - b/bmax
        if (passing .neqv. class_passing(i)) error stop 'mirror class mismatch'
        expected = ((1 - zstart(5, i)**2)*bmax/b - 1)*bmin/(bmax - bmin)
        if (abs(trap_par(i) - expected) > 2.0e-6_dp) &
            error stop 'trapping parameter mismatch'
    end do
    print *, 'Finite boundary and independent surface mirror criterion verified'

contains

    function field_strength(position) result(strength)
        real(dp), intent(in) :: position(3)
        real(dp) :: strength, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie(position, strength, sqrtg, bder, hcovar, hctrvr, hcurl)
    end function field_strength

end program test_classifier_boundary
