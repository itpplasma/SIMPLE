program sheet_accuracy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use candidate_field, only: model_t, local_t, geometry, frame, canonical_pphi
    use candidate_orbits, only: options_t, orbit_t, initialize, advance, particle
    implicit none
    integer, parameter :: nm = 256
    type(model_t) :: pm
    real(dp) :: errors(3, nm), speed, pi, rms(3), valid(nm), active(nm)
    integer :: case_id, candidate, level, im, unit
    character(8) :: names(6)
    names = ['C1L0    ', 'C1L1    ', 'C2fixed ', 'C2pphi  ', 'C4      ', 'C5      ']
    pm%torus%a_in = 0.25_dp
    pm%torus%a_out = 2.0_dp
    pi = acos(-1.0_dp)
    open (newunit=unit, &
          file='concept/interface_candidates/results/sheet_accuracy.csv', &
          status='replace')
    write (unit, '(a)') 'case,method,speed,valid,position_rms,velocity_rms,' &
        // 'pphi_change_rms,switches_mean'
    do case_id = 1, 2
        pm%delta = 0.1_dp*real(case_id - 1, dp)
        do candidate = 1, 6
            if (case_id == 2 .and. candidate == 4) cycle
            do level = 1, 3
                speed = 0.01_dp*2**(level - 1)
                !$omp parallel do default(shared) private(im) schedule(dynamic)
                do im = 1, nm
                    call segment(im, errors(:, im), valid(im), active(im))
                end do
                !$omp end parallel do
                rms = sqrt(sum(errors**2, dim=2)/max(sum(valid), 1.0_dp))
                write (unit, '(i0,3a,2(es24.16,a),4(es24.16,a))') &
                    case_id, ',', trim(names(candidate)), ',', speed, ',', &
                    sum(valid), ',', rms(1), ',', rms(2), ',', rms(3), ',', &
                    sum(active)/nm, ''
            end do
        end do
    end do
    close (unit)
contains
    subroutine segment(im, result, good_marker, switches)
        integer, intent(in) :: im
        real(dp), intent(out) :: result(3), good_marker, switches
        type(options_t) :: opt, refopt
        type(orbit_t) :: ref, orbit
        type(local_t) :: lf
        type(model_t) :: control_pm
        real(dp) :: y0(6), yf(6), yg(6), theta, phase, pitch, rr, r, e1(3), e2(3)
        integer :: j
        logical :: ok
        theta = 2*pi*real(im, dp)/nm
        phase = 2*pi*modulo(0.618033988749895_dp*im, 1.0_dp)
        pitch = 1.6_dp*modulo(0.414213562373095_dp*im, 1.0_dp) - 0.8_dp
        r = 0.5_dp - 0.35_dp*speed
        if (candidate == 3 .or. candidate == 4) r = 0.8_dp
        rr = 10.0_dp + r*cos(theta)
        y0(1:3) = [rr, 0.0_dp, r*sin(theta)]
        call geometry(pm, y0(1:3), merge(1, 2, r < 0.5_dp), lf)
        call frame(lf%b, e1, e2)
        y0(4:6) = speed*(pitch*lf%b + sqrt(1-pitch*pitch)* &
            (cos(phase)*e1 + sin(phase)*e2))
        refopt%tol = 1.0e-11_dp
        opt%tol = 1.0e-11_dp
        select case (candidate)
        case (1, 2)
            opt%method = 1
            opt%old_level = candidate - 1
        case (3, 4)
            opt%method = 2
            opt%phase_rule = candidate - 3
        case (5)
            opt%method = 4
            opt%order = 1
            opt%patch_step = 0.25_dp
        case (6)
            opt%method = 5
            opt%boris_step = 0.025_dp
        end select
        control_pm = pm
        if (candidate == 6) control_pm%torus%width = 0.2_dp*speed
        call initialize(pm, refopt, y0, ref)
        call initialize(control_pm, opt, y0, orbit)
        do j = 1, 100
            call advance(pm, refopt, ref, 0.5_dp)
            call advance(control_pm, opt, orbit, 0.5_dp)
        end do
        call particle(pm, refopt, ref, yf, ok)
        call particle(control_pm, opt, orbit, yg, ok)
        good_marker = merge(1.0_dp, 0.0_dp, ok .and. .not. orbit%failed)
        result = 0.0_dp
        if (good_marker > 0) then
            result(1) = norm2(yg(1:3) - yf(1:3))/speed
            result(2) = norm2(yg(4:6) - yf(4:6))/speed
            result(3) = abs(canonical_pphi(pm, yg) - canonical_pphi(pm, y0))/(10*speed)
        end if
        switches = orbit%switches
    end subroutine
end program
