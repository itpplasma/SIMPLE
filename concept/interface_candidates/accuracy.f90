program accuracy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use candidate_field, only: model_t, local_t, geometry, frame, side_of
    use transforms, only: lift, inverse
    use candidate_orbits, only: options_t, orbit_t, initialize, advance, &
        particle, dp_step
    implicit none
    integer, parameter :: nm = 256
    type(model_t) :: pm
    real(dp) :: errors(4, nm), speed, phase, pi, rms(4)
    integer :: case_id, order, level, im, unit
    pm%torus%a_in = 0.25_dp
    pm%torus%a_out = 2.0_dp
    pi = acos(-1.0_dp)
    open (newunit=unit, file='concept/interface_candidates/results/accuracy.csv', &
          status='replace')
    write (unit, '(a)') 'case,method,order,speed,position_rms,velocity_rms,' &
        // 'reverse_rms,energy_rms'
    do case_id = 1, 2
        pm%delta = 0.1_dp*real(case_id - 1, dp)
        do order = 0, 1
            do level = 1, 3
                speed = 0.01_dp*2**(level - 1)
                !$omp parallel do default(shared) private(im, phase) schedule(dynamic)
                do im = 1, nm
                    phase = 2*pi*real(im - 1, dp)/nm
                    call segment(im, phase, errors(:, im))
                end do
                !$omp end parallel do
                rms = sqrt(sum(errors**2, dim=2)/nm)
                write (unit, '(i0,a,i0,5(a,es24.16))') case_id, ',GC,', order, &
                    ',', speed, (',', rms(im), im=1, 4)
            end do
        end do
    end do
    close (unit)
contains
    subroutine segment(im, phase, result)
        integer, intent(in) :: im
        real(dp), intent(in) :: phase
        real(dp), intent(out) :: result(4)
        type(options_t) :: opt
        type(orbit_t) :: ref
        real(dp) :: g(6), g0(6), y0(6), yf(6), yg(6), residual, energy
        logical :: ok
        g = [10.8_dp, 0.0_dp, 0.1_dp, &
             speed*(0.2_dp + 0.6_dp*real(im, dp)/nm), 0.0_dp, phase]
        g(5) = (speed**2 - g(4)**2)/2.0_dp
        call lift(pm, g, 0, y0, ok)
        call inverse(pm, y0, order, g0, ok, residual)
        if (.not. ok) error stop 'accuracy inverse'
        g = g0
        opt%tol = 1.0e-11_dp
        call initialize(pm, opt, y0, ref)
        call advance(pm, opt, ref, 50.0_dp)
        call particle(pm, opt, ref, yf, ok)
        call gc_solve(g, 50.0_dp)
        call lift(pm, g, order, yg, ok)
        result(1) = norm2(yg(1:3) - yf(1:3))/speed
        result(2) = norm2(yg(4:6) - yf(4:6))/speed
        energy = 0.5_dp*sum(y0(4:6)**2)
        result(4) = 0.5_dp*sum(yg(4:6)**2)/energy - 1
        call gc_solve(g, -50.0_dp)
        call lift(pm, g, order, yg, ok)
        result(3) = norm2(yg - y0)/speed
    end subroutine
    subroutine gc_solve(g, span)
        real(dp), intent(inout) :: g(6)
        real(dp), intent(in) :: span
        real(dp) :: t, h, err, gn(6)
        t = 0.0_dp
        h = 0.2_dp
        do while (t < abs(span) - 1.0e-12_dp)
            h = min(h, abs(span) - t)
            call dp_step(pm, g, 2, .false., order, sign(h, span), speed, gn, err)
            err = err/1.0e-11_dp
            if (err > 1.0_dp) then
                h = h*max(0.2_dp, 0.9_dp*err**(-0.2_dp))
                cycle
            end if
            g = gn
            t = t + h
            h = h*min(5.0_dp, 0.9_dp*max(err, 1.0e-12_dp)**(-0.2_dp))
        end do
    end subroutine
end program
