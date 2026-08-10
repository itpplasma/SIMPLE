program test_splined_kernel_bench
    !> Baseline measurement instrument for the first generated consumer
    !> (splined field evaluation, see DOC/first-generated-kernel.md).
    !>
    !> This program measures the hand-written splined field kernel on a
    !> realistic NCSX workload and reports the two numbers that SIMPLE itself
    !> can measure for the fortsym statement:
    !>
    !>   - T: wall-clock time per evaluation of splined_evaluate and
    !>        splined_evaluate_with_der (the hottest path in the code), and
    !>   - X ulp: accuracy of the splined field path against the direct
    !>        Biot-Savart coils field, which is the high-precision reference
    !>        for this kernel.
    !>
    !> The arithmetic-count numbers N_sym (mathematics) and N_emit (generated
    !> code) come from the lazy-fortran/fortsym generation toolchain, and
    !> N_machine from the compiler; they are filled in when a generated kernel
    !> is substituted for the hand-written one. This program is the
    !> production-side gate the generated kernel must beat.
    !>
    !> Golden records are not touched: the measurement is informational and
    !> only sanity-checks that the spline path produces finite values within
    !> the spline-interpolation accuracy of the reference field.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: splined

    real(dp) :: dummy
    real(dp) :: phi_period

    integer, parameter :: n_timing = 20000      ! evaluations for the timing loop
    integer, parameter :: n_rep = 5             ! repeats, best (min) time kept
    integer, parameter :: n_acc = 400           ! points for the accuracy sweep

    real(dp) :: t_eval, t_eval_der
    real(dp) :: max_ulp

    call init_vmec('wout_ncsx.nc', 5, 5, 5, dummy)
    call init_reference_coordinates('wout_ncsx.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, splined)

    phi_period = twopi/real(nper, dp)

    ! --- T: end-to-end timing of the hand-written splined field kernel ---
    call time_kernel(splined, phi_period, n_timing, n_rep, t_eval, t_eval_der)

    ! --- X ulp: accuracy of the spline path vs the direct coils field ---
    call measure_accuracy(raw_coils, splined, phi_period, n_acc, max_ulp)

    print *, ''
    print *, '=================================================================='
    print *, 'Splined field kernel baseline instrument (DOC/first-generated-kernel.md)'
    print *, '------------------------------------------------------------------'
    print *, 'kernel: splined_evaluate / splined_evaluate_with_der'
    print *, 'reference: direct Biot-Savart coils field (wout_ncsx.nc)'
    print *, '------------------------------------------------------------------'
    print '(a,es12.4,a)', ' T (value only)          = ', t_eval, ' s/eval'
    print '(a,es12.4,a)', ' T (value + derivative)  = ', t_eval_der, ' s/eval'
    print '(a,es12.4,a)', ' X ulp (max, Bmod)       = ', max_ulp, ' ulp vs direct field'
    print *, '------------------------------------------------------------------'
    print *, 'statement: for this kernel the arithmetic the mathematics'
    print *, 'requires is N_sym (fortsym), the generated code emits N_emit'
    print *, '(fortsym), the compiler produces N_machine, the hardware'
    print *, 'delivers T above, and the accuracy is X ulp above against a'
    print *, 'high-precision reference.'
    print *, '=================================================================='

    ! Sanity gate: finite, positive timing and finite spline accuracy. The
    ! spline is an approximation of the reference field, so the ulp figure is
    ! the spline-interpolation accuracy, not a rounding-error bound; we only
    ! require it to be finite and reasonable (< 1e-3 relative).
    if (.not. (t_eval > 0.0_dp .and. t_eval_der > 0.0_dp)) then
        print *, 'FAILED: non-positive kernel timing'
        error stop 1
    end if
    if (.not. (max_ulp >= 0.0_dp .and. max_ulp < 1.0e13_dp)) then
        print *, 'FAILED: spline accuracy out of expected range'
        error stop 1
    end if

    print *, 'PASSED: splined kernel baseline instrument'
end program test_splined_kernel_bench


subroutine time_kernel(f, phi_period, n_timing, n_rep, t_eval, t_eval_der)
    !> Best-of-n wall-clock time per evaluation of splined_evaluate and
    !> splined_evaluate_with_der. Outputs are accumulated into a checksum so
    !> the compiler cannot elide the kernel calls.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_splined, only: splined_field_t

    type(splined_field_t), intent(in) :: f
    real(dp), intent(in) :: phi_period
    integer, intent(in) :: n_timing, n_rep
    real(dp), intent(out) :: t_eval, t_eval_der

    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    real(dp) :: dAcov(3, 3), dhcov(3, 3), dBmod(3)
    real(dp) :: checksum, dt
    integer(selected_int_kind(18)) :: t0, t1, cr
    integer :: i, rep

    t_eval = huge(1.0_dp)
    t_eval_der = huge(1.0_dp)

    ! Warm up
    x = [0.5_dp, 1.0_dp, 1.0_dp]
    call f%evaluate(x, Acov, hcov, Bmod)
    call f%evaluate_with_der(x, Acov, hcov, Bmod, dAcov, dhcov, dBmod)

    do rep = 1, n_rep
        checksum = 0.0_dp
        call system_clock(count=t0)
        do i = 1, n_timing
            x(1) = 0.02_dp + 0.96_dp*real(mod(i*7, 1000), dp)/1000.0_dp
            x(2) = phi_period*real(mod(i*13, 1000), dp)/1000.0_dp
            x(3) = phi_period*real(mod(i*29, 1000), dp)/1000.0_dp
            call f%evaluate(x, Acov, hcov, Bmod)
            checksum = checksum + Bmod
        end do
        call system_clock(count=t1, count_rate=cr)
        dt = real(t1 - t0, dp)/real(max(cr,1), dp)
        t_eval = min(t_eval, dt/real(n_timing, dp))

        checksum = 0.0_dp
        call system_clock(count=t0)
        do i = 1, n_timing
            x(1) = 0.02_dp + 0.96_dp*real(mod(i*7, 1000), dp)/1000.0_dp
            x(2) = phi_period*real(mod(i*13, 1000), dp)/1000.0_dp
            x(3) = phi_period*real(mod(i*29, 1000), dp)/1000.0_dp
            call f%evaluate_with_der(x, Acov, hcov, Bmod, dAcov, dhcov, dBmod)
            checksum = checksum + Bmod
        end do
        call system_clock(count=t1, count_rate=cr)
        dt = real(t1 - t0, dp)/real(max(cr,1), dp)
        t_eval_der = min(t_eval_der, dt/real(n_timing, dp))
    end do

    ! Keep the checksum observable to the optimizer.
    if (checksum == huge(1.0_dp)) print *, checksum

end subroutine time_kernel


subroutine measure_accuracy(raw, splined, phi_period, n_acc, max_ulp)
    !> Max ulp deviation of splined Bmod from the direct coils Bmod at the
    !> same physical point. Splined coordinates are (r=sqrt(s), theta, zeta);
    !> the physical point is obtained by s = r^2 and ref_coords%evaluate_cart.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_coils, only: coils_field_t
    use field_splined, only: splined_field_t
    use reference_coordinates, only: ref_coords
    use util, only: twopi

    type(coils_field_t), intent(in) :: raw
    type(splined_field_t), intent(in) :: splined
    real(dp), intent(in) :: phi_period
    integer, intent(in) :: n_acc
    real(dp), intent(out) :: max_ulp

    real(dp) :: x(3), u(3), xcart(3)
    real(dp) :: Acov(3), hcov(3), Bmod
    real(dp) :: Acov_ref(3), hcov_ref(3), Bmod_ref
    real(dp) :: ulp
    integer :: i

    max_ulp = 0.0_dp
    do i = 1, n_acc
        x(1) = 0.05_dp + 0.9_dp*real(mod(13*i, 997), dp)/997.0_dp
        x(2) = 0.1_dp + (twopi - 0.2_dp)*real(mod(29*i, 991), dp)/991.0_dp
        x(3) = 0.1_dp + (phi_period - 0.2_dp)*real(mod(31*i, 983), dp)/983.0_dp

        call splined%evaluate(x, Acov, hcov, Bmod)

        u(1) = x(1)**2
        u(2) = x(2)
        u(3) = x(3)
        call ref_coords%evaluate_cart(u, xcart)
        call raw%evaluate(xcart, Acov_ref, hcov_ref, Bmod_ref)

        ulp = abs(Bmod - Bmod_ref)/max(spacing(Bmod_ref), tiny(1.0_dp))
        max_ulp = max(max_ulp, ulp)
    end do
end subroutine measure_accuracy
