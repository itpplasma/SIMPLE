program diag_darboux_axis_fields
    !> Near-axis behaviour of p_theta = A_theta/ro0 + h_theta v_par for a
    !> SIMPLE field_can source on a VMEC file: the Darboux chart
    !> (X, Y) = sqrt(2 |p_theta|) (cos, sin) theta needs A_theta -> 0 and
    !> h_theta -> 0 at the axis with a definite sign of dp_theta/dr.
    !> Usage: diag_darboux_axis_fields.x <wout.nc> <isw_field_type>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: field_input, coord_input
    use velo_mod, only: isw_field_type
    use field_can_mod, only: field_can_t, eval_field => evaluate

    implicit none

    type(tracer_t) :: norb
    type(field_can_t) :: f
    character(len=256) :: wout, arg
    real(dp) :: r, th, athmax, hthmax, dathmin, dathmax, cone, b0
    real(dp) :: bmin, bmax, hmin, hmax
    integer :: ir, ith, ftype

    call get_command_argument(1, wout)
    call get_command_argument(2, arg)
    read (arg, *) ftype
    isw_field_type = ftype
    field_input = trim(wout)
    coord_input = trim(wout)
    call init_field(norb, trim(wout), 5, 5, 5, 1)

    print '(a)', '      r      max|Ath|/r  max|hth|/r  min dAth/dr  max dAth/dr'// &
        '  <dB/dr>*2rho/B0  dB_th/B0  dhph_th/hph'
    do ir = 1, 7
        r = 10.0_dp**(-ir - 1)
        athmax = 0.0_dp
        hthmax = 0.0_dp
        dathmin = huge(1.0_dp)
        dathmax = -huge(1.0_dp)
        cone = 0.0_dp
        b0 = 0.0_dp
        bmin = huge(1.0_dp)
        bmax = -huge(1.0_dp)
        hmin = huge(1.0_dp)
        hmax = -huge(1.0_dp)
        do ith = 0, 15
            th = 6.283185307179586_dp*ith/16.0_dp
            call eval_field(f, r, th, 0.3_dp, 0)
            athmax = max(athmax, abs(f%Ath)/r)
            hthmax = max(hthmax, abs(f%hth)/r)
            dathmin = min(dathmin, f%dAth(1))
            dathmax = max(dathmax, f%dAth(1))
            cone = cone + f%dBmod(1)/16.0_dp
            b0 = b0 + f%Bmod/16.0_dp
            bmin = min(bmin, f%Bmod)
            bmax = max(bmax, f%Bmod)
            hmin = min(hmin, f%hph)
            hmax = max(hmax, f%hph)
        end do
        ! m = 0 cone: theta-averaged dB/dr times 2 sqrt(r) (r = s charts) is
        ! the coefficient of a term linear in rho, which is not smooth in
        ! (X, Y); a regular field has it -> 0 like rho.
        ! Poloidal variation at fixed r: a regular field has it ~ rho; a
        ! nonzero limit is a theta-dependent axis value (a 1/rho force).
        print '(es9.1,7es12.4)', r, athmax, hthmax, dathmin, dathmax, &
            cone*2.0_dp*sqrt(r)/b0, (bmax - bmin)/b0, (hmax - hmin)/abs(hmax)
    end do
end program diag_darboux_axis_fields
