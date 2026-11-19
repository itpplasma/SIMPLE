program test_nbi_field_direction_axis
    !> The birth field direction stays finite and correct within one
    !> finite-difference step of the magnetic axis.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use simple_api, only: simple_reset_config, simple_init
    use params, only: field_input, coord_input, startmode, isw_field_type, integmode, &
                      nper, npoiper, ntimstep, deterministic, n_d, n_e, &
                      particle_energy_eV, ntestpart, contr_pp
    use reference_coordinates, only: ref_coords
    use nbi_birth, only: field_direction
    use timing, only: init_timer
    use spline_vmec_sub, only: vmec_field
    implicit none
    real(dp), parameter :: s_values(4) = [0.0_dp, 1.0e-7_dp, 4.0e-6_dp, 1.0e-3_dp]
    real(dp) :: u(3), bhat(3), bmod, e_cov(3, 3), bcart(3), worst
    real(dp) :: at, ap, dat, dap, iota, sqg, alam, dl_ds, dl_dt, dl_dp, bctr_t, bctr_p
    real(dp) :: bco_s, bco_t, bco_p
    character(1024) :: filename
    integer :: i, j

    call get_command_argument(1, filename)
    call init_timer
    call simple_reset_config
    field_input = trim(filename)
    coord_input = trim(filename)
    isw_field_type = 1
    integmode = 0
    startmode = 1
    deterministic = .true.
    nper = 1
    npoiper = 4
    ntimstep = 2
    n_d = 2
    n_e = 1
    particle_energy_eV = 1.0e5_dp
    ntestpart = 4
    contr_pp = -1.0e10_dp
    call simple_init

    worst = 0.0_dp
    do i = 2, size(s_values)
        do j = 0, 3
            u = [s_values(i), 1.5_dp*j + 0.2_dp, 0.3_dp]
            call field_direction(u, bhat, bmod)
            if (.not. all(ieee_is_finite(bhat))) error stop 'non-finite field direction'
            ! Independent direction from VMEC contravariant components on
            ! libneo's analytic covariant basis.
            call vmec_field(u(1), u(2), u(3), at, ap, dat, dap, iota, sqg, alam, &
                            dl_ds, dl_dt, dl_dp, bctr_t, bctr_p, bco_s, bco_t, bco_p)
            call ref_coords%covariant_basis(u, e_cov)
            bcart = (bctr_t - bctr_p*dl_dp)/(1.0_dp + dl_dt)*e_cov(:, 2) + &
                    bctr_p*e_cov(:, 3)
            worst = max(worst, norm2(bhat - bcart/norm2(bcart)))
        end do
    end do
    print '(a,es10.3)', 'worst direction difference near the axis ', worst
    if (worst > 1.0e-3_dp) error stop 'field direction differs near the axis'
    print '(a)', 'PASS'
end program test_nbi_field_direction_axis
