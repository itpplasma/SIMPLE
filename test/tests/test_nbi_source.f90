module nbi_slab
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nbi_source, only: nbi_plasma_t
    implicit none
contains
    function slab(x) result(p)
        !> Uniform plasma for 0.5 <= x <= 2.5 m.
        real(dp), intent(in) :: x(3)
        type(nbi_plasma_t) :: p
        p%inside = x(1) >= 0.5_dp .and. x(1) <= 2.5_dp
        p%ne = 3.0e19_dp
        p%te = 1.0e3_dp
        p%zeff = 1.0_dp
    end function slab

    function peaked(x) result(p)
        !> n_e = 5e19 (1 - 0.96 u^2), T_e = 3 keV (1 - 0.67 u^2), u = x - 1.5,
        !> carbon at Z_eff = 1.6, for 0.5 <= x <= 2.5 m.
        real(dp), intent(in) :: x(3)
        type(nbi_plasma_t) :: p
        p%inside = x(1) >= 0.5_dp .and. x(1) <= 2.5_dp
        p%ne = 5.0e19_dp*(1.0_dp - 0.96_dp*(x(1) - 1.5_dp)**2)
        p%te = 3.0e3_dp*(1.0_dp - 0.67_dp*(x(1) - 1.5_dp)**2)
        p%zeff = 1.6_dp
        p%z_impurity = 6
    end function peaked
end module nbi_slab

program test_nbi_source
    !> Slab attenuation against exp(-n_e sigma L), attenuation in a peaked
    !> profile against an independent quadrature of the optical depth, and the
    !> guiding-centre shift of a gyrating ion in a uniform field.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nbi_source
    use nbi_slab, only: slab, peaked
    use beam_stopping, only: suzuki_stopping
    implicit none
    type(nbi_beamlet_t) :: b
    type(nbi_segment_t), allocatable :: seg(:)
  real(dp) :: injected, shine, sigma, lambda, deposited, l, expected, x(3), v(3), xgc(3)
    real(dp) :: tau
    integer :: i, j
    character(16) :: mode

    call get_command_argument(1, mode)
    b%mass = 2
    b%energy_ev = 60.0e3_dp
    b%power_w = 1.0e6_dp
    b%fraction = [0.6_dp, 0.3_dp, 0.1_dp]
    b%origin = 0.0_dp
    b%direction = [1.0_dp, 0.0_dp, 0.0_dp]
    if (trim(mode) == 'truncated') then
        ! Must stop: the path ends 0.5 m before the far side of the slab.
        call beamlet_deposition(b, 1, slab, 1.0e-3_dp, 2.0_dp, seg, injected, shine)
        print '(a)', 'truncated path accepted'
        stop
    end if
    call beamlet_deposition(b, 1, slab, 1.0e-3_dp, 3.0_dp, seg, injected, shine)
    sigma = suzuki_stopping(30.0_dp, 3.0e19_dp, 1.0e3_dp, 2, 1.0_dp, 0)
    lambda = 3.0e19_dp*sigma
    if (abs(injected/(0.6_dp*1.0e6_dp/(60.0e3_dp*ev_joule)) - 1.0_dp) > 1.0e-12_dp) &
        error stop 'injected neutral rate'
    if (abs(shine/(injected*exp(-2.0_dp*lambda)) - 1.0_dp) > 1.0e-10_dp) &
        error stop 'shine-through differs from exp(-n_e sigma L)'
    deposited = sum(seg%rate)
    if (abs((deposited + shine)/injected - 1.0_dp) > 1.0e-12_dp) &
        error stop 'ionized plus shine-through differs from injected'
    ! Ionization in the first 0.5 m of plasma.
    l = 0.0_dp
    do i = 1, size(seg)
        if (seg(i)%l1 <= 1.0_dp + 1.0e-12_dp) l = l + seg(i)%rate
    end do
    expected = injected*(1.0_dp - exp(-0.5_dp*lambda))
    if (abs(l/expected - 1.0_dp) > 1.0e-10_dp) error stop 'deposition profile not exponential'
    print '(a,f8.5,a,f8.5)', 'shine fraction ', shine/injected, ', analytic ', exp(-2*lambda)

    ! Peaked profile, half-energy component: composite Simpson rule on 2e5
    ! intervals for tau = int n_e sigma dx. The midpoint error of the 1 mm
    ! path step is below 1e-7 of tau.
    do j = 2, 3
        call beamlet_deposition(b, j, peaked, 1.0e-3_dp, 3.0_dp, seg, injected, shine)
        tau = simpson_tau(b%energy_ev/j*1.0e-3_dp/b%mass, 200000)
        if (abs(log(injected/shine)/tau - 1.0_dp) > 1.0e-6_dp) &
            error stop 'optical depth in a peaked profile differs from quadrature'
        if (abs((sum(seg%rate) + shine)/injected - 1.0_dp) > 1.0e-12_dp) &
           error stop 'peaked profile: ionized plus shine-through differs from injected'
    end do

    ! Positive ion moving along -y in B along +z gyrates about x - rho e_x.
    x = [1.0_dp, 2.0_dp, 3.0_dp]
    v = [0.0_dp, -2.0_dp, 0.0_dp]
    xgc = guiding_centre_shift(x, v, [0.0_dp, 0.0_dp, 1.0_dp], 4.0_dp)
    if (norm2(xgc - [0.5_dp, 2.0_dp, 3.0_dp]) > 1.0e-14_dp) error stop 'guiding-centre shift'
    print '(a)', 'PASS'

contains

    function simpson_tau(e_kev_amu, n) result(tau)
        real(dp), intent(in) :: e_kev_amu
        integer, intent(in) :: n
        real(dp) :: tau, h, w
        type(nbi_plasma_t) :: p
        integer :: k
        h = 2.0_dp/n
        tau = 0.0_dp
        do k = 0, n
            p = peaked([0.5_dp + k*h, 0.0_dp, 0.0_dp])
            w = merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(k, 2) == 1), k == 0 .or. k == n)
      tau = tau + w*p%ne*suzuki_stopping(e_kev_amu, p%ne, p%te, 2, p%zeff, p%z_impurity)
        end do
        tau = tau*h/3.0_dp
    end function simpson_tau
end program test_nbi_source
