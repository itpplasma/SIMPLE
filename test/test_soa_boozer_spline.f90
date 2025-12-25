program test_soa_boozer_spline
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, &
                          splint_boozer_coord_many
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm
    use simple_main, only: init_field
    use simple, only: tracer_t

    implicit none

    real(dp), parameter :: RTOL = 1.0d-12
    type(tracer_t) :: norb
    character(len=256) :: config_file
    integer :: i, j, npts, nerrors
    real(dp), allocatable :: r(:), th(:), ph(:)
    real(dp), allocatable :: A_theta(:), A_phi(:), dA_theta_dr(:), dA_phi_dr(:)
    real(dp), allocatable :: d2A_phi_dr2(:), d3A_phi_dr3(:)
    real(dp), allocatable :: B_vartheta_B(:), dB_vartheta_B(:), d2B_vartheta_B(:)
    real(dp), allocatable :: B_varphi_B(:), dB_varphi_B(:), d2B_varphi_B(:)
    real(dp), allocatable :: Bmod_B(:), dBmod_B(:,:), d2Bmod_B(:,:)
    real(dp) :: r_s, th_s, ph_s
    real(dp) :: A_theta_s, A_phi_s, dA_theta_dr_s, dA_phi_dr_s
    real(dp) :: d2A_phi_dr2_s, d3A_phi_dr3_s
    real(dp) :: B_vartheta_B_s, dB_vartheta_B_s, d2B_vartheta_B_s
    real(dp) :: B_varphi_B_s, dB_varphi_B_s, d2B_varphi_B_s
    real(dp) :: Bmod_B_s, dBmod_B_s(3), d2Bmod_B_s(6)
    real(dp) :: B_r_s, dB_r_s(3), d2B_r_s(6)

    print *, '=========================================='
    print *, 'Testing splint_boozer_coord_many'
    print *, '=========================================='

    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, 0)
    call get_boozer_coordinates()

    npts = 100
    allocate(r(npts), th(npts), ph(npts))
    allocate(A_theta(npts), A_phi(npts), dA_theta_dr(npts), dA_phi_dr(npts))
    allocate(d2A_phi_dr2(npts), d3A_phi_dr3(npts))
    allocate(B_vartheta_B(npts), dB_vartheta_B(npts), d2B_vartheta_B(npts))
    allocate(B_varphi_B(npts), dB_varphi_B(npts), d2B_varphi_B(npts))
    allocate(Bmod_B(npts), dBmod_B(3, npts), d2Bmod_B(6, npts))

    do i = 1, npts
        r(i) = 0.1d0 + 0.8d0 * real(i-1, dp) / real(npts-1, dp)
        th(i) = 6.28318530718d0 * real(i, dp) / real(npts, dp)
        ph(i) = 3.14159265359d0 * real(i, dp) / real(npts, dp)
    end do

    call splint_boozer_coord_many(npts, r, th, ph, &
        A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
        d2A_phi_dr2, d3A_phi_dr3, &
        B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
        B_varphi_B, dB_varphi_B, d2B_varphi_B, &
        Bmod_B, dBmod_B, d2Bmod_B)

    nerrors = 0
    do i = 1, npts
        r_s = 0.1d0 + 0.8d0 * real(i-1, dp) / real(npts-1, dp)
        th_s = 6.28318530718d0 * real(i, dp) / real(npts, dp)
        ph_s = 3.14159265359d0 * real(i, dp) / real(npts, dp)

        call splint_boozer_coord(r_s, th_s, ph_s, &
            A_theta_s, A_phi_s, dA_theta_dr_s, dA_phi_dr_s, &
            d2A_phi_dr2_s, d3A_phi_dr3_s, &
            B_vartheta_B_s, dB_vartheta_B_s, d2B_vartheta_B_s, &
            B_varphi_B_s, dB_varphi_B_s, d2B_varphi_B_s, &
            Bmod_B_s, dBmod_B_s, d2Bmod_B_s, B_r_s, dB_r_s, d2B_r_s)

        if (reldiff(A_theta(i), A_theta_s) > RTOL) then
            print *, 'Mismatch A_theta at i=', i, A_theta(i), A_theta_s
            nerrors = nerrors + 1
        end if
        if (reldiff(A_phi(i), A_phi_s) > RTOL) then
            print *, 'Mismatch A_phi at i=', i, A_phi(i), A_phi_s
            nerrors = nerrors + 1
        end if
        if (reldiff(Bmod_B(i), Bmod_B_s) > RTOL) then
            print *, 'Mismatch Bmod_B at i=', i, Bmod_B(i), Bmod_B_s
            nerrors = nerrors + 1
        end if
        if (reldiff(B_vartheta_B(i), B_vartheta_B_s) > RTOL) then
            print *, 'Mismatch B_vartheta_B at i=', i, B_vartheta_B(i), B_vartheta_B_s
            nerrors = nerrors + 1
        end if
        if (reldiff(B_varphi_B(i), B_varphi_B_s) > RTOL) then
            print *, 'Mismatch B_varphi_B at i=', i, B_varphi_B(i), B_varphi_B_s
            nerrors = nerrors + 1
        end if
        do j = 1, 3
            if (reldiff(dBmod_B(j,i), dBmod_B_s(j)) > RTOL) then
                print *, 'Mismatch dBmod_B at i=', i, 'j=', j
                nerrors = nerrors + 1
            end if
        end do
        do j = 1, 6
            if (reldiff(d2Bmod_B(j,i), d2Bmod_B_s(j)) > RTOL) then
                print *, 'Mismatch d2Bmod_B at i=', i, 'j=', j
                nerrors = nerrors + 1
            end if
        end do
    end do

    if (nerrors > 0) then
        print *, 'FAILED: ', nerrors, ' errors found'
        stop 1
    end if

    print *, 'PASSED: splint_boozer_coord_many matches single-point version'

contains

    pure function reldiff(a, b) result(rd)
        real(dp), intent(in) :: a, b
        real(dp) :: rd, denom
        denom = max(abs(a), abs(b), 1.0d-30)
        rd = abs(a - b) / denom
    end function reldiff

end program test_soa_boozer_spline
