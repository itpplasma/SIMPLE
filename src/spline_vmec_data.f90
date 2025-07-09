module spline_vmec_sub
    use iso_fortran_env, only: dp => real64
    use spl_three_to_five_sub
    implicit none
contains
    subroutine spline_vmec_data
        real(dp), dimension(:, :), allocatable :: almnc_rho, rmnc_rho, zmnc_rho
        real(dp), dimension(:, :), allocatable :: almns_rho, rmns_rho, zmns_rho

        call initialize_vmec_data_and_arrays
        call perform_axis_healing(almnc_rho, rmnc_rho, zmnc_rho, almns_rho, rmns_rho, zmns_rho)
        call setup_poloidal_flux_splines
        call setup_angular_grid_and_fourier_synthesis(rmnc_rho, zmnc_rho, almnc_rho, &
                                                      rmns_rho, zmns_rho, almns_rho)
        call spline_over_phi
        call spline_over_theta
        call spline_over_s

        deallocate (almnc_rho, rmnc_rho, zmnc_rho, almns_rho, rmns_rho, zmns_rho)
    end subroutine spline_vmec_data

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine initialize_vmec_data_and_arrays
        use new_vmec_stuff_mod, only: ns_A, ns_s, ns_tp, rmnc, zmns, almns, rmns, zmnc, almnc, &
                                      aiota, phi, sps, axm, axn, s, nsurfm, nstrm, kpar
        use vmecin_sub, only: vmecin
        use vmec_alloc_sub, only: new_allocate_vmec_stuff
        use vector_potentail_mod, only: ns, hs, torflux

        print *, 'Splining VMEC data: ns_A = ', ns_A, '  ns_s = ', ns_s, '  ns_tp = ', ns_tp

        call new_allocate_vmec_stuff

        call vmecin(rmnc, zmns, almns, rmns, zmnc, almnc, aiota, phi, sps, axm, axn, s, &
                    nsurfm, nstrm, kpar, torflux)

        ns = kpar + 1
        hs = s(2) - s(1)
    end subroutine initialize_vmec_data_and_arrays

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine perform_axis_healing(almnc_rho, rmnc_rho, zmnc_rho, almns_rho, rmns_rho, zmns_rho)
        use new_vmec_stuff_mod, only: rmnc, zmns, almns, rmns, zmnc, almnc, &
                                      axm, axn, nstrm, old_axis_healing_boundary
        use vector_potentail_mod, only: ns

        real(dp), dimension(:, :), allocatable, intent(out) :: almnc_rho, rmnc_rho, zmnc_rho
        real(dp), dimension(:, :), allocatable, intent(out) :: almns_rho, rmns_rho, zmns_rho
        integer :: i, m, nrho, nheal, iunit_hs

        nrho = ns
        allocate (almnc_rho(nstrm, 0:nrho - 1), rmnc_rho(nstrm, 0:nrho - 1), zmnc_rho(nstrm, 0:nrho - 1))
        allocate (almns_rho(nstrm, 0:nrho - 1), rmns_rho(nstrm, 0:nrho - 1), zmns_rho(nstrm, 0:nrho - 1))

        iunit_hs = 1357
        open (iunit_hs, file='healaxis.dat')

        do i = 1, nstrm
            m = nint(abs(axm(i)))

            if (old_axis_healing_boundary) then
                nheal = min(m, 4)
            else
                call determine_nheal_for_axis(m, ns, rmnc(i, :), nheal)
                write (iunit_hs, *) 'm = ', m, ' n = ', nint(abs(axn(i))), ' skipped ', nheal, ' / ', ns
            end if

            call s_to_rho_healaxis(m, ns, nrho, nheal, rmnc(i, :), rmnc_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, zmnc(i, :), zmnc_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, almnc(i, :), almnc_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, rmns(i, :), rmns_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, zmns(i, :), zmns_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, almns(i, :), almns_rho(i, :))
        end do

        close (iunit_hs)
    end subroutine perform_axis_healing


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine setup_poloidal_flux_splines
        use new_vmec_stuff_mod, only: ns_A, aiota
        use vector_potentail_mod, only: ns, hs, torflux, sA_phi
        use spl_three_to_five_sub, only: spl_reg

        integer :: i, is, k
        real(dp), dimension(:, :), allocatable :: splcoe

        allocate (splcoe(0:ns_A, ns))

        splcoe(0, :) = aiota

        call spl_reg(ns_A - 1, ns, hs, splcoe(0:ns_A - 1, :))

        do i = ns_A, 1, -1
            splcoe(i, :) = splcoe(i - 1, :)/dble(i)
        end do

        splcoe(0, 1) = 0.d0
        do is = 1, ns - 1
            splcoe(0, is + 1) = splcoe(ns_A, is)
            do k = ns_A - 1, 0, -1
                splcoe(0, is + 1) = splcoe(k, is) + hs*splcoe(0, is + 1)
            end do
        end do

        if (allocated(sA_phi)) deallocate (sA_phi)
        allocate (sA_phi(ns_A + 1, ns))
        do k = 0, ns_A
            sA_phi(k + 1, :) = -torflux*splcoe(k, :)
        end do

        deallocate (splcoe)
    end subroutine setup_poloidal_flux_splines

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine setup_angular_grid_and_fourier_synthesis(rmnc_rho, zmnc_rho, almnc_rho, &
                                                        rmns_rho, zmns_rho, almns_rho)
        use new_vmec_stuff_mod, only: axm, axn, multharm, nper, nstrm, &
                                      n_theta, n_phi, h_theta, h_phi, sR, sZ, slam, ns_s, ns_tp
        use vector_potentail_mod, only: ns
        use vmec_alloc_sub, only: new_deallocate_vmec_stuff

        real(dp), dimension(nstrm, 0:ns - 1), intent(in) :: rmnc_rho, zmnc_rho, almnc_rho
        real(dp), dimension(nstrm, 0:ns - 1), intent(in) :: rmns_rho, zmns_rho, almns_rho

        integer :: i, m, n, is, i_theta, i_phi, m_max, n_max
        integer :: nsize_exp_imt, nsize_exp_inp, iexpt, iexpp
        real(dp) :: twopi, cosphase, sinphase
        complex(8) :: base_exp_imt, base_exp_inp, base_exp_inp_inv, expphase
        complex(8), dimension(:), allocatable :: exp_imt, exp_inp

        ! Setup angular grid
        m_max = nint(maxval(axm))
        n_max = nint(maxval(axn))

        print *, 'VMEC ns = ', ns, ' m_max = ', m_max, ' n_max = ', n_max

        n_theta = m_max*multharm + 1
        n_phi = n_max*multharm + 1
        twopi = 8.d0*atan2(1.d0, 1.d0)
        h_theta = twopi/dble(n_theta - 1)
        h_phi = twopi/dble((n_phi - 1)*nper)

        ! Setup exponential arrays
        nsize_exp_imt = (n_theta - 1)*m_max
        nsize_exp_inp = (n_phi - 1)*n_max

        allocate (exp_imt(0:nsize_exp_imt), exp_inp(-nsize_exp_inp:nsize_exp_inp))

        base_exp_imt = exp(cmplx(0.d0, h_theta, kind=kind(0d0)))
        base_exp_inp = exp(cmplx(0.d0, h_phi, kind=kind(0d0)))
        base_exp_inp_inv = (1.d0, 0.d0)/base_exp_inp
        exp_imt(0) = (1.d0, 0.d0)
        exp_inp(0) = (1.d0, 0.d0)

        do i = 1, nsize_exp_imt
            exp_imt(i) = exp_imt(i - 1)*base_exp_imt
        end do

        do i = 1, nsize_exp_inp
            exp_inp(i) = exp_inp(i - 1)*base_exp_inp
            exp_inp(-i) = exp_inp(1 - i)*base_exp_inp_inv
        end do

        ! Allocate spline arrays
        if (allocated(sR)) deallocate (sR)
        allocate (sR(ns_s + 1, ns_tp + 1, ns_tp + 1, ns, n_theta, n_phi))
        if (allocated(sZ)) deallocate (sZ)
        allocate (sZ(ns_s + 1, ns_tp + 1, ns_tp + 1, ns, n_theta, n_phi))
        if (allocated(slam)) deallocate (slam)
        allocate (slam(ns_s + 1, ns_tp + 1, ns_tp + 1, ns, n_theta, n_phi))

        sR(1, 1, 1, :, :, :) = 0.d0
        sZ(1, 1, 1, :, :, :) = 0.d0
        slam(1, 1, 1, :, :, :) = 0.d0

        ! Fourier synthesis to real space
!$omp parallel private(m, n, i_theta, i_phi, i, is, iexpt, iexpp, expphase, cosphase, sinphase)
!$omp do
        do i_theta = 1, n_theta
            do i_phi = 1, n_phi
                do i = 1, nstrm
                    m = nint(axm(i))
                    n = nint(axn(i))
                    iexpt = m*(i_theta - 1)
                    iexpp = n*(i_phi - 1)
                    expphase = exp_imt(iexpt)*exp_inp(-iexpp)
                    cosphase = dble(expphase)
                    sinphase = aimag(expphase)
                    do is = 1, ns
                        sR(1, 1, 1, is, i_theta, i_phi) = sR(1, 1, 1, is, i_theta, i_phi) &
                                                          + rmnc_rho(i, is - 1)*cosphase &
                                                          + rmns_rho(i, is - 1)*sinphase
                        sZ(1, 1, 1, is, i_theta, i_phi) = sZ(1, 1, 1, is, i_theta, i_phi) &
                                                          + zmnc_rho(i, is - 1)*cosphase &
                                                          + zmns_rho(i, is - 1)*sinphase
                        slam(1, 1, 1, is, i_theta, i_phi) = slam(1, 1, 1, is, i_theta, i_phi) &
                                                            + almnc_rho(i, is - 1)*cosphase &
                                                            + almns_rho(i, is - 1)*sinphase
                    end do
                end do
            end do
        end do
!$omp end do
!$omp barrier
!$omp single
        call new_deallocate_vmec_stuff
!$omp end single
!$omp end parallel

        deallocate (exp_imt, exp_inp)
    end subroutine setup_angular_grid_and_fourier_synthesis

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine spline_over_phi
        use new_vmec_stuff_mod, only: ns_s, ns_tp, n_theta, n_phi, h_phi, sR, sZ, slam
        use vector_potentail_mod, only: ns
        use spl_three_to_five_sub, only: spl_per

        integer :: is, i_theta, k
        real(dp), dimension(:, :), allocatable :: splcoe

        if (n_phi == 1) then
            print *, 'Spline not supported for a Phi period of 1, exiting...'
            call exit(-1)
        end if

!$omp parallel private(is, i_theta, k, splcoe)
        allocate (splcoe(0:ns_tp, n_phi))
!$omp do
        do is = 1, ns
            do i_theta = 1, n_theta
                splcoe(0, :) = sR(1, 1, 1, is, i_theta, :)
                call spl_per(ns_tp, n_phi, h_phi, splcoe)
                do k = 1, ns_tp
                    sR(1, 1, k + 1, is, i_theta, :) = splcoe(k, :)
                end do

                splcoe(0, :) = sZ(1, 1, 1, is, i_theta, :)
                call spl_per(ns_tp, n_phi, h_phi, splcoe)
                do k = 1, ns_tp
                    sZ(1, 1, k + 1, is, i_theta, :) = splcoe(k, :)
                end do

                splcoe(0, :) = slam(1, 1, 1, is, i_theta, :)
                call spl_per(ns_tp, n_phi, h_phi, splcoe)
                do k = 1, ns_tp
                    slam(1, 1, k + 1, is, i_theta, :) = splcoe(k, :)
                end do
            end do
        end do
!$omp end do
        deallocate (splcoe)
!$omp end parallel
    end subroutine spline_over_phi

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine spline_over_theta
        use new_vmec_stuff_mod, only: ns_s, ns_tp, n_theta, n_phi, h_theta, sR, sZ, slam
        use vector_potentail_mod, only: ns
        use spl_three_to_five_sub, only: spl_per

        integer :: is, i_phi, isp, k
        real(dp), dimension(:, :), allocatable :: splcoe

!$omp parallel private(is, i_phi, isp, k, splcoe)
        allocate (splcoe(0:ns_tp, n_theta))
!$omp do
        do is = 1, ns
            do i_phi = 1, n_phi
                do isp = 1, ns_tp + 1
                    splcoe(0, :) = sR(1, 1, isp, is, :, i_phi)
                    call spl_per(ns_tp, n_theta, h_theta, splcoe)
                    do k = 1, ns_tp
                        sR(1, k + 1, isp, is, :, i_phi) = splcoe(k, :)
                    end do

                    splcoe(0, :) = sZ(1, 1, isp, is, :, i_phi)
                    call spl_per(ns_tp, n_theta, h_theta, splcoe)
                    do k = 1, ns_tp
                        sZ(1, k + 1, isp, is, :, i_phi) = splcoe(k, :)
                    end do

                    splcoe(0, :) = slam(1, 1, isp, is, :, i_phi)
                    call spl_per(ns_tp, n_theta, h_theta, splcoe)
                    do k = 1, ns_tp
                        slam(1, k + 1, isp, is, :, i_phi) = splcoe(k, :)
                    end do
                end do
            end do
        end do
!$omp end do
        deallocate (splcoe)
!$omp end parallel
    end subroutine spline_over_theta

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine spline_over_s
        use new_vmec_stuff_mod, only: ns_s, ns_tp, n_theta, n_phi, sR, sZ, slam
        use vector_potentail_mod, only: ns, hs
        use spl_three_to_five_sub, only: spl_reg

        integer :: i_theta, i_phi, ist, isp, k
        real(dp), dimension(:, :), allocatable :: splcoe

!$omp parallel private(i_theta, i_phi, ist, isp, k, splcoe)
        allocate (splcoe(0:ns_s, ns))
!$omp do
        do i_theta = 1, n_theta
            do i_phi = 1, n_phi
                do ist = 1, ns_tp + 1
                    do isp = 1, ns_tp + 1
                        splcoe(0, :) = sR(1, ist, isp, :, i_theta, i_phi)
                        call spl_reg(ns_s, ns, hs, splcoe)
                        do k = 1, ns_s
                            sR(k + 1, ist, isp, :, i_theta, i_phi) = splcoe(k, :)
                        end do

                        splcoe(0, :) = sZ(1, ist, isp, :, i_theta, i_phi)
                        call spl_reg(ns_s, ns, hs, splcoe)
                        do k = 1, ns_s
                            sZ(k + 1, ist, isp, :, i_theta, i_phi) = splcoe(k, :)
                        end do

                        splcoe(0, :) = slam(1, ist, isp, :, i_theta, i_phi)
                        call spl_reg(ns_s, ns, hs, splcoe)
                        do k = 1, ns_s
                            slam(k + 1, ist, isp, :, i_theta, i_phi) = splcoe(k, :)
                        end do
                    end do
                end do
            end do
        end do
!$omp end do
        deallocate (splcoe)
!$omp end parallel
    end subroutine spline_over_s

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine deallocate_vmec_spline(mode)

        use new_vmec_stuff_mod, only: sR, sZ, slam

        integer :: mode

        if (mode .eq. 0) then
            deallocate (sR, sZ, slam)
        elseif (mode .eq. 1) then
            deallocate (sR, sZ)
        elseif (mode .eq. 2) then
            deallocate (slam)
        else
            print *, 'deallocate_vmec_spline: unknown mode'
        end if

    end subroutine deallocate_vmec_spline

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine normalize_coordinates(s, theta, varphi, ds, is, dtheta, i_theta, dphi, i_phi)
        use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, nper
        use vector_potentail_mod, only: ns, hs

        real(dp), parameter :: twopi = 2.d0*3.14159265358979d0

        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: ds, dtheta, dphi
        integer, intent(out) :: is, i_theta, i_phi

        ds = s/hs
        is = max(0, min(ns - 1, int(ds)))
        ds = (ds - dble(is))*hs
        is = is + 1

        dtheta = modulo(theta, twopi)/h_theta
        i_theta = max(0, min(n_theta - 1, int(dtheta)))
        dtheta = (dtheta - dble(i_theta))*h_theta
        i_theta = i_theta + 1

        dphi = modulo(varphi, twopi/dble(nper))/h_phi
        i_phi = max(0, min(n_phi - 1, int(dphi)))
        dphi = (dphi - dble(i_phi))*h_phi
        i_phi = i_phi + 1

    end subroutine normalize_coordinates

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine interpolate_vector_potential(s, ds, is, A_phi, A_theta, dA_phi_ds, dA_theta_ds)
        use vector_potentail_mod, only: ns, hs, torflux, sA_phi
        use new_vmec_stuff_mod, only: ns_A

        real(dp), intent(in) :: s, ds
        integer, intent(in) :: is
        real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds
        integer :: k

        A_theta = torflux*s
        dA_theta_ds = torflux

        if (.not. allocated(sA_phi)) call spline_vmec_data
        A_phi = sA_phi(ns_A + 1, is)
        dA_phi_ds = 0.d0

        do k = ns_A, 1, -1
            A_phi = sA_phi(k, is) + ds*A_phi
            dA_phi_ds = sA_phi(k + 1, is)*dble(k) + ds*dA_phi_ds
        end do

    end subroutine interpolate_vector_potential

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine get_spline_coefficients_3d(is, i_theta, i_phi, stp_R, stp_Z, stp_lam, &
                                          dstp_R_ds, dstp_Z_ds, dstp_lam_ds)
        use new_vmec_stuff_mod, only: sR, sZ, slam, ns_s, ns_tp

        integer, intent(in) :: is, i_theta, i_phi
        real(dp), dimension(:, :), intent(out) :: stp_R, stp_Z, stp_lam
        real(dp), dimension(:, :), intent(out) :: dstp_R_ds, dstp_Z_ds, dstp_lam_ds
        integer :: nstp

        nstp = ns_tp + 1

        stp_R(1:nstp, 1:nstp) = sR(ns_s + 1, :, :, is, i_theta, i_phi)
        dstp_R_ds(1:nstp, 1:nstp) = 0.d0
        stp_Z(1:nstp, 1:nstp) = sZ(ns_s + 1, :, :, is, i_theta, i_phi)
        dstp_Z_ds(1:nstp, 1:nstp) = 0.d0
        stp_lam(1:nstp, 1:nstp) = slam(ns_s + 1, :, :, is, i_theta, i_phi)
        dstp_lam_ds(1:nstp, 1:nstp) = 0.d0

    end subroutine get_spline_coefficients_3d

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine interpolate_s_direction(ds, is, i_theta, i_phi, stp_R, stp_Z, stp_lam, &
                                       dstp_R_ds, dstp_Z_ds, dstp_lam_ds)
        use new_vmec_stuff_mod, only: sR, sZ, slam, ns_s, ns_tp

        real(dp), intent(in) :: ds
        integer, intent(in) :: is, i_theta, i_phi
        real(dp), dimension(:, :), intent(inout) :: stp_R, stp_Z, stp_lam
        real(dp), dimension(:, :), intent(inout) :: dstp_R_ds, dstp_Z_ds, dstp_lam_ds
        integer :: k, nstp

        nstp = ns_tp + 1

        do k = ns_s, 1, -1
            stp_R(1:nstp, 1:nstp) = sR(k, :, :, is, i_theta, i_phi) + ds*stp_R(1:nstp, 1:nstp)
            dstp_R_ds(1:nstp, 1:nstp) = sR(k + 1, :, :, is, i_theta, i_phi)*dble(k) + ds*dstp_R_ds(1:nstp, 1:nstp)
            stp_Z(1:nstp, 1:nstp) = sZ(k, :, :, is, i_theta, i_phi) + ds*stp_Z(1:nstp, 1:nstp)
            dstp_Z_ds(1:nstp, 1:nstp) = sZ(k + 1, :, :, is, i_theta, i_phi)*dble(k) + ds*dstp_Z_ds(1:nstp, 1:nstp)
            stp_lam(1:nstp, 1:nstp) = slam(k, :, :, is, i_theta, i_phi) + ds*stp_lam(1:nstp, 1:nstp)
            dstp_lam_ds(1:nstp, 1:nstp) = slam(k + 1, :, :, is, i_theta, i_phi)*dble(k) + ds*dstp_lam_ds(1:nstp, 1:nstp)
        end do

    end subroutine interpolate_s_direction

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine interpolate_theta_direction(dtheta, stp_R, stp_Z, stp_lam, dstp_R_ds, dstp_Z_ds, dstp_lam_ds, &
                                           sp_R, sp_Z, sp_lam, dsp_R_ds, dsp_Z_ds, dsp_lam_ds, &
                                           dsp_R_dt, dsp_Z_dt, dsp_lam_dt)
        use new_vmec_stuff_mod, only: ns_tp

        real(dp), intent(in) :: dtheta
        real(dp), dimension(:, :), intent(in) :: stp_R, stp_Z, stp_lam
        real(dp), dimension(:, :), intent(in) :: dstp_R_ds, dstp_Z_ds, dstp_lam_ds
        real(dp), dimension(:), intent(out) :: sp_R, sp_Z, sp_lam
        real(dp), dimension(:), intent(out) :: dsp_R_ds, dsp_Z_ds, dsp_lam_ds
        real(dp), dimension(:), intent(out) :: dsp_R_dt, dsp_Z_dt, dsp_lam_dt
        integer :: k, nstp

        nstp = ns_tp + 1

        sp_R(1:nstp) = stp_R(nstp, 1:nstp)
        dsp_R_ds(1:nstp) = dstp_R_ds(nstp, 1:nstp)
        dsp_R_dt(1:nstp) = 0.d0
        sp_Z(1:nstp) = stp_Z(nstp, 1:nstp)
        dsp_Z_ds(1:nstp) = dstp_Z_ds(nstp, 1:nstp)
        dsp_Z_dt(1:nstp) = 0.d0
        sp_lam(1:nstp) = stp_lam(nstp, 1:nstp)
        dsp_lam_ds(1:nstp) = dstp_lam_ds(nstp, 1:nstp)
        dsp_lam_dt(1:nstp) = 0.d0

        do k = ns_tp, 1, -1
            sp_R(1:nstp) = stp_R(k, 1:nstp) + dtheta*sp_R(1:nstp)
            dsp_R_ds(1:nstp) = dstp_R_ds(k, 1:nstp) + dtheta*dsp_R_ds(1:nstp)
            dsp_R_dt(1:nstp) = stp_R(k + 1, 1:nstp)*dble(k) + dtheta*dsp_R_dt(1:nstp)

            sp_Z(1:nstp) = stp_Z(k, 1:nstp) + dtheta*sp_Z(1:nstp)
            dsp_Z_ds(1:nstp) = dstp_Z_ds(k, 1:nstp) + dtheta*dsp_Z_ds(1:nstp)
            dsp_Z_dt(1:nstp) = stp_Z(k + 1, 1:nstp)*dble(k) + dtheta*dsp_Z_dt(1:nstp)

            sp_lam(1:nstp) = stp_lam(k, 1:nstp) + dtheta*sp_lam(1:nstp)
            dsp_lam_ds(1:nstp) = dstp_lam_ds(k, 1:nstp) + dtheta*dsp_lam_ds(1:nstp)
            dsp_lam_dt(1:nstp) = stp_lam(k + 1, 1:nstp)*dble(k) + dtheta*dsp_lam_dt(1:nstp)
        end do

    end subroutine interpolate_theta_direction

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine interpolate_phi_direction(dphi, sp_R, sp_Z, sp_lam, dsp_R_ds, dsp_Z_ds, dsp_lam_ds, &
                                         dsp_R_dt, dsp_Z_dt, dsp_lam_dt, &
                                         R, Z, alam, dR_ds, dZ_ds, dl_ds, dR_dt, dZ_dt, dl_dt, &
                                         dR_dp, dZ_dp, dl_dp)
        use new_vmec_stuff_mod, only: ns_tp

        real(dp), intent(in) :: dphi
        real(dp), dimension(:), intent(in) :: sp_R, sp_Z, sp_lam
        real(dp), dimension(:), intent(in) :: dsp_R_ds, dsp_Z_ds, dsp_lam_ds
        real(dp), dimension(:), intent(in) :: dsp_R_dt, dsp_Z_dt, dsp_lam_dt
        real(dp), intent(out) :: R, Z, alam, dR_ds, dZ_ds, dl_ds
        real(dp), intent(out) :: dR_dt, dZ_dt, dl_dt, dR_dp, dZ_dp, dl_dp
        integer :: k, nstp

        nstp = ns_tp + 1

        R = sp_R(nstp)
        dR_ds = dsp_R_ds(nstp)
        dR_dt = dsp_R_dt(nstp)
        dR_dp = 0.d0
        Z = sp_Z(nstp)
        dZ_ds = dsp_Z_ds(nstp)
        dZ_dt = dsp_Z_dt(nstp)
        dZ_dp = 0.d0
        alam = sp_lam(nstp)
        dl_ds = dsp_lam_ds(nstp)
        dl_dt = dsp_lam_dt(nstp)
        dl_dp = 0.d0

        do k = ns_tp, 1, -1
            R = sp_R(k) + dphi*R
            dR_ds = dsp_R_ds(k) + dphi*dR_ds
            dR_dt = dsp_R_dt(k) + dphi*dR_dt
            dR_dp = sp_R(k + 1)*dble(k) + dphi*dR_dp

            Z = sp_Z(k) + dphi*Z
            dZ_ds = dsp_Z_ds(k) + dphi*dZ_ds
            dZ_dt = dsp_Z_dt(k) + dphi*dZ_dt
            dZ_dp = sp_Z(k + 1)*dble(k) + dphi*dZ_dp

            alam = sp_lam(k) + dphi*alam
            dl_ds = dsp_lam_ds(k) + dphi*dl_ds
            dl_dt = dsp_lam_dt(k) + dphi*dl_dt
            dl_dp = sp_lam(k + 1)*dble(k) + dphi*dl_dp
        end do

    end subroutine interpolate_phi_direction

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
        use vector_potentail_mod, only: ns, hs
        use new_vmec_stuff_mod, only: ns_tp

        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
        real(dp), intent(out) :: R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
        real(dp), intent(out) :: dl_ds, dl_dt, dl_dp

        integer, parameter :: ns_max = 6
        integer :: is, i_theta, i_phi, is_rho
        real(dp) :: ds, dtheta, dphi, ds_rho, rho_tor
        real(dp), dimension(ns_max) :: sp_R, sp_Z, sp_lam
        real(dp), dimension(ns_max) :: dsp_R_ds, dsp_Z_ds, dsp_lam_ds
        real(dp), dimension(ns_max) :: dsp_R_dt, dsp_Z_dt, dsp_lam_dt
        real(dp), dimension(ns_max, ns_max) :: stp_R, stp_Z, stp_lam
        real(dp), dimension(ns_max, ns_max) :: dstp_R_ds, dstp_Z_ds, dstp_lam_ds

        call normalize_coordinates(s, theta, varphi, ds, is, dtheta, i_theta, dphi, i_phi)
        call interpolate_vector_potential(s, ds, is, A_phi, A_theta, dA_phi_ds, dA_theta_ds)

        aiota = -dA_phi_ds/dA_theta_ds
        rho_tor = sqrt(s)
        ds_rho = rho_tor/hs
        is_rho = max(0, min(ns - 2, int(ds_rho)))
        ds_rho = (ds_rho - dble(is_rho))*hs
        is_rho = is_rho + 1

        call get_spline_coefficients_3d(is_rho, i_theta, i_phi, stp_R, stp_Z, stp_lam, &
                                        dstp_R_ds, dstp_Z_ds, dstp_lam_ds)

        call interpolate_s_direction(ds_rho, is_rho, i_theta, i_phi, stp_R, stp_Z, stp_lam, &
                                     dstp_R_ds, dstp_Z_ds, dstp_lam_ds)

        call interpolate_theta_direction(dtheta, stp_R, stp_Z, stp_lam, &
                                         dstp_R_ds, dstp_Z_ds, dstp_lam_ds, &
                                         sp_R, sp_Z, sp_lam, dsp_R_ds, dsp_Z_ds, dsp_lam_ds, &
                                         dsp_R_dt, dsp_Z_dt, dsp_lam_dt)

        call interpolate_phi_direction(dphi, sp_R, sp_Z, sp_lam, &
                                       dsp_R_ds, dsp_Z_ds, dsp_lam_ds, &
                                       dsp_R_dt, dsp_Z_dt, dsp_lam_dt, &
                                       R, Z, alam, dR_ds, dZ_ds, dl_ds, &
                                       dR_dt, dZ_dt, dl_dt, dR_dp, dZ_dp, dl_dp)

        dR_ds = 0.5d0*dR_ds/rho_tor
        dZ_ds = 0.5d0*dZ_ds/rho_tor
        dl_ds = 0.5d0*dl_ds/rho_tor

    end subroutine splint_vmec_data

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                          sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                          Bcovar_r, Bcovar_vartheta, Bcovar_varphi)

        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                                 sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                                 Bcovar_r, Bcovar_vartheta, Bcovar_varphi
        real(dp) :: R, Z, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp

        call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                              R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

        call compute_field_components(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                      dA_theta_ds, dA_phi_ds, dl_ds, dl_dt, dl_dp, &
                                      sqg, Bctrvr_vartheta, Bctrvr_varphi, Bcovar_r, Bcovar_vartheta, Bcovar_varphi)

    end subroutine vmec_field

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine compute_field_components(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                        dA_theta_ds, dA_phi_ds, dl_ds, dl_dt, dl_dp, &
                                        sqg, Bctrvr_vartheta, Bctrvr_varphi, Bcovar_r, Bcovar_vartheta, Bcovar_varphi)

        real(dp), intent(in) :: R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                dA_theta_ds, dA_phi_ds, dl_ds, dl_dt, dl_dp
        real(dp), intent(out) :: sqg, &
                                 Bctrvr_vartheta, Bctrvr_varphi, Bcovar_vartheta, Bcovar_varphi, Bcovar_r

        real(dp) :: g(3, 3)

        call compute_metric_tensor(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                   dl_ds, dl_dt, dl_dp, g, sqg)

        Bctrvr_vartheta = -dA_phi_ds/sqg
        Bctrvr_varphi = dA_theta_ds/sqg

        Bcovar_r = g(1, 2)*Bctrvr_vartheta + g(1, 3)*Bctrvr_varphi
        Bcovar_vartheta = g(2, 2)*Bctrvr_vartheta + g(2, 3)*Bctrvr_varphi
        Bcovar_varphi = g(3, 2)*Bctrvr_vartheta + g(3, 3)*Bctrvr_varphi
    end subroutine compute_field_components

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine compute_metric_tensor(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                     dl_ds, dl_dt, dl_dp, g, sqg)

        real(dp), intent(in) :: R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                dl_ds, dl_dt, dl_dp
        real(dp), intent(out) :: g(3, 3), sqg

        real(dp), dimension(3, 3) :: cmat, gV
        real(dp) :: cjac, sqgV

        gV(1, 1) = dR_ds**2 + dZ_ds**2
        gV(1, 2) = dR_ds*dR_dt + dZ_ds*dZ_dt
        gV(1, 3) = dR_ds*dR_dp + dZ_ds*dZ_dp
        gV(2, 1) = gV(1, 2)
        gV(2, 2) = dR_dt**2 + dZ_dt**2
        gV(2, 3) = dR_dt*dR_dp + dZ_dt*dZ_dp
        gV(3, 1) = gV(1, 3)
        gV(3, 2) = gV(2, 3)
        gV(3, 3) = R**2 + dR_dp**2 + dZ_dp**2
        sqgV = R*(dR_dt*dZ_ds - dR_ds*dZ_dt)

        cjac = 1.d0/(1.d0 + dl_dt)
        sqg = sqgV*cjac

        cmat(1, 2:3) = 0.d0
        cmat(3, 1:2) = 0.d0
        cmat(1, 1) = 1.d0
        cmat(3, 3) = 1.d0
        cmat(2, 1) = -dl_ds*cjac
        cmat(2, 2) = cjac
        cmat(2, 3) = -dl_dp*cjac

        g = matmul(transpose(cmat), matmul(gV, cmat))
    end subroutine compute_metric_tensor

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine splint_iota(s, aiota, daiota_ds)

        use vector_potentail_mod, only: ns, hs, torflux, sA_phi
        use new_vmec_stuff_mod, only: ns_A

        real(dp), parameter :: twopi = 2.d0*3.14159265358979d0

        integer :: is, i_theta, i_phi, k
        real(dp) :: ds, dtheta, dphi
        real(dp) :: s, dA_phi_ds, dA_theta_ds, d2A_phi_ds2, aiota, daiota_ds

        integer, parameter :: ns_max = 6

        integer :: nstp

        dA_theta_ds = torflux

        ds = s/hs
        is = max(0, min(ns - 1, int(ds)))
        ds = (ds - dble(is))*hs
        is = is + 1

        dA_phi_ds = 0.d0

        do k = ns_A, 1, -1
            dA_phi_ds = sA_phi(k + 1, is)*dble(k) + ds*dA_phi_ds
        end do

        d2A_phi_ds2 = 0.d0

        do k = ns_A, 2, -1
            d2A_phi_ds2 = sA_phi(k + 1, is)*dble(k)*dble(k - 1) + ds*d2A_phi_ds2
        end do

        aiota = -dA_phi_ds/dA_theta_ds
        daiota_ds = -d2A_phi_ds2/dA_theta_ds

    end subroutine splint_iota

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine splint_lambda(s, theta, varphi, alam, dl_dt)

        use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, slam, nper, ns_s, ns_tp
        use vector_potentail_mod, only: ns, hs

        real(dp), parameter :: twopi = 2.d0*3.14159265358979d0

        integer :: is, i_theta, i_phi, k
        real(dp) :: ds, dtheta, dphi
        real(dp) :: s, theta, varphi, alam, dl_dt

        integer, parameter :: ns_max = 6

        integer :: nstp

        real(dp), dimension(ns_max) :: sp_lam
        real(dp), dimension(ns_max) :: dsp_lam_dt
        real(dp), dimension(ns_max, ns_max) :: stp_lam

        nstp = ns_tp + 1

        ds = sqrt(s)/hs
        is = max(0, min(ns - 1, int(ds)))
        ds = (ds - dble(is))*hs
        is = is + 1

        dtheta = modulo(theta, twopi)/h_theta
        i_theta = max(0, min(n_theta - 1, int(dtheta)))
        dtheta = (dtheta - dble(i_theta))*h_theta
        i_theta = i_theta + 1

        dphi = modulo(varphi, twopi/dble(nper))/h_phi
        i_phi = max(0, min(n_phi - 1, int(dphi)))
        dphi = (dphi - dble(i_phi))*h_phi
        i_phi = i_phi + 1

! Begin interpolation over $s$

        stp_lam(1:nstp, 1:nstp) = slam(ns_s + 1, :, :, is, i_theta, i_phi)

        do k = ns_s, 1, -1
            stp_lam(1:nstp, 1:nstp) = slam(k, :, :, is, i_theta, i_phi) + ds*stp_lam(1:nstp, 1:nstp)
        end do

! End interpolation over $s$
!----------------------------

! Begin interpolation over $\theta$

        sp_lam(1:nstp) = stp_lam(nstp, 1:nstp)
        dsp_lam_dt(1:nstp) = 0.d0

        do k = ns_tp, 1, -1
            sp_lam(1:nstp) = stp_lam(k, 1:nstp) + dtheta*sp_lam(1:nstp)
            dsp_lam_dt(1:nstp) = stp_lam(k + 1, 1:nstp)*dble(k) + dtheta*dsp_lam_dt(1:nstp)
        end do

! End interpolation over $\theta$
!--------------------------------

! Begin interpolation over $\varphi$

        alam = sp_lam(nstp)
        dl_dt = dsp_lam_dt(nstp)

        do k = ns_tp, 1, -1
            alam = sp_lam(k) + dphi*alam
            dl_dt = dsp_lam_dt(k) + dphi*dl_dt
        end do

    end subroutine splint_lambda

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Go from s to rho grid, with special treatment of the axis.

! Interpolate values from s to rho grid. It is assumed that the
! innermost points of the input grid are not valid/noisy and thus need
! special treatment.
! This is done by extrapolating from values outside of this region to
! the axis.
! Extrapolation is done linear (more robust).
! An intermediate rescaling with rho can be used (might be useful to
! enforce behaviour near the axis). This will be in effect for
! extrapolating to the axis and for the interpolation to the new grid.

! input:
! ------
! m: integer, exponent, values <= 0 are ignored. Intermediate scaling of
!   values is done with rho**m.
! ns: integer, size of input array.
! nrho: integer, size of output array.
! nheal: integer,
! arr_in: real(dp) 1d array, with ns elements.

! output:
! -------
! arr_out: real(dp) 1d array, with nrho elements.

! sideeffects:
! ------------
! none
    subroutine s_to_rho_healaxis(m, ns, nrho, nheal, arr_in, arr_out)

        use new_vmec_stuff_mod, only: ns_s, old_axis_healing

        integer, intent(in) :: m, ns, nrho, nheal
        real(dp), dimension(ns), intent(in) :: arr_in
        real(dp), dimension(nrho), intent(out) :: arr_out

        integer :: irho, is, k, nhe
        real(dp) :: hs, hrho, s, ds, rho, a, b, c
        real(dp), dimension(:, :), allocatable :: splcoe

        hs = 1.d0/dble(ns - 1)
        hrho = 1.d0/dble(nrho - 1)

        nhe = max(1, nheal) + 1

        ! Rescale
        do is = nhe, ns
            if (m .gt. 0) then
                rho = sqrt(hs*dble(is - 1))
                arr_out(is) = arr_in(is)/rho**m
            else
                arr_out(is) = arr_in(is)
            end if
        end do

        if (old_axis_healing) then
            ! parabolic extrapolation:
            a = arr_out(nhe)
            b = 0.5d0*(4.d0*arr_out(nhe + 1) - 3.d0*arr_out(nhe) - arr_out(nhe + 2))
            c = 0.5d0*(arr_out(nhe) + arr_out(nhe + 2) - 2.d0*arr_out(nhe + 1))

            do is = 1, nhe - 1
                arr_out(is) = a + b*dble(is - nhe) + c*dble(is - nhe)**2
            end do

        else
            ! linear extrapolation ("less accurate" but more robust):
            a = arr_out(nhe)
            b = arr_out(nhe + 1) - arr_out(nhe)

            do is = 1, nhe - 1
                arr_out(is) = a + b*dble(is - nhe)
            end do

        end if

        allocate (splcoe(0:ns_s, ns))

        splcoe(0, :) = arr_out

        call spl_reg(ns_s, ns, hs, splcoe)

        do irho = 1, nrho
            rho = hrho*dble(irho - 1)
            s = rho**2

            ds = s/hs
            is = max(0, min(ns - 1, int(ds)))
            ds = (ds - dble(is))*hs
            is = is + 1

            arr_out(irho) = splcoe(ns_s, is)

            do k = ns_s - 1, 0, -1
                arr_out(irho) = splcoe(k, is) + ds*arr_out(irho)
            end do

            ! Undo rescaling
            if (m .gt. 0) arr_out(irho) = arr_out(irho)*rho**m
        end do

        deallocate (splcoe)

    end subroutine s_to_rho_healaxis

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine determine_nheal_for_axis(m, ns, arr_in, nheal)
        !> Determines the number of first radial points, nheal, where data is
        !> replaced by extrapolation.
        !>
        !> Takes cosine harmonic amplitude of arr_in and checks the difference
        !> between the data value at point i and the value obtained by 3-rd
        !> order Lagrange polynomial extrapolation from the points i+1,i+2,i+3
        !> and i+4. If relative value exceeds given tolerance (parameter set
        !> to 30%) this point is regarded as "bad" so that amplitude values
        !> for all points with indices smaller than i are replaced by linear
        !> extrapolation from points i and i+1. Note that harmonic function at
        !> 60 points per period is extrapolated with relative error about
        !> 1e-4, therefore 30% is quite a large tolerance which is exceeded if
        !> there is noise in the data. Even with this high tolerance, some
        !> harmonics have to be extrapolated almost from the very edge. Those,
        !> however, have amplitudes about 8 orders of magnitude smaller than
        !> main harmonics and, therefore, play no role.
        !>
        !> input:
        !> ------
        !> m:
        !> ns:
        !> arr_in: real(dp) array (ns entries), data from which to
        !>   determine the number of points to extrapolate at the axis.
        !>
        !> output:
        !> -------
        !> nheal: integer, number of points to extrapolate at the axis.

        ! Lagrange polynomial stencil size for checking the data by extraplation:
        integer, parameter :: nplag = 4
        ! tolerance for Lagrange polynomial extrapolation by one point (to check if data is noisy):
        real(dp), parameter :: tol = 3.d-1
        real(dp), parameter :: tiny = 1.d-200
        ! 3-rd order Lagrange polynomial extrapolation coefficients from points (1,2,3,4) to point 0:
        real(dp), parameter, dimension(nplag) :: weight = (/4.d0, -6.d0, 4.d0, -1.d0/)

        integer, intent(in) :: m, ns
        integer, intent(out) :: nheal
        real(dp), dimension(ns), intent(in) :: arr_in

        integer :: is, k, nhe, ncheck

        real(dp) :: hs, s, ds, rho, rho_nonzero, errmax

        real(dp), dimension(:), allocatable :: arr

        ! We check points which are away by more than 3 stencils from the edge:
        ncheck = ns - 3*nplag

        hs = 1.d0/dble(ns - 1)
        allocate (arr(ns))

        do is = 2, ns
            if (m > 0) then
                rho = sqrt(hs*dble(is - 1))
                rho_nonzero = max(rho**m, tiny)
                arr(is) = arr_in(is)/rho_nonzero
            else
                arr(is) = arr_in(is)
            end if
        end do

        nheal = 1
        do is = ncheck, 2, -1
            nheal = is
            errmax = maxval(abs(arr(is:is + nplag)))*tol
            if (abs(arr(is) - sum(arr(is + 1:is + nplag)*weight)) > errmax) then
                exit
            end if
        end do

        deallocate (arr)

    end subroutine determine_nheal_for_axis

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine volume_and_B00(volume, B00)

        use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, nper

        integer :: is, i_theta, i_phi, k
        real(dp) :: volume, B00
        real(dp) :: B3, B2, bmod2
        real(dp) :: s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                    R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp, &
                    sqg, Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi

        s = 0.9999999999d0
        volume = 0.d0

        do i_theta = 0, n_theta - 2
            theta = h_theta*dble(i_theta)
            do i_phi = 0, n_phi - 2
                varphi = h_phi*dble(i_phi)

                call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                      R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

                volume = volume + R**2*dZ_dt
            end do
        end do

        volume = 0.5d0*abs(volume)*h_theta*h_phi*dble(nper)

        s = 1d-8
        theta = 0.d0
        B2 = 0.d0
        B3 = 0.d0
        do i_phi = 0, n_phi - 2
            varphi = h_phi*dble(i_phi)

            call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                            sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                            Bcovar_r, Bcovar_vartheta, Bcovar_varphi)

            bmod2 = Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi
            B2 = B2 + bmod2/Bctrvr_varphi
            B3 = B3 + bmod2*sqrt(bmod2)/Bctrvr_varphi
        end do

        B00 = B3/B2

    end subroutine volume_and_B00
end module spline_vmec_sub
