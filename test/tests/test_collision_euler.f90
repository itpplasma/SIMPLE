!> collis_alp::stost_euler and the counter-based noise shared with ASCOT5.
!> The noise fixture is the one of itpplasma/ascot5
!> src/unit_tests/test_collision_noise.c: both codes must draw these numbers.
program test_collision_euler
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_alp, only: stost_euler, collision_noise_normal, coleff, efcolf, &
                          velrat, enrat, set_uniform_coeffs, compute_collision_coeffs, &
                          clog_ref, clog_bratio
    implicit none

    call check_shared_fixture()
    call setup_cold_plasma()
    call check_zero_noise_step()
    call check_drag_only_ignores_noise()
    call check_kernel_moments()
    call check_reflections()
    print *, 'PASS: euler collision step with shared noise'
contains

    subroutine check_shared_fixture()
        integer(8), parameter :: key(4, 4) = reshape([0_8, 1_8, 0_8, 3_8, &
                       20261110_8, 1_8, 0_8, 4_8, 20261110_8, 16384_8, 1000000_8, 3_8, &
                                                   7_8, 42_8, 123456789_8, 0_8], [4, 4])
        real(dp), parameter :: expected(4) = [1.0115071454929696_dp, &
                -0.70803017369489951_dp, 1.1794644411780411_dp, -0.65541438399747554_dp]
        integer :: i
        real(dp) :: x
        do i = 1, 4
            x = collision_noise_normal(key(1, i), key(2, i), key(3, i), key(4, i))
            if (abs(x - expected(i)) > 4.0e-15_dp) then
                print *, i, x, expected(i)
                error stop 'Counter noise differs from the ASCOT5 fixture'
            end if
        end do
    end subroutine check_shared_fixture

    subroutine setup_cold_plasma()
        ! 55 keV H test ion in a 20 eV H plasma with 2.5 % C6+, n_e = 2e19 m^-3.
        real(dp) :: ni1, ni2, v0
        ni1 = 2.0e13_dp/(1.0_dp + 6.0_dp*0.025_dp)
        ni2 = 0.025_dp*ni1
        call compute_collision_coeffs(1.0_dp, 12.0_dp, 1.0_dp, 6.0_dp, ni1, ni2, &
                                      20.0_dp, 20.0_dp, 20.0_dp, 55.0e3_dp, v0, &
                                      efcolf, velrat, enrat, &
                                      clog_out=clog_ref, bratio_out=clog_bratio)
        call set_uniform_coeffs(efcolf, velrat, enrat, clog_ref, clog_bratio)
    end subroutine setup_cold_plasma

    subroutine check_zero_noise_step()
        real(dp) :: z(5), p, lam, h, dpp, dhh, fpeff
        integer :: ierr
        p = 0.3_dp
        lam = 0.4_dp
        call coleff(p, dpp, dhh, fpeff)
        h = 1.0e-3_dp/dhh
        z = [0.5_dp, 0.2_dp, 0.3_dp, p, lam]
        call stost_euler(z, h, 1, [0.0_dp, 0.0_dp], ierr)
        if (ierr /= 0) error stop 'Zero-noise step failed'
        if (abs(z(4) - (p + fpeff*h)) > 1.0e-15_dp .or. &
            abs(z(5) - (lam - 2.0_dp*dhh*h*lam)) > 1.0e-15_dp) &
            error stop 'Zero-noise step is not the Euler drift'
    end subroutine check_zero_noise_step

    subroutine check_drag_only_ignores_noise()
        real(dp) :: z(5), y(5), h, dpp, dhh, fpeff
        integer :: ierr
        call coleff(0.3_dp, dpp, dhh, fpeff)
        h = 1.0e-3_dp/dhh
        z = [0.5_dp, 0.2_dp, 0.3_dp, 0.3_dp, 0.4_dp]
        y = z
        call stost_euler(z, h, 3, [1.3_dp, -0.7_dp], ierr)
        call stost_euler(y, h, 3, [0.0_dp, 0.0_dp], ierr)
        if (any(z /= y) .or. z(5) /= 0.4_dp) &
            error stop 'drag_only step depends on the noise or moves the pitch'
    end subroutine check_drag_only_ignores_noise

    subroutine check_kernel_moments()
        ! Over the counter stream of one marker, increments must have the
        ! drift and variance of the Euler-Maruyama kernel.
        integer, parameter :: n = 400000
        real(dp), parameter :: p = 0.3_dp, lam = 0.4_dp
        real(dp) :: z(5), w(2), h, dpp, dhh, fpeff, sp, spp, sl, sll
        real(dp) :: var_p, var_l
        integer :: k, ierr
        call coleff(p, dpp, dhh, fpeff)
        h = 1.0e-3_dp/dhh
        sp = 0.0_dp; spp = 0.0_dp; sl = 0.0_dp; sll = 0.0_dp
        do k = 1, n
            w = [collision_noise_normal(5_8, 9_8, int(k, 8), 3_8), &
                 collision_noise_normal(5_8, 9_8, int(k, 8), 4_8)]
            if (collision_noise_normal(5_8, 9_8, int(k, 8), 3_8) /= w(1)) &
                error stop 'Counter noise is not reproducible'
            z = [0.5_dp, 0.2_dp, 0.3_dp, p, lam]
            call stost_euler(z, h, 1, w, ierr)
            if (ierr /= 0) error stop 'Euler step failed'
            sp = sp + (z(4) - p); spp = spp + (z(4) - p)**2
            sl = sl + (z(5) - lam); sll = sll + (z(5) - lam)**2
        end do
        var_p = 2.0_dp*dpp*h
        var_l = 2.0_dp*dhh*h*(1.0_dp - lam**2)
        if (abs(sp/n - fpeff*h) > 5.0_dp*sqrt(var_p/n) .or. &
            abs(sl/n + 2.0_dp*dhh*h*lam) > 5.0_dp*sqrt(var_l/n) .or. &
            abs(spp/n - (sp/n)**2 - var_p) > 0.02_dp*var_p .or. &
            abs(sll/n - (sl/n)**2 - var_l) > 0.02_dp*var_l) &
            error stop 'Euler increments do not have the kernel moments'
    end subroutine check_kernel_moments

    subroutine check_reflections()
        ! A large pitch kick at lambda near 1 and a speed below ASCOT5's
        ! cutoff 0.1*sqrt(T_e/m) both reflect.
        real(dp) :: z(5), h, dpp, dhh, fpeff, cutoff
        integer :: ierr
        call coleff(0.3_dp, dpp, dhh, fpeff)
        h = 0.05_dp/dhh
        z = [0.5_dp, 0.2_dp, 0.3_dp, 0.3_dp, 0.99_dp]
        call stost_euler(z, h, 4, [0.0_dp, 3.0_dp], ierr)
        if (ierr /= 0 .or. abs(z(5)) > 1.0_dp) error stop 'Pitch not reflected at 1'
        cutoff = 0.1_dp/sqrt(2.0_dp*enrat(3))
        z = [0.5_dp, 0.2_dp, 0.3_dp, 0.5_dp*cutoff, 0.0_dp]
        call stost_euler(z, 1.0e-30_dp, 3, [0.0_dp, 0.0_dp], ierr)
        if (ierr /= 0 .or. abs(z(4) - 1.5_dp*cutoff) > 1.0e-12_dp*cutoff) &
            error stop 'Speed not reflected at the cutoff'
    end subroutine check_reflections
end program test_collision_euler
