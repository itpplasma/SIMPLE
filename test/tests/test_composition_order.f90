program test_composition_order
    ! Composition schemes alternate the Euler step with its adjoint; their order
    ! of convergence on a smooth orbit shows whether each stage really runs
    ! its own Euler variant. The state compared is (p_theta, theta, phi, p_phi):
    ! after an implicit-explicit Euler stage the stored radius still belongs to
    ! the previous p_theta and is only first-order accurate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_symplectic, only: symplectic_integrator_t, multistage_integrator_t, &
                   orbit_sympl_init, orbit_timestep_sympl, orbit_timestep_sympl_multi, &
                                orbit_sympl_init_verlet, orbit_sympl_init_order4, &
                                orbit_sympl_init_mclachlan4, orbit_sympl_init_blanes4, &
                                orbit_sympl_init_kahan6, orbit_sympl_init_kahan8, GAUSS4
    use orbit_symplectic_quasi, only: f_quasi => f, orbit_timestep_multi_quasi
    use field_can_mod, only: field_can_t, field_can_init, eval_field => evaluate, &
                             field_can_from_name

    implicit none

    real(dp), parameter :: t_end = 1950.0_dp  ! a quarter bounce period
    real(dp), parameter :: rtol = 1.0e-13_dp
    integer, parameter :: nscheme = 6
    character(len=10), parameter :: names(nscheme) = &
                    [character(len=10) :: 'verlet', 'order4', 'mclachlan4', 'blanes4', &
                                                           'kahan6', 'kahan8']
    integer, parameter :: expected(nscheme) = [2, 4, 4, 4, 6, 8]
    integer, parameter :: nsteps(nscheme) = [16, 8, 8, 8, 4, 4]
    real(dp), parameter :: order_tol = 0.5_dp

    type(field_can_t) :: f
    real(dp) :: z0(4), zref(4), err_coarse, err_fine, order
    integer :: k, iq, nfail
    logical :: quasi

    call field_can_from_name('test')
    call field_can_init(f, 1.0e-5_dp, 1.0_dp, 0.0_dp)
    z0 = [0.1_dp, 1.5_dp, 0.0_dp, 0.0_dp]
    call eval_field(f, z0(1), z0(2), z0(3), 0)
    z0(4) = f%vpar*f%hph + f%Aph/f%ro0

    ! Eighth-order Gauss-Legendre, far below the composed errors.
    zref = reference_state(256)

    nfail = 0
    do iq = 0, 1
        quasi = iq == 1
        do k = 1, nscheme
            err_coarse = maxval(abs(composed_state(k, nsteps(k), quasi) - zref))
            err_fine = maxval(abs(composed_state(k, 2*nsteps(k), quasi) - zref))
            order = log(err_coarse/err_fine)/log(2.0_dp)
            write (*, '(a10, a6, 2es11.3, a, f6.2)') names(k), &
                merge(' quasi', '      ', quasi), err_coarse, err_fine, &
                '  order ', order
            if (abs(order - expected(k)) > order_tol) nfail = nfail + 1
        end do
    end do
    if (nfail > 0) error stop 'composition scheme misses its order'

contains

    function reference_state(nt) result(z)
        integer, intent(in) :: nt
        real(dp) :: z(4)
        type(symplectic_integrator_t) :: si
        integer :: kt, ierr

        call orbit_sympl_init(si, f, z0, t_end/nt, 1, rtol, GAUSS4)
        do kt = 1, nt
            call orbit_timestep_sympl(si, f, ierr)
            if (ierr /= 0) error stop 'reference orbit failed'
        end do
        z = [f%pth, si%z(2:4)]
    end function reference_state

    function composed_state(scheme, nt, quasi) result(z)
        integer, intent(in) :: scheme, nt
        logical, intent(in) :: quasi
        real(dp) :: z(4)
        type(multistage_integrator_t) :: mi
        integer :: kt, ierr
        real(dp) :: dt

        dt = t_end/nt
        select case (scheme)
        case (1)
            call orbit_sympl_init_verlet(mi, f, z0, dt, 1, rtol)
        case (2)
            call orbit_sympl_init_order4(mi, f, z0, dt, 1, rtol)
        case (3)
            call orbit_sympl_init_mclachlan4(mi, f, z0, dt, 1, rtol)
        case (4)
            call orbit_sympl_init_blanes4(mi, f, z0, dt, 1, rtol)
        case (5)
            call orbit_sympl_init_kahan6(mi, f, z0, dt, 1, rtol)
        case (6)
            call orbit_sympl_init_kahan8(mi, f, z0, dt, 1, rtol)
        end select
        if (quasi) f_quasi = f
        do kt = 1, nt
            if (quasi) then
                call orbit_timestep_multi_quasi(mi, ierr)
            else
                call orbit_timestep_sympl_multi(mi, f, ierr)
            end if
            if (ierr /= 0) error stop 'composed orbit failed'
        end do
        if (quasi) then
            z = [f_quasi%pth, mi%stages(1)%z(2:4)]
        else
            z = [f%pth, mi%stages(1)%z(2:4)]
        end if
    end function composed_state

end program test_composition_order
