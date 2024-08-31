module field_can_albert

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: SplineData3D, construct_splines_3d, &
    evaluate_splines_3d, evaluate_splines_3d_der, evaluate_splines_3d_der2
use field_can_base, only: FieldCan, n_field_evaluations
use field_can_meiss, only: xmin, xmax, n_r, n_th, n_phi, order, periodic, twopi, &
    get_grid_point, generate_regular_grid, spl_Ath, spl_Aph, spl_hth, spl_hph, spl_Bmod, &
    init_albert => init_meiss, init_transformation, spline_transformation, &
    init_canonical_field_components
use psi_transform, only: grid_r_to_psi

implicit none

! For splining psi
real(dp) :: psi_inner, psi_outer
real(dp), dimension(:,:,:), allocatable :: psi_of_x
real(dp), dimension(:), allocatable :: psi_grid
logical :: dpsi_dr_positive

! For splining field components over canonical coordinates
real(dp), dimension(:,:,:), allocatable :: r_of_xc, &
Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc

type(SplineData3D) :: spl_r_of_xc, &
spl_Aphi_of_xc, spl_hth_of_xc, spl_hph_of_xc, spl_Bmod_of_xc

real(8) :: Ath_norm

contains

subroutine evaluate(f, r, th_c, ph_c, mode_secders)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    real(dp) :: x(3)

    n_field_evaluations = n_field_evaluations + 1

    x = [r, th_c, ph_c]

    f%Ath = Ath_norm*x(1)
    f%dAth = [Ath_norm, 0d0, 0d0]

    if (mode_secders > 0) then
        f%d2Ath = 0d0
        call evaluate_splines_3d_der2(spl_Aphi_of_xc, x, f%Aph, f%dAph, f%d2Aph)

        call evaluate_splines_3d_der2(spl_hth_of_xc, x, f%hth, f%dhth, f%d2hth)
        call evaluate_splines_3d_der2(spl_hph_of_xc, x, f%hph, f%dhph, f%d2hph)

        call evaluate_splines_3d_der2(spl_Bmod_of_xc, x, f%Bmod, f%dBmod, f%d2Bmod)

        return
    end if

    call evaluate_splines_3d_der(spl_Aphi_of_xc, x, f%Aph, f%dAph)

    call evaluate_splines_3d_der(spl_hth_of_xc, x, f%hth, f%dhth)
    call evaluate_splines_3d_der(spl_hph_of_xc, x, f%hph, f%dhph)

    call evaluate_splines_3d_der(spl_Bmod_of_xc, x, f%Bmod, f%dBmod)
end subroutine evaluate


subroutine get_albert_coordinates
    print *, 'field_can_meiss.init_transformation'
    call init_transformation

    print *, 'field_can_meiss.spline_transformation'
    call spline_transformation

    print *, 'field_can_meiss.init_canonical_field_components'
    call init_canonical_field_components

    print *, 'field_can_albert.init_splines_with_psi'
    call init_splines_with_psi
end subroutine get_albert_coordinates


subroutine init_splines_with_psi
    use psi_transform, only: grid_r_to_psi
    real(dp) :: x(3)
    integer :: i_r, i_th, i_phi

    allocate( &
        r_of_xc(n_r, n_th, n_phi), &
        Aph_of_xc(n_r, n_th, n_phi), &
        hth_of_xc(n_r, n_th, n_phi), &
        hph_of_xc(n_r, n_th, n_phi), &
        Bmod_of_xc(n_r, n_th, n_phi) &
    )

    call init_psi_grid()

    associate(Aphi => spl_Aph%coeff(order(1), 0, 0, :, :, :), &
        hth => spl_hth%coeff(order(1), 0, 0, :, :, :), &
        hph => spl_hph%coeff(order(1), 0, 0, :, :, :), &
        Bmod => spl_Bmod%coeff(order(1), 0, 0, :, :, :))

        call grid_r_to_psi(xmin(1), xmax(1), psi_inner, psi_outer, psi_of_x, &
            Aphi, hth, hph, Bmod, r_of_xc, Aph_of_xc, hth_of_xc, hph_of_xc, Bmod_of_xc)
    end associate

    call construct_splines_3d([psi_inner, 0.0d0, xmin(3)], [psi_outer, twopi, xmax(3)],&
        r_of_xc, order, periodic, spl_r_of_xc)

    !$omp parallel private(i_r, i_th, i_phi, x)
    !$omp do
    do i_phi = 1, n_phi
        do i_th = 1, n_th
            do i_r = 1, n_r
                x = get_grid_point(i_r, i_th, i_phi)
                x(1) = r_of_xc(i_r, i_th, i_phi)
                call evaluate_splines_3d(spl_Aph, x, Aph_of_xc(i_r, i_th, i_phi))
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    call construct_splines_3d([psi_inner, xmin(2), xmin(3)], &
        [psi_outer, xmax(2), xmax(3)], Aph_of_xc, order, periodic, spl_Aphi_of_xc)
    call construct_splines_3d([psi_inner, xmin(2), xmin(3)], &
        [psi_outer, xmax(2), xmax(3)], hth_of_xc, order, periodic, spl_hth_of_xc)
    call construct_splines_3d([psi_inner, xmin(2), xmin(3)], &
        [psi_outer, xmax(2), xmax(3)], hph_of_xc, order, periodic, spl_hph_of_xc)
    call construct_splines_3d([psi_inner, xmin(2), xmin(3)], &
        [psi_outer, xmax(2), xmax(3)], Bmod_of_xc, order, periodic, spl_Bmod_of_xc)
end subroutine init_splines_with_psi


subroutine init_psi_grid
    real(dp) :: x(3)
    integer :: i_r, i_th, i_phi

    allocate(psi_of_x(n_r, n_th, n_phi), psi_grid(n_r))

    psi_of_x(:, :, :) = spl_Ath%coeff(order(1), 0, 0, :, :, :)
    Ath_norm = maxval(psi_of_x)
    psi_of_x = psi_of_x / Ath_norm

    ! Here we use the "safe side" approach (new grid is fully within the old grid).
    ! For the risky approach (old grid within the new grid) exchange
    ! "minval" and "maxval".
    if(psi_of_x(n_r, n_th/2, n_phi/2) > psi_of_x(1, n_th/2, n_phi/2)) then
        dpsi_dr_positive = .true.
        psi_inner = maxval(psi_of_x(1,:,:))
        psi_outer = minval(psi_of_x(n_r,:,:))
    else
        dpsi_dr_positive = .false.
        psi_inner = maxval(psi_of_x(n_r,:,:))
        psi_outer = minval(psi_of_x(1,:,:))
    endif

    do i_r = 1, n_r
        psi_grid(i_r) = psi_inner + (psi_outer - psi_inner) * (i_r - 1) / (n_r - 1)
    end do
end subroutine init_psi_grid

end module field_can_albert
