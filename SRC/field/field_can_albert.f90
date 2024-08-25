module field_can_albert

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: SplineData3D, construct_splines_3d, &
    evaluate_splines_3d, evaluate_splines_3d_der2
use field_can_base, only: FieldCan
use field_can_meiss, only: xmin, xmax, n_r, n_phi, n_th, order, periodic, twopi, &
    get_grid_point, generate_regular_grid, spl_A2, spl_A3, init
use psi_transform, only: grid_r_to_psi

implicit none

! For splining psi
real(dp) :: psi_inner, psi_outer
real(dp), dimension(:,:,:), allocatable :: psi_of_x, Bmod
real(dp), dimension(:,:,:,:), allocatable :: hcan, Acan
real(dp), dimension(:), allocatable :: psi_grid
logical :: dpsi_dr_positive

! For splining field components over canonical coordinates
real(dp), dimension(:,:,:), allocatable :: r_of_xc, &
Aph_of_xc, hph_of_xc, hth_of_xc, Bmod_of_xc

type(SplineData3D) :: spl_r_of_xc, &
spl_Aphi_of_xc, spl_hph_of_xc, spl_hth_of_xc, spl_Bmod_of_xc

contains

subroutine evaluate(f, r, th_c, ph_c, mode_secders)

    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    error stop 'field_can_albert.evaluate: not implemented.'

end subroutine evaluate


subroutine init_splines_with_psi
    use psi_transform, only: grid_r_to_psi
    real(dp) :: x(3)
    integer :: i_r, i_phi, i_th

    allocate( &
        r_of_xc(n_r, n_phi, n_th), &
        Aph_of_xc(n_r, n_phi, n_th), &
        hph_of_xc(n_r, n_phi, n_th), &
        hth_of_xc(n_r, n_phi, n_th), &
        Bmod_of_xc(n_r, n_phi, n_th) &
    )

    call init_psi_grid()

    call grid_r_to_psi(n_r, n_phi, n_th, xmin(1), xmax(1), psi_inner, psi_outer, &
        psi_of_x, Acan(2,:,:,:), hcan(2,:,:,:), hcan(3,:,:,:), Bmod, &
        r_of_xc, Aph_of_xc, hph_of_xc, hth_of_xc, Bmod_of_xc)

    call construct_splines_3d([psi_inner, 0.0d0, xmax(3)], [psi_outer, twopi, xmax(3)],&
        r_of_xc, order, periodic, spl_r_of_xc)

    !$omp parallel private(i_r, i_phi, i_th, x)
    !$omp do
    do i_th = 1, n_th
        do i_phi = 1, n_phi
            do i_r = 1, n_r
                x = get_grid_point(i_r, i_phi, i_th)
                x(1) = r_of_xc(i_r, i_phi, i_th)
                call evaluate_splines_3d(spl_A2, x, Aph_of_xc(i_r, i_phi, i_th))
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    call construct_splines_3d([psi_inner, 0.0d0, xmin(3)], [psi_outer, twopi, xmax(3)],&
        Aph_of_xc, order, periodic, spl_Aphi_of_xc)
    call construct_splines_3d([psi_inner, 0.0d0, xmin(3)], [psi_outer, twopi, xmax(3)],&
        hph_of_xc, order, periodic, spl_hph_of_xc)
    call construct_splines_3d([psi_inner, 0.0d0, xmin(3)], [psi_outer, twopi, xmax(3)],&
        hth_of_xc, order, periodic, spl_hth_of_xc)
    call construct_splines_3d([psi_inner, 0.0d0, xmin(3)], [psi_outer, twopi, xmax(3)],&
        Bmod_of_xc, order, periodic, spl_Bmod_of_xc)
end subroutine init_splines_with_psi


subroutine init_psi_grid
    real(dp) :: x(3)
    integer :: i_r, i_phi, i_th

    allocate(psi_of_x(n_r, n_phi, n_th), psi_grid(n_r))

    !$omp parallel private(i_r, i_phi, i_th, x)
    !$omp do
    do i_th = 1, n_th
        do i_phi = 1, n_phi
            do i_r = 1, n_r
                x = get_grid_point(i_r, i_phi, i_th)
                call evaluate_splines_3d(spl_A3, x, psi_of_x(i_r, i_phi, i_th))
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    ! Here we use the "safe side" approach (new grid is fully within the old grid).
    ! For the risky approach (old grid within the new grid) exchange
    ! "minval" and "maxval".
    if(psi_of_x(n_r, n_phi/2, n_th/2) > psi_of_x(1, n_phi/2, n_th/2)) then
        dpsi_dR_positive = .true.
        psi_inner = maxval(psi_of_x(1,:,:))
        psi_outer = minval(psi_of_x(n_r,:,:))
    else
        dpsi_dR_positive = .false.
        psi_inner = maxval(psi_of_x(n_r,:,:))
        psi_outer = minval(psi_of_x(1,:,:))
    endif

    do i_r = 1, n_r
        psi_grid(i_r) = psi_inner + (psi_outer - psi_inner) * (i_r - 1) / (n_r - 1)
    end do
end subroutine init_psi_grid

end module field_can_albert
