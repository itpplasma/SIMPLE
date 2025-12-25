module field_can_albert
    !> Albert canonical coordinates extending Meiss coordinates.
    !>
    !> This module builds on field_can_meiss to provide Albert-type canonical
    !> coordinates (psi, theta_c, phi_c) where psi is proportional to the
    !> poloidal magnetic flux A_theta, giving the simple form A = (0, psi, Aphi).
    !>
    !> The Albert coordinates differ from Meiss in the radial coordinate:
    !>   - Meiss: r = sqrt(s), fixed radial coordinate
    !>   - Albert: psi ~ A_theta, flux-based radial coordinate
    !>
    !> The transformation is built by:
    !>   1. Computing Meiss coordinates via field_can_meiss
    !>   2. Regridding onto a uniform psi grid via psi_transform
    !>
    !> Key subroutines:
    !>   - get_albert_coordinates: Build full transformation (calls Meiss first)
    !>   - evaluate_albert: Fast splined field evaluation in Albert coords
    !>   - integ_to_ref_albert: Convert integrator coords to reference (s,th,ph)
    !>
    !> The Albert form simplifies the symplectic integrator since dA_theta/dr = 0.

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: &
    BatchSplineData3D, construct_batch_splines_3d, &
    evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
    evaluate_batch_splines_3d_der2
use field_can_base, only: field_can_t, n_field_evaluations
use field_can_meiss, only: xmin, xmax, n_r, n_th, n_phi, order, periodic, twopi, &
    get_grid_point, &
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

! Batch spline for r_of_xc transformation (1 component: r)
type(BatchSplineData3D) :: spl_r_batch

! Batch spline for optimized field evaluation (4 components: Aphi, hth, hph, Bmod)
type(BatchSplineData3D) :: spl_albert_batch

real(dp) :: Ath_norm

contains

subroutine evaluate_albert(f, r, th_c, ph_c, mode_secders)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    real(dp) :: x(3)

    n_field_evaluations = n_field_evaluations + 1

    x = [r, th_c, ph_c]

    f%Ath = Ath_norm*x(1)
    f%dAth = [Ath_norm, 0d0, 0d0]

    if (mode_secders > 0) then
        f%d2Ath = 0d0
        call evaluate_albert_batch_der2(f, x)
        return
    end if

    call evaluate_albert_batch_der(f, x)
end subroutine evaluate_albert


subroutine integ_to_ref_albert(xinteg, xref)
    use field_can_meiss, only: integ_to_ref_meiss

    real(dp), intent(in) :: xinteg(3)
    real(dp), intent(out) :: xref(3)
    real(dp) :: xmeiss(3), y_batch(1), x_spl(3)

    ! Swap coordinates for spline: physics [psi, th, phi] -> spline [phi, th, psi]
    x_spl = [xinteg(3), xinteg(2), xinteg(1)]
    call evaluate_batch_splines_3d(spl_r_batch, x_spl, y_batch)
    xmeiss(1) = y_batch(1)  ! r component
    xmeiss(2:3) = xinteg(2:3)
    call integ_to_ref_meiss(xmeiss, xref)
end subroutine integ_to_ref_albert


subroutine ref_to_integ_albert(xref, xinteg)
    use field_can_meiss, only: ref_to_integ_meiss, spl_field_batch

    real(dp), intent(in) :: xref(3)
    real(dp), intent(out) :: xinteg(3)

    real(dp) :: Ath, xmeiss(3), x_spl(3), y_batch_local(5)

    call ref_to_integ_meiss(xref, xmeiss)
    ! Swap coordinates for Meiss spline: physics [r, th, phi] -> spline [phi, th, r]
    x_spl = [xmeiss(3), xmeiss(2), xmeiss(1)]
    call evaluate_batch_splines_3d(spl_field_batch, x_spl, y_batch_local)
    Ath = y_batch_local(1)  ! Extract Ath component
    xinteg(1) = Ath/Ath_norm
    xinteg(2:3) = xmeiss(2:3)
end subroutine ref_to_integ_albert


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
    use field_can_meiss, only: spl_field_batch
    real(dp), dimension(:,:,:,:), allocatable :: y_batch
    real(dp), dimension(:,:,:), allocatable :: Aphi_grid, hth_grid, hph_grid, Bmod_grid
    real(dp) :: x_grid(3), x_spl(3), y_batch_temp(5)
    real(dp) :: xmin_spl(3), xmax_spl(3)
    logical :: periodic_spl(3)
    integer :: i_r, i_th, i_phi

    allocate( &
        r_of_xc(n_r, n_th, n_phi), &
        Aph_of_xc(n_r, n_th, n_phi), &
        hth_of_xc(n_r, n_th, n_phi), &
        hph_of_xc(n_r, n_th, n_phi), &
        Bmod_of_xc(n_r, n_th, n_phi) &
    )

    call init_psi_grid

    allocate(Aphi_grid(n_r, n_th, n_phi), hth_grid(n_r, n_th, n_phi), &
             hph_grid(n_r, n_th, n_phi), Bmod_grid(n_r, n_th, n_phi))

    ! Evaluate Meiss batch spline on grid to get field components
    ! spl_field_batch uses swapped coordinates: [phi, th, r]
    do i_phi = 1, n_phi
        do i_th = 1, n_th
            do i_r = 1, n_r
                x_grid = [xmin(1) + (i_r - 1)*(xmax(1) - xmin(1))/(n_r - 1), &
                          xmin(2) + (i_th - 1)*(xmax(2) - xmin(2))/(n_th - 1), &
                          xmin(3) + (i_phi - 1)*(xmax(3) - xmin(3))/(n_phi - 1)]
                ! Swap for Meiss spline: [r, th, phi] -> [phi, th, r]
                x_spl = [x_grid(3), x_grid(2), x_grid(1)]
                call evaluate_batch_splines_3d(spl_field_batch, x_spl, y_batch_temp)
                Aphi_grid(i_r, i_th, i_phi) = y_batch_temp(2)  ! Aph component
                hth_grid(i_r, i_th, i_phi) = y_batch_temp(3)   ! hth component
                hph_grid(i_r, i_th, i_phi) = y_batch_temp(4)   ! hph component
                Bmod_grid(i_r, i_th, i_phi) = y_batch_temp(5)  ! Bmod component
            end do
        end do
    end do

    ! Center Aphi around zero
    Aphi_grid = Aphi_grid - 0.5d0*sum(Aphi_grid)/real(n_r*n_th*n_phi, dp)

    call grid_r_to_psi(xmin(1), xmax(1), psi_inner, psi_outer, psi_of_x, &
        Aphi_grid, hth_grid, hph_grid, Bmod_grid, r_of_xc, Aph_of_xc, &
        hth_of_xc, hph_of_xc, Bmod_of_xc)

    ! Swap coordinate bounds and periodic flags: [psi, th, phi] -> [phi, th, psi]
    xmin_spl = [xmin(3), xmin(2), psi_inner]
    xmax_spl = [xmax(3), xmax(2), psi_outer]
    periodic_spl = [periodic(3), periodic(2), periodic(1)]

    ! Construct batch spline for r_of_xc (1 component: r)
    ! Transpose: (n_r, n_th, n_phi) -> (n_phi, n_th, n_r) for [phi, th, psi] layout
    block
        real(dp), dimension(:, :, :, :), allocatable :: y_r_batch
        allocate(y_r_batch(n_phi, n_th, n_r, 1))
        do i_phi = 1, n_phi
            do i_th = 1, n_th
                do i_r = 1, n_r
                    y_r_batch(i_phi, i_th, i_r, 1) = r_of_xc(i_r, i_th, i_phi)
                end do
            end do
        end do
        call construct_batch_splines_3d(xmin_spl, xmax_spl, y_r_batch, &
            [order(3), order(2), order(1)], periodic_spl, spl_r_batch)
    end block

    ! Construct batch spline for 4 Albert field components: [Aphi, hth, hph, Bmod]
    ! Transpose: (n_r, n_th, n_phi) -> (n_phi, n_th, n_r)
    allocate(y_batch(n_phi, n_th, n_r, 4))

    do i_phi = 1, n_phi
        do i_th = 1, n_th
            do i_r = 1, n_r
                y_batch(i_phi, i_th, i_r, 1) = Aph_of_xc(i_r, i_th, i_phi)
                y_batch(i_phi, i_th, i_r, 2) = hth_of_xc(i_r, i_th, i_phi)
                y_batch(i_phi, i_th, i_r, 3) = hph_of_xc(i_r, i_th, i_phi)
                y_batch(i_phi, i_th, i_r, 4) = Bmod_of_xc(i_r, i_th, i_phi)
            end do
        end do
    end do

    call construct_batch_splines_3d(xmin_spl, xmax_spl, y_batch, &
        [order(3), order(2), order(1)], periodic_spl, spl_albert_batch)
end subroutine init_splines_with_psi


subroutine init_psi_grid
    use field_can_meiss, only: spl_field_batch
    real(dp) :: x(3), x_spl(3), y_batch_local(5)
    integer :: i_r, i_th, i_phi

    allocate(psi_of_x(n_r, n_th, n_phi), psi_grid(n_r))

    ! Evaluate Meiss batch spline to get Ath (component 1) on grid
    ! spl_field_batch uses swapped coordinates: [phi, th, r]
    do i_phi = 1, n_phi
        do i_th = 1, n_th
            do i_r = 1, n_r
                x = [xmin(1) + (i_r - 1)*(xmax(1) - xmin(1))/(n_r - 1), &
                     xmin(2) + (i_th - 1)*(xmax(2) - xmin(2))/(n_th - 1), &
                     xmin(3) + (i_phi - 1)*(xmax(3) - xmin(3))/(n_phi - 1)]
                ! Swap for Meiss spline: [r, th, phi] -> [phi, th, r]
                x_spl = [x(3), x(2), x(1)]
                call evaluate_batch_splines_3d(spl_field_batch, x_spl, y_batch_local)
                psi_of_x(i_r, i_th, i_phi) = y_batch_local(1)  ! Ath component
            end do
        end do
    end do

    Ath_norm = sign(maxval(abs(psi_of_x)), psi_of_x(n_r, n_th/2, n_phi/2))
    psi_of_x = psi_of_x/Ath_norm

    ! Here we use the safe side approach (new grid is fully within the old grid).
    if (psi_of_x(n_r, n_th/2, n_phi/2) > psi_of_x(1, n_th/2, n_phi/2)) then
        dpsi_dr_positive = .true.
        psi_inner = maxval(psi_of_x(1, :, :))
        psi_outer = minval(psi_of_x(n_r, :, :))
    else
        dpsi_dr_positive = .false.
        psi_inner = maxval(psi_of_x(n_r, :, :))
        psi_outer = minval(psi_of_x(1, :, :))
    end if

    do i_r = 1, n_r
        psi_grid(i_r) = psi_inner + (psi_outer - psi_inner)*(i_r - 1)/(n_r - 1)
    end do
end subroutine init_psi_grid


subroutine magfie_albert(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!  Computes magnetic field and derivatives with bmod in units of the magnetic code
!
!  Input parameters:
!            formal:  x               - array of canonicalized coordinates r, th, ph
!  Output parameters:
!            formal:  bmod            - magnetic field module
!                     sqrtg           - metric determinant
!                     bder            - covariant components of (grad B)/B
!                     hcovar          - covariant components of \bB/B
!                     hctrvr          - contravariant components of \bB/B
!                     hcurl           - contravariant components of curl (\bB/B)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

    type(field_can_t) :: f
    real(dp) :: sqrtg_bmod

    call evaluate_albert(f, x(1), x(2), x(3), 0)

    bmod = f%Bmod

    sqrtg_bmod = f%hph*Ath_norm - f%hth*f%dAph(1)
    sqrtg = sqrtg_bmod/bmod
    bder = f%dBmod/bmod

    hcovar(1) = 0.d0
    hcovar(2) = f%hth
    hcovar(3) = f%hph

    hctrvr(1) = f%dAph(2)/sqrtg_bmod
    hctrvr(2) = -f%dAph(1)/sqrtg_bmod
    hctrvr(3) = Ath_norm/sqrtg_bmod

    hcurl(1) = (f%dhph(2) - f%dhth(3))/sqrtg
    hcurl(2) = -f%dhph(1)/sqrtg
    hcurl(3) = f%dhth(1)/sqrtg
end subroutine magfie_albert


! Batch evaluation helper for first derivatives
subroutine evaluate_albert_batch_der(f, x)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: x(3)

    real(dp) :: x_spl(3)
    real(dp) :: y_batch(4), dy_batch(3, 4)

    ! Swap coordinates for spline: physics [psi, th, phi] -> spline [phi, th, psi]
    x_spl = [x(3), x(2), x(1)]
    call evaluate_batch_splines_3d_der(spl_albert_batch, x_spl, y_batch, dy_batch)

    ! Unpack results: order is [Aphi, hth, hph, Bmod]
    f%Aph = y_batch(1)
    f%hth = y_batch(2)
    f%hph = y_batch(3)
    f%Bmod = y_batch(4)

    ! Unswap derivatives: spline [phi, th, psi] -> physics [psi, th, phi]
    ! dy_batch(1,:) = d/dphi, dy_batch(2,:) = d/dth, dy_batch(3,:) = d/dpsi
    f%dAph = [dy_batch(3, 1), dy_batch(2, 1), dy_batch(1, 1)]
    f%dhth = [dy_batch(3, 2), dy_batch(2, 2), dy_batch(1, 2)]
    f%dhph = [dy_batch(3, 3), dy_batch(2, 3), dy_batch(1, 3)]
    f%dBmod = [dy_batch(3, 4), dy_batch(2, 4), dy_batch(1, 4)]
end subroutine evaluate_albert_batch_der


! Batch evaluation helper for second derivatives
subroutine evaluate_albert_batch_der2(f, x)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: x(3)

    real(dp) :: x_spl(3)
    real(dp) :: y_batch(4), dy_batch(3, 4), d2y_batch(6, 4)

    ! Swap coordinates for spline: physics [psi, th, phi] -> spline [phi, th, psi]
    x_spl = [x(3), x(2), x(1)]
    call evaluate_batch_splines_3d_der2(spl_albert_batch, x_spl, y_batch, dy_batch, &
                                        d2y_batch)

    ! Unpack results: order is [Aphi, hth, hph, Bmod]
    f%Aph = y_batch(1)
    f%hth = y_batch(2)
    f%hph = y_batch(3)
    f%Bmod = y_batch(4)

    ! Unswap first derivatives: spline [phi, th, psi] -> physics [psi, th, phi]
    f%dAph = [dy_batch(3, 1), dy_batch(2, 1), dy_batch(1, 1)]
    f%dhth = [dy_batch(3, 2), dy_batch(2, 2), dy_batch(1, 2)]
    f%dhph = [dy_batch(3, 3), dy_batch(2, 3), dy_batch(1, 3)]
    f%dBmod = [dy_batch(3, 4), dy_batch(2, 4), dy_batch(1, 4)]

    ! Unswap second derivatives: spline [phi,th,psi] order to physics [psi,th,phi]
    ! Spline: (1)=phi-phi, (2)=phi-th, (3)=phi-psi, (4)=th-th, (5)=th-psi, (6)=psi-psi
    ! Physics: (1)=psi-psi, (2)=psi-th, (3)=psi-phi, (4)=th-th, (5)=th-phi, (6)=phi-phi
    f%d2Aph = [d2y_batch(6, 1), d2y_batch(5, 1), d2y_batch(3, 1), &
               d2y_batch(4, 1), d2y_batch(2, 1), d2y_batch(1, 1)]
    f%d2hth = [d2y_batch(6, 2), d2y_batch(5, 2), d2y_batch(3, 2), &
               d2y_batch(4, 2), d2y_batch(2, 2), d2y_batch(1, 2)]
    f%d2hph = [d2y_batch(6, 3), d2y_batch(5, 3), d2y_batch(3, 3), &
               d2y_batch(4, 3), d2y_batch(2, 3), d2y_batch(1, 3)]
    f%d2Bmod = [d2y_batch(6, 4), d2y_batch(5, 4), d2y_batch(3, 4), &
                d2y_batch(4, 4), d2y_batch(2, 4), d2y_batch(1, 4)]
end subroutine evaluate_albert_batch_der2

end module field_can_albert
