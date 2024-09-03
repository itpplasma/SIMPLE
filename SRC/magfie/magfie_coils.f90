module simple_magfie_coils

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: SplineData3D, construct_splines_3d, &
    evaluate_splines_3d, evaluate_splines_3d_der, evaluate_splines_3d_der2
use util, only: twopi
use neo_biotsavart, only: coils_t, load_coils_from_file, &
    compute_vector_potential, compute_magnetic_field
use simple_magfie_base, only: MagneticField
use simple_coordinates, only: transform_vmec_to_cart

implicit none

type, extends(MagneticField) :: CoilsField
    type(coils_t) :: coils
    type(SplineData3D) :: spl_x, spl_y, spl_z
    type(SplineData3D) :: spl_Ar, spl_Ath, spl_Aphi, spl_hr, spl_hth, spl_hphi, spl_Bmod
contains
    procedure :: evaluate
    procedure :: evaluate_direct
    procedure :: init_splines
end type CoilsField

contains

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    class(CoilsField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    call evaluate_splines_3d(self%spl_Ar, x, Acov(1))
    call evaluate_splines_3d(self%spl_Ath, x, Acov(2))
    call evaluate_splines_3d(self%spl_Aphi, x, Acov(3))

    call evaluate_splines_3d(self%spl_hr, x, hcov(1))
    call evaluate_splines_3d(self%spl_hth, x, hcov(2))
    call evaluate_splines_3d(self%spl_hphi, x, hcov(3))

    call evaluate_splines_3d(self%spl_Bmod, x, Bmod)

    if (present(sqgBctr)) then
        error stop 'sqgBctr not implemented'
    end if
end subroutine evaluate


subroutine evaluate_direct(self, x, Acov, hcov, Bmod, sqgBctr)
    class(CoilsField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    real(dp) :: xcart(3), dxcart_dxvmec(3,3)
    real(dp), dimension(3) :: Acart, Bcart

    call transform_vmec_to_cart(x, xcart, dxcart_dxvmec)

    Acart = compute_vector_potential(self%coils, xcart)
    Bcart = compute_magnetic_field(self%coils, xcart)

    Acov = matmul(Acart, dxcart_dxvmec)

    Bmod = sqrt(Bcart(1)**2 + Bcart(2)**2 + Bcart(3)**2)

    hcov = matmul(Bcart, dxcart_dxvmec)/Bmod

    if (present(sqgBctr)) then
        error stop 'sqgBctr not implemented'
    end if
end subroutine evaluate_direct


function create_coils_field(coils_file, should_spline) result(field)
    class(CoilsField), allocatable :: field
    character(*), intent(in) :: coils_file
    logical, intent(in), optional :: should_spline

    real(dp), parameter :: M_TO_CM = 100.0d0
    real(dp), parameter :: A_TO_STATA = 2997924536.8431d0

    logical :: should_spline_ = .True.

    allocate(CoilsField :: field)
    call load_coils_from_file(coils_file, field%coils)

    field%coils%x = field%coils%x * M_TO_CM
    field%coils%y = field%coils%y * M_TO_CM
    field%coils%z = field%coils%z * M_TO_CM
    field%coils%current = field%coils%current * A_TO_STATA

    if(present(should_spline)) should_spline_ = should_spline
    if(should_spline_) call field%init_splines
end function create_coils_field


subroutine init_splines(self)
    use new_vmec_stuff_mod, only : nper

    class(CoilsField), intent(inout) :: self

    integer :: n_r=62, n_th=63, n_phi=64
    real(dp) :: xmin(3) = [1d-12, 0d0, 0d0]
    real(dp) :: xmax(3) = [1d0, twopi, twopi]

    real(dp) :: h_r, h_th, h_phi

    integer, parameter :: order(3) = [5, 5, 5]
    logical, parameter :: periodic(3) = [.False., .True., .True.]

    real(dp) :: x(3)
    real(dp), dimension(:,:,:), allocatable :: Ar, Ath, Aphi, hr, hth, hphi, Bmod
    real(dp) :: Acov(3), hcov(3), lam, dlam(3), chi, dchi(3)
    integer :: i_r, i_th, i_phi, i_ctr

    xmax(3) = twopi/nper

    h_r = (xmax(1)-xmin(1))/(n_r-1)
    h_th = (xmax(2)-xmin(2))/(n_th-1)
    h_phi = (xmax(3)-xmin(3))/(n_phi-1)

    allocate(Ar(n_r,n_th,n_phi), Ath(n_r,n_th,n_phi), Aphi(n_r,n_th,n_phi))
    allocate(hr(n_r,n_th,n_phi), hth(n_r,n_th,n_phi), hphi(n_r,n_th,n_phi))
    allocate(Bmod(n_r,n_th,n_phi))

    i_ctr = 0d0
    !$omp parallel private(i_r, i_th, i_phi)
    !$omp do
    do i_phi=1,n_phi
        !$omp atomic
        i_ctr = i_ctr + 1
        call print_progress(i_ctr)
        do i_th=1,n_th
            do i_r=1,n_r
                x = get_grid_point(i_r, i_th, i_phi)

                call self%evaluate_direct(x, Acov, hcov, Bmod(i_r, i_th, i_phi))

                Ar(i_r, i_th, i_phi) = Acov(1)
                Ath(i_r, i_th, i_phi) = Acov(2)
                Aphi(i_r, i_th, i_phi) = Acov(3)

                hr(i_r, i_th, i_phi) = hcov(1)
                hth(i_r, i_th, i_phi) = hcov(2)
                hphi(i_r, i_th, i_phi) = hcov(3)
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    call construct_splines_3d(xmin, xmax, Ar, order, periodic, self%spl_Ar)
    call construct_splines_3d(xmin, xmax, Ath, order, periodic, self%spl_Ath)
    call construct_splines_3d(xmin, xmax, Aphi, order, periodic, self%spl_Aphi)
    call construct_splines_3d(xmin, xmax, hr, order, periodic, self%spl_hr)
    call construct_splines_3d(xmin, xmax, hth, order, periodic, self%spl_hth)
    call construct_splines_3d(xmin, xmax, hphi, order, periodic, self%spl_hphi)
    call construct_splines_3d(xmin, xmax, Bmod, order, periodic, self%spl_Bmod)

    contains

    pure function get_grid_point(i_r_, i_th_, i_phi_)
    real(dp) :: get_grid_point(3)
    integer, intent(in) :: i_r_, i_th_, i_phi_

    get_grid_point = [ &
        xmin(1) + h_r*dble(i_r_-1), &
        xmin(2) + h_th*dble(i_th_-1), &
        xmin(3) + h_phi*dble(i_phi_-1) &
    ]
    end function get_grid_point

    subroutine print_progress(i_ctr_)
        integer, intent(in) :: i_ctr_
        !$omp critical
        write(*,'(A, I4, A, I4)',advance='no') 'Biot-Savart: ', i_ctr_, ' of ', n_phi
        if (i_ctr_ < n_phi) then
            write(*, '(A)', advance="no") char(13)
        else
            write(*, *)
        end if
        !$omp end critical
    end subroutine print_progress
end subroutine init_splines

end module simple_magfie_coils
