module field_coils

use, intrinsic :: iso_fortran_env, only: dp => real64
use interpolate, only: &
    BatchSplineData3D, construct_batch_splines_3d, &
    evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
    evaluate_batch_splines_3d_der2
use util, only: twopi
use neo_biotsavart, only: coils_t, load_coils_from_file, &
    compute_vector_potential, compute_magnetic_field
use field_base, only: MagneticField
use simple_coordinates, only: transform_vmec_to_cart

implicit none

type, extends(MagneticField) :: CoilsField
    type(coils_t) :: coils
    
    ! Batch spline for optimized field evaluation (7 components: Ar, Ath, Aphi, hr, hth, hphi, Bmod)
    type(BatchSplineData3D) :: spl_coils_batch
    logical :: splines_initialized = .false.
contains
    procedure :: evaluate
    procedure :: evaluate_direct
    procedure :: init_splines
    procedure :: evaluate_coils_batch
    final :: coils_field_cleanup
end type CoilsField

contains

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    class(CoilsField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    call self%evaluate_coils_batch(x, Acov, hcov, Bmod)

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

    real(dp) :: s, ds_dr

    s = x(1)**2
    ds_dr = 2d0*x(1)

    call transform_vmec_to_cart([s, x(2), x(3)], xcart, dxcart_dxvmec)

    Acart = compute_vector_potential(self%coils, xcart)
    Bcart = compute_magnetic_field(self%coils, xcart)

    Acov = matmul(Acart, dxcart_dxvmec)
    Acov(1) = Acov(1)*ds_dr

    Bmod = sqrt(Bcart(1)**2 + Bcart(2)**2 + Bcart(3)**2)

    hcov = matmul(Bcart, dxcart_dxvmec)/Bmod
    hcov(1) = hcov(1)*ds_dr

    if (present(sqgBctr)) then
        error stop 'sqgBctr not implemented'
    end if
end subroutine evaluate_direct


subroutine create_coils_field(coils_file, coils_field, should_spline)
    character(*), intent(in) :: coils_file
    class(CoilsField), allocatable, intent(out) :: coils_field
    logical, intent(in), optional :: should_spline

    real(dp), parameter :: M_TO_CM = 100.0d0
    real(dp), parameter :: A_TO_STATA = 2997924536.8431d0

    logical :: should_spline_ = .True.

    allocate(CoilsField :: coils_field)
    call load_coils_from_file(coils_file, coils_field%coils)

    coils_field%coils%x = coils_field%coils%x * M_TO_CM
    coils_field%coils%y = coils_field%coils%y * M_TO_CM
    coils_field%coils%z = coils_field%coils%z * M_TO_CM
    coils_field%coils%current = coils_field%coils%current * A_TO_STATA

    if(present(should_spline)) should_spline_ = should_spline
    if(should_spline_) call coils_field%init_splines
end subroutine create_coils_field


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
    real(dp), dimension(:,:,:,:), allocatable :: y_batch
    real(dp) :: Acov(3), hcov(3)
    integer :: i_r, i_th, i_phi, i_ctr, dims(3)

    xmax(3) = twopi/nper

    h_r = (xmax(1)-xmin(1))/(n_r-1)
    h_th = (xmax(2)-xmin(2))/(n_th-1)
    h_phi = (xmax(3)-xmin(3))/(n_phi-1)

    allocate(Ar(n_r,n_th,n_phi), Ath(n_r,n_th,n_phi), Aphi(n_r,n_th,n_phi))
    allocate(hr(n_r,n_th,n_phi), hth(n_r,n_th,n_phi), hphi(n_r,n_th,n_phi))
    allocate(Bmod(n_r,n_th,n_phi))

    i_ctr = 0d0
    !$omp parallel private(i_r, i_th, i_phi, x, Acov, hcov)
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

    Ar = Ar - Ar(1,1,1)
    Ath = Ath - Ath(1,1,1)
    Aphi = Aphi - Aphi(1,1,1)

    ! Construct batch spline for all 7 field components: [Ar, Ath, Aphi, hr, hth, hphi, Bmod]
    dims = shape(Ar)
    allocate(y_batch(dims(1), dims(2), dims(3), 7))
    
    y_batch(:,:,:,1) = Ar
    y_batch(:,:,:,2) = Ath
    y_batch(:,:,:,3) = Aphi
    y_batch(:,:,:,4) = hr
    y_batch(:,:,:,5) = hth
    y_batch(:,:,:,6) = hphi
    y_batch(:,:,:,7) = Bmod
    
    call construct_batch_splines_3d(xmin, xmax, y_batch, order, periodic, self%spl_coils_batch)
    self%splines_initialized = .true.

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

    subroutine print_progress(i)
        integer, intent(in) :: i
        !$omp critical
        if (mod(i, n_phi / 10) == 0) then
            write(*, '(a,f6.2,a)', advance='no') char(13), 100.0*i/n_phi, ' %'
            call flush(6)
        else
            write(*, *)
        end if
        !$omp end critical
    end subroutine print_progress
end subroutine init_splines


subroutine evaluate_coils_batch(self, x, Acov, hcov, Bmod)
    class(CoilsField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    
    real(dp) :: y_batch(7)
    
    call evaluate_batch_splines_3d(self%spl_coils_batch, x, y_batch)
    
    ! Unpack results: order is [Ar, Ath, Aphi, hr, hth, hphi, Bmod]
    Acov(1) = y_batch(1)
    Acov(2) = y_batch(2)
    Acov(3) = y_batch(3)
    
    hcov(1) = y_batch(4)
    hcov(2) = y_batch(5)
    hcov(3) = y_batch(6)
    
    Bmod = y_batch(7)
end subroutine evaluate_coils_batch


subroutine coils_field_cleanup(self)
    ! Clean up batch splines when CoilsField is destroyed
    use interpolate, only: destroy_batch_splines_3d
    
    type(CoilsField), intent(inout) :: self
    
    if (self%splines_initialized) then
        call destroy_batch_splines_3d(self%spl_coils_batch)
        self%splines_initialized = .false.
    end if
end subroutine coils_field_cleanup

end module field_coils