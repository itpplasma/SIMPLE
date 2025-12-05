module field_coils

use, intrinsic :: iso_fortran_env, only : dp => real64
use interpolate, only : &
    BatchSplineData3D, &
    evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
    evaluate_batch_splines_3d_der2
use util, only : twopi
use neo_biotsavart, only : coils_t, load_coils_from_file, &
    compute_vector_potential, compute_magnetic_field
use field_base, only : magnetic_field_t
use simple_coordinates, only : transform_vmec_to_cart

implicit none

type, extends(magnetic_field_t) :: CoilsField
    type(coils_t) :: coils

    ! Batch spline for optimized field evaluation
    ! 7 components: Ar, Ath, Aphi, hr, hth, hphi, Bmod
    type(BatchSplineData3D) :: spl_coils_batch
    logical :: splines_initialized = .false.
contains
    procedure :: evaluate
    procedure :: evaluate_direct
    procedure :: init_splines
    procedure :: evaluate_coils_batch
    final :: coils_field_cleanup
end type CoilsField

! Module variable for evaluate_direct_wrapper callback
class(CoilsField), pointer :: current_coils_field => null()

contains

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    class(CoilsField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    call self%evaluate_coils_batch(x, Acov, hcov, Bmod)

    if (present(sqgBctr)) then
        error stop 'sqgBctr not implemented for CoilsField'
    end if
end subroutine evaluate


subroutine evaluate_direct(self, x, Acov, hcov, Bmod, sqgBctr)
    !> Direct Biot-Savart evaluation without spline interpolation.
    !> Used during spline construction.
    class(CoilsField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    real(dp) :: xcart(3), dxcart_dxvmec(3,3)
    real(dp), dimension(3) :: Acart, Bcart
    real(dp) :: s, ds_dr

    s = x(1)**2
    ds_dr = 2d0 * x(1)

    call transform_vmec_to_cart([s, x(2), x(3)], xcart, dxcart_dxvmec)

    Acart = compute_vector_potential(self%coils, xcart)
    Bcart = compute_magnetic_field(self%coils, xcart)

    Acov = matmul(Acart, dxcart_dxvmec)
    Acov(1) = Acov(1) * ds_dr

    Bmod = sqrt(Bcart(1)**2 + Bcart(2)**2 + Bcart(3)**2)

    hcov = matmul(Bcart, dxcart_dxvmec) / Bmod
    hcov(1) = hcov(1) * ds_dr

    if (present(sqgBctr)) then
        error stop 'sqgBctr not implemented for CoilsField'
    end if
end subroutine evaluate_direct


subroutine create_coils_field(coils_file, coils_field, should_spline)
    character(*), intent(in) :: coils_file
    class(CoilsField), allocatable, intent(out) :: coils_field
    logical, intent(in), optional :: should_spline

    real(dp), parameter :: M_TO_CM = 100.0_dp
    real(dp), parameter :: A_TO_STATA = 2997924536.8431_dp

    logical :: should_spline_ = .true.

    allocate(CoilsField :: coils_field)
    call load_coils_from_file(coils_file, coils_field%coils)

    coils_field%coils%x = coils_field%coils%x * M_TO_CM
    coils_field%coils%y = coils_field%coils%y * M_TO_CM
    coils_field%coils%z = coils_field%coils%z * M_TO_CM
    coils_field%coils%current = coils_field%coils%current * A_TO_STATA

    if (present(should_spline)) should_spline_ = should_spline
    if (should_spline_) call coils_field%init_splines
end subroutine create_coils_field


subroutine init_splines(self)
    use new_vmec_stuff_mod, only : nper
    use field_spline, only : build_field_splines

    class(CoilsField), intent(inout), target :: self

    integer, parameter :: n_r = 62, n_th = 63, n_phi = 64
    real(dp) :: xmin(3), xmax(3)

    xmin = [1d-12, 0d0, 0d0]
    xmax = [1d0, twopi, twopi / nper]

    ! Set module pointer for callback
    current_coils_field => self

    call build_field_splines(self, xmin, xmax, n_r, n_th, n_phi, &
                              self%spl_coils_batch, evaluate_direct_wrapper)

    current_coils_field => null()
    self%splines_initialized = .true.
end subroutine init_splines


subroutine evaluate_direct_wrapper(field, x, Acov, hcov, Bmod)
    !> Wrapper to pass evaluate_direct to build_field_splines.
    !> Uses module pointer since we need to call the type-bound procedure.
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: Acov(3), hcov(3), Bmod

    if (.not. associated(current_coils_field)) then
        error stop 'evaluate_direct_wrapper: current_coils_field not set'
    end if
    call current_coils_field%evaluate_direct(x, Acov, hcov, Bmod)
end subroutine evaluate_direct_wrapper


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
    use interpolate, only : destroy_batch_splines_3d

    type(CoilsField), intent(inout) :: self

    if (self%splines_initialized) then
        call destroy_batch_splines_3d(self%spl_coils_batch)
        self%splines_initialized = .false.
    end if
end subroutine coils_field_cleanup

end module field_coils
