module field_coils

use, intrinsic :: iso_fortran_env, only : dp => real64
use interpolate, only : &
    BatchSplineData3D, construct_batch_splines_3d, &
    evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
    evaluate_batch_splines_3d_der2
use util, only : twopi
use neo_biotsavart, only : coils_t, load_coils_from_file, &
    compute_vector_potential, compute_magnetic_field
use field_base, only : magnetic_field_t
use simple_coordinates, only : transform_cyl_to_cart
use reference_coordinates, only : ref_coords

implicit none

type, extends(magnetic_field_t) :: CoilsField
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

    real(dp) :: xcart(3), dxcart_dcyl(3,3)
    real(dp), dimension(3) :: Acart, Bcart
    real(dp) :: s, ds_dr
    real(dp) :: u_ref(3), xcyl(3)
    real(dp) :: e_cov(3,3)
    real(dp) :: Acov_u(3), Bcov_u(3)

    if (.not. allocated(ref_coords)) then
        print *, 'field_coils.evaluate_direct: ref_coords not initialized'
        error stop
    end if

    s = x(1)**2
    ds_dr = 2.0_dp*x(1)

    u_ref = [s, x(2), x(3)]

    call ref_coords%evaluate_point(u_ref, xcyl)
    call ref_coords%covariant_basis(u_ref, e_cov)

    call transform_cyl_to_cart(xcyl, xcart, dxcart_dcyl)

    Acart = compute_vector_potential(self%coils, xcart)
    Bcart = compute_magnetic_field(self%coils, xcart)

    Acov_u = matmul(transpose(e_cov), Acart)
    Acov(1) = Acov_u(1)*ds_dr
    Acov(2) = Acov_u(2)
    Acov(3) = Acov_u(3)

    Bmod = sqrt(Bcart(1)**2 + Bcart(2)**2 + Bcart(3)**2)

    Bcov_u = matmul(transpose(e_cov), Bcart)

    hcov(1) = Bcov_u(1)*ds_dr/Bmod
    hcov(2) = Bcov_u(2)/Bmod
    hcov(3) = Bcov_u(3)/Bmod

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
    use field_spline, only : build_field_splines

    class(CoilsField), intent(inout) :: self

    integer :: n_r, n_th, n_phi
    real(dp) :: xmin(3), xmax(3)

    n_r = 62
    n_th = 63
    n_phi = 64

    xmin = [1.0e-12_dp, 0.0_dp, 0.0_dp]
    xmax = [1.0_dp, twopi, twopi/nper]

    call build_field_splines(self, xmin, xmax, n_r, n_th, n_phi, self%spl_coils_batch)
    self%splines_initialized = .true.

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
