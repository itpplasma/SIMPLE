module field_coils
    !> Biot-Savart coil field evaluation in Cartesian coordinates.
    !> Single responsibility: evaluate A and B from coil geometry.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart, only: coils_t, load_coils_from_file, &
        compute_vector_potential, compute_magnetic_field
    use field_base, only: magnetic_field_t
    use cartesian_coordinates, only: cartesian_coordinate_system_t

    implicit none

    type, extends(magnetic_field_t) :: coils_field_t
        type(coils_t) :: coils
    contains
        procedure :: evaluate => coils_evaluate
    end type coils_field_t

contains

    subroutine coils_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        !> Evaluate Biot-Savart field at Cartesian position x.
        !> Returns A and h = B/|B| as Cartesian vectors (covariant = contravariant).
        class(coils_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: Bcart(3)

        Acov = compute_vector_potential(self%coils, x)
        Bcart = compute_magnetic_field(self%coils, x)

        Bmod = sqrt(Bcart(1)**2 + Bcart(2)**2 + Bcart(3)**2)
        hcov = Bcart / Bmod

        if (present(sqgBctr)) then
            sqgBctr = Bcart
        end if
    end subroutine coils_evaluate


    subroutine create_coils_field(coils_file, field)
        !> Load coils from file and create field with Cartesian coordinates.
        character(*), intent(in) :: coils_file
        type(coils_field_t), intent(out) :: field

        real(dp), parameter :: M_TO_CM = 100.0d0
        real(dp), parameter :: A_TO_STATA = 2997924536.8431d0

        call load_coils_from_file(coils_file, field%coils)

        field%coils%x = field%coils%x * M_TO_CM
        field%coils%y = field%coils%y * M_TO_CM
        field%coils%z = field%coils%z * M_TO_CM
        field%coils%current = field%coils%current * A_TO_STATA

        allocate(cartesian_coordinate_system_t :: field%coords)
    end subroutine create_coils_field

end module field_coils
