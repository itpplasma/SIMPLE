program test_gvec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_gvec, only: gvec_field_t, create_gvec_field

    implicit none

    class(gvec_field_t), allocatable :: field
    real(dp) :: x(3)
    real(dp) :: Acov(3)
    real(dp) :: hcov(3)
    real(dp) :: Bmod
    real(dp) :: sqgBctr(3)

    call create_gvec_field('wout.gvec_export.nc', field)

    x = [0.35_dp, 1.1_dp, 0.2_dp]
    call field%evaluate(x, Acov, hcov, Bmod, sqgBctr)

    if (.not. all(abs(Acov) > 0.0_dp)) error stop 'test_gvec: Acov must be finite'
    if (.not. all(abs(hcov) > 0.0_dp)) error stop 'test_gvec: hcov must be finite'
    if (Bmod <= 0.0_dp) error stop 'test_gvec: Bmod must be positive'
    if (abs(sqgBctr(2)) <= 0.0_dp) error stop 'test_gvec: sqgBctr(2) must be nonzero'
    if (abs(sqgBctr(3)) <= 0.0_dp) error stop 'test_gvec: sqgBctr(3) must be nonzero'

    print *, 'test_gvec passed'
end program test_gvec
