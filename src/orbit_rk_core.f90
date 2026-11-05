module orbit_rk_core
  ! Device-callable shared core for the implicit Gauss-collocation step used by
  ! the guiding-center (GC) integrator. No procedure pointers and no class()
  ! polymorphism in the hot path: the residual/Jacobian/Newton kernels are
  ! inlinable !$acc routine seq routines, the way orbit_symplectic_euler1
  ! already shares the symplectic-Euler residual/Jacobian/Newton between the
  ! CPU integrator and the OpenACC GPU kernel.
  !
    ! The residual and Jacobian algebra is generated (src/generated/
    ! orbit_systems_generated.f90). coeff_rk_gauss stays in
  ! orbit_symplectic_base (plain data) and is called inside the kernels. The
  ! Newton shell here uses the LU solve rk_solve, which the CPU symplectic
  ! integrators in orbit_symplectic also use (no LAPACK call per step).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: twopi
  use field_can_mod, only: field_can_t, get_derivatives2, eval_field => evaluate
  use orbit_symplectic_base, only: symplectic_integrator_t, coeff_rk_gauss, sympl_rmax
    use orbit_systems_generated, only: gauss1_residual, gauss2_residual, &
                   gauss3_residual, gauss4_residual, gauss1_jacobian, gauss2_jacobian, &
                                       gauss3_jacobian, gauss4_jacobian

  implicit none
  private

  ! Model code kept as a named constant for the residual/Jacobian/Newton
  ! signatures so the GC caller in orbit_symplectic reads explicitly. Only the
  ! GC kernel exists here; the CPP and symplectic full-orbit models were shelved.
  integer, parameter, public :: MODEL_GC = 0

  public :: rk_gauss_residual, rk_gauss_jacobian, rk_solve

contains

  ! GC Gauss residual: the field-canonical degenerate-Lagrangian residual.
  subroutine rk_gauss_residual(model, si, fs, s, x, fvec)
    !$acc routine seq
    integer, intent(in) :: model
    type(symplectic_integrator_t), intent(inout) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(inout) :: fs(:)
    real(dp), intent(in) :: x(4*s)
    real(dp), intent(out) :: fvec(4*s)

    call gauss_canfield_residual(si, fs, s, x, fvec)
  end subroutine rk_gauss_residual

  ! GC Gauss Jacobian, analytic from the 2nd derivatives of field_can_t.
  subroutine rk_gauss_jacobian(model, si, fs, s, jac)
    !$acc routine seq
    integer, intent(in) :: model
    type(symplectic_integrator_t), intent(in) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(in) :: fs(:)
    real(dp), intent(out) :: jac(4*s, 4*s)

    call gauss_canfield_jacobian(si, fs, s, jac)
  end subroutine rk_gauss_jacobian

    ! Residual: layout x(4*k-3:4*k) = (r,theta,phi,pphi) per stage. The
    ! algebra is generated from gen/generate_orbit_systems.f90.
  subroutine gauss_canfield_residual(si, fs, s, x, fvec)
    !$acc routine seq
    type(symplectic_integrator_t), intent(inout) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(inout) :: fs(:)
    real(dp), intent(in) :: x(4*s)
    real(dp), intent(out) :: fvec(4*s)

        real(dp) :: a(s, s), b(s), c(s)
        integer :: k

    call coeff_rk_gauss(s, a, b, c)

    do k = 1, s
      call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 2)
      call get_derivatives2(fs(k), x(4*k))
    end do

        select case (s)
        case (1)
            call gauss1_residual(si%dt, si%pthold, si%z, a, x, fs, fvec)
        case (2)
            call gauss2_residual(si%dt, si%pthold, si%z, a, x, fs, fvec)
        case (3)
            call gauss3_residual(si%dt, si%pthold, si%z, a, x, fs, fvec)
        case (4)
            call gauss4_residual(si%dt, si%pthold, si%z, a, x, fs, fvec)
        end select
  end subroutine gauss_canfield_residual

    ! Jacobian, analytic from the 2nd derivatives of field_can_t.
  subroutine gauss_canfield_jacobian(si, fs, s, jac)
    !$acc routine seq
    type(symplectic_integrator_t), intent(in) :: si
    integer, intent(in) :: s
    type(field_can_t), intent(in) :: fs(:)
    real(dp), intent(out) :: jac(4*s, 4*s)

        real(dp) :: a(s, s), b(s), c(s)

    call coeff_rk_gauss(s, a, b, c)

        select case (s)
        case (1)
            call gauss1_jacobian(si%dt, a, fs, jac)
        case (2)
            call gauss2_jacobian(si%dt, a, fs, jac)
        case (3)
            call gauss3_jacobian(si%dt, a, fs, jac)
        case (4)
            call gauss4_jacobian(si%dt, a, fs, jac)
        end select
  end subroutine gauss_canfield_jacobian

  ! Device-portable dense LU solve A x = rhs with partial pivoting, in place on
  ! rhs. n <= 4*S_MAX (n <= 16 for s <= 4). info = 0 on success, else the failing
  ! pivot column. Replaces dgesv on CPU and device: thread-safe, no allocation.
    pure subroutine rk_solve(n, A, rhs, info)
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n, n), rhs(n)
        integer, intent(out) :: info
        integer :: i, j, k, ipiv
        real(dp) :: amax, akj, tmp

        ! Column-oriented loops: the innermost index runs down a column, which is
        ! contiguous in Fortran storage.
        info = 0
        do k = 1, n
            ipiv = k
            amax = abs(A(k, k))
            do i = k + 1, n
                if (abs(A(i, k)) > amax) then
                    amax = abs(A(i, k))
                    ipiv = i
                end if
            end do
            ! Singular column: amax = max|A(i,k)| >= 0, so .not.(amax > 0) means the
            ! whole pivot column is zero. Phrased this way to avoid a real-equality
            ! comparison (-Wcompare-reals).
            if (.not. (amax > 0d0)) then
                info = k
                return
            end if
            if (ipiv /= k) then
                do j = k, n
                    tmp = A(k, j); A(k, j) = A(ipiv, j); A(ipiv, j) = tmp
                end do
                tmp = rhs(k); rhs(k) = rhs(ipiv); rhs(ipiv) = tmp
            end if
            tmp = 1d0/A(k, k)
            do i = k + 1, n
                A(i, k) = A(i, k)*tmp
            end do
            do j = k + 1, n
                akj = A(k, j)
                do i = k + 1, n
                    A(i, j) = A(i, j) - A(i, k)*akj
                end do
            end do
            do i = k + 1, n
                rhs(i) = rhs(i) - A(i, k)*rhs(k)
            end do
        end do

        do j = n, 1, -1
            rhs(j) = rhs(j)/A(j, j)
            do i = 1, j - 1
                rhs(i) = rhs(i) - A(i, j)*rhs(j)
            end do
        end do
    end subroutine rk_solve

end module orbit_rk_core
