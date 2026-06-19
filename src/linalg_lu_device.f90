module linalg_lu_device
  ! Leaf device-portable dense LU solve, no field/spline dependencies. Shared by
  ! the GC/CPP/full-orbit Newton shells (orbit_rk_core) and the 6D canonical
  ! integrator (orbit_cpp_canonical). Keeping it a leaf module lets the COORD_TOK
  ! device kernel chain link without pulling the field-canonical / boozer stack.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: rk_solve

contains

  ! Device-portable dense LU solve A x = rhs with partial pivoting, in place on
  ! rhs. info = 0 on success, else the failing pivot column. Replaces dgesv on the
  ! device; the GC CPU path keeps dgesv.
  pure subroutine rk_solve(n, A, rhs, info)
    !$acc routine seq
    integer, intent(in) :: n
    real(dp), intent(inout) :: A(n,n), rhs(n)
    integer, intent(out) :: info
    integer :: i, j, k, ipiv
    real(dp) :: piv, amax, factor, tmp

    info = 0
    do k = 1, n
      ipiv = k
      amax = abs(A(k,k))
      do i = k+1, n
        if (abs(A(i,k)) > amax) then
          amax = abs(A(i,k))
          ipiv = i
        end if
      end do
      if (amax == 0d0) then
        info = k
        return
      end if
      if (ipiv /= k) then
        do j = 1, n
          tmp = A(k,j); A(k,j) = A(ipiv,j); A(ipiv,j) = tmp
        end do
        tmp = rhs(k); rhs(k) = rhs(ipiv); rhs(ipiv) = tmp
      end if
      piv = A(k,k)
      do i = k+1, n
        factor = A(i,k)/piv
        A(i,k) = factor
        do j = k+1, n
          A(i,j) = A(i,j) - factor*A(k,j)
        end do
        rhs(i) = rhs(i) - factor*rhs(k)
      end do
    end do

    do i = n, 1, -1
      tmp = rhs(i)
      do j = i+1, n
        tmp = tmp - A(i,j)*rhs(j)
      end do
      rhs(i) = tmp/A(i,i)
    end do
  end subroutine rk_solve

end module linalg_lu_device
