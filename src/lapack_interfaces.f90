module lapack_interfaces
  implicit none

  interface
    ! DGESV solves the system of linear equations A*X = B
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: ipiv(n), info
      real(8), intent(inout) :: a(lda,n), b(ldb,nrhs)
    end subroutine dgesv
  end interface

end module lapack_interfaces
