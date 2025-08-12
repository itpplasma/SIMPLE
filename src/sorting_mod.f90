module sorting_mod
  implicit none

  interface
    subroutine sortin(a, ipoi, n)
      integer, intent(in) :: n
      double precision, intent(in) :: a(n)
      integer, intent(inout) :: ipoi(n)
    end subroutine sortin
  end interface

end module sorting_mod
