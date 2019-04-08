module test_matlab

implicit none

contains

subroutine test(y, x)
  real(8), intent(out) :: y
  real(8), intent(in)  :: x
  
  y = 2d0*x
end subroutine test

end module test_matlab