

module test_example
   use funit
   implicit none
   
contains

   @test
   subroutine test_example_1()
      integer :: i, j, k

      i = 3
      j = 2*i
      k = 6
      
      @assertEqual(k, j, message = "test example 1")

   end subroutine test_example_1
   
end module test_example
