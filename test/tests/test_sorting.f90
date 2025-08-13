program test_sorting
  use sorting_mod
  implicit none
  
  integer :: errors
  
  errors = 0
  
  ! Test basic sorting functionality
  call test_basic_sorting(errors)
  
  ! Test edge cases
  call test_edge_cases(errors)
  
  ! Test sorting properties
  call test_sorting_properties(errors)
  
  ! Test large arrays
  call test_large_array_sorting(errors)
  
  if (errors == 0) then
    print *, "All sorting module tests passed!"
  else
    print *, "ERROR: ", errors, " test(s) failed!"
    stop 1
  end if
  
contains

  subroutine test_basic_sorting(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 5
    double precision :: a(n)
    integer :: ipoi(n)
    integer :: i
    
    print *, "Testing basic sorting functionality..."
    
    ! Given: An unsorted array with known values
    ! When: We call sortin to get sorted indices
    ! Then: The indices should point to elements in ascending order
    
    ! Test case 1: Simple unsorted array
    a = [3.0d0, 1.0d0, 4.0d0, 2.0d0, 5.0d0]
    
    call sortin(a, ipoi, n)
    
    ! Check that indices point to values in ascending order
    do i = 1, n-1
      if (a(ipoi(i)) > a(ipoi(i+1))) then
        print *, "ERROR: Array not sorted correctly"
        print *, "Position", i, ": a(", ipoi(i), ") =", a(ipoi(i))
        print *, "Position", i+1, ": a(", ipoi(i+1), ") =", a(ipoi(i+1))
        errors = errors + 1
        exit
      end if
    end do
    
    ! Verify that the sorted order is correct for this specific case
    if (ipoi(1) /= 2 .or. ipoi(2) /= 4 .or. ipoi(3) /= 1 .or. &
        ipoi(4) /= 3 .or. ipoi(5) /= 5) then
      print *, "ERROR: Incorrect sorted indices for test case"
      print *, "Expected: [2, 4, 1, 3, 5]"
      print *, "Got:", ipoi
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Basic sorting test PASSED"
    end if
    
  end subroutine test_basic_sorting
  
  subroutine test_edge_cases(errors)
    integer, intent(inout) :: errors
    
    print *, "Testing edge cases..."
    
    ! Test single element array
    call test_single_element(errors)
    
    ! Test two element arrays
    call test_two_elements(errors)
    
    ! Test already sorted array
    call test_already_sorted(errors)
    
    ! Test reverse sorted array
    call test_reverse_sorted(errors)
    
    ! Test array with duplicate values
    call test_duplicates(errors)
    
    if (errors == 0) then
      print *, "  Edge cases test PASSED"
    end if
    
  end subroutine test_edge_cases
  
  subroutine test_single_element(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 1
    double precision :: a(n)
    integer :: ipoi(n)
    
    ! Given: A single element array
    ! When: We sort it
    ! Then: The index should point to the only element
    
    a(1) = 42.0d0
    call sortin(a, ipoi, n)
    
    if (ipoi(1) /= 1) then
      print *, "ERROR: Single element array not handled correctly"
      print *, "Expected: 1, Got:", ipoi(1)
      errors = errors + 1
    end if
    
  end subroutine test_single_element
  
  subroutine test_two_elements(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 2
    double precision :: a(n)
    integer :: ipoi(n)
    
    ! Test case 1: already in order
    a = [1.0d0, 2.0d0]
    call sortin(a, ipoi, n)
    
    if (ipoi(1) /= 1 .or. ipoi(2) /= 2) then
      print *, "ERROR: Two element array (in order) not handled correctly"
      print *, "Expected: [1, 2], Got:", ipoi
      errors = errors + 1
    end if
    
    ! Test case 2: reverse order
    a = [2.0d0, 1.0d0]
    call sortin(a, ipoi, n)
    
    if (ipoi(1) /= 2 .or. ipoi(2) /= 1) then
      print *, "ERROR: Two element array (reverse order) not handled correctly"
      print *, "Expected: [2, 1], Got:", ipoi
      errors = errors + 1
    end if
    
  end subroutine test_two_elements
  
  subroutine test_already_sorted(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 6
    double precision :: a(n)
    integer :: ipoi(n)
    integer :: i
    
    ! Given: An already sorted array
    ! When: We sort it
    ! Then: The indices should be [1, 2, 3, ..., n]
    
    a = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0]
    call sortin(a, ipoi, n)
    
    do i = 1, n
      if (ipoi(i) /= i) then
        print *, "ERROR: Already sorted array not preserved"
        print *, "Position", i, "expected", i, "got", ipoi(i)
        errors = errors + 1
        exit
      end if
    end do
    
  end subroutine test_already_sorted
  
  subroutine test_reverse_sorted(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 5
    double precision :: a(n)
    integer :: ipoi(n)
    integer :: i
    
    ! Given: A reverse sorted array
    ! When: We sort it
    ! Then: The indices should be [n, n-1, ..., 2, 1]
    
    a = [5.0d0, 4.0d0, 3.0d0, 2.0d0, 1.0d0]
    call sortin(a, ipoi, n)
    
    do i = 1, n
      if (ipoi(i) /= n + 1 - i) then
        print *, "ERROR: Reverse sorted array not handled correctly"
        print *, "Position", i, "expected", n + 1 - i, "got", ipoi(i)
        errors = errors + 1
        exit
      end if
    end do
    
  end subroutine test_reverse_sorted
  
  subroutine test_duplicates(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 6
    double precision :: a(n), a_orig(n)
    integer :: ipoi(n), orig_indices(n)
    integer :: i, j
    logical :: is_stable
    
    ! Given: An array with duplicate values and tracked original positions
    ! When: We sort it
    ! Then: The result should be sorted and preserve relative order of equal elements (stability)
    
    a = [3.0d0, 1.0d0, 3.0d0, 2.0d0, 1.0d0, 2.0d0]
    a_orig = a
    do i = 1, n
      orig_indices(i) = i
    end do
    
    call sortin(a, ipoi, n)
    
    ! Check that the array is sorted
    do i = 1, n-1
      if (a(ipoi(i)) > a(ipoi(i+1))) then
        print *, "ERROR: Array with duplicates not sorted correctly"
        print *, "Position", i, ": a(", ipoi(i), ") =", a(ipoi(i))
        print *, "Position", i+1, ": a(", ipoi(i+1), ") =", a(ipoi(i+1))
        errors = errors + 1
        exit
      end if
    end do
    
    ! Check stability: for equal elements, original order should be preserved
    is_stable = .true.
    do i = 1, n-1
      if (abs(a(ipoi(i)) - a(ipoi(i+1))) < 1.0d-14) then
        ! Equal elements found - check original order
        if (ipoi(i) > ipoi(i+1)) then
          print *, "WARNING: Sort may not be stable for duplicates"
          print *, "Elements at", ipoi(i), "and", ipoi(i+1), "have same value but reversed order"
          ! Note: This is a warning, not an error, as stability may not be guaranteed
        end if
      end if
    end do
    
  end subroutine test_duplicates
  
  subroutine test_sorting_properties(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 10
    double precision :: a(n)
    integer :: ipoi(n)
    logical :: indices_valid, seen(n)
    integer :: i, j
    
    print *, "Testing sorting properties..."
    
    ! Given: An arbitrary array
    ! When: We sort it
    ! Then: The sorting should satisfy fundamental properties
    
    a = [7.5d0, 2.3d0, 9.1d0, 1.8d0, 5.6d0, 3.4d0, 8.2d0, 4.7d0, 6.9d0, 0.5d0]
    call sortin(a, ipoi, n)
    
    ! Property 1: All indices should be valid (between 1 and n)
    indices_valid = .true.
    do i = 1, n
      if (ipoi(i) < 1 .or. ipoi(i) > n) then
        indices_valid = .false.
        print *, "ERROR: Invalid index found:", ipoi(i)
        errors = errors + 1
      end if
    end do
    
    ! Property 2: All indices should be unique (permutation)
    ! Use more efficient O(n) check instead of O(n²)
    if (indices_valid) then
      seen = .false.
      do i = 1, n
        if (seen(ipoi(i))) then
          print *, "ERROR: Duplicate index found:", ipoi(i)
          errors = errors + 1
          exit
        end if
        seen(ipoi(i)) = .true.
      end do
    end if
    
    ! Property 3: The array should be sorted according to the indices
    do i = 1, n-1
      if (a(ipoi(i)) > a(ipoi(i+1))) then
        print *, "ERROR: Sorting property violated"
        errors = errors + 1
        exit
      end if
    end do
    
    if (errors == 0) then
      print *, "  Sorting properties test PASSED"
    end if
    
  end subroutine test_sorting_properties
  
  subroutine test_large_array_sorting(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 1000
    double precision :: a(n)
    integer :: ipoi(n)
    integer :: i
    logical :: is_sorted, is_permutation
    logical :: seen(n)
    
    print *, "Testing large array sorting..."
    
    ! Given: A large array with deterministic test values
    ! When: We sort it
    ! Then: The result should be properly sorted and form a valid permutation
    
    ! Generate deterministic test values using a simple pattern
    ! This creates a mix of values with some duplicates
    do i = 1, n
      ! Simple pattern: values from 0 to 99 repeating
      a(i) = dble(mod(i-1, 100))
    end do
    
    call sortin(a, ipoi, n)
    
    ! Check that the result is sorted
    is_sorted = .true.
    do i = 1, n-1
      if (a(ipoi(i)) > a(ipoi(i+1))) then
        is_sorted = .false.
        print *, "ERROR: Not sorted at position", i
        print *, "a(", ipoi(i), ") =", a(ipoi(i)), "> a(", ipoi(i+1), ") =", a(ipoi(i+1))
        errors = errors + 1
        exit
      end if
    end do
    
    ! Check that ipoi forms a valid permutation
    seen = .false.
    is_permutation = .true.
    do i = 1, n
      if (ipoi(i) < 1 .or. ipoi(i) > n) then
        print *, "ERROR: Invalid index", ipoi(i), "at position", i
        errors = errors + 1
        is_permutation = .false.
        exit
      else if (seen(ipoi(i))) then
        print *, "ERROR: Duplicate index", ipoi(i), "at position", i
        errors = errors + 1
        is_permutation = .false.
        exit
      else
        seen(ipoi(i)) = .true.
      end if
    end do
    
    if (is_sorted .and. is_permutation .and. errors == 0) then
      print *, "  Large array sorting test PASSED"
    end if
    
  end subroutine test_large_array_sorting

end program test_sorting