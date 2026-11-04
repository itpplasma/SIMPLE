program test_cgal_openmp_link
    implicit none
    integer :: value

    value = 0
    !$omp target map(tofrom:value)
    value = 7
    !$omp end target
    if (value /= 7) error stop 'OpenMP target fallback did not execute'
end program test_cgal_openmp_link
