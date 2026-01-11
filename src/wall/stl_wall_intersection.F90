module stl_wall_intersection

    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int
    use, intrinsic :: iso_c_binding, only: c_double, c_char, c_null_char
    implicit none

    private
    public :: stl_wall_t
    public :: stl_wall_init
    public :: stl_wall_finalize
    public :: stl_wall_first_hit_segment

    type :: stl_wall_t
        type(c_ptr) :: handle = c_null_ptr
    end type stl_wall_t

#ifdef SIMPLE_ENABLE_CGAL
    interface
        function stl_wall_create(filename_c, scale_to_m) &
            bind(C, name="stl_wall_create") result(h)
            import :: c_ptr, c_char, c_double
            character(kind=c_char), intent(in) :: filename_c(*)
            real(c_double), value, intent(in) :: scale_to_m
            type(c_ptr) :: h
        end function stl_wall_create

        subroutine stl_wall_destroy(h) bind(C, name="stl_wall_destroy")
            import :: c_ptr
            type(c_ptr), value, intent(in) :: h
        end subroutine stl_wall_destroy

        function stl_wall_first_hit_segment_c(h, p0_m, p1_m, hit_m) &
            bind(C, name="stl_wall_first_hit_segment") result(hit)
            import :: c_ptr, c_double, c_int
            type(c_ptr), value, intent(in) :: h
            real(c_double), intent(in) :: p0_m(3)
            real(c_double), intent(in) :: p1_m(3)
            real(c_double), intent(out) :: hit_m(3)
            integer(c_int) :: hit
        end function stl_wall_first_hit_segment_c
    end interface
#endif

contains

    subroutine stl_wall_init(wall, filename, scale_to_m)
        type(stl_wall_t), intent(inout) :: wall
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: scale_to_m

        character(kind=c_char), allocatable :: filename_c(:)
        integer :: n, i

        if (len_trim(filename) == 0) then
            error stop "stl_wall_init: filename must be non-empty"
        end if
        if (.not. (scale_to_m > 0.0_dp)) then
            error stop "stl_wall_init: scale_to_m must be > 0"
        end if

        n = len_trim(filename)
        allocate (filename_c(n + 1))
        do i = 1, n
            filename_c(i) = filename(i:i)
        end do
        filename_c(n + 1) = c_null_char

#ifdef SIMPLE_ENABLE_CGAL
        wall%handle = stl_wall_create(filename_c, real(scale_to_m, c_double))
        if (.not. c_associated(wall%handle)) then
            error stop "stl_wall_init: failed to create STL wall (CGAL)"
        end if
#else
        print *, "stl_wall_init: SIMPLE built without CGAL support."
        print *, "Rebuild with -DSIMPLE_ENABLE_CGAL=ON to use wall_input."
        error stop "stl_wall_init: CGAL disabled"
#endif
    end subroutine stl_wall_init

    subroutine stl_wall_finalize(wall)
        type(stl_wall_t), intent(inout) :: wall

#ifdef SIMPLE_ENABLE_CGAL
        if (c_associated(wall%handle)) then
            call stl_wall_destroy(wall%handle)
        end if
#endif
        wall%handle = c_null_ptr
    end subroutine stl_wall_finalize

    subroutine stl_wall_first_hit_segment(wall, p0_m, p1_m, hit, hit_m)
        type(stl_wall_t), intent(in) :: wall
        real(dp), intent(in) :: p0_m(3)
        real(dp), intent(in) :: p1_m(3)
        logical, intent(out) :: hit
        real(dp), intent(out) :: hit_m(3)

        integer(c_int) :: hit_i

#ifdef SIMPLE_ENABLE_CGAL
        if (.not. c_associated(wall%handle)) then
            error stop "stl_wall_first_hit_segment: wall not initialized"
        end if
        hit_i = stl_wall_first_hit_segment_c(wall%handle, &
                                             real(p0_m, c_double), &
                                             real(p1_m, c_double), &
                                             hit_m)
        hit = (hit_i /= 0)
#else
        hit = .false.
        hit_m = 0.0_dp
        error stop "stl_wall_first_hit_segment: CGAL disabled"
#endif
    end subroutine stl_wall_first_hit_segment

end module stl_wall_intersection
