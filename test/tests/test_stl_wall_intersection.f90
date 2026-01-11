program test_stl_wall_intersection
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use stl_wall_intersection, only: stl_wall_t, stl_wall_init, stl_wall_finalize, &
                                     stl_wall_first_hit_segment_with_normal
    implicit none

    type(stl_wall_t) :: wall
    logical :: hit
    real(dp) :: p0(3), p1(3), hit_p(3)
    real(dp) :: normal(3), dotx

    call write_tetra_stl_m("tetra_m.stl")
    call stl_wall_init(wall, "tetra_m.stl", 1.0_dp)

    p0 = [0.25_dp, 0.25_dp, 0.25_dp]
    p1 = [2.0_dp, 0.25_dp, 0.25_dp]
    call stl_wall_first_hit_segment_with_normal(wall, p0, p1, hit, hit_p, normal)
    call assert_true(hit, "expected intersection for meters STL")
    call assert_close_vec(hit_p, [0.5_dp, 0.25_dp, 0.25_dp], 1.0e-12_dp, &
                          "meters STL hit point mismatch")
    call assert_close_scalar(norm2(normal), 1.0_dp, 1.0e-12_dp, &
                             "meters STL normal not unit length")
    dotx = abs(normal(1))
    call assert_close_scalar(dotx, 1.0_dp/sqrt(3.0_dp), 1.0e-12_dp, &
                             "meters STL normal mismatch")
    call stl_wall_finalize(wall)

    call write_tetra_stl_mm("tetra_mm.stl")
    call stl_wall_init(wall, "tetra_mm.stl", 1.0e-3_dp)
    call stl_wall_first_hit_segment_with_normal(wall, p0, p1, hit, hit_p, normal)
    call assert_true(hit, "expected intersection for millimeters STL (scaled)")
    call assert_close_vec(hit_p, [0.5_dp, 0.25_dp, 0.25_dp], 1.0e-12_dp, &
                          "mm STL scaled hit point mismatch")
    call assert_close_scalar(norm2(normal), 1.0_dp, 1.0e-12_dp, &
                             "mm STL normal not unit length")
    call stl_wall_finalize(wall)
contains

    subroutine write_tetra_stl_m(path)
        character(len=*), intent(in) :: path
        integer :: u

        open (newunit=u, file=path, status="replace", action="write")
        write (u, "(A)") "solid tetra"
        call write_tri(u, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp], &
                       [0.0_dp, 1.0_dp, 0.0_dp])
        call write_tri(u, [0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 1.0_dp], &
                       [1.0_dp, 0.0_dp, 0.0_dp])
        call write_tri(u, [0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], &
                       [0.0_dp, 0.0_dp, 1.0_dp])
        call write_tri(u, [1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 1.0_dp], &
                       [0.0_dp, 1.0_dp, 0.0_dp])
        write (u, "(A)") "endsolid tetra"
        close (u)
    end subroutine write_tetra_stl_m

    subroutine write_tetra_stl_mm(path)
        character(len=*), intent(in) :: path
        integer :: u
        real(dp), parameter :: s = 1000.0_dp

        open (newunit=u, file=path, status="replace", action="write")
        write (u, "(A)") "solid tetra_mm"
        call write_tri(u, s*[0.0_dp, 0.0_dp, 0.0_dp], s*[1.0_dp, 0.0_dp, 0.0_dp], &
                       s*[0.0_dp, 1.0_dp, 0.0_dp])
        call write_tri(u, s*[0.0_dp, 0.0_dp, 0.0_dp], s*[0.0_dp, 0.0_dp, 1.0_dp], &
                       s*[1.0_dp, 0.0_dp, 0.0_dp])
        call write_tri(u, s*[0.0_dp, 0.0_dp, 0.0_dp], s*[0.0_dp, 1.0_dp, 0.0_dp], &
                       s*[0.0_dp, 0.0_dp, 1.0_dp])
        call write_tri(u, s*[1.0_dp, 0.0_dp, 0.0_dp], s*[0.0_dp, 0.0_dp, 1.0_dp], &
                       s*[0.0_dp, 1.0_dp, 0.0_dp])
        write (u, "(A)") "endsolid tetra_mm"
        close (u)
    end subroutine write_tetra_stl_mm

    subroutine write_tri(u, a, b, c)
        integer, intent(in) :: u
        real(dp), intent(in) :: a(3), b(3), c(3)

        write (u, "(A)") "  facet normal 0 0 0"
        write (u, "(A)") "    outer loop"
        write (u, "(A,3ES24.16)") "      vertex ", a
        write (u, "(A,3ES24.16)") "      vertex ", b
        write (u, "(A,3ES24.16)") "      vertex ", c
        write (u, "(A)") "    endloop"
        write (u, "(A)") "  endfacet"
    end subroutine write_tri

    subroutine assert_true(cond, msg)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: msg

        if (.not. cond) then
            print *, "ASSERTION FAILED: ", trim(msg)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close_vec(a, b, tol, msg)
        real(dp), intent(in) :: a(3), b(3)
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: msg

        if (maxval(abs(a - b)) > tol) then
            print *, "ASSERTION FAILED: ", trim(msg)
            print *, "got:      ", a
            print *, "expected: ", b
            error stop 1
        end if
    end subroutine assert_close_vec

    function norm2(v) result(n)
        real(dp), intent(in) :: v(3)
        real(dp) :: n
        n = sqrt(sum(v*v))
    end function norm2

    subroutine assert_close_scalar(a, b, tol, msg)
        real(dp), intent(in) :: a, b
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: msg

        if (abs(a - b) > tol) then
            print *, "ASSERTION FAILED: ", trim(msg)
            print *, "got:      ", a
            print *, "expected: ", b
            error stop 1
        end if
    end subroutine assert_close_scalar

end program test_stl_wall_intersection
