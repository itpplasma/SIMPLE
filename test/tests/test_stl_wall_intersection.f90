program test_stl_wall_intersection
    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use stl_wall_intersection, only: stl_wall_init, stl_wall_finalize, &
                                     stl_wall_first_hit_segment_with_normal
    use test_utils, only: check, check_close, check_close_vec
    use simple_main, only: main_wall => wall, chartmap_cart_scale_to_m, &
                           locate_wall_segment
    use params, only: wall_hit, wall_hit_cart, wall_hit_normal_cart, &
                      wall_hit_cos_incidence, wall_hit_angle_rad, &
                      wall_query_rho_lcfs
    use reference_coordinates, only: ref_coords
    use cartesian_coordinates, only: cartesian_coordinate_system_t
    implicit none

    integer :: errors
    logical :: hit
    real(dp) :: p0(3), p1(3), hit_p(3)
    real(dp) :: normal(3), dotx
    real(dp) :: z(5), z_start(5), z_end(5), hit_step
    real(dp) :: u_start(3), u_end(3)

    errors = 0
    call write_tetra_stl_m("tetra_m.stl")
    call stl_wall_init(main_wall, "tetra_m.stl", 1.0_dp)

    p0 = [0.25_dp, 0.25_dp, 0.25_dp]
    p1 = [2.0_dp, 0.25_dp, 0.25_dp]
    call stl_wall_first_hit_segment_with_normal( &
        main_wall, p0, p1, hit, hit_p, normal)
    call check(hit, "expected intersection for meters STL", errors)
    call check_close_vec(hit_p, [0.5_dp, 0.25_dp, 0.25_dp], 1.0e-12_dp, &
        "meters STL hit point mismatch", errors)
    call check_close(norm2(normal), 1.0_dp, 1.0e-12_dp, &
        "meters STL normal not unit length", errors)
    dotx = abs(normal(1))
    call check_close(dotx, 1.0_dp/sqrt(3.0_dp), 1.0e-12_dp, &
        "meters STL normal mismatch", errors)

    allocate (cartesian_coordinate_system_t :: ref_coords)
    allocate (wall_hit(1), wall_hit_cart(3, 1), wall_hit_normal_cart(3, 1), &
              wall_hit_cos_incidence(1), wall_hit_angle_rad(1))
    wall_hit = 0_int8
    wall_hit_cart = 0.0_dp
    wall_hit_normal_cart = 0.0_dp
    wall_hit_cos_incidence = 0.0_dp
    wall_hit_angle_rad = 0.0_dp
    chartmap_cart_scale_to_m = 1.0_dp
    z_start = [p0, 4.0_dp, 5.0_dp]
    z_end = [p1, 11.0_dp, 19.0_dp]

    ! A Cartesian chord can intersect a non-convex external wall even though
    ! its accepted endpoint remains inside the LCFS. Such a chord is not a
    ! physical wall loss and must not be queried.
    wall_query_rho_lcfs = 1.0_dp
    u_start = [0.25_dp, 0.25_dp, 0.25_dp]
    u_end = [0.50_dp, 0.25_dp, 0.25_dp]
    z = z_end
    call locate_wall_segment(z, z_start, z_end, 1, p0, u_start, p1, u_end, &
        6.0_dp, 8.0_dp, 7.0_dp, hit, hit_step)
    call check(.not. hit, "interior-LCFS chord was queried against wall", errors)
    call check(wall_hit(1) == 0_int8, "interior chord set wall hit flag", errors)
    call check_close_vec(z, z_end, 1.0e-12_dp, &
        "interior chord changed orbit state", errors)
    call check_close(hit_step, 15.0_dp, 1.0e-12_dp, &
        "interior chord changed event time", errors)

    u_start = p0
    u_end = p1
    z = z_end
    call locate_wall_segment(z, z_start, z_end, 1, p0, u_start, p1, u_end, 6.0_dp, &
        8.0_dp, 7.0_dp, hit, hit_step)
    call check(hit, "wall event locator missed crossing", errors)
    call check_close_vec(z(1:3), hit_p, 1.0e-12_dp, &
        "wall event position mismatch", errors)
    call check_close_vec(z(4:5), [5.0_dp, 7.0_dp], 1.0e-12_dp, &
        "wall event phase-space interpolation mismatch", errors)
    call check_close(hit_step, 9.0_dp, 1.0e-12_dp, &
        "wall event time mismatch", errors)
    call check(wall_hit(1) == 1_int8, "wall hit flag missing", errors)
    call check_close_vec(wall_hit_cart(:, 1), hit_p, 1.0e-12_dp, &
        "stored wall hit position mismatch", errors)
    call check_close_vec(wall_hit_normal_cart(:, 1), normal, 1.0e-12_dp, &
        "stored wall normal mismatch", errors)
    call check_close(wall_hit_cos_incidence(1), 1.0_dp/sqrt(3.0_dp), &
        1.0e-12_dp, "wall incidence cosine mismatch", errors)
    call check_close(wall_hit_angle_rad(1), acos(1.0_dp/sqrt(3.0_dp)), &
        1.0e-12_dp, "wall incidence angle mismatch", errors)

    u_start = [0.25_dp, 6.2_dp, 3.0_dp]
    u_end = [2.0_dp, 0.1_dp, -3.0_dp]
    z = z_end
    call locate_wall_segment(z, z_start, z_end, 1, p0, u_start, p1, u_end, &
        6.0_dp, 8.0_dp, 7.0_dp, hit, hit_step)
    call check_close(z(2), 6.2_dp + (0.1_dp - 6.2_dp + &
        2.0_dp*acos(-1.0_dp))/7.0_dp, 1.0e-12_dp, &
        "wall event poloidal seam mismatch", errors)
    call check_close(z(3), 3.0_dp, 1.0e-12_dp, &
        "wall event field-period seam mismatch", errors)
    call stl_wall_finalize(main_wall)
    wall_query_rho_lcfs = -1.0_dp

    call write_tetra_stl_mm("tetra_mm.stl")
    call stl_wall_init(main_wall, "tetra_mm.stl", 1.0e-3_dp)
    call stl_wall_first_hit_segment_with_normal( &
        main_wall, p0, p1, hit, hit_p, normal)
    call check(hit, "expected intersection for millimeters STL (scaled)", errors)
    call check_close_vec(hit_p, [0.5_dp, 0.25_dp, 0.25_dp], 1.0e-12_dp, &
        "mm STL scaled hit point mismatch", errors)
    call check_close(norm2(normal), 1.0_dp, 1.0e-12_dp, &
        "mm STL normal not unit length", errors)
    call stl_wall_finalize(main_wall)
    if (errors /= 0) error stop 1
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

    function norm2(v) result(n)
        real(dp), intent(in) :: v(3)
        real(dp) :: n
        n = sqrt(sum(v*v))
    end function norm2

end program test_stl_wall_intersection
