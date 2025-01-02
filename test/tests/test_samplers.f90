
program test_samplers
    use samplers
    use params, only: isw_field_type

    implicit none

    isw_field_type = 1

    call test_sample_load_save
    call test_sample_volume
    call test_sample_surface

    contains

    subroutine test_sample_load_save()
        integer, parameter :: npart = 5
        double precision, dimension(5,npart):: zstart
        double precision, dimension(5,npart):: zreadout
        integer :: i

        zstart = reshape([(1.0*i, i = 1, 25)], shape(zstart))
        call save_starting_points(zstart)
        call load_starting_points(zreadout)

        if (any(zstart /= zreadout)) then
            print *, "Error: read and write do not match"
            error stop
        end if

    end subroutine test_sample_load_save

    subroutine test_sample_volume()
        double precision, parameter :: s_in = 0.3
        double precision, parameter :: s_out = 0.5
        integer, parameter :: npart = 10
        double precision, dimension(5,npart):: zstart
        double precision, dimension(5,npart):: zreadout

        call sample_volume_single(zstart, s_in, s_out)
        call load_starting_points(zreadout)

        if (any(zreadout(1,:) < s_in) .or. any(zreadout(1,:) > s_out)) then
            print *, "Error: points are outside the specified volume"
            error stop
        end if

    end subroutine test_sample_volume

    subroutine test_sample_surface()
        use simple_main, only: init_starting_surf
        use params, only: sbeg
        ! check if all points are on the specified surface
        integer, parameter :: npart = 100
        double precision :: s_in
        double precision, dimension(5,npart):: zstart
        double precision, dimension(5,npart):: zreadout

        print *, "test_sample_surface not implemented"
        return

        s_in = sbeg(1)
        call init_starting_surf
        call sample_surface_fieldline(zstart)
        call load_starting_points(zreadout)

        if (any(zreadout(1,:) /= s_in)) then
            print *, "Error: points are outside the specified volume"
            error stop
        end if

    end subroutine test_sample_surface
end program test_samplers
