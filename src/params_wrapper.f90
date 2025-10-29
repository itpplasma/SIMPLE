module params_wrapper
    ! Wrapper module for params arrays that f90wrap struggles with.
    ! Provides simple subroutines to get/set array elements safely.
    use params, only: dp, sbeg, zstart, zend, ntestpart
    implicit none

contains

    subroutine get_sbeg(index, val)
        integer, intent(in) :: index
        real(dp), intent(out) :: val
        val = sbeg(index)
    end subroutine get_sbeg

    subroutine set_sbeg(index, val)
        integer, intent(in) :: index
        real(dp), intent(in) :: val
        sbeg(index) = val
    end subroutine set_sbeg

    ! Wrapper for zstart to avoid allocatable array access issues
    subroutine get_zstart_column(ipart, column)
        integer, intent(in) :: ipart
        real(dp), dimension(:), intent(out) :: column
        column = zstart(:, ipart)
    end subroutine get_zstart_column

    subroutine set_zstart_column(ipart, column)
        integer, intent(in) :: ipart
        real(dp), dimension(:), intent(in) :: column
        zstart(:, ipart) = column
    end subroutine set_zstart_column

    subroutine get_zstart_shape(dim1, dim2)
        integer, intent(out) :: dim1, dim2
        dim1 = size(zstart, 1)
        dim2 = size(zstart, 2)
    end subroutine get_zstart_shape

    subroutine get_zend_column(ipart, column)
        integer, intent(in) :: ipart
        real(dp), dimension(:), intent(out) :: column
        column = zend(:, ipart)
    end subroutine get_zend_column

    subroutine get_zend_shape(dim1, dim2)
        integer, intent(out) :: dim1, dim2
        dim1 = size(zend, 1)
        dim2 = size(zend, 2)
    end subroutine get_zend_shape

    ! Bulk array access for performance
    subroutine get_zstart_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        real(dp), dimension(:,:), intent(out) :: data
        data(:, 1:n_particles) = zstart(:, 1:n_particles)
    end subroutine get_zstart_bulk

    subroutine set_zstart_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        real(dp), dimension(:,:), intent(in) :: data
        zstart(:, 1:n_particles) = data(:, 1:n_particles)
    end subroutine set_zstart_bulk

    subroutine get_zend_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        real(dp), dimension(:,:), intent(out) :: data
        data(:, 1:n_particles) = zend(:, 1:n_particles)
    end subroutine get_zend_bulk

end module params_wrapper
