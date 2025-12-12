module params_wrapper
    ! Wrapper module for params arrays that f90wrap struggles with.
    ! Provides simple subroutines to get/set array elements safely.
    use params, only: dp, sbeg, zstart, zend, times_lost, trap_par, perp_inv, &
                      iclass, class_passing, class_lost, ntestpart
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

    subroutine get_times_lost_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        real(dp), dimension(:), intent(out) :: data
        data(1:n_particles) = times_lost(1:n_particles)
    end subroutine get_times_lost_bulk

    subroutine get_trap_par_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        real(dp), dimension(:), intent(out) :: data
        data(1:n_particles) = trap_par(1:n_particles)
    end subroutine get_trap_par_bulk

    subroutine get_perp_inv_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        real(dp), dimension(:), intent(out) :: data
        data(1:n_particles) = perp_inv(1:n_particles)
    end subroutine get_perp_inv_bulk

    subroutine get_iclass_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        integer, dimension(:,:), intent(out) :: data
        data(:, 1:n_particles) = iclass(:, 1:n_particles)
    end subroutine get_iclass_bulk

    subroutine get_class_passing_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        logical, dimension(:), intent(out) :: data
        data(1:n_particles) = class_passing(1:n_particles)
    end subroutine get_class_passing_bulk

    subroutine get_class_lost_bulk(n_particles, data)
        integer, intent(in) :: n_particles
        logical, dimension(:), intent(out) :: data
        data(1:n_particles) = class_lost(1:n_particles)
    end subroutine get_class_lost_bulk

    ! Convert a canonical trajectory to reference coordinates in bulk.
    ! traj_can and traj_ref are (5, n_steps) arrays.
    subroutine integ_traj_to_ref(traj_can, traj_ref)
        use field_can_mod, only: integ_to_ref
        real(dp), intent(in) :: traj_can(:,:)
        real(dp), intent(out) :: traj_ref(:,:)
        integer :: it, n_steps
        real(dp) :: xref(3)

        if (size(traj_can, 1) < 5) error stop 'integ_traj_to_ref: traj_can dim1 < 5'
        if (size(traj_ref, 1) < 5) error stop 'integ_traj_to_ref: traj_ref dim1 < 5'

        n_steps = size(traj_can, 2)
        if (size(traj_ref, 2) /= n_steps) then
            error stop 'integ_traj_to_ref: dim2 mismatch'
        end if

        do it = 1, n_steps
            call integ_to_ref(traj_can(1:3, it), xref)
            traj_ref(1:3, it) = xref
            traj_ref(4:5, it) = traj_can(4:5, it)
        end do
    end subroutine integ_traj_to_ref

end module params_wrapper
