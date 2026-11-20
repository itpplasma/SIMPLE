module coarse_coils
    ! Biot-Savart tables cost (coil segments) x (nodes). Tests that compare
    ! routes for the same coils keep every stride-th vertex of each coil and
    ! its zero-current closing vertex, which is still a closed coil set.
    use neo_biotsavart, only: coils_t, load_coils_from_file, save_coils_to_file
    implicit none
    private

    public :: coarsen_coils, write_coarse_coils

contains

    subroutine coarsen_coils(coils, stride)
        type(coils_t), intent(inout) :: coils
        integer, intent(in) :: stride
        logical, allocatable :: keep(:)
        integer :: i, first

        allocate (keep(size(coils%x)))
        keep = .false.
        first = 1
        do i = 1, size(coils%x)
            if (coils%current(i) /= 0.0d0) cycle
            keep(first:i - 1:stride) = .true.
            keep(i) = .true.
            first = i + 1
        end do
        coils%x = pack(coils%x, keep)
        coils%y = pack(coils%y, keep)
        coils%z = pack(coils%z, keep)
        coils%current = pack(coils%current, keep)
    end subroutine coarsen_coils

    subroutine write_coarse_coils(source, path, stride)
        character(*), intent(in) :: source, path
        integer, intent(in) :: stride
        type(coils_t) :: coils

        call load_coils_from_file(source, coils)
        call coarsen_coils(coils, stride)
        call save_coils_to_file(path, coils)
    end subroutine write_coarse_coils

end module coarse_coils
