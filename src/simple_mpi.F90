module simple_mpi
    ! Optional MPI layer (SIMPLE_ENABLE_MPI). Particle i belongs to rank
    ! mod(i - 1, mpi_ranks); rank 0 collects every per-particle result and is
    ! the only rank that writes files. Without MPI this is a single rank.
    use, intrinsic :: iso_fortran_env, only: dp => real64
#ifdef SIMPLE_ENABLE_MPI
    use mpi
#endif
    implicit none
    private

    integer, public, protected :: mpi_rank = 0, mpi_ranks = 1
    ! Rank within the node; ranks on one node start at different GPUs.
    integer, public, protected :: mpi_node_rank = 0

    public :: simple_mpi_init, simple_mpi_finalize, mpi_is_root
    public :: mpi_owned_count, mpi_broadcast_r64, mpi_sum_r64
    public :: mpi_gather_r64, mpi_gather_int, mpi_reduce_sum_r64, mpi_reduce_sum_int

contains

    subroutine simple_mpi_init
#ifdef SIMPLE_ENABLE_MPI
        integer :: provided, ierr, node_comm
        logical :: started

        call MPI_Initialized(started, ierr)
        ! Only the master thread of a rank calls MPI.
        if (.not. started) call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, mpi_ranks, ierr)
        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                                 MPI_INFO_NULL, node_comm, ierr)
        call MPI_Comm_rank(node_comm, mpi_node_rank, ierr)
        call MPI_Comm_free(node_comm, ierr)
#endif
    end subroutine simple_mpi_init

    subroutine simple_mpi_finalize
#ifdef SIMPLE_ENABLE_MPI
        integer :: ierr
        logical :: done

        call MPI_Finalized(done, ierr)
        if (.not. done) call MPI_Finalize(ierr)
#endif
    end subroutine simple_mpi_finalize

    logical function mpi_is_root()
        mpi_is_root = mpi_rank == 0
    end function mpi_is_root

    pure integer function mpi_owned_count(rank, n)
        ! Number of the particles 1..n that belong to rank.
        integer, intent(in) :: rank, n

        mpi_owned_count = max(0, (n - rank - 1)/mpi_ranks + 1)
        if (rank >= n) mpi_owned_count = 0
    end function mpi_owned_count

    subroutine mpi_broadcast_r64(a, n)
        integer, intent(in) :: n
        real(dp), intent(inout) :: a(n)
#ifdef SIMPLE_ENABLE_MPI
        integer :: ierr

        if (mpi_ranks > 1) call MPI_Bcast(a, n, MPI_DOUBLE_PRECISION, 0, &
                                          MPI_COMM_WORLD, ierr)
#endif
    end subroutine mpi_broadcast_r64

    function mpi_sum_r64(x) result(total)
        real(dp), intent(in) :: x
        real(dp) :: total
#ifdef SIMPLE_ENABLE_MPI
        integer :: ierr

        call MPI_Allreduce(x, total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)
#else
        total = x
#endif
    end function mpi_sum_r64

    subroutine mpi_reduce_sum_r64(a, n)
        ! Sum a over the ranks into rank 0.
        integer, intent(in) :: n
        real(dp), intent(inout) :: a(n)
#ifdef SIMPLE_ENABLE_MPI
        real(dp) :: unused(1)
        integer :: ierr

        if (mpi_ranks == 1) return
        if (mpi_is_root()) then
            call MPI_Reduce(MPI_IN_PLACE, a, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                            0, MPI_COMM_WORLD, ierr)
        else
            call MPI_Reduce(a, unused, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                            MPI_COMM_WORLD, ierr)
        end if
#endif
    end subroutine mpi_reduce_sum_r64

    subroutine mpi_reduce_sum_int(a, n)
        integer, intent(in) :: n
        integer, intent(inout) :: a(n)
#ifdef SIMPLE_ENABLE_MPI
        integer :: unused(1), ierr

        if (mpi_ranks == 1) return
        if (mpi_is_root()) then
            call MPI_Reduce(MPI_IN_PLACE, a, n, MPI_INTEGER, MPI_SUM, 0, &
                            MPI_COMM_WORLD, ierr)
        else
            call MPI_Reduce(a, unused, n, MPI_INTEGER, MPI_SUM, 0, &
                            MPI_COMM_WORLD, ierr)
        end if
#endif
    end subroutine mpi_reduce_sum_int

    subroutine mpi_gather_r64(a, m, n)
        ! Collect the columns a(:, i) that other ranks own into rank 0.
        integer, intent(in) :: m, n
        real(dp), intent(inout) :: a(m, n)
#ifdef SIMPLE_ENABLE_MPI
        real(dp), allocatable :: mine(:, :), all(:, :)
        integer, allocatable :: counts(:), displs(:)
        integer :: r, ierr

        if (mpi_ranks == 1) return
        call gather_layout(m, n, counts, displs)
        allocate (mine(m, mpi_owned_count(mpi_rank, n)))
        mine = a(:, mpi_rank + 1:n:mpi_ranks)
        allocate (all(m, merge(n, 1, mpi_is_root())))
        call MPI_Gatherv(mine, size(mine), MPI_DOUBLE_PRECISION, all, counts, &
                         displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (.not. mpi_is_root()) return
        do r = 0, mpi_ranks - 1
            a(:, r + 1:n:mpi_ranks) = &
                all(:, displs(r + 1)/m + 1:(displs(r + 1) + counts(r + 1))/m)
        end do
#endif
    end subroutine mpi_gather_r64

    subroutine mpi_gather_int(a, m, n)
        integer, intent(in) :: m, n
        integer, intent(inout) :: a(m, n)
#ifdef SIMPLE_ENABLE_MPI
        integer, allocatable :: mine(:, :), all(:, :)
        integer, allocatable :: counts(:), displs(:)
        integer :: r, ierr

        if (mpi_ranks == 1) return
        call gather_layout(m, n, counts, displs)
        allocate (mine(m, mpi_owned_count(mpi_rank, n)))
        mine = a(:, mpi_rank + 1:n:mpi_ranks)
        allocate (all(m, merge(n, 1, mpi_is_root())))
        call MPI_Gatherv(mine, size(mine), MPI_INTEGER, all, counts, displs, &
                         MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (.not. mpi_is_root()) return
        do r = 0, mpi_ranks - 1
            a(:, r + 1:n:mpi_ranks) = &
                all(:, displs(r + 1)/m + 1:(displs(r + 1) + counts(r + 1))/m)
        end do
#endif
    end subroutine mpi_gather_int

#ifdef SIMPLE_ENABLE_MPI
    subroutine gather_layout(m, n, counts, displs)
        integer, intent(in) :: m, n
        integer, allocatable, intent(out) :: counts(:), displs(:)
        integer :: r

        allocate (counts(mpi_ranks), displs(mpi_ranks))
        displs(1) = 0
        do r = 0, mpi_ranks - 1
            counts(r + 1) = m*mpi_owned_count(r, n)
            if (r > 0) displs(r + 1) = displs(r) + counts(r)
        end do
    end subroutine gather_layout
#endif

end module simple_mpi
