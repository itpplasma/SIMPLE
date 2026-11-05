subroutine simple_entry_main()
    use simple_main, only: main
    use simple_mpi, only: simple_mpi_init, simple_mpi_finalize

    implicit none

    call simple_mpi_init
    call main
    call simple_mpi_finalize
end subroutine simple_entry_main
