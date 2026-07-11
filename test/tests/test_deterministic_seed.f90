program test_deterministic_seed
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use params, only: read_config

    implicit none

    real(dp) :: first, repeated, different

    call write_config("seed_12345.in", 12345)
    call write_config("seed_22345.in", 22345)

    call read_config("seed_12345.in")
    call random_number(first)
    call read_config("seed_12345.in")
    call random_number(repeated)

    if (first /= repeated) then
        error stop "same ran_seed does not reproduce the random stream"
    end if

    call read_config("seed_22345.in")
    call random_number(different)

    if (first == different) then
        error stop "different ran_seed values produce the same random stream"
    end if

    call delete_file("seed_12345.in")
    call delete_file("seed_22345.in")

contains

    subroutine write_config(path, seed)
        character(len=*), intent(in) :: path
        integer, intent(in) :: seed
        integer :: unit

        open (newunit=unit, file=path, status="replace", action="write")
        write (unit, '(A)') "&config"
        write (unit, '(A)') "  deterministic = .true."
        write (unit, '(A,I0)') "  ran_seed = ", seed
        write (unit, '(A)') "/"
        close (unit)
    end subroutine write_config

    subroutine delete_file(path)
        character(len=*), intent(in) :: path
        integer :: unit

        open (newunit=unit, file=path, status="old")
        close (unit, status="delete")
    end subroutine delete_file
end program test_deterministic_seed
