
program test_against_legacy_behaviour
! Compare outputs from old and current versions.
!  Assumes the files were renamed correctly.
    integer file_size_old
    integer file_size_new
    integer idx
    character(len=1024) :: old_line, new_line


    inquire(FILE='times_lost_old.dat', SIZE=file_size_old)
    inquire(FILE='times_lost_new.dat', SIZE=file_size_new)

    if (file_size_old /= file_size_new) error stop

    open(1,file='times_lost_old.dat',recl=1024)
    open(2,file='times_lost_new.dat',recl=1024)
        do idx=1,file_size_old
            read(1,*) old_line
            read(2,*) new_line
            if (old_line /= new_line) error stop
        enddo
    close(1)
    close(2)
end program test_against_legacy_behaviour
