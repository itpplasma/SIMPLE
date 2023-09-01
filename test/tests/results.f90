! Add anything wanted for result parsing of tests here.

module test_results
    implicit none
    save
    
    public result_handling
    interface result_handling
        module procedure :: result_printer
    end interface result_handling

    contains
    
subroutine result_printer(res, num_res, info)
    logical :: res(num_res)
    integer num_res, i
    Character (len = *), optional :: info
    
    write (*,'(A)', advance="no") trim(info)
    write (*, '(A)', advance="no") ":    "
    do i=1,num_res
        if (res(i)) then
            write (*,fmt='(A)', advance="no") "o"
        else
            write (*,fmt='(A)', advance="no") "."
        end if
    end do
    write(*,*) ""
    
end subroutine result_printer
    
    
end module test_results
