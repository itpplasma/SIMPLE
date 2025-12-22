program simple
    implicit none

    interface
        subroutine simple_entry_main()
            implicit none
        end subroutine simple_entry_main
    end interface

    call simple_entry_main()
end program simple
