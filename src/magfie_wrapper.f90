module magfie_wrapper
    ! Wrapper for magfie_sub with explicit subroutines that f90wrap can handle.
    ! F90wrap cannot wrap procedure pointers or re-exported items.
    use magfie_sub
    implicit none

    ! Export field type constants as module parameters
    integer, parameter :: field_type_test = TEST
    integer, parameter :: field_type_canflux = CANFLUX
    integer, parameter :: field_type_vmec = VMEC
    integer, parameter :: field_type_boozer = BOOZER
    integer, parameter :: field_type_meiss = MEISS
    integer, parameter :: field_type_albert = ALBERT

contains

    subroutine wrapper_init_magfie(field_type)
        ! Explicit wrapper for init_magfie
        integer, intent(in) :: field_type
        call init_magfie(field_type)
    end subroutine wrapper_init_magfie

    function get_field_type_vmec() result(ftype)
        ! Return VMEC field type constant
        integer :: ftype
        ftype = VMEC
    end function get_field_type_vmec

    function get_field_type_boozer() result(ftype)
        ! Return BOOZER field type constant
        integer :: ftype
        ftype = BOOZER
    end function get_field_type_boozer

    function get_field_type_canflux() result(ftype)
        ! Return CANFLUX field type constant
        integer :: ftype
        ftype = CANFLUX
    end function get_field_type_canflux

    function get_field_type_meiss() result(ftype)
        ! Return MEISS field type constant
        integer :: ftype
        ftype = MEISS
    end function get_field_type_meiss

    function get_field_type_albert() result(ftype)
        ! Return ALBERT field type constant
        integer :: ftype
        ftype = ALBERT
    end function get_field_type_albert

    function get_field_type_test() result(ftype)
        ! Return TEST field type constant
        integer :: ftype
        ftype = TEST
    end function get_field_type_test

end module magfie_wrapper
