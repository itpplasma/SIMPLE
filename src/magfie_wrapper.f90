module magfie_wrapper
    ! Thin wrapper for magfie_sub that f90wrap can handle.
    ! F90wrap cannot wrap procedure pointers, so we expose only what Python needs:
    ! - init_magfie subroutine
    ! - Field type constants
    use magfie_sub, only: init_magfie, TEST, CANFLUX, VMEC, BOOZER, MEISS, ALBERT
    implicit none
end module magfie_wrapper
