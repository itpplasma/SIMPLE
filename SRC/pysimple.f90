module pysimple
    use neo_orb, only: NeoOrb
    implicit none

    type(NeoOrb) :: norb
end module pysimple

subroutine init(vmec_file, ns_s, ns_tp, multharm, integmode)
    use neo_orb, only: init_field
    use pysimple, only: norb
    implicit none

    character(len=*), intent(in) :: vmec_file
    integer, intent(in) :: ns_s, ns_tp, multharm, integmode

    call init_field(norb, vmec_file, ns_s, ns_tp, multharm, integmode)
end subroutine init
