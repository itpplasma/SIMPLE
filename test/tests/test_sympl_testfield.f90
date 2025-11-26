program test_sympl_testfield
  use simple_main, only : init_field
  use simple, only : Tracer
  use params, only : isw_field_type, field_input, integmode
  use magfie_sub, only : TEST
  use field_can_mod, only : evaluate

  implicit none

  type(Tracer) :: norb
  character(*), parameter :: vmec_file = WOUT_FILE
  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5

  ! Configure symplectic GC with TEST field to ensure initialization succeeds
  isw_field_type = TEST
  field_input = ''
  integmode = 1

  call init_field(norb, vmec_file, ans_s, ans_tp, amultharm, integmode)

  if (.not. associated(evaluate)) then
    print *, 'evaluate pointer not associated for TEST field'
    stop 1
  end if

  print *, 'TEST field initialized for symplectic mode'
end program test_sympl_testfield
