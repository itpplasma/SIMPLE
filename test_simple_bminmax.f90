program test_simple_bminmax
  use find_bminmax_sub
  use simple, only: Tracer
  use simple_main, only: init_field
  implicit none

  type(Tracer) :: norb
  real(8) :: s, bmin, bmax
  
  ! Initialize VMEC field (required for magfie calls)
  call init_field(norb, 'wout.nc', 5, 5, 5, -1)
  
  print *, "Testing simple bminmax call"
  
  s = 0.5d0
  call find_bminmax(s, bmin, bmax)
  
  print *, "s=", s, " bmin=", bmin, " bmax=", bmax
  
  if (bmin > bmax) then
    error stop "bmin > bmax"
  end if
  
  print *, "Test passed"
  
end program test_simple_bminmax