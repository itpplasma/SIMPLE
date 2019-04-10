program main
  use neo_orb, only: init_field, init_params, timestep

  implicit none

  real(8) :: s, th, ph, lam
  integer :: ierr

  call init_field(5, 5, 3, -1)
  call init_params(2, 4, 3.5d6, 1d-2, 1d-2, 1d-8)  ! fusion alphas

  s = 0.5d0
  th = 0d0
  ph = 0d0
  lam = 0.5d0

  call timestep(s, th, ph, lam, ierr)

  print *, s, th, ph, lam, ierr
end program main
