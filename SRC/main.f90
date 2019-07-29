program main
#ifdef _OPENMP
  use omp_lib
#endif
  use common, only: twopi
  use neo_orb, only: init_field, init_params, timestep
  use cut_detector, only: init_cut => init, trace_to_cut
  use new_vmec_stuff_mod, only: rmajor

  implicit none

  real(8) :: s, th, ph, lam
  integer :: ierr

  integer :: npoiper2, ntimstep, ncut, i, cut_type
  real(8) :: rbig, dtau, dtaumax, t

  real(8) :: starttime, endtime
  real(8) :: var_tip(6)

  call init_field(5, 5, 3, 0)

  npoiper2 = 64
  rbig = rmajor*1.0d2
  dtaumax = twopi*rbig/npoiper2
  dtau = dtaumax

  call init_params(2, 4, 3.5d6, dtau, dtaumax, 1d-8)  ! fusion alphas

  s = 0.5d0
  th = 0d0
  ph = 0.314d0
  lam = 0.22d0

  ntimstep = 10000

  starttime = omp_get_wtime()
  do i = 1,ntimstep
    call timestep(s, th, ph, lam, ierr)
  end do
  endtime = omp_get_wtime()

  print *, s, th, ph, lam, ierr
  write(*,*) 'Time elapsed: ', endtime - starttime,' s'

  ncut = 1000
  call init_cut

  starttime = omp_get_wtime()
  do i = 1,ncut
    call trace_to_cut(t, var_tip, cut_type, ierr)
  end do
  endtime = omp_get_wtime()
  write(*,*) 'Time elapsed: ', endtime - starttime,' s'

end program main
