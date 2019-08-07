program main
#ifdef _OPENMP
  use omp_lib
#endif
  use common, only: twopi
  use neo_orb, only: NeoOrb, init_field, init_params, init_integrator, timestep_sympl_z
  use field_can_mod, only: FieldCan
  use orbit_symplectic, only: SymplecticIntegrator, orbit_sympl_init
  use cut_detector, only: CutDetector, init_cut => init, trace_to_cut
  use new_vmec_stuff_mod, only: rmajor

  implicit none

  real(8) :: s, th, ph, lam
  integer :: ierr

  integer :: npoiper2, ntimstep, ncut, i, cut_type
  real(8) :: rbig, dtau, dtaumax, t

  real(8) :: starttime, endtime
  real(8) :: var_tip(6)

  real(8) :: z(5)

  type(NeoOrb) :: norb
  type(FieldCan) :: f
  type(SymplecticIntegrator) :: si
  type(CutDetector) :: cutter

  call init_field(norb, 5, 5, 3, 2)

  npoiper2 = 64
  rbig = rmajor*1.0d2
  dtaumax = twopi*rbig/npoiper2
  dtau = dtaumax

  call init_params(norb, 2, 4, 3.5d6, dtau, dtaumax, 1d-8)  ! fusion alphas

  s = 0.5d0
  th = 0d0
  ph = 0.314d0
  lam = 0.22d0

  z(1) = s
  z(2) = th
  z(3) = ph
  z(4) = 1d0
  z(5) = lam

  call init_integrator(norb, z)

  ntimstep = 10000

  starttime = omp_get_wtime()
  do i = 1,ntimstep
    call timestep_sympl_z(norb, z, ierr)
  end do
  endtime = omp_get_wtime()

  print *, z
  write(*,*) 'Time elapsed: ', endtime - starttime,' s'

  ! ncut = 1000
  ! call init_cut(cutter, norb, z)

  ! starttime = omp_get_wtime()
  ! do i = 1, ncut
  !   call trace_to_cut(cutter, z, var_tip, cut_type, ierr)
  ! end do
  ! endtime = omp_get_wtime()
  ! write(*,*) 'Time elapsed: ', endtime - starttime,' s'
  ! print *, z

end program main
