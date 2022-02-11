program test_collis
  use collis_alp
  use common, only: pi, twopi, c, e_charge, e_mass, p_mass, ev

  implicit none

  double precision, dimension(5) :: z
  double precision :: am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha
  double precision :: v0,dchichi,slowrate,dchichi_norm,slowrate_norm,dtau
  double precision :: n_d=4.0d0, n_e=2.0d0  ! Test particle mass and charge no.
  double precision :: trace_time = 1.0d-1   ! Total trace time in seconds

  integer :: i,j,k,ntimstep,ierr

  namelist /config/ n_d, n_e, ealpha, trace_time, ntimstep
  namelist /collisions/ am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe

  open(1, file='collis.in')
  read(nml=config, unit=1)
  read(nml=collisions, unit=1)
  close(1)

  call loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha, &
                    v0,dchichi,slowrate,dchichi_norm,slowrate_norm)

  print *, 'v0 = ', v0, ' = ', sqrt(2.d0*ealpha*ev/(n_d*p_mass))

  dtau = trace_time*v0/dble(ntimstep-1)

  z = 0.0d0
  z(4) = 1.0d0
  k = 10000
  open(1, file='energyslow.dat')
  write (1, *) 0.d0, z(4)**2, z(5)
  do i = 2, ntimstep
    do j = 1, k
      call stost(z, dtau/dble(k), 1, ierr)
    enddo
    write (1,*) dble(i-1)*dtau/v0, z(4)**2, z(5)
  enddo
  close(1)
end program test_collis
