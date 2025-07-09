program test_collis
  use collis_alp
  use util, only: p_mass, ev

  implicit none

  double precision, dimension(5) :: z
  double precision :: am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha
  double precision :: v0,dchichi,slowrate,dchichi_norm,slowrate_norm,dtau
  double precision :: n_d=4.0d0, n_e=2.0d0  ! Test particle mass and charge no.
  double precision :: trace_time = 1.0d-1   ! Total trace time in seconds

  integer, parameter :: ntimstep=100000 
  integer :: i,j,k,ierr,npart,ipart
  double precision, dimension(ntimstep) :: enoftim

!  namelist /config/ n_d, n_e, ealpha, trace_time, ntimstep
!  namelist /collisions/ am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe

!  open(1, file='collis.in')
!  read(nml=config, unit=1)
!  read(nml=collisions, unit=1)
!  close(1)

  ealpha=3.5d6
  trace_time=1.d0

  am1=2.0d0 
  am2=3.0d0 
  Z1=1.0d0 
  Z2=1.0d0
  densi1=0.5d14 
  densi2=0.5d14
  tempi1=1.0d4
  tempi2=1.0d4
  tempe=1.0d4

  call loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha, &
                    v0,dchichi,slowrate,dchichi_norm,slowrate_norm)

  print *, 'v0 = ', v0, ' = ', sqrt(2.d0*ealpha*ev/(n_d*p_mass))

  dtau = trace_time*v0/dble(ntimstep-1)

  print *, 'dtau = ', dtau

enoftim=0.d0
npart=10000
do ipart=1,npart
  z = 0.0d0
  z(4) = 1.0d0
  enoftim(1)=enoftim(1)+1.d0
  do i = 2, ntimstep
    call stost(z, dtau, 1, ierr)
!    if(ierr.ne.0) print *,ierr
    enoftim(i)=enoftim(i)+z(4)**2
  enddo
print *,ipart
enddo
enoftim=enoftim/dble(npart)
  open(1, file='energyslow_aver.dat')
do i = 1, ntimstep
    write (1,*) dble(i-1)*dtau/v0, enoftim(i)
enddo
  close(1)
end program test_collis
