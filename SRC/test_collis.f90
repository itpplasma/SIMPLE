program test_collis
  use collis_alp

  implicit none

  double precision, dimension(5) :: z
  double precision :: am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha
  double precision :: v0,dchichi,slowrate,dchichi_norm,slowrate_norm,dtau

  integer :: i,j,k,ntimstep,ierr

  am1 = 1.0d0
  am2 = 1.0d0
  Z1 = 1.0d0
  Z2 = 1.0d0
  densi1 = 1.0d20
  densi2 = 1.0d20
  tempi1 = 1.0d4
  tempi2 = 1.0d4
  tempi2 = 1.0d4
  ealpha = 3.5d6

  dtau = 1.0d0 ! TODO
  ntimstep = 3

  call loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha, &
                    v0,dchichi,slowrate,dchichi_norm,slowrate_norm)

  z=0.d0
  z(4)=1.d0
  k=10
  open(1,file='energyslow.dat')
  write (1,*) 0.d0,z(4)**2
  do i=2,ntimstep
    do j=1,k

      call stost(z,dtau/dble(k),3,ierr)

    enddo
    write (1,*) dble(i-1)*dtau/v0,z(4)**2
  enddo
  close(1)
end program test_collis
