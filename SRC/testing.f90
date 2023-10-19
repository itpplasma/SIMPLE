  subroutine testing
  use alpha_lifetime_sub, only : velo_axis
!
  implicit real(8) (a-h,o-z),integer(i-n)
!
  double precision, dimension(:), allocatable :: dummy,dummy1
  double precision, dimension(:,:), allocatable :: dummy2d
  double precision :: tau
  double precision, dimension(5) :: z,vz
!
  pi=4.d0*atan2(1.d0,1.d0)
!
  xrange=0.01d0
  yrange=0.01d0
!
  nx=100
  ny=100
!
  hx=xrange/dble(nx)
  hy=yrange/dble(ny)
!
  allocate(dummy2d(-ny:ny,5))
!
  z(4)=1.d0
  z(5)=0.1d0
!
  do ix=-nx,nx
    z(1)=hx*dble(ix)
    do iy=-ny,ny
      z(2)=hy*dble(iy)+epsilon(1.d0)
!
      call velo_axis(tau,z,vz)
!
      dummy2d(iy,:)=vz
    enddo
    do k=1,5
      write(200+k,*) dummy2d(:,k)
    enddo
  enddo
!
  stop
  end
