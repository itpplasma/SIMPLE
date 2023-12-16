!
  implicit none
  !
    integer, parameter :: ns=50,nperp=100
    integer :: npart,i,ipart,j,k
    double precision :: dummy,pmax,pmin,hs,hp
    double precision, dimension(nperp,ns,2)     :: promp,regul,stoch
    integer,          dimension(:), allocatable :: icl_jpa,icl_rec
    double precision, dimension(:), allocatable :: s,perp_inv
  !
    npart=0
    open(1,file='class_parts.dat')
    do
      read (1,*,end=1) j
      npart=npart+1
    enddo
  1 close(1)
    print *,'npart = ',npart
  !
    allocate(icl_jpa(npart),icl_rec(npart),s(npart),perp_inv(npart))
  !
    open(1,file='class_parts.dat')
    do ipart=1,npart
      read (1,*) j,s(ipart),perp_inv(ipart),icl_jpa(ipart),icl_rec(ipart)
    enddo
    close(1)
  !
    hs=1.d0/dble(ns)
  !
    pmin=0.d0
    pmax=maxval(perp_inv)
    promp=0.d0
    regul=0.d0
    stoch=0.d0
  !
    hp=(pmax-pmin)/dble(nperp)
  !
    do ipart=1,npart
      i=ceiling(s(ipart)/hs)
      i=min(ns,max(1,i))
  !
      k=ceiling(perp_inv(ipart)/hp)
      k=min(nperp,max(1,k))
  !
      select case(icl_jpa(ipart))
      case(0)
        promp(k,i,1)=promp(k,i,1)+1.d0
      case(1)
        regul(k,i,1)=regul(k,i,1)+1.d0
      case(2)
        stoch(k,i,1)=stoch(k,i,1)+1.d0
      end select
  !
      select case(icl_rec(ipart))
      case(0)
        promp(k,i,2)=promp(k,i,2)+1.d0
      case(1)
        regul(k,i,2)=regul(k,i,2)+1.d0
      case(2)
        stoch(k,i,2)=stoch(k,i,2)+1.d0
      end select
    enddo
  !
    open(1,file='prompt1.dat')
    do i=1,ns
      write(1,*) promp(:,i,1)
    enddo
    close(1)
  !
    open(1,file='prompt2.dat')
    do i=1,ns
      write(1,*) promp(:,i,2)
    enddo
    close(1)
  !
    open(1,file='regular1.dat')
    do i=1,ns
      write(1,*) regul(:,i,1)
    enddo
    close(1)
  !
    open(1,file='regular2.dat')
    do i=1,ns
      write(1,*) regul(:,i,2)
    enddo
    close(1)
  !
    open(1,file='stochastic1.dat')
    do i=1,ns
      write(1,*) stoch(:,i,1)
    enddo
    close(1)
  !
    open(1,file='stochastic2.dat')
    do i=1,ns
      write(1,*) stoch(:,i,2)
    enddo
    close(1)
  !
    end
