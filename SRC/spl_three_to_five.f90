!
  subroutine spl_five_reg(n,h,a,b,c,d,e,f)
!
  implicit none
!
  integer :: n,i,ip1,ip2
  double precision :: h,rhop,rhom,fac,fpl31,fpl40,fmn31,fmn40          ,x
  double precision :: a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,det
  double precision :: abeg,bbeg,cbeg,dbeg,ebeg,fbeg
  double precision :: aend,bend,cend,dend,eend,fend
  double precision, dimension(n) :: a,b,c,d,e,f
  double precision, dimension(:), allocatable :: alp,bet,gam
!
  rhop=13.d0+sqrt(105.d0)
  rhom=13.d0-sqrt(105.d0)
!
  a11=1.d0
  a12=1.d0/4.d0
  a13=1.d0/16.d0
  a21=3.d0
  a22=27.d0/4.d0
  a23=9.d0*27.d0/16.d0
  a31=5.d0
  a32=125.d0/4.d0
  a33=5.d0**5/16.d0
  det=a11*a22*a33+a12*a23*a31+a13*a21*a32-a12*a21*a33-a13*a22*a31-a11*a23*a32
  b1=a(4)-a(3)
  b2=a(5)-a(2)
  b3=a(6)-a(1)
  bbeg=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  bbeg=bbeg/det
  dbeg=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  dbeg=dbeg/det
  fbeg=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  fbeg=fbeg/det
  b1=a(n-2)-a(n-3)
  b2=a(n-1)-a(n-4)
  b3=a(n)-a(n-5)
  bend=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  bend=bend/det
  dend=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  dend=dend/det
  fend=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  fend=fend/det
  a11=2.d0
  a12=1.d0/2.d0
  a13=1.d0/8.d0
  a21=2.d0
  a22=9.d0/2.d0
  a23=81.d0/8.d0
  a31=2.d0
  a32=25.d0/2.d0
  a33=625.d0/8.d0
  det=a11*a22*a33+a12*a23*a31+a13*a21*a32-a12*a21*a33-a13*a22*a31-a11*a23*a32
  b1=a(4)+a(3)
  b2=a(5)+a(2)
  b3=a(6)+a(1)
  abeg=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  abeg=abeg/det
  cbeg=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  cbeg=cbeg/det
  ebeg=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  ebeg=ebeg/det
  b1=a(n-2)+a(n-3)
  b2=a(n-1)+a(n-4)
  b3=a(n)+a(n-5)
  aend=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  aend=aend/det
  cend=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  cend=cend/det
  eend=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  eend=eend/det
!
  allocate(alp(n),bet(n),gam(n))
!
  alp(1)=0.0d0
  bet(1)=ebeg*(2.d0+rhom)-5.d0*fbeg*(3.d0+1.5d0*rhom) !gamma1
!
  do i=1,n-4
    ip1=i+1
    alp(ip1)=-1.d0/(rhop+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)- &
             5.d0*(a(i+4)-4.d0*a(i+3)+6.d0*a(i+2)-4.d0*a(ip1)+a(i)))
  enddo
!
  gam(n-2)=eend*(2.d0+rhom)+5.d0*fend*(3.d0+1.5d0*rhom) !gamma
  do i=n-3,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  enddo
!
  alp(1)=0.0d0
  bet(1)=ebeg-2.5d0*5.d0*fbeg !e1
!
  do i=1,n-2
    ip1=i+1
    alp(ip1)=-1.d0/(rhom+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-gam(i))
  enddo
!
  e(n)=eend+2.5d0*5.d0*fend
  e(n-1)=e(n)*alp(n-1)+bet(n-1)
  f(n-1)=(e(n)-e(n-1))/5.d0
  e(n-2)=e(n-1)*alp(n-2)+bet(n-2)
  f(n-2)=(e(n-1)-e(n-2))/5.d0
  d(n-2)=dend+1.5d0*4.d0*eend+1.5d0**2*10.d0*fend
!
  do i=n-3,1,-1
    e(i)=e(i+1)*alp(i)+bet(i)
    f(i)=(e(i+1)-e(i))/5.d0
    d(i)=(a(i+3)-3.d0*a(i+2)+3.d0*a(i+1)-a(i))/6.d0 &
        -(e(i+3)+27.d0*e(i+2)+93.d0*e(i+1)+59.d0*e(i))/30.d0
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.5d0*d(i+1)-2.5d0*d(i) &
        -0.1d0*(e(i+2)+18.d0*e(i+1)+31.d0*e(i))
    b(i)=a(i+1)-a(i)-c(i)-d(i)-0.2d0*(4.d0*e(i)+e(i+1))
  enddo
!
  do i=n-3,n
    b(i)=b(i-1)+2.d0*c(i-1)+3.d0*d(i-1)+4.d0*e(i-1)+5.d0*f(i-1)
    c(i)=c(i-1)+3.d0*d(i-1)+6.d0*e(i-1)+10.d0*f(i-1)
    d(i)=d(i-1)+4.d0*e(i-1)+10.d0*f(i-1)
    if(i.ne.n) f(i)= a(i+1)-a(i)-b(i)-c(i)-d(i)-e(i)
  enddo
  f(n)=f(n-1)
!
  fac=1.d0/h
  b=b*fac
  fac=fac/h
  c=c*fac
  fac=fac/h
  d=d*fac
  fac=fac/h
  e=e*fac
  fac=fac/h
  f=f*fac
!
  deallocate(alp,bet,gam)
!
  return
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spl_five_per(n,h,a,b,c,d,e,f)
!
  implicit none
!
  integer :: n,i,ip1,ip2
  double precision :: h,rhop,rhom,fac,xplu,xmin,gammao_m,gammao_p,gammao_m_redef
  double precision :: dummy
  double precision, dimension(n) :: a,b,c,d,e,f
  double precision, dimension(:), allocatable :: alp,bet,gam
!
  rhop=13.d0+sqrt(105.d0)
  rhom=13.d0-sqrt(105.d0)
!
  allocate(alp(n),bet(n),gam(n))
!
  alp(1)=0.0d0
  bet(1)=0.0d0
!
  do i=1,n-4
    ip1=i+1
    alp(ip1)=-1.d0/(rhop+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)- &
             5.d0*(a(i+4)-4.d0*a(i+3)+6.d0*a(i+2)-4.d0*a(ip1)+a(i)))
  enddo
  alp(n-2)=-1.d0/(rhop+alp(n-3))
  bet(n-2)=alp(n-2)*(bet(n-3)- &
           5.d0*(a(2)-4.d0*a(1)+6.d0*a(n-1)-4.d0*a(n-2)+a(n-3)))
  alp(n-1)=-1.d0/(rhop+alp(n-2))
  bet(n-1)=alp(n-1)*(bet(n-2)- &
           5.d0*(a(3)-4.d0*a(2)+6.d0*a(1)-4.d0*a(n-1)+a(n-2)))
  alp(n)=-1.d0/(rhop+alp(n-1))
  bet(n)=alp(n)*(bet(n-1)- &
           5.d0*(a(4)-4.d0*a(3)+6.d0*a(2)-4.d0*a(1)+a(n-1)))
!
  gam(n)=bet(n)
  do i=n-1,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  enddo
!
  xplu=sqrt(0.25d0*rhop**2-1.d0)-0.5d0*rhop
  xmin=-sqrt(0.25d0*rhop**2-1.d0)-0.5d0*rhop
!  gammao_m=(gam(2)+xplu*gam(n))/(xmin**(n-1)-1.d0)/(xmin-xplu)
  dummy=(1.d0/xmin)**(n-1)
  gammao_m_redef=(gam(2)+xplu*gam(n))/(1.d0-dummy)/(xmin-xplu)
  gammao_p=(gam(2)+xmin*gam(n))/(xplu**(n-1)-1.d0)/(xplu-xmin)
!  gam(1)=gam(1)+gammao_m+gammao_p
  gam(1)=gam(1)+gammao_m_redef*dummy+gammao_p
  do i=2,n
!    gam(i)=gam(i)+gammao_m*xmin**(i-1)+gammao_p*xplu**(i-1)
    gam(i)=gam(i)+gammao_m_redef*(1.d0/xmin)**(n-i)+gammao_p*xplu**(i-1)
  enddo
!
  alp(1)=0.0d0
  bet(1)=0.d0
!
  do i=1,n-1
    ip1=i+1
    alp(ip1)=-1.d0/(rhom+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-gam(i))
  enddo
!
  e(n)=bet(n)
  do i=n-1,1,-1
    e(i)=e(i+1)*alp(i)+bet(i)
  enddo
!
  xplu=sqrt(0.25d0*rhom**2-1.d0)-0.5d0*rhom
  xmin=-sqrt(0.25d0*rhom**2-1.d0)-0.5d0*rhom
!  gammao_m=(e(2)+xplu*e(n))/(xmin**(n-1)-1.d0)/(xmin-xplu)
  dummy=(1.d0/xmin)**(n-1)
  gammao_m_redef=(e(2)+xplu*e(n))/(1.d0-dummy)/(xmin-xplu)
  gammao_p=(e(2)+xmin*e(n))/(xplu**(n-1)-1.d0)/(xplu-xmin)
!  e(1)=e(1)+gammao_m+gammao_p
  e(1)=e(1)+gammao_m_redef*dummy+gammao_p
  do i=2,n
!    e(i)=e(i)+gammao_m*xmin**(i-1)+gammao_p*xplu**(i-1)
    e(i)=e(i)+gammao_m_redef*(1.d0/xmin)**(n-i)+gammao_p*xplu**(i-1)
  enddo
!
  do i=n-1,1,-1
    f(i)=(e(i+1)-e(i))/5.d0
  enddo
  f(n)=f(1)
!
  d(n-1)=(a(3)-3.d0*a(2)+3.d0*a(1)-a(n-1))/6.d0 &
      -(e(3)+27.d0*e(2)+93.d0*e(1)+59.d0*e(n-1))/30.d0
  d(n-2)=(a(2)-3.d0*a(1)+3.d0*a(n-1)-a(n-2))/6.d0 &
      -(e(2)+27.d0*e(1)+93.d0*e(n-1)+59.d0*e(n-2))/30.d0
  do i=n-3,1,-1
    d(i)=(a(i+3)-3.d0*a(i+2)+3.d0*a(i+1)-a(i))/6.d0 &
        -(e(i+3)+27.d0*e(i+2)+93.d0*e(i+1)+59.d0*e(i))/30.d0
  enddo
  d(n)=d(1)
  c(n-1)=0.5d0*(a(2)+a(n-1))-a(1)-0.5d0*d(1)-2.5d0*d(n-1) &
      -0.1d0*(e(2)+18.d0*e(1)+31.d0*e(n-1))
  b(n-1)=a(1)-a(n-1)-c(n-1)-d(n-1)-0.2d0*(4.d0*e(n-1)+e(1))
!
  do i=n-2,1,-1
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.5d0*d(i+1)-2.5d0*d(i) &
        -0.1d0*(e(i+2)+18.d0*e(i+1)+31.d0*e(i))
    b(i)=a(i+1)-a(i)-c(i)-d(i)-0.2d0*(4.d0*e(i)+e(i+1))
  enddo
  b(n)=b(1)
  c(n)=c(1)
!
  fac=1.d0/h
  b=b*fac
  fac=fac/h
  c=c*fac
  fac=fac/h
  d=d*fac
  fac=fac/h
  e=e*fac
  fac=fac/h
  f=f*fac
!
  deallocate(alp,bet,gam)
!
  return
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spl_four_reg(n,h,a,b,c,d,e)
!
  implicit none
!
  integer :: n,i,ip1,ip2
  double precision :: h,fac,fpl31,fpl40,fmn31,fmn40
  double precision, dimension(n) :: a,b,c,d,e
  double precision, dimension(:), allocatable :: alp,bet,gam
!
  allocate(alp(n),bet(n),gam(n))
!
  fpl31=.5d0*(a(2)+a(4))-a(3)
  fpl40=.5d0*(a(1)+a(5))-a(3)
  fmn31=.5d0*(a(4)-a(2))
  fmn40=.5d0*(a(5)-a(1))
  d(3)=(fmn40-2.d0*fmn31)/6.d0
  e(3)=(fpl40-4.d0*fpl31)/12.d0
  d(2)=d(3)-4.d0*e(3)
  d(1)=d(3)-8.d0*e(3)
!
  alp(1)=0.0d0
  bet(1)=d(1)+d(2)
!
  do i=1,n-3
    ip1=i+1
    alp(ip1)=-1.d0/(10.d0+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-4.d0*(a(i+3)-3.d0*(a(i+2)-a(ip1))-a(i)))
  enddo
!
  fpl31=.5d0*(a(n-3)+a(n-1))-a(n-2)
  fpl40=.5d0*(a(n-4)+a(n))-a(n-2)
  fmn31=.5d0*(a(n-1)-a(n-3))
  fmn40=.5d0*(a(n)-a(n-4))
  d(n-2)=(fmn40-2.d0*fmn31)/6.d0
  e(n-2)=(fpl40-4.d0*fpl31)/12.d0
  d(n-1)=d(n-2)+4.d0*e(n-2)
  d(n)=d(n-2)+8.d0*e(n-2)
!
  gam(n-1)=d(n)+d(n-1)
!
  do i=n-2,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
    d(i)=gam(i)-d(i+1)
    e(i)=(d(i+1)-d(i))/4.d0
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.125d0*(d(i+2)+12.d0*d(i+1)+11.d0*d(i))
    b(i)=a(i+1)-a(i)-c(i)-(3.d0*d(i)+d(i+1))/4.d0
  enddo
!
  b(n-1)=b(n-2)+2.d0*c(n-2)+3.d0*d(n-2)+4.d0*e(n-2)
  c(n-1)=c(n-2)+3.d0*d(n-2)+6.d0*e(n-2)
  e(n-1)= a(n)-a(n-1)-b(n-1)-c(n-1)-d(n-1)
  b(n)=b(n-1)+2.d0*c(n-1)+3.d0*d(n-1)+4.d0*e(n-1)
  c(n)=c(n-1)+3.d0*d(n-1)+6.d0*e(n-1)
  e(n)=e(n-1) 
!
  fac=1.d0/h
  b=b*fac
  fac=fac/h
  c=c*fac
  fac=fac/h
  d=d*fac
  fac=fac/h
  e=e*fac
!
  deallocate(alp,bet,gam)
!
  return
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spl_four_per(n,h,a,b,c,d,e)
!
  implicit none
!
  integer :: n,i,ip1,ip2
  double precision :: h,fac,base1,base2,phi1,phi2,phi
  double precision, dimension(n) :: a,b,c,d,e
  double precision, dimension(:), allocatable :: alp,bet,gam
!
  allocate(alp(n),bet(n),gam(n))
!
  base1=-5.d0+2.d0*sqrt(6.d0)
  base2=-5.d0-2.d0*sqrt(6.d0)
!
  alp(1)=0.0d0
  bet(1)=0.d0
!
  do i=1,n-3
    ip1=i+1
    alp(ip1)=-1.d0/(10.d0+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-4.d0*(a(i+3)-3.d0*(a(i+2)-a(ip1))-a(i)))
  enddo
  alp(n-1)=-1.d0/(10.d0+alp(n-2))
  bet(n-1)=alp(n-1)*(bet(n-2)-4.d0*(a(2)-3.d0*(a(n)-a(n-1))-a(n-2)))
  alp(n)=-1.d0/(10.d0+alp(n-1))
  bet(n)=alp(n)*(bet(n-1)-4.d0*(a(3)-3.d0*(a(2)-a(n))-a(n-1)))
!
  gam(n)=bet(n)
!
  do i=n-1,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  enddo
!
  phi1=(gam(n)*base2+gam(2))/(base2-base1)/(1.d0-base1**(n-1))        !first
  phi2=(gam(n)*base1+gam(2))/(base2-base1)/(1.d0-(1.d0/base2)**(n-1)) !last
!
  do i=n,1,-1
    gam(i)=gam(i)+phi2
    phi2=phi2/base2
  enddo
!
  do i=1,n
    gam(i)=gam(i)+phi1
    phi1=phi1*base1
  enddo
!
  d(n)=0.d0
  do i=n-1,1,-1
    d(i)=gam(i)-d(i+1)
  enddo
!
  phi=-.5d0*d(1)
  do i=1,n
    d(i)=d(i)+phi
    phi=-phi
  enddo
!
  e(n)=(d(2)-d(n))/4.d0
  c(n)=0.5d0*(a(3)+a(n))-a(2)-0.125d0*(d(3)+12.d0*d(2)+11.d0*d(n))
  b(n)=a(2)-a(n)-c(n)-(3.d0*d(n)+d(2))/4.d0
  e(n-1)=(d(1)-d(n-1))/4.d0
  c(n-1)=0.5d0*(a(2)+a(n-1))-a(1)-0.125d0*(d(2)+12.d0*d(1)+11.d0*d(n-1))
  b(n-1)=a(1)-a(n-1)-c(n-1)-(3.d0*d(n-1)+d(1))/4.d0
!
  do i=n-2,1,-1
    e(i)=(d(i+1)-d(i))/4.d0
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.125d0*(d(i+2)+12.d0*d(i+1)+11.d0*d(i))
    b(i)=a(i+1)-a(i)-c(i)-(3.d0*d(i)+d(i+1))/4.d0
  enddo
!
  fac=1.d0/h
  b=b*fac
  fac=fac/h
  c=c*fac
  fac=fac/h
  d=d*fac
  fac=fac/h
  e=e*fac
!
  deallocate(alp,bet,gam)
!
  return
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE splreg(n,h,y,bi,ci,di)

! Makes a cubic spline of function y(x)
!
! Input:  n                   number of values in y  
!         h                   step size in x (equidistant)
!         y(n)                y-values
! Output: bi(n),ci(n),di(n)   Spline parameters
  
  IMPLICIT NONE
  
  INTEGER,                     INTENT(in)  :: n
  REAL(kind=kind(1.d0)),               INTENT(in)  :: h
  REAL(kind=kind(1.d0)), DIMENSION(n), INTENT(in)  :: y
  REAL(kind=kind(1.d0)), DIMENSION(n), INTENT(out) :: bi, ci, di

  REAL(kind=kind(1.d0))                         :: ak1, ak2, am1, am2, c, e, c1
  REAL(kind=kind(1.d0)), DIMENSION(:), ALLOCATABLE :: al, bt 
  INTEGER                                  :: k, n2, i, i5
  
  ALLOCATE ( al(n), bt(n) )
  
  ak1 = 0.d0
  ak2 = 0.d0
  am1 = 0.d0
  am2 = 0.d0
  k = n-1 
  al(1) = ak1
  bt(1) = am1
  n2 = n-2
  c = -4.d0*h
  DO i = 1,n2
     e = -3.d0*((y(i+2)-y(i+1))-(y(i+1)-y(i)))/h
     c1 = c-al(i)*h
     al(i+1) = h/c1
     bt(i+1) = (h*bt(i)+e)/c1
  END DO
  ci(n) = (am2+ak2*bt(k))/(1.d0-al(k)*ak2)
  DO i = 1,k
     i5 = n-i 
     ci(i5) = al(i5)*ci(i5+1)+bt(i5)
  END DO
  n2 = n-1
  DO i = 1,n2
     bi(i) = (y(i+1)-y(i))/h-h*(ci(i+1)+2.d0*ci(i))/3.d0
     di(i) = (ci(i+1)-ci(i))/h/3.d0
  END DO
  bi(n)=0.d0
  di(n)=0.d0
  DEALLOCATE ( al, bt )
  
  RETURN
END SUBROUTINE splreg
!=====================================================
SUBROUTINE splper(n,h,y,bi,ci,di)

! Makes a cubic spline of periodic function y(x)
!
! Input:  n                   number of values in y
!         h                   step size in x (equidistant)
!         y(n)                y-values
! Output: bi(n),ci(n),di(n)   Spline parameters

  IMPLICIT NONE

  INTEGER,                     INTENT(in)  :: n
  REAL(kind=kind(1.d0)),               INTENT(in)  :: h
  REAL(kind=kind(1.d0)), DIMENSION(n), INTENT(in)  :: y
  REAL(kind=kind(1.d0)), DIMENSION(n), INTENT(out) :: bi, ci, di

  REAL(kind=kind(1.d0))                            :: psi, ss
  REAL(kind=kind(1.d0)), DIMENSION(:), ALLOCATABLE :: bmx, yl
  REAL(kind=kind(1.d0)), DIMENSION(:), ALLOCATABLE :: amx1, amx2, amx3
  INTEGER                                  :: nmx, n1, n2, i, i1

  ALLOCATE ( bmx(n), yl(n), amx1(n), amx2(n), amx3(n) )

  bmx(1) = 1.d30

  nmx=n-1
  n1=nmx-1
  n2=nmx-2
  psi=3.d0/h/h

  CALL spfper(n,amx1,amx2,amx3)

  bmx(nmx) = (y(nmx+1)-2.d0*y(nmx)+y(nmx-1))*psi
  bmx(1)   =(y(2)-y(1)-y(nmx+1)+y(nmx))*psi
  DO i = 3,nmx
     bmx(i-1) = (y(i)-2.d0*y(i-1)+y(i-2))*psi
  END DO
  yl(1) = bmx(1)/amx1(1)
  DO i = 2,n1
     i1 = i-1
     yl(i) = (bmx(i)-yl(i1)*amx2(i1))/amx1(i)
  END DO
  ss = 0.d0
  DO i = 1,n1
     ss = ss+yl(i)*amx3(i)
  END DO
  yl(nmx) = (bmx(nmx)-ss)/amx1(nmx)
  bmx(nmx) = yl(nmx)/amx1(nmx)
  bmx(n1) = (yl(n1)-amx2(n1)*bmx(nmx))/amx1(n1)
  DO i = n2,1,-1
     bmx(i) = (yl(i)-amx3(i)*bmx(nmx)-amx2(i)*bmx(i+1))/amx1(i)
  END DO
  DO i = 1,nmx
     ci(i) = bmx(i)
  END DO

  DO i = 1,n1
     bi(i) = (y(i+1)-y(i))/h-h*(ci(i+1)+2.d0*ci(i))/3.d0
     di(i) = (ci(i+1)-ci(i))/h/3.d0
  END DO
  bi(nmx) = (y(n)-y(n-1))/h-h*(ci(1)+2.d0*ci(nmx))/3.d0
  di(nmx) = (ci(1)-ci(nmx))/h/3.d0
!
! Fix of problems at upper periodicity boundary
!
  bi(n) = bi(1)
  ci(n) = ci(1)
  di(n) = di(1)

  DEALLOCATE ( bmx, yl, amx1, amx2, amx3 )

  RETURN
END SUBROUTINE splper
!=====================================================
SUBROUTINE spfper(np1,amx1,amx2,amx3)

! Helper routine for splfi

  IMPLICIT NONE

  INTEGER,                       INTENT(in)  :: np1
  REAL(kind=kind(1.d0)), DIMENSION(np1), INTENT(out) :: amx1, amx2, amx3
  REAL(kind=kind(1.d0))                              :: beta, ss
  INTEGER                                    :: n, n1, i, i1

  n = np1-1

  n1 = n-1
  amx1(1) = 2.d0
  amx2(1) = 0.5d0
  amx3(1) = 0.5d0
  amx1(2) = SQRT(15.d0)/2.d0
  amx2(2) = 1.d0/amx1(2)
  amx3(2) = -.25d0/amx1(2)
  beta = 3.75d0
  DO i = 3,n1
     i1 = i-1
     beta = 4.d0-1.d0/beta
     amx1(i) = SQRT(beta)
     amx2(i) = 1.d0/amx1(i)
     amx3(i) = -amx3(i1)/amx1(i)/amx1(i1)
  END DO
  amx3(n1) = amx3(n1)+1.d0/amx1(n1)
  amx2(n1) = amx3(n1)
  ss = 0.0d0
  DO i = 1,n1
     ss = ss+amx3(i)*amx3(i)
  END DO
  amx1(n) = SQRT(4.d0-ss)

  RETURN
END SUBROUTINE spfper
!=====================================================
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spl_reg(ns,n,h,splcoe)
!
! Generic routine for equidistant regular spline of order 3,4,5
!
! Input/output parameters:
! ns              - spline order (input)
! n               - number of data points
! h               - abscissa step size
! splcoe(0:ns,n)  - spline coefficients; coefficient splcoe(0,:) are the input (function to spline)
!
  implicit none
!
  integer :: ns,n
  double precision :: h
  double precision, dimension(0:ns,n) :: splcoe
  double precision, dimension(:), allocatable :: a,b,c,d,e,f
!
  if(ns.eq.3) then
    allocate(a(n),b(n),c(n),d(n))
    a=splcoe(0,:)
!
    call splreg(n,h,a,b,c,d)
!
    splcoe(1,:)=b
    splcoe(2,:)=c
    splcoe(3,:)=d
    deallocate(a,b,c,d)
  elseif(ns.eq.4) then
    allocate(a(n),b(n),c(n),d(n),e(n))
    a=splcoe(0,:)
!
    call spl_four_reg(n,h,a,b,c,d,e)
!
    splcoe(1,:)=b
    splcoe(2,:)=c
    splcoe(3,:)=d
    splcoe(4,:)=e
    deallocate(a,b,c,d,e)
  elseif(ns.eq.5) then
    allocate(a(n),b(n),c(n),d(n),e(n),f(n))
    a=splcoe(0,:)
!
    call spl_five_reg(n,h,a,b,c,d,e,f)
!
    splcoe(1,:)=b
    splcoe(2,:)=c
    splcoe(3,:)=d
    splcoe(4,:)=e
    splcoe(5,:)=f
    deallocate(a,b,c,d,e,f)
  else
    print *,'wrong spline order'
  endif
!
  end subroutine spl_reg
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spl_per(ns,n,h,splcoe)
!
! Generic routine for equidistant periodic spline of order 3,4,5
!
! Input/output parameters:
! ns              - spline order (input)
! n               - number of data points
! h               - abscissa step size
! splcoe(0:ns,n)  - spline coefficients; coefficient splcoe(0,:) are the input (function to spline)
!
  implicit none
!
  integer :: ns,n
  double precision :: h
  double precision, dimension(0:ns,n) :: splcoe
  double precision, dimension(:), allocatable :: a,b,c,d,e,f
!
  if(ns.eq.3) then
    allocate(a(n),b(n),c(n),d(n))
    a=splcoe(0,:)
!
    call splper(n,h,a,b,c,d)
!
    splcoe(1,:)=b
    splcoe(2,:)=c
    splcoe(3,:)=d
    deallocate(a,b,c,d)
  elseif(ns.eq.4) then
    allocate(a(n),b(n),c(n),d(n),e(n))
    a=splcoe(0,:)
!
    call spl_four_per(n,h,a,b,c,d,e)
!
    splcoe(1,:)=b
    splcoe(2,:)=c
    splcoe(3,:)=d
    splcoe(4,:)=e
    deallocate(a,b,c,d,e)
  elseif(ns.eq.5) then
    allocate(a(n),b(n),c(n),d(n),e(n),f(n))
    a=splcoe(0,:)
!
    call spl_five_per(n,h,a,b,c,d,e,f)
!
    splcoe(1,:)=b
    splcoe(2,:)=c
    splcoe(3,:)=d
    splcoe(4,:)=e
    splcoe(5,:)=f
    deallocate(a,b,c,d,e,f)
  else
    print *,'wrong spline order'
  endif
!
  end subroutine spl_per
