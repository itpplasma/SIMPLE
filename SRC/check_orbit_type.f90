!
  module detect_oneline_mod
  !$omp threadprivate(iper, igroup, iret, fprs) 
    logical :: prop=.true.
    integer, parameter :: nfp_dim=3,ipermin=10
    integer :: nfp_max,iper,igroup
    double precision :: pi,twopi
    integer,          dimension(:),   allocatable :: iret
    double precision, dimension(:,:), allocatable :: fprs
  end module detect_oneline_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine check_orbit_type(nturns,nfp,fpr_in,ideal,ijpar,ierr)
!
  use detect_oneline_mod, only : prop,nfp_dim,nfp_max,iper, &
                                 pi,twopi, &
                                 fprs,igroup, &
                                 ipermin,iret
!
  implicit none
!
  logical :: did_recur
!
  integer :: nturns,ideal,ijpar
!
  integer :: nfp,iunit,i,j,k,ngroup,info,ierr,noline
  integer :: ibeg,iend,npl,nturnm1
!
  double precision :: dummy,dt1,dt2,derfix,per,ds1,ds2
  double precision :: theta_1,theta_2,theta_test,errint
  double precision :: perpinv_beg,tol_perpinv
  double precision, dimension(nfp_dim) :: fpr_in
  integer, dimension(:), allocatable :: ipoi
!
  tol_perpinv=15.d0
!
  if(prop) then
    prop=.false.
    pi=atan(1.d0)*4.d0
    twopi=2.d0*pi
!
!
    nfp_max=30000
    allocate(fprs(nfp_dim,nfp_max))
    allocate(iret(nfp_max))
!
  endif
!
  ierr=0
!
  if(nfp.eq.0) then
    iper=0
    igroup=0
    ideal=0
    ijpar=0
  endif
print *,igroup
!
  if(igroup.ge.nturns) then
    ierr=2
    return
  endif
!
  nfp=nfp+1
  if(nfp.gt.nfp_max) then
    print *,'detect_oneline: nfp > nfp_max, nturns not reached'
    ierr=1
    return
  endif
  fprs(:,nfp)=fpr_in
!
  if(nfp.lt.3) return
!
  theta_1=fprs(2,1)
  theta_2=modulo(fprs(2,2)-theta_1+pi,twopi)+theta_1-pi
  theta_test=modulo(fprs(2,nfp)-theta_1+pi,twopi)+theta_1-pi
  dt1=theta_test-theta_1
  dt2=theta_test-theta_2
  did_recur=dt1*dt2.lt.0.d0
  if(did_recur) then
    if(iper.eq.0) then
      iper=nfp
      if(iper.le.ipermin) then
        print *,'detect_oneline: stochastic, iper < ipermin'
        ideal=2
        ijpar=2
        igroup=nturns
        return
      endif
      print *,'detect_oneline: return period',iper-1
      igroup=1
    else
      igroup=igroup+1
    endif
    iret(igroup)=nfp-1
  endif
!
  if(igroup.eq.nturns) then
!check all recurrencties:
!
    call check_recurs(igroup,iret(1:igroup),noline)
!
    ideal=noline+1
!
! check sequencies in all intervals:
    if(ideal.eq.1) then
      nturnm1=nturns-1
      allocate(ipoi(nturnm1))
!
      call sortin(fprs(2,iret(1:nturnm1)+1),ipoi,nturnm1)
!
      do k=1,iret(nturns)-iret(nturnm1)
        do j=2,nturnm1
          theta_1=fprs(2,iret(ipoi(j-1))+k+1)
          theta_2=modulo(fprs(2,iret(ipoi(j))+k+1)-theta_1+pi,twopi)+theta_1-pi
          if(theta_2.lt.theta_1) then
            ideal=2
            exit
          endif
        enddo
        if(ideal.eq.2) exit
      enddo
      deallocate(ipoi)
    endif
!
! check J_perp change
    ijpar=1
    perpinv_beg=fprs(3,2)
    do k=2,iret(nturns)
      if(abs(fprs(3,k)-perpinv_beg).gt.tol_perpinv) then
        ijpar=2
        exit
      endif
    enddo
  endif
!
  end subroutine check_orbit_type
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine check_recurs(n,nrec,noline)
!
  implicit none
!
  integer :: n,noline,k,kmax,mult,nfirst,nnext,maxmult
  integer, dimension(n) :: nrec
!
  noline=0
!
  maxmult=n/2
!
  do mult=1,maxmult
    kmax=n/mult
    if(kmax.lt.2) exit
    nfirst=nrec(mult)
    do k=2,kmax
      nnext=nrec(mult*k)-nrec(mult*(k-1))
      if(nnext.gt.nfirst.or.nnext.lt.nfirst-1) then
        noline=1
        return
      endif
    enddo
  enddo
!
  end subroutine check_recurs
!
