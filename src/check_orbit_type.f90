!
  module detect_oneline_mod
    logical :: prop=.true.
    integer, parameter :: nfp_dim=3, ipermin=10
    integer :: nfp_max,iper,igroup
    double precision :: pi,twopi
    integer,          dimension(:),   allocatable :: iret
    double precision, dimension(:,:), allocatable :: fprs
    !> Continuous margins behind the integer classes. The J_parallel class
    !> thresholds jpar_spread against tol_perpinv, and the ideal-orbit class
    !> is the sign test on topology_margin. Both are the quantities the
    !> classifier already forms and then discards, kept so that a smooth
    !> confinement score can be built without re-tracing.
    !> score_status: 0 unresolved; 1 resolved with a valid topology margin;
    !> 2 early stochastic exit; 3 resolved but the ideal-orbit class came from
    !> the recurrence test, which forms no margin. jpar_spread is valid for
    !> status 1 and 3, topology_margin only for status 1.
    double precision :: jpar_spread, jpar_ref, topology_margin
    integer :: score_status
  !$omp threadprivate(prop, nfp_max, iper, igroup, iret, fprs, &
  !$omp               jpar_spread, jpar_ref, topology_margin, score_status)
  end module detect_oneline_mod
!
module check_orbit_type_sub

implicit none

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine check_orbit_type(nturns,nfp,fpr_in,ideal,ijpar,ierr)
!
  use detect_oneline_mod, only : prop,nfp_dim,nfp_max,iper, &
                                 pi,twopi, &
                                 fprs,igroup, &
                                 ipermin,iret, &
                                 jpar_spread,jpar_ref,topology_margin, &
                                 score_status
  use sorting_mod, only : sortin
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
  double precision :: perpinv_beg,tol_perpinv,drift,margin
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
    jpar_spread=0.d0
    jpar_ref=0.d0
    topology_margin=0.d0
    score_status=0
  endif
!
  if(igroup.ge.nturns) then
    ierr=2
    return
  endif
!
  nfp=nfp+1
  if(nfp.gt.nfp_max) then
#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
    print *,'detect_oneline: nfp > nfp_max, nturns not reached'
#endif
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
#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
        print *,'detect_oneline: stochastic, iper < ipermin'
#endif
        ideal=2
        ijpar=2
        igroup=nturns
        score_status=2
        return
      endif
#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
      print *,'detect_oneline: return period',iper-1
#endif
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
    margin=huge(1.d0)
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
          margin=min(margin,theta_2-theta_1)
          if(theta_2.lt.theta_1) ideal=2
        enddo
      enddo
      deallocate(ipoi)
    endif
    !
    ! Status 1 means topology_margin is a real margin. The recurrence test
    ! decides without forming one, and the monotonicity loop is empty when
    ! the last two return intervals coincide.
    if(margin.lt.huge(1.d0)) then
      topology_margin=margin
      score_status=1
    else
      score_status=3
    endif
!
! check J_perp change
    ijpar=1
    perpinv_beg=fprs(3,2)
    jpar_ref=abs(perpinv_beg)
    jpar_spread=0.d0
    do k=2,iret(nturns)
      drift=abs(fprs(3,k)-perpinv_beg)
      jpar_spread=max(jpar_spread,drift)
      if(drift.gt.tol_perpinv) ijpar=2
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
end module check_orbit_type_sub
