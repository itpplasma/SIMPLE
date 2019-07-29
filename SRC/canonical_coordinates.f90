module canonical_coordinates_main

  use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
  use parmot_mod, only : rmu,ro0
  use velo_mod,   only : isw_field_type
  use diag_mod, only : icounter
  use orbit_symplectic, only : orbit_sympl_init, orbit_timestep_sympl
  use common, only: twopi,c,e_charge,e_mass,p_mass,ev,newunit
  use cut_detector, only : fract_dimension

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  double precision,parameter  :: snear_axis=0.05d0

  logical          :: near_axis
  integer          :: npoi,ierr,L1i,nper,npoiper,ntimstep,ntestpart
  integer          :: ipart,notrace_passing,loopskip,iskip,ilost,it
  real             :: zzg
  double precision :: dphi,rbeg,phibeg,zbeg,bmod00,rcham,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: RT0,R0i,cbfi,bz0i,bf0,trap_par
  double precision :: sbeg,thetabeg
  double precision :: rbig,z1,z2
  double precision, dimension(5) :: z
  integer          :: npoiper2
  double precision :: contr_pp
  double precision :: facE_al
  integer          :: ibins
  integer          :: n_e,n_d,n_b
  integer, parameter :: npart = 960 !100 !984 !960
  double precision :: r,vartheta_c(npart),varphi_c(npart),theta_vmec,varphi_vmec,alam0(npart)

  integer :: i_ctr ! for nice counting in parallel
  integer, parameter :: mode_sympl = 1 ! 0 = RK canonical, 1 = Euler1, 2 = Euler2, 3 = Verlet

  double precision :: relerr ! relative error for RK

  !---------------------------------------------------------------------------
  ! buffer for Poincare plot:

  !---------------------------------------------------------------------------
  ! Prepare calculation of orbit tip by interpolation and buffer for Poincare plot:

  integer                                       :: nplagr,nder,npl_half
  integer                                       :: npass_regular(2),npass_chaotic(2),ntr_regular(2),ntr_chaotic(2)
  integer(8)                                    :: nstep_tot
  integer(4)                                    :: norbper,ifp_tip,ifp_per
  integer                                       :: nfp,nfp_tip,nfp_per
  double precision                              :: zerolam,fraction,fper
  integer, parameter                            :: n_tip_vars = 6  ! variables to evaluate at tip: z(1..5), par_inv
  double precision, dimension(n_tip_vars)       :: var_tip
  integer,          dimension(:),   allocatable :: ipoi
  double precision, dimension(:),   allocatable :: xp
  double precision, dimension(:,:), allocatable :: coef,orb_sten
  double precision, dimension(:,:), allocatable :: zpoipl_tip,zpoipl_per,dummy2d
  integer :: stat

  !$omp threadprivate(z,ifp_tip,ifp_per,ierr,orb_sten,xp,var_tip, &
  !$omp& zpoipl_tip, zpoipl_per, dummy2d, nfp, nfp_tip, nfp_per, ipoi, coef)

contains

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine trace_orbit(ipart, mode, orb_kind)

    implicit none

    integer, intent(in) :: ipart
    integer, intent(in) :: mode
    integer, intent(out) :: orb_kind

    double precision :: phiper, alam_prev, par_inv

    integer(8) :: i
    integer iper, itip, kper

    z(1)=r
    z(2)=vartheta_c(ipart)
    z(3)=varphi_c(ipart)
    z(4)=1.d0
    z(5)=alam0(ipart)

    allocate(zpoipl_tip(2,nfp_tip),zpoipl_per(2,nfp_per))

    if(mode == 0) then
        open(unit=10000+ipart, iostat=stat, status='old')
        if (stat == 0) close(10000+ipart, status='delete')
        open(unit=10000+ipart, recl=1024, position='append')

        open(unit=20000+ipart, iostat=stat, status='old')
        if (stat == 0) close(20000+ipart, status='delete')
        open(unit=20000+ipart, recl=1024, position='append')
    else
        open(unit=11000+ipart, iostat=stat, status='old')
        if (stat == 0) close(11000+ipart, status='delete')
        open(unit=11000+ipart, recl=1024, position='append')

        open(unit=21000+ipart, iostat=stat, status='old')
        if (stat == 0) close(21000+ipart, status='delete')
        open(unit=21000+ipart, recl=1024, position='append')
    endif

    ifp_tip=0               !<= initialize footprint counter on tips
    ifp_per=0               !<= initialize footprint counter on periods

    icounter=0
    phiper=0.0d0
    if (mode>0) call orbit_sympl_init(z, dtau, dtaumin, 1d-12, mode_sympl)

      !--------------------------------
      ! Initialize tip detector

      itip=npl_half+1
      alam_prev=z(5)

      ! End initialize tip detector
      !--------------------------------
      ! Initialize period crossing detector

      iper=npl_half+1
      kper=int(z(3)/fper)

      ! End initialize period crossing detector
      !--------------------------------

      par_inv = 0d0
      do i=1,nstep_tot

        if (mode==0) call orbit_timestep_axis(z, dtau, dtaumin,relerr,ierr)
        if (mode>0)  call orbit_timestep_sympl(z, ierr)

        if(ierr.ne.0) exit

        par_inv = par_inv+z(5)**2*dtau ! parallel adiabatic invariant

        if(i.le.nplagr) then          !<=first nplagr points to initialize stencil
          orb_sten(1:5,i)=z
          orb_sten(6,i)=par_inv
        else                          !<=normal case, shift stencil
          orb_sten(1:5,ipoi(1))=z
          orb_sten(6,ipoi(1))=par_inv
          ipoi=cshift(ipoi,1)
        endif

        !-------------------------------------------------------------------------
        ! Tip detection and interpolation

        if(alam_prev.lt.0.d0.and.z(5).gt.0.d0) itip=0   !<=tip has been passed
        itip=itip+1
        alam_prev=z(5)
        if(i.gt.nplagr) then          !<=use only initialized stencil
          if(itip.eq.npl_half) then   !<=stencil around tip is complete, interpolate
            xp=orb_sten(5,ipoi)

            call plag_coeff(nplagr,nder,zerolam,xp,coef)

            var_tip=matmul(orb_sten(:,ipoi),coef(0,:))
            var_tip(2)=modulo(var_tip(2),twopi)
            var_tip(3)=modulo(var_tip(3),twopi)

            if (mode==0) write(10000+ipart,*) var_tip
            if (mode>0) write(11000+ipart,*) var_tip

            ifp_tip=ifp_tip+1
            if(ifp_tip.gt.nfp_tip) then   !<=increase the buffer for banana tips
              allocate(dummy2d(2,ifp_tip-1))
              dummy2d=zpoipl_tip(:,1:ifp_tip-1)
              deallocate(zpoipl_tip)
              nfp_tip=nfp_tip+nfp
              allocate(zpoipl_tip(2,nfp_tip))
              zpoipl_tip(:,1:ifp_tip-1)=dummy2d
              deallocate(dummy2d)
            endif
            zpoipl_tip(:,ifp_tip)=var_tip(1:2)
            par_inv = par_inv - var_tip(6)
          endif
        endif

    ! End tip detection and interpolation
    !-------------------------------------------------------------------------
    ! Periodic boundary footprint detection and interpolation

        if(z(3).gt.dble(kper+1)*fper) then
          iper=0   !<=periodic boundary has been passed
          phiper=dble(kper+1)*fper
          kper=kper+1
        elseif(z(3).lt.dble(kper)*fper) then
          iper=0   !<=periodic boundary has been passed
          phiper=dble(kper)*fper
          kper=kper-1
        endif
        iper=iper+1
        if(i.gt.nplagr) then          !<=use only initialized stencil
          if(iper.eq.npl_half) then   !<=stencil around periodic boundary is complete, interpolate
            xp=orb_sten(3,ipoi)-phiper

            call plag_coeff(nplagr,nder,zerolam,xp,coef)

            var_tip=matmul(orb_sten(:,ipoi),coef(0,:))
            var_tip(2)=modulo(var_tip(2),twopi)
            var_tip(3)=modulo(var_tip(3),twopi)
    if (mode==0) write(20000+ipart,*) var_tip
    if (mode>0) write(21000+ipart,*) var_tip
            ifp_per=ifp_per+1
            if(ifp_per.gt.nfp_per) then   !<=increase the buffer for periodic boundary footprints
              allocate(dummy2d(2,ifp_per-1))
              dummy2d=zpoipl_per(:,1:ifp_per-1)
              deallocate(zpoipl_per)
              nfp_per=nfp_per+nfp
              allocate(zpoipl_per(2,nfp_per))
              zpoipl_per(:,1:ifp_per-1)=dummy2d
              deallocate(dummy2d)
            endif
            zpoipl_per(:,ifp_per)=var_tip(1:2)
          endif
        endif

    ! End periodic boundary footprint detection and interpolation
    !-------------------------------------------------------------------------

      enddo

      if(ifp_tip.eq.0) then

        call fract_dimension(ifp_per,zpoipl_per(:,1:ifp_per),fraction)

        if(fraction.gt.0.2d0) then
    !print *,'chaotic passing orbit'
          orb_kind = 2
        else
    !print *,'regular passing orbit'
          orb_kind = 1
        endif
      else

        call fract_dimension(ifp_tip,zpoipl_tip(:,1:ifp_tip),fraction)

        if(fraction.gt.0.2d0) then
    !print *,'chaotic trapped orbit'
          orb_kind = 4
        else
    !print *,'regular trapped orbit'
          orb_kind = 3
        endif
      endif

      deallocate(zpoipl_tip,zpoipl_per)
      if(mode == 0) close(10000+ipart)
      if(mode == 0) close(20000+ipart)
      if(mode > 0) close(11000+ipart)
      if(mode > 0) close(21000+ipart)
    end subroutine trace_orbit


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


end module canonical_coordinates_main

program canonical_coordinates
  use canonical_coordinates_main

  implicit none

  integer(8) :: calls_rk(npart), calls_sympl(npart)
  integer :: i, funit
  integer :: orb_kind_rk(npart), orb_kind_sympl(npart)


  ! run with fixed random seed
  integer :: seedsize
  integer, allocatable :: seed(:)

  call random_seed(size = seedsize)
  allocate(seed(seedsize))
  seed = 0
  call random_seed(put=seed)

  zerolam=0.d0
  nplagr=6
  nder=0
  npl_half=nplagr/2

  ! End prepare calculation of orbit tip by interpolation
  !--------------------------------------------------------------------------

  !

  open(1,file='alpha_lifetime_m.inp')
  read (1,*) notrace_passing   !skip tracing passing prts if notrace_passing=1
  read (1,*) nper              !number of periods for initial field line
  read (1,*) npoiper           !number of points per period on this field line
  read (1,*) ntimstep          !number of time steps per slowing down time
  read (1,*) ntestpart         !number of test particles
  read (1,*) bmod_ref          !reference field, G, for Boozer $B_{00}$
  read (1,*) trace_time        !slowing down time, s
  read (1,*) sbeg              !starting s for field line                       !<=2017
  read (1,*) phibeg            !starting phi for field line                     !<=2017
  read (1,*) thetabeg          !starting theta for field line                   !<=2017
  read (1,*) loopskip          !how many loops to skip to shift random numbers
  read (1,*) contr_pp          !control of passing particle fraction
  read (1,*) facE_al           !facE_al test particle energy reduction factor
  read (1,*) npoiper2          !additional integration step split factor
  read (1,*) n_e               !test particle charge number (the same as Z)
  read (1,*) n_d               !test particle mass number (the same as A)
  read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
  read (1,*) relerr            !relative error for RK integrator
  close(1)

  ! inverse relativistic temperature
  rmu=1d8

  ! alpha particle energy, eV:
  E_alpha=3.5d6/facE_al
  ! alpha particle velocity, cm/s
  v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))

  ! Larmor radius:
  rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)
  ! normalized slowing down time:
  tau=trace_time*v0

  bmod00=281679.46317784750d0
  ! Larmor radius corresponds to the field stregth egual to $B_{00}$ harmonic
  ! in Boozer coordinates:
  ! 14.11.2011  bmod00=bmod_ref  !<=deactivated, use value from the 'alpha_lifetime.inp'
  ro0=rlarm*bmod00  ! 23.09.2013

  multharm=3 !3 !7
  ns_A=5
  ns_s=5
  ns_tp=5

  call spline_vmec_data
  !call testing

  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)         !<=2017

  rbig=rt0
  ! field line integration step step over phi (to check chamber wall crossing)
  dphi=twopi/(L1i*npoiper)
  ! orbit integration time step (to check chamber wall crossing)
  dtaumin=twopi*rbig/npoiper2
  nstep_tot = ceiling(tau/dtaumin)
  dtau=tau/nstep_tot
  dtaumin = dtau

  print *, 'tau: ', tau, dtau, dtaumin, min(dabs(mod(tau, dtaumin)), dabs(mod(tau, dtaumin)-dtaumin))/dtaumin
  print *, 'total steps: ', nstep_tot

  fper=twopi/dble(L1i)   !<= field period

  call get_canonical_coordinates

  npass_regular=0
  npass_chaotic=0
  ntr_regular=0
  ntr_chaotic=0

  print *, 'generating random initial conditions'
  r=0.5d0

  do ipart=1,npart

    call random_number(zzg)

    vartheta_c(ipart)=twopi*zzg

    call random_number(zzg)

    varphi_c(ipart)=twopi*zzg

    call random_number(zzg)

    alam0(ipart)=2.d0*zzg-1.d0
  enddo

  calls_rk = 0
  calls_sympl = 0

  open(unit=newunit(funit), file='orbit_kinds.out', recl=1024)
  write(funit,*) '# ipart', ' r vartheta_c varphi_c p alam0', 'orb_kind_rk', &
      ' orb_kind_sympl', ' calls_rk', ' calls_sympl'

  !$omp parallel
  isw_field_type=0
  i_ctr=0

  norbper=ceiling(1d0*nstep_tot/(L1i*npoiper2))
  nfp=L1i*norbper         !<= guess for footprint number
  nfp_tip=nfp             !<= initial array dimension for tips
  nfp_per=nfp             !<= initial array dimension for periods

  print *, 'guessed orbit periods: ', norbper

  allocate(ipoi(nplagr),coef(0:nder,nplagr),orb_sten(6,nplagr),xp(nplagr))
  do i=1,nplagr
    ipoi(i)=i
  enddo
  !$omp do
  do ipart=1,npart
    icounter = 0
    call trace_orbit(ipart, 0, orb_kind_rk(ipart))
    calls_rk(ipart) = icounter
    icounter = 0
    call trace_orbit(ipart, 1, orb_kind_sympl(ipart))
    calls_sympl(ipart) = icounter
    if (orb_kind_rk(ipart) == 1) npass_regular(1) = npass_regular(1) + 1
    if (orb_kind_rk(ipart) == 2) npass_chaotic(1) = npass_chaotic(1) + 1
    if (orb_kind_rk(ipart) == 3) ntr_regular(1) = ntr_regular(1) + 1
    if (orb_kind_rk(ipart) == 4) ntr_chaotic(1) = ntr_chaotic(1) + 1

    if (orb_kind_sympl(ipart) == 1) npass_regular(2) = npass_regular(2) + 1
    if (orb_kind_sympl(ipart) == 2) npass_chaotic(2) = npass_chaotic(2) + 1
    if (orb_kind_sympl(ipart) == 3) ntr_regular(2) = ntr_regular(2) + 1
    if (orb_kind_sympl(ipart) == 4) ntr_chaotic(2) = ntr_chaotic(2) + 1
    !$omp critical
    i_ctr = i_ctr+1
    print *, 'ipart', ipart, '(',i_ctr,'/', npart,')',' RK calls:',calls_rk(ipart), &
      ' kind:', orb_kind_rk(ipart), ' Sympl calls:', calls_sympl(ipart), ' kind:', orb_kind_sympl(ipart)
    if (orb_kind_rk(ipart) /= orb_kind_sympl(ipart)) then
      print *, 'difference in classification - ipart', ipart, ': RK kind:',orb_kind_rk(ipart),&
      ', Sympl kind:', orb_kind_sympl(ipart), 'z0=', r, vartheta_c(ipart), varphi_c(ipart), 1.d0, alam0(ipart)
    !read(*,*)
    end if

    print *,npass_regular,' passing regular ',npass_chaotic,' passing chaotic ', &
          ntr_regular,' trapped regular ',ntr_chaotic,' trapped chaotic'

    write(funit,*) ipart, r, vartheta_c(ipart), varphi_c(ipart), 1.d0, &
    alam0(ipart), orb_kind_rk(ipart), orb_kind_sympl(ipart), calls_rk(ipart), calls_sympl(ipart)
    close(funit)
    open(funit,file='orbit_kinds.out', position='append', recl=1024)
    !$omp end critical
  enddo
  !$omp end do
  !$omp end parallel

  call deallocate_can_coord

end program canonical_coordinates
