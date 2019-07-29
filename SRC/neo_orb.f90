module neo_orb
  use common, only: c, e_charge, p_mass, ev, twopi
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp

  use parmot_mod, only : rmu, ro0
  use velo_mod,   only : isw_field_type
  use orbit_symplectic, only : orbit_sympl_init, orbit_timestep_sympl
  use diag_mod, only : icounter

  implicit none

  interface tstep
      module procedure timestep
      module procedure timestep_z
      module procedure timestep_global
   end interface tstep

  double precision :: fper  ! field period
  double precision :: dtau, dtaumax, v0
  double precision, dimension(5) :: z
  integer          :: n_e, n_d

  integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet
  double precision :: relerr

  logical :: firstrun = .True.
  logical :: debug = .False.

contains

  subroutine init_field(ans_s, ans_tp, amultharm, aintegmode)
    ! initialize field geometry
    ! character*32, intent(in) :: vmec_file
    integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
    integer             :: ierr
    integer             :: L1i
    double precision    :: RT0, R0i, cbfi, bz0i, bf0

    netcdffile = 'wout.nc'  ! TODO: don't hardcode this
    ns_s = ans_s
    ns_tp = ans_tp
    multharm = amultharm
    integmode = aintegmode

    call spline_vmec_data ! initialize splines for VMEC field
    call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius
    fper = twopi/dble(L1i)   !<= field period
    print *, 'R0 = ', RT0, ' cm'
    isw_field_type = 1 ! evaluate fields in VMEC coords (0 = CAN, 1 = VMEC)
    if (integmode>=0) then
      call get_canonical_coordinates ! pre-compute transformation to canonical coords
      isw_field_type = 0 ! evaluate fields in canonical coords (0 = CAN, 1 = VMEC)
    end if

    ! initialize position and do first check if z is inside vacuum chamber
    z = 0.0d0
    call chamb_can(z(1:2), z(3), ierr)
    if(ierr.ne.0) stop
    z = 1.0d0
  end subroutine init_field

  subroutine init_params(Z_charge, m_mass, E_kin, adtau, adtaumax, arelerr)
    ! Initializes normalization for velocity and Larmor radius based on kinetic energy
    ! of plasma particles (= temperature for thermal particles).

    integer, intent(in) :: Z_charge, m_mass
    real(8), intent(in) :: E_kin, adtau, adtaumax
    double precision :: bmod_ref=5d4 ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    double precision :: bmod00, rlarm ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    double precision :: pi=3.14159265359 ! added by johanna 14.05.2019 in order to correct orbit calculation (missing 2pi in vmec in phitor)
    real(8), intent(in) :: arelerr

    n_e = Z_charge
    n_d = m_mass
    relerr = arelerr

    ! Neglect relativistic effects by large inverse relativistic temperature
    rmu=1d8

    ! Reference velocity and normalized Larmor radius
    v0 = sqrt(2.d0*E_kin*ev/(n_d*p_mass))
    !ro0 = v0*n_d*p_mass*c/(n_e*e_charge) !commented by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    ! Larmor radius:
    rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref) ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    bmod00=281679.46317784750d0 ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    ro0=rlarm*bmod00  ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    !ro0=v0*n_d*p_mass*c*2*pi/(n_e*e_charge) ! added by johanna 14.05.2019 in order to correct orbit calculation (missing 2pi in vmec in phitor)
    dtau = adtau ! timestep where to get results
    dtaumax = adtaumax ! maximum timestep for adaptive integration

  end subroutine init_params

  subroutine timestep(s, th, ph, lam, ierr)
    real(8), intent(inout) :: s, th, ph, lam
    integer, intent(out) :: ierr

    z(1) = s
    z(2) = th
    z(3) = ph
    z(4) = 1d0
    z(5) = lam

    if (integmode <= 0) then
      call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
    else
      if (firstrun) then
        call orbit_sympl_init(z, dtau, dtaumax, relerr, integmode)
        firstrun = .False.
      endif
      call orbit_timestep_sympl(z, ierr)
    endif

    s = z(1)
    th = z(2)
    ph = z(3)
    lam = z(5)
  end subroutine timestep

  subroutine timestep_z(az, ierr)
    real(8), intent(inout) :: az(:)
    integer, intent(out) :: ierr

    call timestep(az(1), az(2), az(3), az(4), ierr)
  end subroutine timestep_z

  subroutine timestep_global(ierr)
    integer, intent(out) :: ierr

    call timestep(z(1), z(2), z(3), z(5), ierr)
  end subroutine timestep_global

end module neo_orb

module cut_detector
  use common, only: twopi
  use neo_orb, only: z, fper, debug, dtau, tstep

  implicit none

  ! for Poincare cuts
  integer                                       :: nplagr,nder,npl_half
  integer,          dimension(:),   allocatable :: ipoi
  double precision, dimension(:),   allocatable :: xp
  double precision, dimension(:,:), allocatable :: coef,orb_sten
  double precision, dimension(:,:), allocatable :: zpoipl_tip,zpoipl_per,dummy2d
  double precision :: alam_prev, par_inv, zerolam
  integer iper, itip, kper
  integer, parameter :: n_tip_vars = 6

contains

  subroutine init()
    integer :: i

    !---------------------------------------------------------------------------
    ! Prepare calculation of orbit tip by interpolation and buffer for Poincare plot:

    nplagr=6
    nder=0
    npl_half=nplagr/2

    ! End prepare calculation of orbit tip by interpolation
    !--------------------------------------------------------------------------


    allocate(ipoi(nplagr),coef(0:nder,nplagr),orb_sten(6,nplagr),xp(nplagr))
    do i=1,nplagr
      ipoi(i)=i
    enddo

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
    par_inv = 0.0d0
    zerolam = 0.0d0
    !
  end subroutine init

  subroutine trace_to_cut(t, var_cut, cut_type, ierr)
    double precision, intent(inout) :: t
    ! variables to evaluate at tip: z(1..5), par_inv
    double precision, dimension(:), intent(inout) :: var_cut
    integer, intent(out) :: cut_type
    integer, intent(out) :: ierr

    integer, parameter :: nstep_max = 1000000000
    integer :: i
    double precision :: phiper

    phiper = 0.0d0

    do i=1,nstep_max
      call tstep(ierr)
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

          var_cut=matmul(orb_sten(:,ipoi),coef(0,:))
          var_cut(2)=modulo(var_cut(2),twopi)
          var_cut(3)=modulo(var_cut(3),twopi)

          par_inv = par_inv - var_cut(6)
          cut_type = 0
          return
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

            var_cut=matmul(orb_sten(:,ipoi),coef(0,:))
            var_cut(2)=modulo(var_cut(2),twopi)
            var_cut(3)=modulo(var_cut(3),twopi)

            cut_type = 1
            return
          endif   
        endif

    ! End periodic boundary footprint detection and interpolation
    !-------------------------------------------------------------------------


    end do
  end subroutine trace_to_cut

  subroutine fract_dimension(ntr,rt,fraction)

    integer, parameter :: iunit=1003
    integer :: itr,ntr,ngrid,nrefine,irefine,kr,kt,nboxes
    double precision :: fraction,rmax,rmin,tmax,tmin,hr,ht
    double precision, dimension(2,ntr)              :: rt
    logical,          dimension(:,:),   allocatable :: free

    rmin=minval(rt(1,:))
    rmax=maxval(rt(1,:))
    tmin=minval(rt(2,:))
    tmax=maxval(rt(2,:))

    nrefine=int(log(dble(ntr))/log(4.d0))

    ngrid=1
    nrefine=nrefine+3       !<=add 3 for curiousity
    do irefine=1,nrefine
      ngrid=ngrid*2
      allocate(free(0:ngrid,0:ngrid))
      free=.true.
      hr=(rmax-rmin)/dble(ngrid)
      ht=(tmax-tmin)/dble(ngrid)
      nboxes=0
      do itr=1,ntr
        kr=int((rt(1,itr)-rmin)/hr)
        kr=min(ngrid-1,max(0,kr))
        kt=int((rt(2,itr)-tmin)/ht)
        kt=min(ngrid-1,max(0,kt))
        if(free(kr,kt)) then
          free(kr,kt)=.false.
          nboxes=nboxes+1
        endif
      enddo
      deallocate(free)
      if(debug) write(iunit,*) dble(irefine),dble(nboxes)/dble(ngrid**2)
      if(irefine.eq.nrefine-3) fraction=dble(nboxes)/dble(ngrid**2)
    enddo
    close(iunit)

  end subroutine fract_dimension

end module cut_detector
