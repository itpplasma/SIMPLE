module neo_orb
  use common, only: c, e_charge, p_mass, ev, twopi
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp

  use parmot_mod, only : rmu, ro0
  use velo_mod,   only : isw_field_type
  use field_can_mod, only : FieldCan
  use orbit_symplectic, only : SymplecticIntegrator, MultistageIntegrator, &
    orbit_sympl_init, orbit_timestep_sympl
  use field_can_mod, only : FieldCan, eval_field
  use diag_mod, only : icounter

implicit none
save

  logical :: debug = .False.

public

  type :: NeoOrb
    double precision :: fper  ! field period
    double precision :: dtau, dtaumin, v0
    integer          :: n_e, n_d

    integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet
    double precision :: relerr

    type(FieldCan) :: f
    type(SymplecticIntegrator) :: si
    type(MultistageIntegrator) :: mi
  end type NeoOrb

  interface tstep
      module procedure timestep
      module procedure timestep_z
   end interface tstep

contains

  subroutine init_field(self, ans_s, ans_tp, amultharm, aintegmode)
    ! initialize field geometry
    ! character*32, intent(in) :: vmec_file
    type(NeoOrb), intent(inout) :: self
    integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
    integer             :: ierr
    integer             :: L1i
    double precision    :: RT0, R0i, cbfi, bz0i, bf0
    double precision    :: z(5)

    netcdffile = 'wout.nc'  ! TODO: don't hardcode this
    ns_s = ans_s
    ns_tp = ans_tp
    multharm = amultharm
    self%integmode = aintegmode

    call spline_vmec_data ! initialize splines for VMEC field
    call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius
    self%fper = twopi/dble(L1i)   !<= field period
    print *, 'R0 = ', RT0, ' cm, fper = ', self%fper
    isw_field_type = 1 ! evaluate fields in VMEC coords (0 = CAN, 1 = VMEC)
    if (self%integmode>=0) then
      call get_canonical_coordinates ! pre-compute transformation to canonical coords
      isw_field_type = 0 ! evaluate fields in canonical coords (0 = CAN, 1 = VMEC)
    end if

    ! initialize position and do first check if z is inside vacuum chamber
    z = 0.0d0
    call chamb_can(z(1:2), z(3), ierr)
    if(ierr.ne.0) stop
    z = 1.0d0
  end subroutine init_field

  subroutine init_params(self, Z_charge, m_mass, E_kin, adtau, adtaumin, arelerr)
    ! Initializes normalization for velocity and Larmor radius based on kinetic energy
    ! of plasma particles (= temperature for thermal particles).

    type(NeoOrb) :: self
    integer, intent(in) :: Z_charge, m_mass
    real(8), intent(in) :: E_kin, adtau, adtaumin
    double precision :: bmod_ref=5d4 ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    double precision :: bmod00, rlarm ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    real(8), intent(in) :: arelerr

    self%n_e = Z_charge
    self%n_d = m_mass
    self%relerr = arelerr

    ! Neglect relativistic effects by large inverse relativistic temperature
    rmu=1d8

    ! Reference velocity and normalized Larmor radius
    self%v0 = sqrt(2.d0*E_kin*ev/(self%n_d*p_mass))
    !ro0 = v0*n_d*p_mass*c/(n_e*e_charge) !commented by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    ! Larmor radius:
    rlarm=self%v0*self%n_d*p_mass*c/(self%n_e*e_charge*bmod_ref) ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    bmod00=281679.46317784750d0 ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    ro0=rlarm*bmod00  ! added by johanna 10.05.2019 in order to correct orbit calculation (in analogy to test_orbits_vmec)
    !ro0=v0*n_d*p_mass*c*2*pi/(n_e*e_charge) ! added by johanna 14.05.2019 in order to correct orbit calculation (missing 2pi in vmec in phitor)
    self%dtau = adtau ! timestep where to get results
    self%dtaumin = adtaumin ! minimum timestep for adaptive integration

  end subroutine init_params

  subroutine init_sympl(si, f, z0, dtau, dtaumin, rtol_init, mode_init)
    !
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    double precision, intent(in) :: z0(:)
    double precision, intent(in) :: dtau, dtaumin
    double precision, intent(in) :: rtol_init
    integer, intent(in) :: mode_init ! 1 = expl.-impl. Euler, 2 = impl.-expl. Euler

    double precision :: z(4)

    if(min(dabs(mod(dtau, dtaumin)), dabs(mod(dtau, dtaumin)-dtaumin)) > 1d-9*dtaumin) then
      stop 'orbit_sympl_init - error: dtau/dtaumin not integer'
    endif

    ! Initialize symplectic integrator
    call eval_field(f, z0(1), z0(2), z0(3), 0)

    si%pabs = z0(4)
  
    f%mu = .5d0*z0(4)**2*(1.d0-z0(5)**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
    f%ro0 = ro0/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
    f%vpar = z0(4)*z0(5)*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules
      
    z(1:3) = z0(1:3)  ! s, th, ph
    z(4) = f%vpar*f%hph + f%Aph/f%ro0 ! pphi

    ! factor 1/sqrt(2) due to velocity normalisation different from other modules
    call orbit_sympl_init(si, f, z, dtaumin/dsqrt(2d0), nint(dtau/dtaumin), &
                          rtol_init, mode_init, 1) 
  end subroutine init_sympl

  subroutine init_integrator(self, z0)
    type(NeoOrb), intent(inout) :: self
    double precision, intent(in) :: z0(:)
    
    call init_sympl(self%si, self%f, z0, self%dtau, self%dtaumin, &
      self%relerr, self%integmode)
  end subroutine init_integrator

  subroutine timestep(self, s, th, ph, lam, ierr)
    type(NeoOrb), intent(inout) :: self
    real(8), intent(inout) :: s, th, ph, lam
    integer, intent(out) :: ierr

    double precision, dimension(5) :: z

    z(1) = s
    z(2) = th
    z(3) = ph
    z(4) = 1d0
    z(5) = lam

    call timestep_z(self, z, ierr)

    s = z(1)
    th = z(2)
    ph = z(3)
    lam = z(5)
  end subroutine timestep

  subroutine timestep_z(self, z, ierr)
    type(NeoOrb), intent(inout) :: self
    real(8), intent(inout) :: z(:)
    integer, intent(out) :: ierr

    call orbit_timestep_axis(z, self%dtau, self%dtau, self%relerr, ierr)
  end subroutine timestep_z

  subroutine timestep_sympl_z(norb, z, ierr)
    type(NeoOrb), intent(inout) :: norb
    real(8), intent(inout) :: z(:)
    integer, intent(out) :: ierr

    if (z(1) < 0.0 .or. z(1) > 1.0) then
      ierr = 1
      return
    end if

    call orbit_timestep_sympl(norb%si, norb%f, ierr)

    z(1:3) = norb%si%z(1:3)
    z(5) = norb%f%vpar/(norb%si%pabs*dsqrt(2d0))  ! alambda

  end subroutine timestep_sympl_z

end module neo_orb

module cut_detector
  use common, only: twopi
  use neo_orb, only: NeoOrb, debug, tstep

  implicit none
  save
  
  integer, parameter :: n_tip_vars = 6
  integer, parameter :: nplagr = 6
  integer, parameter :: nder = 0

  public

  type :: CutDetector
    ! for Poincare cuts
    integer          :: npl_half
    double precision :: alam_prev, par_inv
    integer          :: iper, itip, kper

    double precision :: orb_sten(6,nplagr), coef(0:nder,nplagr)
    integer :: ipoi(nplagr)

    type(NeoOrb), pointer :: norb
  end type CutDetector

contains

  subroutine init(self, norb, z)
    type(CutDetector) :: self
    type(NeoOrb), target :: norb
    double precision, intent(in) :: z(:)
    integer :: i

    self%norb => norb

    !---------------------------------------------------------------------------
    ! Prepare calculation of orbit tip by interpolation and buffer for Poincare plot:

    self%npl_half=nplagr/2

    do i=1,nplagr
      self%ipoi(i)=i
    enddo

    !--------------------------------
    ! Initialize tip detector

    self%itip=self%npl_half+1
    self%alam_prev=z(5)

    ! End initialize tip detector
    !--------------------------------
    ! Initialize period crossing detector

    self%iper=self%npl_half+1
    self%kper=int(z(3)/self%norb%fper)

    ! End initialize period crossing detector
    !--------------------------------
    self%par_inv = 0.0d0
    !

    ! End prepare calculation of orbit tip by interpolation
    !--------------------------------------------------------------------------
  end subroutine init

  subroutine trace_to_cut(self, z, var_cut, cut_type, ierr)
    type(CutDetector) :: self
    double precision, intent(inout) :: z(:)
    ! variables to evaluate at tip: z(1..5), par_inv
    double precision, dimension(:), intent(inout) :: var_cut
    integer, intent(out) :: cut_type
    integer, intent(out) :: ierr

    integer, parameter :: nstep_max = 1000000000
    integer :: i
    double precision :: phiper = 0.0d0

    do i=1, nstep_max
      call tstep(self%norb, z, ierr)
      if(ierr.ne.0) exit

      self%par_inv = self%par_inv+z(5)**2*self%norb%dtau ! parallel adiabatic invariant

      if(i.le.nplagr) then          !<=first nplagr points to initialize stencil
        self%orb_sten(1:5,i)=z
        self%orb_sten(6,i)=self%par_inv
      else                          !<=normal case, shift stencil
        self%orb_sten(1:5,self%ipoi(1))=z
        self%orb_sten(6,self%ipoi(1))=self%par_inv
        self%ipoi=cshift(self%ipoi,1)
      endif

      !-------------------------------------------------------------------------
      ! Tip detection and interpolation

      if(self%alam_prev.lt.0.d0.and.z(5).gt.0.d0) self%itip=0   !<=tip has been passed
      self%itip=self%itip+1
      self%alam_prev=z(5)

      !<=use only initialized stencil
      if(i.gt.nplagr .and. self%itip.eq.self%npl_half) then
        cut_type = 0
        call plag_coeff(nplagr, nder, 0d0, self%orb_sten(5, self%ipoi), self%coef)
        var_cut = matmul(self%orb_sten(:, self%ipoi), self%coef(0,:))
        var_cut(2) = modulo(var_cut(2), twopi)
        var_cut(3) = modulo(var_cut(3), twopi)
        self%par_inv = self%par_inv - var_cut(6)
        return
      end if

      ! End tip detection and interpolation

      !-------------------------------------------------------------------------
      ! Periodic boundary footprint detection and interpolation

      if(z(3).gt.dble(self%kper+1)*self%norb%fper) then
        self%iper=0   !<=periodic boundary has been passed
        phiper=dble(self%kper+1)*self%norb%fper
        self%kper=self%kper+1
      elseif(z(3).lt.dble(self%kper)*self%norb%fper) then
        self%iper=0   !<=periodic boundary has been passed
        phiper=dble(self%kper)*self%norb%fper
        self%kper=self%kper-1
      endif
      self%iper=self%iper+1

      !<=use only initialized stencil
      if(i.gt.nplagr .and. self%iper.eq.self%npl_half) then
        cut_type = 1
        !<=stencil around periodic boundary is complete, interpolate
        call plag_coeff(nplagr, nder, 0d0, self%orb_sten(3, self%ipoi) - phiper, self%coef)
        var_cut = matmul(self%orb_sten(:, self%ipoi), self%coef(0,:))
        var_cut(2) = modulo(var_cut(2), twopi)
        var_cut(3) = modulo(var_cut(3), twopi)
        return
      end if

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
