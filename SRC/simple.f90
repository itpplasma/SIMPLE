module simple
  use util, only: c, e_charge, p_mass, ev, twopi
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp, &
                                 vmec_B_scale, vmec_RZ_scale

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

  type :: Tracer
    double precision :: fper
    double precision :: dtau, dtaumin, v0
    integer          :: n_e, n_d

    integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet
    double precision :: relerr

    type(FieldCan) :: f
    type(SymplecticIntegrator) :: si
    type(MultistageIntegrator) :: mi
  end type Tracer

  interface tstep
      module procedure timestep
      module procedure timestep_z
      module procedure timestep_sympl_z
   end interface tstep

contains

  subroutine init_field(self, vmec_file, ans_s, ans_tp, amultharm, aintegmode)
    ! initialize field geometry
    character(len=*), intent(in) :: vmec_file
    type(Tracer), intent(inout) :: self
    integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
    integer             :: ierr
    integer             :: L1i
    double precision    :: RT0, R0i, cbfi, bz0i, bf0, volume, B00
    double precision    :: z(5)

    netcdffile = vmec_file
    ns_s = ans_s
    ns_tp = ans_tp
    multharm = amultharm
    self%integmode = aintegmode

    call spline_vmec_data ! initialize splines for VMEC field
    call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius
    self%fper = twopi/dble(L1i)   !<= field period
    print *, 'R0 = ', RT0, ' cm, fper = ', self%fper
    call volume_and_B00(volume,B00)
    print *,'volume = ',volume,' cm^3,  B_00 = ',B00,' G'

    if (self%integmode>=0) then
      if (isw_field_type == 0) then
        call get_canonical_coordinates ! pre-compute transformation to canonical coords
      elseif (isw_field_type == 2) then
        print *, 'Boozer field'
        call boozer_converter
      else
        print *, 'Unknown field type ', isw_field_type
      endif

    end if

    ! initialize position and do first check if z is inside vacuum chamber
    z = 0.0d0
    call chamb_can(z(1:2), z(3), ierr)
    if(ierr.ne.0) stop
    z = 1.0d0
  end subroutine init_field

  subroutine init_params(self, Z_charge, m_mass, E_kin, npoints, store_step, relerr)
    ! Initializes normalization for velocity and Larmor radius based on kinetic energy
    ! of plasma particles (= temperature for thermal particles).
    use new_vmec_stuff_mod, only : rmajor

    type(Tracer) :: self
    integer, intent(in), optional :: Z_charge, m_mass
    double precision, intent(in), optional :: E_kin
    integer, intent(in), optional :: npoints  ! Integrator resolution. Number of
    ! integration points for strongly passing particles around the torus

    integer, intent(in), optional :: store_step       ! Store every X timesteps
    double precision, intent(in), optional :: relerr  ! Relative error

    if (present(Z_charge)) self%n_e = Z_charge
    if (present(m_mass)) self%n_d = m_mass
    if (present(relerr)) self%relerr = relerr

    ! Neglect relativistic effects by large inverse relativistic temperature
    rmu=1d8

    ! Reference velocity and normalized Larmor radius
    if (present(E_kin)) then
      self%v0 = sqrt(2.d0*E_kin*ev/(self%n_d*p_mass))
    else
      self%v0 = sqrt(2.d0*3.5d6*ev/(self%n_d*p_mass))
    end if

    ! Larmor radius times magnetic field:
    ro0=self%v0*self%n_d*p_mass*c/(self%n_e*e_charge)

    if (present(npoints)) then
      self%dtaumin=twopi*rmajor*1.0d2/npoints
    else
      self%dtaumin=twopi*rmajor*1.0d2/256.0d0
    end if

   ! timestep where to get results
    if (present(store_step)) then
      self%dtau = store_step*self%dtaumin
    else
      self%dtau = self%dtaumin
    end if

  end subroutine init_params

  subroutine init_sympl(si, f, z0, dtau, dtaumin, rtol_init, mode_init)
    !
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    double precision, intent(in) :: z0(:)
    double precision, intent(in) :: dtau, dtaumin
    double precision, intent(in) :: rtol_init
    integer, intent(in) :: mode_init ! 1 = expl.-impl. Euler, 2 = impl.-expl. Euler, 3 = Midpoint

    double precision :: z(4)

    if(min(dabs(mod(dtau, dtaumin)), dabs(mod(dtau, dtaumin)-dtaumin)) > 1d-9*dtaumin) then
      stop 'orbit_sympl_init - error: dtau/dtaumin not integer'
    endif

    ! Initialize symplectic integrator
    f%field_type = isw_field_type
    call eval_field(f, z0(1), z0(2), z0(3), 0)

    si%pabs = z0(4)

    f%mu = .5d0*z0(4)**2*(1.d0-z0(5)**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
    f%ro0 = ro0/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
    f%vpar = z0(4)*z0(5)*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

    z(1:3) = z0(1:3)  ! s, th, ph
    z(4) = f%vpar*f%hph + f%Aph/f%ro0 ! pphi

    ! factor 1/sqrt(2) due to velocity normalisation different from other modules
    call orbit_sympl_init(si, f, z, dtaumin/dsqrt(2d0), nint(dtau/dtaumin), &
                          rtol_init, mode_init, 0)
  end subroutine init_sympl

  subroutine init_integrator(self, z0)
    type(Tracer), intent(inout) :: self
    double precision, intent(in) :: z0(:)

    call init_sympl(self%si, self%f, z0, self%dtau, self%dtaumin, &
      self%relerr, self%integmode)
  end subroutine init_integrator

  subroutine timestep(self, s, th, ph, lam, ierr)
    type(Tracer), intent(inout) :: self
    double precision, intent(inout) :: s, th, ph, lam
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
    type(Tracer), intent(inout) :: self
    double precision, intent(inout) :: z(:)
    integer, intent(out) :: ierr

    call orbit_timestep_axis(z, self%dtau, self%dtau, self%relerr, ierr)
  end subroutine timestep_z

  subroutine timestep_sympl_z(si, f, z, ierr)
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    double precision, intent(inout) :: z(:)
    integer, intent(out) :: ierr

    if (z(1) < 0.0 .or. z(1) > 1.0) then
      ierr = 1
      return
    end if

    call orbit_timestep_sympl(si, f, ierr)

    z(1:3) = si%z(1:3)
    z(5) = f%vpar/(si%pabs*dsqrt(2d0))  ! alambda

  end subroutine timestep_sympl_z

end module simple

module cut_detector
  use util, only: twopi
  use simple, only: debug, tstep
  use orbit_symplectic, only: SymplecticIntegrator
  use field_can_mod, only: FieldCan

  implicit none
  save

  integer, parameter :: n_tip_vars = 6
  integer, parameter :: nplagr = 6
  integer, parameter :: nder = 0

  public

  type :: CutDetector
    double precision :: fper  ! field period

    ! for Poincare cuts
    double precision :: alam_prev, par_inv
    integer          :: iper, itip, kper

    double precision :: orb_sten(6,nplagr), coef(0:nder,nplagr)
    integer :: ipoi(nplagr)
  end type CutDetector

contains

  subroutine init(self, fper, z)
    type(CutDetector) :: self
    double precision, intent(in) :: fper
    double precision, intent(in) :: z(:)
    integer :: i

    self%fper = fper

    !---------------------------------------------------------------------------
    ! Prepare calculation of orbit tip by interpolation and buffer for Poincare plot:

    do i=1,nplagr
      self%ipoi(i)=i
    enddo

    !--------------------------------
    ! Initialize tip detector

    self%itip=nplagr/2+1
    self%alam_prev=z(5)

    ! End initialize tip detector
    !--------------------------------
    ! Initialize period crossing detector

    self%iper=nplagr/2+1
    self%kper=int(z(3)/self%fper)

    ! End initialize period crossing detector
    !--------------------------------
    self%par_inv = 0.0d0
    !

    ! End prepare calculation of orbit tip by interpolation
    !--------------------------------------------------------------------------
  end subroutine init

  subroutine trace_to_cut(self, si, f, z, var_cut, cut_type, ierr)
    type(CutDetector) :: self
    type(SymplecticIntegrator) :: si
    type(FieldCan) :: f

    double precision, intent(inout) :: z(:)
    ! variables to evaluate at tip: z(1..5), par_inv
    double precision, dimension(:), intent(inout) :: var_cut
    integer, intent(out) :: cut_type
    integer, intent(out) :: ierr

    integer, parameter :: nstep_max = 1000000000
    integer :: i
    double precision :: phiper = 0.0d0

    do i=1, nstep_max
      call tstep(si, f, z, ierr)
      if(ierr.ne.0) exit

      self%par_inv = self%par_inv+z(5)**2 ! parallel adiabatic invariant

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
      if(i.gt.nplagr .and. self%itip.eq.nplagr/2) then
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

      if(z(3).gt.dble(self%kper+1)*self%fper) then
        self%iper=0   !<=periodic boundary has been passed
        phiper=dble(self%kper+1)*self%fper
        self%kper=self%kper+1
      elseif(z(3).lt.dble(self%kper)*self%fper) then
        self%iper=0   !<=periodic boundary has been passed
        phiper=dble(self%kper)*self%fper
        self%kper=self%kper-1
      endif
      self%iper=self%iper+1

      !<=use only initialized stencil
      if(i.gt.nplagr .and. self%iper.eq.nplagr/2) then
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
    double precision, dimension(2,ntr)              :: rt!0, rt
    logical,          dimension(:,:),   allocatable :: free

! TODO: check if this works better on tips only
!    rt(1,:) = rt0(1,:)*cos(rt0(2,:))
!    rt(2,:) = rt0(1,:)*sin(rt0(2,:))

    rmin=minval(rt(1,:))
    rmax=maxval(rt(1,:))
    tmin=minval(rt(2,:))
    tmax=maxval(rt(2,:))

    nrefine=int(log(dble(ntr))/log(4.d0))

    ngrid=1
    nrefine=nrefine+3       !<=add 3 for curiousity
    do irefine=1,nrefine
      ngrid=ngrid*2
!$omp critical
      allocate(free(0:ngrid,0:ngrid))
!$omp end critical
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
!$omp critical
      deallocate(free)
!$omp end critical
      if(debug) then
!$omp critical
        !
        ! Right now criterion for regular is at nboxes/ngrid**2 < 0.2 .
        ! For fractal dimension d = log(nboxes)/log(ngrid) this means
        ! d_thresh = 2 + log(0.2)/log(ngrid)
        !
        write(iunit,*) irefine, nboxes, ngrid, dble(nboxes)/dble(ngrid**2), 0.2d0,&
                       log(1d0*nboxes)/log(1d0*ngrid), 2d0 + log(0.2d0)/log(1d0*ngrid)
!$omp end critical
      end if
      if(irefine.eq.nrefine-3) fraction=dble(nboxes)/dble(ngrid**2)
    enddo
    close(iunit)

  end subroutine fract_dimension

end module cut_detector
