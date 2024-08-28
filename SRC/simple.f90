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
  use params, only : Tracer, idx
  use boozer_sub, only : boozer_converter
  use chamb_sub, only : chamb_can

implicit none
save

public

  interface tstep
      module procedure timestep
      module procedure timestep_z
      module procedure timestep_sympl_z
   end interface tstep

contains

  subroutine init_field(self, vmec_file, ans_s, ans_tp, amultharm, aintegmode)
    use spline_vmec_sub, only : spline_vmec_data, volume_and_B00
    use get_can_sub, only : get_canonical_coordinates
    use magfie_sub, only : init_magfie
    use vmecin_sub, only : stevvo

    ! initialize field geometry
    character(*), intent(in) :: vmec_file
    type(Tracer), intent(inout) :: self
    integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
    integer             :: ierr
    integer             :: L1i
    double precision    :: RT0, R0i, cbfi, bz0i, bf0, volume, B00
    double precision    :: z(5)

    call init_magfie(isw_field_type)

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
      print *, 'dtau = ', dtau, ' dtaumin = ', dtaumin
      error stop 'orbit_sympl_init - error: dtau/dtaumin not integer'
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
    use alpha_lifetime_sub, only : orbit_timestep_axis

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
  use simple, only: tstep
  use orbit_symplectic, only: SymplecticIntegrator
  use field_can_mod, only: FieldCan
  use params, only: debug

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
    use plag_coeff_sub, only : plag_coeff

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

module simple_main
  use omp_lib
  use util, only: pi, twopi, c, e_charge, e_mass, p_mass, ev, sqrt2
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp, &
    vmec_B_scale, vmec_RZ_scale

  use velo_mod,   only : isw_field_type
  use orbit_symplectic, only : orbit_timestep_sympl, get_val
  use simple, only : init_field, init_sympl, eval_field
  use cut_detector, only : fract_dimension
  use diag_mod, only : icounter
  use collis_alp, only : loacol_alpha, stost
  use params
  use binsrc_sub, only : binsrc
  use boozer_sub, only : vmec_to_boozer, boozer_to_vmec
  use check_orbit_type_sub, only : check_orbit_type

  implicit none

  public

contains

subroutine run(norb)
  type(Tracer), intent(inout) :: norb
  integer :: i

  call print_parameters
  if (swcoll) call init_collisions

  ! pre-compute starting flux surface
  call init_starting_surf

  ! local?
  if(1 == num_surf) then
    call init_starting_points
  else
    call init_starting_points_global
  endif
  if (startmode == 0) stop 'startmode == 0, stopping after generating start.dat'

  call init_counters

  !$omp parallel firstprivate(norb)
  !$omp do
  do i = 1, ntestpart
    !$omp critical
    kpart = kpart+1
    print *, kpart, ' / ', ntestpart, 'particle: ', i, 'thread: ', omp_get_thread_num()
    !$omp end critical
    call trace_orbit(norb, i)
  end do
  !$omp end do
  !$omp end parallel

  confpart_pass=confpart_pass/ntestpart
  confpart_trap=confpart_trap/ntestpart

end subroutine run

subroutine print_parameters
  print *, 'tau: ', dtau, dtaumin, min(dabs(mod(dtau, dtaumin)), &
                    dabs(mod(dtau, dtaumin)-dtaumin))/dtaumin, ntau
  print *, 'v0 = ', v0
end subroutine print_parameters

subroutine init_collisions
  real(8) :: v0_coll

  call loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe, &
    3.5d6/facE_al,v0_coll,dchichi,slowrate,dchichi_norm,slowrate_norm)

  if (abs(v0_coll - v0) > 1d-6) then
    error stop 'simple_main.init_collisions: v0_coll != v0'
  end if
end subroutine init_collisions

subroutine init_counters
  icounter=0 ! evaluation counter
  kpart=0

  ! initialize array of confined particle percentage
  confpart_trap=0.d0
  confpart_pass=0.d0

  ! initialize lost times when particles get lost
  times_lost = -1.d0
end subroutine init_counters


subroutine finalize
  use get_can_sub, only : deallocate_can_coord

  if (integmode >= 0) call deallocate_can_coord

  deallocate(times_lost, confpart_trap, confpart_pass, trap_par, &
    perp_inv, iclass)
  deallocate(xstart, bstart, volstart, zstart)
end subroutine finalize

subroutine write_output

  integer :: i, num_lost
  double precision :: inverse_times_lost_sum

  open(1,file='times_lost.dat',recl=1024)
  num_lost = 0
  inverse_times_lost_sum = 0.0d0
  do i=1,ntestpart
    write(1,*) i, times_lost(i), trap_par(i), zstart(1,i), perp_inv(i), zend(:,i)
    if (times_lost(i) > 0.0d0 .and. times_lost(i) < trace_time) then
      num_lost = num_lost + 1
      inverse_times_lost_sum = inverse_times_lost_sum + 1.0/times_lost(i)
    end if
  enddo
  close(1)
  open(1,file='avg_inverse_t_lost.dat',recl=1024) ! Write average loss time
  write(1,*) inverse_times_lost_sum/num_lost
  close(1)

  open(1,file='confined_fraction.dat',recl=1024)
  do i=1,ntimstep
    write(1,*) dble(i-1)*dtau/v0,confpart_pass(i),confpart_trap(i),ntestpart
  enddo
  close(1)

  open(1,file='class_parts.dat',recl=1024)
  do i=1,ntestpart
    write(1,*) i, zstart(1,i), perp_inv(i), iclass(:,i)
  enddo
  close(1)

end subroutine write_output

subroutine init_starting_surf
  use alpha_lifetime_sub, only : integrate_mfl_can

  xstart=0.d0
  bstart=0.d0
  volstart=0.d0 !ToDo add loop for all sbeg, check dimension

  call integrate_mfl_can( &
    npoiper*nper,dphi,sbeg(1),phibeg,thetabeg, &
    xstart,bstart,volstart,bmod00,ierr)

  if(ierr.ne.0) then
    print *,'starting field line has points outside the chamber'
    stop
  endif

! maximum value of B module:
  bmax=maxval(bstart)
  bmin=minval(bstart)

  print *, 'bmod00 = ', bmod00, 'bmin = ', bmin, 'bmax = ', bmax
end subroutine init_starting_surf

subroutine init_starting_points_ants(unit)
  use parse_ants, only : process_line
  use get_can_sub, only : vmec_to_can

  integer, intent(in) :: unit

  integer, parameter :: maxlen = 4096
  character(len=maxlen) :: line
  real(8) :: v_par, v_perp, u, v, s
  real(8) :: th, ph, th_c, ph_c  ! Canonical flux coordinate angles
  integer :: ipart

  do ipart=1,ntestpart
    read(unit, '(A)') line
    call process_line(line, v_par, v_perp, u, v, s)
    ! In the test case, u runs from 0 to 1 and v from 0 to 4
    th = 2d0*pi*u
    ph = 2d0*pi*v/4d0
    call vmec_to_can(s, th, ph, th_c, ph_c)
    zstart(1, ipart) = s
    zstart(2, ipart) = ph_c
    zstart(3, ipart) = th_c
    zstart(4, ipart) = 1.d0
    zstart(5, ipart) = v_par / sqrt(v_par**2 + v_perp**2)
  enddo
end subroutine

subroutine init_starting_points
  use get_can_sub, only: can_to_vmec

  integer :: i, ipart, iskip
  double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec

  ! skip random numbers according to configuration
    do iskip=1,loopskip
      do ipart=1,ntestpart
        call random_number(xi)
        call random_number(xi)
      enddo
    enddo

  ! files for storing starting coords
  open(1,file='start.dat',recl=1024)
  ! determine the starting point:
  if (startmode == 0 .or. startmode == 1) then
    do ipart=1,ntestpart
      call random_number(xi)
      call binsrc(volstart,1,npoiper*nper,xi,i)
      ibins=i
      ! coordinates: z(1) = r, z(2) = vartheta, z(3) = varphi
      r=xstart(1,i)
      vartheta=xstart(2,i)
      varphi=xstart(3,i)
!
! we store starting points in VMEC coordinates:
      if(isw_field_type.eq.0) then
        call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      elseif(isw_field_type.eq.1) then
        theta_vmec=vartheta
        varphi_vmec=varphi
      elseif(isw_field_type.eq.2) then
        call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      else
        print *,'init_starting_points: unknown field type'
      endif
!
      zstart(1,ipart)=r
      zstart(2,ipart)=theta_vmec
      zstart(3,ipart)=varphi_vmec
      ! normalized velocity module z(4) = v / v_0:
      zstart(4,ipart)=1.d0
      ! starting pitch z(5)=v_\parallel / v:
      call random_number(xi)
      zstart(5,ipart)=2.d0*(xi-0.5d0)
      write(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 2) then
    do ipart=1,ntestpart
      read(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 3) then  ! ANTS input mode
    call init_starting_points_ants(1)
  else if (startmode == 4) then
    !select only the indices from batch and overwrite zstart.
    do ipart=idx(0),idx(ntestpart)
      read(1,*) zstart(:,ipart)
    enddo
    ! indices no longer needed
    deallocate(idx)
  endif

  close(1)
end subroutine init_starting_points

subroutine init_starting_points_global

  use find_bminmax_sub, only : get_bminmax
  use get_can_sub, only: can_to_vmec

  integer, parameter :: ns=1000
  integer :: ipart,iskip,is,s_idx,parts_per_s
  double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec
  double precision :: s,bmin,bmax
!
  open(1,file='bminmax.dat',recl=1024)
  do is=0,ns
    s=dble(is)/dble(ns)
!
    call get_bminmax(s,bmin,bmax)
!
    write(1,*) s,bmin,bmax
  enddo
  close(1)

  ! skip random numbers according to configuration
    do iskip=1,loopskip
      do ipart=1,ntestpart
        call random_number(xi)
        call random_number(xi)
        call random_number(xi)
        call random_number(xi)
      enddo
    enddo

  ! files for storing starting coords
  open(1,file='start.dat',recl=1024)
  ! determine the starting point:
  if (startmode == 0 .or. startmode == 1) then
    print *, "Initialising for", num_surf, "surfaces."
    do ipart=1,ntestpart
      if (0 == num_surf) then
        call random_number(xi)
        r = xi
      else if (1 < num_surf) then
        parts_per_s = int(ntestpart/num_surf)
        s_idx = (ipart/parts_per_s)+1
        r = sbeg(s_idx)
      else ! Should not happen (as we are not in "local mode"), however let's catch it anyway.
        r = sbeg(1)
      endif

      call random_number(xi)
      vartheta=twopi*xi
      call random_number(xi)
      varphi=twopi*xi
!
! we store starting points in VMEC coordinates:
      if(isw_field_type.eq.0) then
        call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      elseif(isw_field_type.eq.1) then
        theta_vmec=vartheta
        varphi_vmec=varphi
      elseif(isw_field_type.eq.2) then
        call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      else
        print *,'init_starting_points: unknown field type'
      endif
!
      zstart(1,ipart)=r
      zstart(2,ipart)=theta_vmec
      zstart(3,ipart)=varphi_vmec
      ! normalized velocity module z(4) = v / v_0:
      zstart(4,ipart)=1.d0
      ! starting pitch z(5)=v_\parallel / v:
      call random_number(xi)
      zstart(5,ipart)=2.d0*(xi-0.5d0)
      write(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 2) then
    do ipart=1,ntestpart
      read(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 3) then  ! ANTS input mode
    call init_starting_points_ants(1)
  else if (startmode == 4) then
    !select only the indices from batch and overwrite zstart.
    do ipart=idx(0),idx(ntestpart)
      read(1,*) zstart(:,ipart)
    enddo
    ! indices no longer needed
    deallocate(idx)
  endif

  close(1)
end subroutine init_starting_points_global

subroutine trace_orbit(anorb, ipart)
  use find_bminmax_sub, only : get_bminmax
  use get_can_sub, only : vmec_to_can
  use magfie_sub, only : magfie
  use plag_coeff_sub, only : plag_coeff
  use alpha_lifetime_sub, only : orbit_timestep_axis

  type(Tracer), intent(inout) :: anorb
  integer, intent(in) :: ipart
  integer :: ierr, ierr_coll
  double precision, dimension(5) :: z
  double precision :: bmod,sqrtg
  double precision, dimension(3) :: bder, hcovar, hctrvr, hcurl
  integer :: it, ktau
  integer(8) :: kt
  logical :: passing

  integer                                       :: ifp_tip,ifp_per
  integer,          dimension(:),   allocatable :: ipoi
  double precision, dimension(:),   allocatable :: xp
  double precision, dimension(:,:), allocatable :: coef,orb_sten
  double precision, dimension(:,:), allocatable :: zpoipl_tip,zpoipl_per,dummy2d
  double precision, dimension(n_tip_vars)       :: var_tip
  double precision :: phiper, alam_prev, par_inv
  integer :: iper, itip, kper, nfp_tip, nfp_per

  double precision :: fraction
  double precision :: r,theta_vmec,varphi_vmec
  logical :: regular

! Variables and settings for classification by J_parallel and ideal orbit condition:
  integer, parameter :: nfp_dim=3, nturns=8
  integer :: nfp_cot,ideal,ijpar,ierr_cot,iangvar
  double precision, dimension(nfp_dim) :: fpr_in

! for run with fixed random seed
  integer :: seedsize
  integer, allocatable :: seed(:)

  zend(:,ipart) = 0d0

  if (deterministic) then
    call random_seed(size = seedsize)
    if (.not. allocated(seed)) allocate(seed(seedsize))
    seed = 0
    call random_seed(put=seed)
  endif
!
  iangvar=2
! End variables and settings for classification by J_parallel and ideal orbit condition
!

!  open(unit=10000+ipart, iostat=stat, status='old')
!  if (stat == 0) close(10000+ipart, status='delete')
!  open(unit=20000+ipart, iostat=stat, status='old')
!  if (stat == 0) close(20000+ipart, status='delete')

! Write out trapped-passing boundary at the classification cut:
  if(class_plot) then
    if(ipart.eq.1) then
      z(1)=zstart(1,ipart)
      z(3)=cut_in_per*fper
      do kt=0,1000
        z(2)=1d-3*twopi*dble(kt)
        call magfie(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
        write(iaaa_bou,*) z(2),sqrt(1.d0-bmod/bmax)
      enddo
    endif
  endif
! End write out trapped-passing boundary at the classification cut
!
  z = zstart(:, ipart)
  r=z(1)
  theta_vmec=z(2)
  varphi_vmec=z(3)
!
  if(isw_field_type.eq.0) then
      call vmec_to_can(r,theta_vmec,varphi_vmec,z(2),z(3))
  elseif(isw_field_type.eq.2) then
      call vmec_to_boozer(r,theta_vmec,varphi_vmec,z(2),z(3))
  endif

! In case of classification plot all starting points are moved to the classification cut:
  if(class_plot) then
    z(3)=cut_in_per*fper
    zstart(2,ipart)=modulo(zstart(2,ipart),twopi)
  endif
! End moving starting points to the classification cut

  if (integmode>0) call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)

  call magfie(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)

!$omp critical
  if(num_surf > 1) then
    call get_bminmax(z(1),bmin,bmax)
  endif
  passing = z(5)**2.gt.1.d0-bmod/bmax
  trap_par(ipart) = ((1.d0-z(5)**2)*bmax/bmod-1.d0)*bmin/(bmax-bmin)
  perp_inv(ipart) = z(4)**2*(1.d0-z(5)**2)/bmod
  iclass(:,ipart) = 0
!$omp end critical

! Forced classification of passing as regular:
  if(passing.and.(notrace_passing.eq.1 .or. trap_par(ipart).le.contr_pp)) then
    ! passing particle
    ! no tracing of passing particles, assume that all are confined
    ! or: strongly passing particles that are certainly confined
    !$omp critical
    confpart_pass=confpart_pass+1.d0
    !$omp end critical
    if(class_plot) then
!$omp critical
      write (iaaa_pnt,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
!$omp end critical
    endif
    iclass(:,ipart) = 1
    return
  endif
! End forced classification of passing as regular

!$omp critical
  if (.not. allocated(ipoi)) &
  allocate(ipoi(nplagr),coef(0:nder,nplagr),orb_sten(6,nplagr),xp(nplagr))
!$omp end critical
  do it=1,nplagr
    ipoi(it)=it
  enddo

  nfp_tip=nfp             !<= initial array dimension for tips
  nfp_per=nfp             !<= initial array dimension for periods
!$omp critical
  if (.not. allocated(zpoipl_tip)) &
  allocate(zpoipl_tip(2,nfp_tip),zpoipl_per(2,nfp_per))
!$omp end critical

!  open(unit=10000+ipart, recl=1024, position='append')
!  open(unit=20000+ipart, recl=1024, position='append')

  ifp_tip=0               !<= initialize footprint counter on tips
  ifp_per=0               !<= initialize footprint counter on periods

  icounter=0
  phiper=0.0d0


  kt = 0
  if (passing) then
    !$omp atomic
    confpart_pass(1)=confpart_pass(1)+1.d0
  else
    !$omp atomic
    confpart_trap(1)=confpart_trap(1)+1.d0
  end if

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
!
! Initialize classification by J_parallel and ideal orbit condition:
  nfp_cot=0
! End Initialize classification by J_parallel and ideal orbit condition
!
  par_inv = 0d0
  regular = .False.
  do it=2,ntimstep
    if (regular) then  ! regular orbit, will not be lost
      if(passing) then
        !$omp atomic
        confpart_pass(it)=confpart_pass(it)+1.d0
      else
        !$omp atomic
        confpart_trap(it)=confpart_trap(it)+1.d0
      endif
      kt = kt+ntau
      cycle
    endif
    do ktau=1,ntau
      if (integmode <= 0) then
        call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
      else
        if (swcoll) then  ! TODO: move this inside modules
          anorb%si%pabs = z(4)
          anorb%f%vpar = z(4)*z(5)*sqrt2
          anorb%f%mu = z(4)**2*(1.d0-z(5)**2)/anorb%f%Bmod
          anorb%si%z(4) = anorb%f%vpar*anorb%f%hph + anorb%f%Aph/anorb%f%ro0
          call get_val(anorb%f, anorb%si%z(4)) ! for pth
          anorb%si%pthold = anorb%f%pth
        endif
        call orbit_timestep_sympl(anorb%si, anorb%f, ierr)
        z(1:3) = anorb%si%z(1:3)
        z(4) = dsqrt(anorb%f%mu*anorb%f%Bmod+0.5d0*anorb%f%vpar**2)
        z(5) = anorb%f%vpar/(z(4)*sqrt2)
      endif

      ! Collisions
      if (swcoll) then
        call stost(z, dtaumin, 1, ierr_coll)
        if (ierr_coll /= 0) then
          print *, 'Error in stost: ', ierr_coll, 'z = ', z, 'dtaumin = ', dtaumin
        endif
      endif

      ! Write starting data for orbits which were lost in case of classification plot
      if(class_plot .and. ierr.ne.0) then
        call output_lost_orbit_starting_data(ipart, passing)
      endif
      ! End write starting data for orbits which were lost in case of classification plot

      if(ierr.ne.0) exit
      kt = kt+1

      par_inv = par_inv+z(5)**2*dtaumin ! parallel adiabatic invariant
      if(kt.le.nplagr) then          !<=first nplagr points to initialize stencil
        orb_sten(1:5,kt)=z
        orb_sten(6,kt)=par_inv
      else                          !<=normal case, shift stencil
        orb_sten(1:5,ipoi(1))=z
        orb_sten(6,ipoi(1))=par_inv
        ipoi=cshift(ipoi,1)
      endif

      ! Tip detection and interpolation
      if(alam_prev.lt.0.d0.and.z(5).gt.0.d0) itip=0   !<=tip has been passed
      itip=itip+1
      alam_prev=z(5)
      if(kt.gt.nplagr) then          !<=use only initialized stencil
        if(itip.eq.npl_half) then   !<=stencil around tip is complete, interpolate
          xp=orb_sten(5,ipoi)

          call plag_coeff(nplagr,nder,zerolam,xp,coef)

          var_tip=matmul(orb_sten(:,ipoi),coef(0,:))
          var_tip(2)=modulo(var_tip(2),twopi)
          var_tip(3)=modulo(var_tip(3),twopi)

!          write(10000+ipart,*) var_tip

          ifp_tip=ifp_tip+1
          if(ifp_tip.gt.nfp_tip) then   !<=increase the buffer for banana tips
            !$omp critical
            allocate(dummy2d(2,ifp_tip-1))
            !$omp end critical
            dummy2d=zpoipl_tip(:,1:ifp_tip-1)
            !$omp critical
            deallocate(zpoipl_tip)
            !$omp end critical
            nfp_tip=nfp_tip+nfp
            !$omp critical
            allocate(zpoipl_tip(2,nfp_tip))
            !$omp end critical
            zpoipl_tip(:,1:ifp_tip-1)=dummy2d
            !$omp critical
            deallocate(dummy2d)
            !$omp end critical
          endif
          zpoipl_tip(:,ifp_tip)=var_tip(1:2)
          par_inv = par_inv - var_tip(6)
!
! Classification by J_parallel and ideal orbit conditions:
          fpr_in(1)=var_tip(1)
          fpr_in(2)=var_tip(iangvar)
          fpr_in(3)=var_tip(6)
!
          call check_orbit_type(nturns,nfp_cot,fpr_in,ideal,ijpar,ierr_cot)
!
          iclass(1,ipart) = ijpar
          iclass(2,ipart) = ideal
          if(fast_class) ierr=ierr_cot
!
! End classification by J_parallel and ideal orbit conditions
        endif
      endif
      ! End tip detection and interpolation

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
      if(kt.gt.nplagr) then          !<=use only initialized stencil
        if(iper.eq.npl_half) then   !<=stencil around periodic boundary is complete, interpolate
          xp=orb_sten(3,ipoi)-phiper

          call plag_coeff(nplagr,nder,zerolam,xp,coef)

          var_tip=matmul(orb_sten(:,ipoi),coef(0,:))
          var_tip(2)=modulo(var_tip(2),twopi)
          var_tip(3)=modulo(var_tip(3),twopi)
! write(20000+ipart,*) var_tip
          ifp_per=ifp_per+1
          if(ifp_per.gt.nfp_per) then   !<=increase the buffer for periodic boundary footprints
            !$omp critical
            allocate(dummy2d(2,ifp_per-1))
            !$omp end critical
            dummy2d=zpoipl_per(:,1:ifp_per-1)
            !$omp critical
            deallocate(zpoipl_per)
            !$omp end critical
            nfp_per=nfp_per+nfp
            !$omp critical
            allocate(zpoipl_per(2,nfp_per))
            !$omp end critical
            zpoipl_per(:,1:ifp_per-1)=dummy2d
            !$omp critical
            deallocate(dummy2d)
            !$omp end critical
          endif
          zpoipl_per(:,ifp_per)=var_tip(1:2)
        endif
      endif
      ! End periodic boundary footprint detection and interpolation

      ! Cut classification into regular or chaotic
      if (kt == ntcut) then
        regular = .True.

        if(ifp_per > 0) then

          call fract_dimension(ifp_per,zpoipl_per(:,1:ifp_per),fraction)

          if(fraction.gt.0.2d0) then
            print *, ipart, ' chaotic per ', ifp_per
            regular = .False.
          else
            print *, ipart, ' regular per', ifp_per
          endif
        endif

        if(ifp_tip > 0) then

          call fract_dimension(ifp_tip,zpoipl_tip(:,1:ifp_tip),fraction)

          if(fraction.gt.0.2d0) then
            print *, ipart, ' chaotic tip ', ifp_tip
            regular = .False.
            iclass(3,ipart) = 2
          else
            print *, ipart, ' regular tip ', ifp_tip
            iclass(3,ipart) = 1
          endif
        endif

        if(class_plot) then
          call output_minkowsky_class(ipart, regular, passing)
          ierr=1
        endif
      endif
!
      if(ierr.ne.0) then
        if(class_plot .and. .not. passing) then
          call output_jpar_class(ipart, ijpar)
          call output_topological_class(ipart, ideal)
          exit
        endif
      endif
!    write(999, *) kt*dtaumin/v0, z
    enddo
    if(ierr.ne.0) exit
    if(passing) then
      !$omp atomic
      confpart_pass(it)=confpart_pass(it)+1.d0
    else
      !$omp atomic
      confpart_trap(it)=confpart_trap(it)+1.d0
    endif
  enddo

  !$omp critical
  zend(:,ipart) = z
  if(isw_field_type.eq.0) then
    ! TODO need to add can_to_vmec
    ! call can_to_vmec(z(1),z(2),z(3),zend(2,ipart),zend(3,ipart))
  elseif(isw_field_type.eq.2) then
    call boozer_to_vmec(z(1),z(2),z(3),zend(2,ipart),zend(3,ipart))
  endif
  times_lost(ipart) = kt*dtaumin/v0
  deallocate(zpoipl_tip, zpoipl_per)
  !$omp end critical
!  close(unit=10000+ipart)
!  close(unit=10000+ipart)
end subroutine trace_orbit

subroutine output_lost_orbit_starting_data(ipart, passing)
  use params, only : iaaa_prp, iaaa_prt

  integer, intent(in) :: ipart
  logical, intent(in) :: passing

  !$omp critical
  if(passing) then
    write (iaaa_prp,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
  else
    write (iaaa_prt,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
  endif
  !$omp end critical
end subroutine output_lost_orbit_starting_data

subroutine output_minkowsky_class(ipart, regular, passing)
  use params, only : iaaa_rep, iaaa_ret, iaaa_stp, iaaa_stt

  integer, intent(in) :: ipart
  logical, intent(in) :: regular, passing

  !$omp critical
    if(regular) then
      if(passing) then
        write (iaaa_rep,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
      else
        write (iaaa_ret,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
      endif
    else
      if(passing) then
        write (iaaa_stp,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
      else
        write (iaaa_stt,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
      endif
    endif
  !$omp end critical
end subroutine output_minkowsky_class

subroutine output_jpar_class(ipart, ijpar)
  use params, only : iaaa_jer, iaaa_jre, iaaa_jst

  integer, intent(in) :: ipart, ijpar

  !$omp critical
  select case(ijpar)
  case(0)
    write (iaaa_jer,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
  case(1)
    write (iaaa_jre,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
  case(2)
    write (iaaa_jst,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
  end select
  !$omp end critical
end subroutine output_jpar_class

subroutine output_topological_class(ipart, ideal)
  use params, only : iaaa_ier, iaaa_ire, iaaa_ist

  integer, intent(in) :: ipart, ideal

  select case(ideal)
  case(0)
    write (iaaa_ier,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
  case(1)
    write (iaaa_ire,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
  case(2)
    write (iaaa_ist,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
  end select
end subroutine output_topological_class

end module simple_main
