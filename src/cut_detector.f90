module cut_detector
    use util, only: twopi
    use simple, only: tstep
    use orbit_symplectic, only: symplectic_integrator_t
    use field_can_mod, only: field_can_t
#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
    use params, only: debug
#endif

    implicit none

    ! Define real(dp) kind parameter
    integer, parameter :: dp = kind(1.0d0)
    save

    integer, parameter :: n_tip_vars = 6
    integer, parameter :: nplagr = 6
    integer, parameter :: nder = 0

    public

    type :: cut_detector_t
      real(dp) :: fper  ! field period

      ! for Poincare cuts
      real(dp) :: alam_prev, par_inv
      integer          :: iper, itip, kper

      real(dp) :: orb_sten(6,nplagr), coef(0:nder,nplagr)
      integer :: ipoi(nplagr)
    end type cut_detector_t

  contains

    subroutine init(self, fper, z)
      type(cut_detector_t) :: self
      real(dp), intent(in) :: fper
      real(dp), intent(in) :: z(:)
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

      type(cut_detector_t) :: self
      type(symplectic_integrator_t) :: si
      type(field_can_t) :: f

      real(dp), intent(inout) :: z(:)
      ! variables to evaluate at tip: z(1..5), par_inv
      real(dp), dimension(:), intent(inout) :: var_cut
      integer, intent(out) :: cut_type
      integer, intent(out) :: ierr

      integer, parameter :: nstep_max = 1000000000
      integer :: i
      real(dp) :: phiper = 0.0d0

      do i=1, nstep_max
        call tstep(si, f, z, ierr)
        if(ierr.ne.0) exit

        self%par_inv = self%par_inv+z(5)**2 ! parallel adiabatic invariant

        if(i.le.nplagr) then          !<=first nplagr points to initialize stencil
          ! Note: i is guaranteed to be <= nplagr here, so array access is safe
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
      real(dp) :: fraction,rmax,rmin,tmax,tmin,hr,ht
      real(dp), dimension(2,ntr)              :: rt!0, rt
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
#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
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
#endif
        if(irefine.eq.nrefine-3) fraction=dble(nboxes)/dble(ngrid**2)
      enddo
      close(iunit)

    end subroutine fract_dimension

  end module cut_detector
