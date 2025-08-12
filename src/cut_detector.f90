module cut_detector
    use util, only: twopi
    use simple, only: tstep
    use orbit_symplectic, only: SymplecticIntegrator
    use field_can_mod, only: FieldCan
    use params, only: debug

    implicit none

    ! Define real(dp) kind parameter
    integer, parameter :: dp = kind(1.0d0)
    save

    integer, parameter :: n_tip_vars = 6
    integer, parameter :: nplagr = 6
    integer, parameter :: nder = 0

    public

    type :: CutDetector
      real(dp) :: fper              ! Field period
      
      ! Poincaré cut detection state
      real(dp) :: alam_prev         ! Previous lambda value for tip detection
      real(dp) :: par_inv           ! Parallel adiabatic invariant
      integer  :: iper              ! Period counter
      integer  :: itip              ! Tip counter
      integer  :: kper              ! Period index
      
      ! Interpolation arrays
      real(dp) :: orb_sten(6,nplagr)      ! Orbit stencil for interpolation
      real(dp) :: coef(0:nder,nplagr)     ! Lagrange coefficients
      integer  :: ipoi(nplagr)            ! Point indices for stencil
    contains
      ! Type-bound procedures for better encapsulation
      procedure :: detect_tip_crossing => detector_detect_tip_crossing
      procedure :: detect_period_crossing => detector_detect_period_crossing
      procedure :: update_stencil => detector_update_stencil
      procedure :: interpolate_cut => detector_interpolate_cut
    end type CutDetector

  contains

    subroutine init(self, fper, z)
      type(CutDetector) :: self
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
      type(CutDetector) :: self
      type(SymplecticIntegrator) :: si
      type(FieldCan) :: f

      real(dp), intent(inout) :: z(:)
      real(dp), dimension(:), intent(inout) :: var_cut  ! variables to evaluate at tip: z(1..5), par_inv
      integer, intent(out) :: cut_type
      integer, intent(out) :: ierr

      integer, parameter :: nstep_max = 1000000000
      integer :: i
      real(dp) :: phiper = 0.0d0
      logical :: cut_found

      do i = 1, nstep_max
        call tstep(si, f, z, ierr)
        if (ierr /= 0) exit

        ! Update parallel adiabatic invariant
        self%par_inv = self%par_inv + z(5)**2

        ! Update orbit stencil for interpolation
        call self%update_stencil(z, i)

        ! Only check for cuts after stencil is initialized
        if (i > nplagr) then
          ! Check for tip crossing
          call self%detect_tip_crossing(z(5), cut_found)
          if (cut_found) then
            cut_type = 0  ! Tip cut
            call self%interpolate_cut(5, 0.0d0, var_cut)
            self%par_inv = self%par_inv - var_cut(6)
            return
          end if

          ! Check for period crossing
          call self%detect_period_crossing(z(3), phiper, cut_found)
          if (cut_found) then
            cut_type = 1  ! Period cut
            call self%interpolate_cut(3, phiper, var_cut)
            return
          end if
        end if
      end do

      ! If we reach here, no cut was found within nstep_max
      ierr = 1
    end subroutine trace_to_cut

    subroutine fract_dimension(ntr, rt, fraction)
      integer, intent(in) :: ntr
      real(dp), dimension(2,ntr), intent(in) :: rt
      real(dp), intent(out) :: fraction

      integer, parameter :: iunit = 1003
      integer :: itr, ngrid, nrefine, irefine, kr, kt, nboxes
      real(dp) :: rmax, rmin, tmax, tmin, hr, ht
      logical, dimension(:,:), allocatable :: free

      ! Compute data bounds
      call compute_data_bounds(rt, ntr, rmin, rmax, tmin, tmax)

      ! Determine refinement levels based on data size
      nrefine = int(log(dble(ntr))/log(4.d0)) + 3

      ngrid = 1
      do irefine = 1, nrefine
        ngrid = ngrid * 2
        
        ! Allocate grid for this refinement level
        !$omp critical
        allocate(free(0:ngrid,0:ngrid))
        !$omp end critical
        
        ! Count occupied boxes using optimized algorithm
        call count_occupied_boxes(rt, ntr, rmin, rmax, tmin, tmax, ngrid, free, nboxes)
        
        !$omp critical
        deallocate(free)
        !$omp end critical
        
        ! Debug output and fraction calculation
        call output_fractal_debug(iunit, irefine, nboxes, ngrid, nrefine, fraction)
      end do
      
      close(iunit)
    end subroutine fract_dimension

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Update orbit stencil for interpolation
  subroutine detector_update_stencil(self, z, step)
    class(CutDetector), intent(inout) :: self
    real(dp), intent(in) :: z(:)
    integer, intent(in) :: step
    
    if (step <= nplagr) then
      ! Initialize stencil with first nplagr points
      self%orb_sten(1:5, step) = z
      self%orb_sten(6, step) = self%par_inv
    else
      ! Shift stencil and add new point
      self%orb_sten(1:5, self%ipoi(1)) = z
      self%orb_sten(6, self%ipoi(1)) = self%par_inv
      self%ipoi = cshift(self%ipoi, 1)
    end if
  end subroutine detector_update_stencil
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Detect tip crossing
  subroutine detector_detect_tip_crossing(self, lambda, cut_found)
    class(CutDetector), intent(inout) :: self
    real(dp), intent(in) :: lambda
    logical, intent(out) :: cut_found
    
    ! Check for sign change indicating tip crossing
    if (self%alam_prev < 0.0d0 .and. lambda > 0.0d0) then
      self%itip = 0  ! Reset tip counter
    end if
    
    self%itip = self%itip + 1
    self%alam_prev = lambda
    
    ! Cut found when stencil is complete around tip
    cut_found = (self%itip == nplagr/2)
  end subroutine detector_detect_tip_crossing
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Detect period crossing
  subroutine detector_detect_period_crossing(self, phi, phiper, cut_found)
    class(CutDetector), intent(inout) :: self
    real(dp), intent(in) :: phi
    real(dp), intent(out) :: phiper
    logical, intent(out) :: cut_found
    logical :: crossed
    
    ! Initialize output variables
    phiper = 0.0d0
    crossed = .false.
    
    ! Check for period boundary crossings
    if (phi > dble(self%kper + 1) * self%fper) then
      self%iper = 0  ! Reset period counter
      phiper = dble(self%kper + 1) * self%fper
      self%kper = self%kper + 1
      crossed = .true.
    else if (phi < dble(self%kper) * self%fper) then
      self%iper = 0  ! Reset period counter
      phiper = dble(self%kper) * self%fper
      self%kper = self%kper - 1
      crossed = .true.
    end if
    
    self%iper = self%iper + 1
    
    ! Cut found when stencil is complete around period boundary AND crossing occurred
    cut_found = crossed .and. (self%iper == nplagr/2)
  end subroutine detector_detect_period_crossing
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Interpolate cut using Lagrange coefficients
  subroutine detector_interpolate_cut(self, coord_index, target_value, var_cut)
    use plag_coeff_sub, only : plag_coeff
    
    class(CutDetector), intent(inout) :: self
    integer, intent(in) :: coord_index
    real(dp), intent(in) :: target_value
    real(dp), dimension(:), intent(out) :: var_cut
    
    ! Compute Lagrange coefficients
    call plag_coeff(nplagr, nder, target_value, &
                   self%orb_sten(coord_index, self%ipoi), self%coef)
    
    ! Interpolate all variables at the cut
    var_cut = matmul(self%orb_sten(:, self%ipoi), self%coef(0,:))
    
    ! Normalize angular coordinates
    var_cut(2) = modulo(var_cut(2), twopi)
    var_cut(3) = modulo(var_cut(3), twopi)
  end subroutine detector_interpolate_cut
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute data bounds for fractal dimension calculation
  subroutine compute_data_bounds(rt, ntr, rmin, rmax, tmin, tmax)
    integer, intent(in) :: ntr
    real(dp), dimension(2,ntr), intent(in) :: rt
    real(dp), intent(out) :: rmin, rmax, tmin, tmax
    
    rmin = minval(rt(1,:))
    rmax = maxval(rt(1,:))
    tmin = minval(rt(2,:))
    tmax = maxval(rt(2,:))
  end subroutine compute_data_bounds
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Count occupied boxes for fractal dimension (optimized)
  subroutine count_occupied_boxes(rt, ntr, rmin, rmax, tmin, tmax, ngrid, free, nboxes)
    integer, intent(in) :: ntr, ngrid
    real(dp), dimension(2,ntr), intent(in) :: rt
    real(dp), intent(in) :: rmin, rmax, tmin, tmax
    logical, dimension(0:ngrid,0:ngrid), intent(inout) :: free
    integer, intent(out) :: nboxes
    
    integer :: itr, kr, kt
    real(dp) :: hr, ht
    
    free = .true.
    hr = (rmax - rmin) / dble(ngrid)
    ht = (tmax - tmin) / dble(ngrid)
    nboxes = 0
    
    do itr = 1, ntr
      kr = int((rt(1,itr) - rmin) / hr)
      kr = min(ngrid-1, max(0, kr))
      kt = int((rt(2,itr) - tmin) / ht)
      kt = min(ngrid-1, max(0, kt))
      
      if (free(kr,kt)) then
        free(kr,kt) = .false.
        nboxes = nboxes + 1
      end if
    end do
  end subroutine count_occupied_boxes
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Output fractal dimension debug information
  subroutine output_fractal_debug(iunit, irefine, nboxes, ngrid, nrefine, fraction)
    integer, intent(in) :: iunit, irefine, nboxes, ngrid, nrefine
    real(dp), intent(inout) :: fraction
    
    if (debug) then
      !$omp critical
      ! Right now criterion for regular is at nboxes/ngrid**2 < 0.2
      ! For fractal dimension d = log(nboxes)/log(ngrid) this means
      ! d_thresh = 2 + log(0.2)/log(ngrid)
      write(iunit,*) irefine, nboxes, ngrid, dble(nboxes)/dble(ngrid**2), 0.2d0, &
                     log(1d0*nboxes)/log(1d0*ngrid), 2d0 + log(0.2d0)/log(1d0*ngrid)
      !$omp end critical
    end if
    
    if (irefine == nrefine - 3) then
      fraction = dble(nboxes) / dble(ngrid**2)
    end if
  end subroutine output_fractal_debug

  end module cut_detector
