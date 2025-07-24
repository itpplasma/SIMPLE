!> Module for reading and evaluating magnetic fields from BOOZXFORM output files
!> This allows direct input of pre-computed Boozer coordinate fields
module field_booz_xform
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_base, only: MagneticField
  implicit none

  private
  public :: BoozXformField
  
  ! Module variable for warning message
  logical :: booz_warning_printed = .false.

  !> Magnetic field representation from BOOZXFORM NetCDF files
  type, extends(MagneticField) :: BoozXformField
    ! BOOZXFORM data arrays
    integer :: ns_b        !< Number of flux surfaces
    integer :: nfp_b       !< Number of field periods
    integer :: mboz_b      !< Maximum poloidal mode number
    integer :: nboz_b      !< Maximum toroidal mode number (divided by nfp)
    integer :: mnboz       !< Total number of Fourier modes
    
    real(dp) :: aspect_b   !< Aspect ratio
    real(dp) :: rmax_b     !< Maximum R
    real(dp) :: rmin_b     !< Minimum R
    real(dp) :: betaxis_b  !< Beta on axis
    
    ! Radial arrays (dimension ns_b)
    real(dp), allocatable :: s_b(:)      !< Normalized toroidal flux
    real(dp), allocatable :: iota_b(:)   !< Rotational transform
    real(dp), allocatable :: buco_b(:)   !< Boozer I (poloidal covariant)
    real(dp), allocatable :: bvco_b(:)   !< Boozer G (toroidal covariant)
    real(dp), allocatable :: beta_b(:)   !< Plasma beta
    real(dp), allocatable :: phip_b(:)   !< d(phi)/ds
    real(dp), allocatable :: chi_b(:)    !< Poloidal flux
    real(dp), allocatable :: pres_b(:)   !< Pressure
    real(dp), allocatable :: phi_b(:)    !< Toroidal flux
    
    ! Fourier mode arrays
    integer, allocatable :: ixm_b(:)     !< Poloidal mode numbers
    integer, allocatable :: ixn_b(:)     !< Toroidal mode numbers
    integer, allocatable :: jlist(:)     !< Surface indices for packed arrays
    
    ! Fourier coefficient arrays (dimension pack_rad x mnboz)
    real(dp), allocatable :: rmnc_b(:,:)  !< R Fourier coefficients (cos)
    real(dp), allocatable :: zmns_b(:,:)  !< Z Fourier coefficients (sin)
    real(dp), allocatable :: pmns_b(:,:)  !< p = theta_VMEC - theta_Boozer (sin)
    real(dp), allocatable :: gmn_b(:,:)   !< g = zeta_VMEC - zeta_Boozer
    real(dp), allocatable :: bmnc_b(:,:)  !< |B| Fourier coefficients (cos)
    
    ! Stellarator symmetry flag
    logical :: lasym_b
    
  contains
    procedure :: load => load_booz_xform
    procedure :: evaluate => evaluate_booz_xform
    procedure :: cleanup => cleanup_booz_xform
  end type BoozXformField

contains

  !> Load BOOZXFORM data from NetCDF file
  subroutine load_booz_xform(this, filename)
    use nctools_module, only: nc_open, nc_close, nc_get
    use netcdf, only: nf90_inquire_dimension, nf90_inq_dimid, nf90_noerr, &
                      nf90_strerror
    class(BoozXformField), intent(inout) :: this
    character(len=*), intent(in) :: filename
    
    integer :: ncid, lasym_int
    integer :: ns_dim, mn_dim, pack_dim
    integer :: dimid, ierr, i
    
    ! Open NetCDF file
    call nc_open(filename, ncid)
    
    ! Read dimensions using NetCDF directly
    ierr = nf90_inq_dimid(ncid, 'radius', dimid)
    if (ierr == nf90_noerr) then
      ierr = nf90_inquire_dimension(ncid, dimid, len=ns_dim)
    else
      print *, 'Error reading radius dimension:', trim(nf90_strerror(ierr))
      error stop
    end if
    
    ierr = nf90_inq_dimid(ncid, 'mn_modes', dimid)
    if (ierr == nf90_noerr) then
      ierr = nf90_inquire_dimension(ncid, dimid, len=mn_dim)
    else
      print *, 'Error reading mn_modes dimension:', trim(nf90_strerror(ierr))
      error stop
    end if
    
    ierr = nf90_inq_dimid(ncid, 'pack_rad', dimid)
    if (ierr == nf90_noerr) then
      ierr = nf90_inquire_dimension(ncid, dimid, len=pack_dim)
    else
      print *, 'Error reading pack_rad dimension:', trim(nf90_strerror(ierr))
      error stop
    end if
    
    this%ns_b = ns_dim
    this%mnboz = mn_dim
    
    ! Read scalar variables
    call nc_get(ncid, 'ns_b', this%ns_b)
    call nc_get(ncid, 'nfp_b', this%nfp_b)
    call nc_get(ncid, 'mboz_b', this%mboz_b)
    call nc_get(ncid, 'nboz_b', this%nboz_b)
    call nc_get(ncid, 'mnboz_b', this%mnboz)
    call nc_get(ncid, 'aspect_b', this%aspect_b)
    call nc_get(ncid, 'rmax_b', this%rmax_b)
    call nc_get(ncid, 'rmin_b', this%rmin_b)
    call nc_get(ncid, 'betaxis_b', this%betaxis_b)
    
    ! Check stellarator symmetry
    call nc_get(ncid, 'lasym__logical__', lasym_int)
    this%lasym_b = (lasym_int == 1)
    
    ! Allocate arrays
    allocate(this%s_b(ns_dim))
    allocate(this%iota_b(ns_dim))
    allocate(this%buco_b(ns_dim))
    allocate(this%bvco_b(ns_dim))
    allocate(this%beta_b(ns_dim))
    allocate(this%phip_b(ns_dim))
    allocate(this%chi_b(ns_dim))
    allocate(this%pres_b(ns_dim))
    allocate(this%phi_b(ns_dim))
    
    allocate(this%ixm_b(mn_dim))
    allocate(this%ixn_b(mn_dim))
    allocate(this%jlist(pack_dim))
    
    allocate(this%rmnc_b(pack_dim, mn_dim))
    allocate(this%zmns_b(pack_dim, mn_dim))
    allocate(this%pmns_b(pack_dim, mn_dim))
    allocate(this%gmn_b(pack_dim, mn_dim))
    allocate(this%bmnc_b(pack_dim, mn_dim))
    
    ! Read radial arrays
    call nc_get(ncid, 'iota_b', this%iota_b)
    call nc_get(ncid, 'buco_b', this%buco_b)
    call nc_get(ncid, 'bvco_b', this%bvco_b)
    call nc_get(ncid, 'beta_b', this%beta_b)
    call nc_get(ncid, 'phip_b', this%phip_b)
    call nc_get(ncid, 'chi_b', this%chi_b)
    call nc_get(ncid, 'pres_b', this%pres_b)
    call nc_get(ncid, 'phi_b', this%phi_b)
    
    ! Create s array (normalized toroidal flux)
    this%s_b = [(real(i-1, dp) / real(ns_dim-1, dp), i = 1, ns_dim)]
    
    ! Read mode arrays
    call nc_get(ncid, 'ixm_b', this%ixm_b)
    call nc_get(ncid, 'ixn_b', this%ixn_b)
    call nc_get(ncid, 'jlist', this%jlist)
    
    ! Read Fourier coefficients
    call nc_get(ncid, 'rmnc_b', this%rmnc_b)
    call nc_get(ncid, 'zmns_b', this%zmns_b)
    call nc_get(ncid, 'pmns_b', this%pmns_b)
    call nc_get(ncid, 'gmn_b', this%gmn_b)
    call nc_get(ncid, 'bmnc_b', this%bmnc_b)
    
    ! Close file
    call nc_close(ncid)
    
    print *, 'Loaded BOOZXFORM file:', trim(filename)
    print *, '  ns =', this%ns_b, ', nfp =', this%nfp_b
    print *, '  mboz =', this%mboz_b, ', nboz =', this%nboz_b
    print *, '  mnboz =', this%mnboz, ' Fourier modes'
    
  end subroutine load_booz_xform

  !> Evaluate magnetic field at given coordinates
  !> Note: This expects VMEC coordinates but we have Boozer data
  !> Need to coordinate transform or integrate differently
  subroutine evaluate_booz_xform(self, x, Acov, hcov, Bmod, sqgBctr)
    class(BoozXformField), intent(in) :: self
    real(dp), intent(in) :: x(3)  ! r=sqrt(s_vmec), theta_vmec, phi_vmec
    real(dp), intent(out) :: Acov(3)
    real(dp), intent(out) :: hcov(3)
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)
    
    ! TODO: This is a placeholder implementation
    ! The challenge is that the base interface expects VMEC coordinates
    ! but we have data in Boozer coordinates
    ! Need to either:
    ! 1. Transform from VMEC to Boozer coordinates first
    ! 2. Store transformation data in the BOOZXFORM file
    ! 3. Use a different integration approach
    
    ! For now, just return dummy values
    Acov = 0.0_dp
    hcov = [0.0_dp, 1.0_dp, 1.0_dp]  ! Dummy covariant B components
    Bmod = 1.0_dp  ! Dummy |B|
    
    if (present(sqgBctr)) then
      sqgBctr = 0.0_dp
    end if
    
    ! Print warning once
    if (.not. booz_warning_printed) then
      print *, 'WARNING: BoozXformField evaluate not fully implemented yet'
      booz_warning_printed = .true.
    end if
    
  end subroutine evaluate_booz_xform

  !> Clean up allocated arrays
  subroutine cleanup_booz_xform(self)
    class(BoozXformField), intent(inout) :: self
    
    if (allocated(self%s_b)) deallocate(self%s_b)
    if (allocated(self%iota_b)) deallocate(self%iota_b)
    if (allocated(self%buco_b)) deallocate(self%buco_b)
    if (allocated(self%bvco_b)) deallocate(self%bvco_b)
    if (allocated(self%beta_b)) deallocate(self%beta_b)
    if (allocated(self%phip_b)) deallocate(self%phip_b)
    if (allocated(self%chi_b)) deallocate(self%chi_b)
    if (allocated(self%pres_b)) deallocate(self%pres_b)
    if (allocated(self%phi_b)) deallocate(self%phi_b)
    if (allocated(self%ixm_b)) deallocate(self%ixm_b)
    if (allocated(self%ixn_b)) deallocate(self%ixn_b)
    if (allocated(self%jlist)) deallocate(self%jlist)
    if (allocated(self%rmnc_b)) deallocate(self%rmnc_b)
    if (allocated(self%zmns_b)) deallocate(self%zmns_b)
    if (allocated(self%pmns_b)) deallocate(self%pmns_b)
    if (allocated(self%gmn_b)) deallocate(self%gmn_b)
    if (allocated(self%bmnc_b)) deallocate(self%bmnc_b)
    
  end subroutine cleanup_booz_xform

end module field_booz_xform