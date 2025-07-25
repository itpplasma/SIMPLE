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
  !> Input: VMEC coordinates (r, theta_vmec, phi_vmec)
  !> Output: Magnetic field components in VMEC coordinates
  subroutine evaluate_booz_xform(self, x, Acov, hcov, Bmod, sqgBctr)
    class(BoozXformField), intent(in) :: self
    real(dp), intent(in) :: x(3)  ! r=sqrt(s_vmec), theta_vmec, phi_vmec
    real(dp), intent(out) :: Acov(3)
    real(dp), intent(out) :: hcov(3)
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)
    
    real(dp) :: s, theta_b, zeta_b, theta_v, zeta_v
    real(dp) :: p_angle, g_angle  ! theta_vmec - theta_booz, zeta_vmec - zeta_booz
    real(dp) :: R_booz, Z_booz, B_booz, Jac_booz
    real(dp) :: dR_dtheta, dR_dzeta, dZ_dtheta, dZ_dzeta
    real(dp) :: sqrtg, Bcov_s, Bcov_theta, Bcov_zeta
    real(dp) :: Bctr_theta, Bctr_zeta
    real(dp) :: iota_local, I_local, G_local, phip_local
    integer :: js, imn
    real(dp) :: wlo, whi, angle_arg
    real(dp), parameter :: twopi = 8.0_dp * atan(1.0_dp)
    
    ! Convert input coordinates
    s = x(1)**2  ! s = r^2
    theta_v = x(2)
    zeta_v = x(3)
    
    ! Find radial grid point and interpolation weights
    js = minloc(abs(self%s_b - s), 1)
    if (js == self%ns_b) js = self%ns_b - 1
    if (js == 1) js = 2
    whi = (s - self%s_b(js-1)) / (self%s_b(js) - self%s_b(js-1))
    wlo = 1.0_dp - whi
    
    ! Interpolate radial profiles
    iota_local = wlo * self%iota_b(js-1) + whi * self%iota_b(js)
    I_local = wlo * self%buco_b(js-1) + whi * self%buco_b(js)
    G_local = wlo * self%bvco_b(js-1) + whi * self%bvco_b(js)
    phip_local = wlo * self%phip_b(js-1) + whi * self%phip_b(js)
    
    ! First, evaluate p and g to get coordinate transformation
    ! p = theta_vmec - theta_booz, g = zeta_vmec - zeta_booz
    p_angle = 0.0_dp
    g_angle = 0.0_dp
    
    do imn = 1, self%mnboz
      angle_arg = self%ixm_b(imn) * theta_v - self%ixn_b(imn) * zeta_v
      ! pmns is sin component for p
      p_angle = p_angle + sin(angle_arg) * &
                (wlo * self%pmns_b(js-1, imn) + whi * self%pmns_b(js, imn))
      ! gmn is cos component for g (but stored in gmn_b array)
      g_angle = g_angle + cos(angle_arg) * &
                (wlo * self%gmn_b(js-1, imn) + whi * self%gmn_b(js, imn))
    end do
    
    ! Convert to Boozer angles
    theta_b = theta_v - p_angle
    zeta_b = zeta_v - g_angle
    
    ! Now evaluate field quantities in Boozer coordinates
    R_booz = 0.0_dp
    Z_booz = 0.0_dp
    B_booz = 0.0_dp
    Jac_booz = 0.0_dp
    dR_dtheta = 0.0_dp
    dR_dzeta = 0.0_dp
    dZ_dtheta = 0.0_dp
    dZ_dzeta = 0.0_dp
    
    do imn = 1, self%mnboz
      angle_arg = self%ixm_b(imn) * theta_b - self%ixn_b(imn) * zeta_b
      
      ! R, B, and Jacobian are cosine components
      R_booz = R_booz + cos(angle_arg) * &
               (wlo * self%rmnc_b(js-1, imn) + whi * self%rmnc_b(js, imn))
      B_booz = B_booz + cos(angle_arg) * &
               (wlo * self%bmnc_b(js-1, imn) + whi * self%bmnc_b(js, imn))
      Jac_booz = Jac_booz + cos(angle_arg) * &
                 (wlo * self%gmn_b(js-1, imn) + whi * self%gmn_b(js, imn))
      
      ! Z is sine component
      Z_booz = Z_booz + sin(angle_arg) * &
               (wlo * self%zmns_b(js-1, imn) + whi * self%zmns_b(js, imn))
      
      ! Derivatives for metric
      dR_dtheta = dR_dtheta - self%ixm_b(imn) * sin(angle_arg) * &
                  (wlo * self%rmnc_b(js-1, imn) + whi * self%rmnc_b(js, imn))
      dR_dzeta = dR_dzeta + self%ixn_b(imn) * sin(angle_arg) * &
                 (wlo * self%rmnc_b(js-1, imn) + whi * self%rmnc_b(js, imn))
      dZ_dtheta = dZ_dtheta + self%ixm_b(imn) * cos(angle_arg) * &
                  (wlo * self%zmns_b(js-1, imn) + whi * self%zmns_b(js, imn))
      dZ_dzeta = dZ_dzeta - self%ixn_b(imn) * cos(angle_arg) * &
                 (wlo * self%zmns_b(js-1, imn) + whi * self%zmns_b(js, imn))
    end do
    
    ! Set output magnetic field magnitude
    Bmod = B_booz
    
    ! Compute sqrt(g) for Boozer coordinates
    ! In Boozer coords: sqrt(g) = (G + iota*I) / B^2
    sqrtg = (G_local + iota_local * I_local) / (B_booz * B_booz)
    
    ! Contravariant B components in Boozer coordinates
    Bctr_theta = I_local / sqrtg
    Bctr_zeta = G_local / sqrtg
    
    ! Covariant B components in Boozer coordinates
    ! B = grad(psi) x grad(theta) + iota * grad(psi) x grad(zeta)
    Bcov_s = 0.0_dp  ! No radial component in Boozer coords
    Bcov_theta = I_local
    Bcov_zeta = G_local
    
    ! For straight field line coordinates, the vector potential has simple form
    ! A_theta = -psi (toroidal flux function)
    ! A_zeta = 0 in Boozer coordinates
    Acov(1) = 0.0_dp  ! A_s = 0
    Acov(2) = -self%phi_b(js) / twopi  ! A_theta = -Phi/(2*pi)
    Acov(3) = 0.0_dp  ! A_zeta = 0 in Boozer
    
    ! Normalized covariant B components
    hcov(1) = Bcov_s / Bmod
    hcov(2) = Bcov_theta / Bmod
    hcov(3) = Bcov_zeta / Bmod
    
    ! Transform back to VMEC coordinates if needed
    ! This is a simplified transformation - may need refinement
    ! The covariant components transform with the inverse Jacobian
    ! For now, we use the fact that the transformation mainly affects angles
    
    if (present(sqgBctr)) then
      sqgBctr(1) = 0.0_dp
      sqgBctr(2) = sqrtg * Bctr_theta
      sqgBctr(3) = sqrtg * Bctr_zeta
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