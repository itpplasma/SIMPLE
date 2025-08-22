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
    use netcdf
    class(BoozXformField), intent(inout) :: this
    character(len=*), intent(in) :: filename
    
    integer :: ncid, lasym_int
    integer :: ns_dim, mn_dim, pack_dim
    integer :: dimid, ierr, i, varid
    
    ! Open NetCDF file using native NetCDF
    ierr = nf90_open(filename, NF90_NOWRITE, ncid)
    if (ierr /= nf90_noerr) then
      print *, 'Error opening file:', trim(filename)
      print *, 'Error:', trim(nf90_strerror(ierr))
      error stop
    end if
    print *, 'Opened NetCDF file:', trim(filename)
    
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
    
    ! Read scalar variables using native NetCDF
    print *, 'Reading scalar variables...'
    
    ierr = nf90_inq_varid(ncid, 'nfp_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%nfp_b)
    
    ierr = nf90_inq_varid(ncid, 'mboz_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%mboz_b)
    
    ierr = nf90_inq_varid(ncid, 'nboz_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%nboz_b)
    
    ierr = nf90_inq_varid(ncid, 'aspect_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%aspect_b)
    
    ierr = nf90_inq_varid(ncid, 'rmax_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%rmax_b)
    
    ierr = nf90_inq_varid(ncid, 'rmin_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%rmin_b)
    
    ierr = nf90_inq_varid(ncid, 'betaxis_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%betaxis_b)
    
    ! Check stellarator symmetry
    ierr = nf90_inq_varid(ncid, 'lasym__logical__', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, lasym_int)
      this%lasym_b = (lasym_int == 1)
    end if
    
    ! Allocate arrays
    print *, 'Allocating arrays: ns_dim =', ns_dim, 'mn_dim =', mn_dim, 'pack_dim =', pack_dim
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
    print *, 'Reading radial arrays...'
    
    ierr = nf90_inq_varid(ncid, 'iota_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%iota_b)
      if (ierr /= nf90_noerr) print *, 'Error reading iota_b:', trim(nf90_strerror(ierr))
    end if
    
    ierr = nf90_inq_varid(ncid, 'buco_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%buco_b)
      if (ierr /= nf90_noerr) print *, 'Error reading buco_b:', trim(nf90_strerror(ierr))
    end if
    
    ierr = nf90_inq_varid(ncid, 'bvco_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%bvco_b)
      if (ierr /= nf90_noerr) print *, 'Error reading bvco_b:', trim(nf90_strerror(ierr))
    end if
    
    ! Read remaining arrays
    ierr = nf90_inq_varid(ncid, 'beta_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%beta_b)
    
    ierr = nf90_inq_varid(ncid, 'phip_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%phip_b)
    
    ierr = nf90_inq_varid(ncid, 'chi_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%chi_b)
    
    ierr = nf90_inq_varid(ncid, 'pres_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%pres_b)
    
    ierr = nf90_inq_varid(ncid, 'phi_b', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, this%phi_b)
    
    ! Create s array (normalized toroidal flux)
    this%s_b = [(real(i-1, dp) / real(ns_dim-1, dp), i = 1, ns_dim)]
    
    ! Read mode arrays and Fourier coefficients
    print *, 'Reading mode arrays and Fourier coefficients...'
    
    ierr = nf90_inq_varid(ncid, 'ixm_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%ixm_b)
      if (ierr /= nf90_noerr) print *, 'Error reading ixm_b:', trim(nf90_strerror(ierr))
    end if
    
    ierr = nf90_inq_varid(ncid, 'ixn_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%ixn_b)
      if (ierr /= nf90_noerr) print *, 'Error reading ixn_b:', trim(nf90_strerror(ierr))
    end if
    
    ierr = nf90_inq_varid(ncid, 'jlist', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%jlist)
      if (ierr /= nf90_noerr) print *, 'Error reading jlist:', trim(nf90_strerror(ierr))
    end if
    
    ! Read Fourier coefficients
    ierr = nf90_inq_varid(ncid, 'rmnc_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%rmnc_b)
      if (ierr /= nf90_noerr) print *, 'Error reading rmnc_b:', trim(nf90_strerror(ierr))
    end if
    
    ierr = nf90_inq_varid(ncid, 'zmns_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%zmns_b)
      if (ierr /= nf90_noerr) print *, 'Error reading zmns_b:', trim(nf90_strerror(ierr))
    end if
    
    ierr = nf90_inq_varid(ncid, 'pmns_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%pmns_b)
      if (ierr /= nf90_noerr) print *, 'Error reading pmns_b:', trim(nf90_strerror(ierr))
    end if
    
    ierr = nf90_inq_varid(ncid, 'gmn_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%gmn_b)
      if (ierr /= nf90_noerr) print *, 'Error reading gmn_b:', trim(nf90_strerror(ierr))
    end if
    
    ierr = nf90_inq_varid(ncid, 'bmnc_b', varid)
    if (ierr == nf90_noerr) then
      ierr = nf90_get_var(ncid, varid, this%bmnc_b)
      if (ierr /= nf90_noerr) print *, 'Error reading bmnc_b:', trim(nf90_strerror(ierr))
    end if
    
    ! Close file
    ierr = nf90_close(ncid)
    if (ierr /= nf90_noerr) then
      print *, 'Error closing file:', trim(nf90_strerror(ierr))
    end if
    
    print *, 'Loaded BOOZXFORM file:', trim(filename)
    print *, '  ns =', this%ns_b, ', nfp =', this%nfp_b
    print *, '  mboz =', this%mboz_b, ', nboz =', this%nboz_b
    print *, '  mnboz =', this%mnboz, ' Fourier modes'
    print *, '  iota_b(1) =', this%iota_b(1), 'iota_b(ns) =', this%iota_b(this%ns_b)
    print *, '  First few mode numbers:'
    print *, '    m =', this%ixm_b(1:min(5,size(this%ixm_b)))
    print *, '    n =', this%ixn_b(1:min(5,size(this%ixn_b)))
    
  end subroutine load_booz_xform

  !> Evaluate magnetic field at given coordinates
  !> Input: VMEC coordinates (r, theta_vmec, phi_vmec)
  !> Output: Magnetic field components
  !> Note: This is a simplified implementation that evaluates B directly in Boozer coordinates
  !> without the full coordinate transformation
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
    integer :: js, imn, js_pack, js_pack_m1
    real(dp) :: wlo, whi, angle_arg
    real(dp), parameter :: twopi = 8.0_dp * atan(1.0_dp)
    
    print *, 'evaluate_booz_xform: Entry, x =', x
    print *, 'evaluate_booz_xform: ns_b =', self%ns_b, 'size(s_b) =', size(self%s_b)
    print *, 'evaluate_booz_xform: size(bmnc_b) =', size(self%bmnc_b, 1), size(self%bmnc_b, 2)
    
    ! Check if data is loaded
    if (.not. allocated(self%s_b)) then
      print *, 'ERROR: s_b not allocated!'
      error stop 'evaluate_booz_xform: BOOZXFORM data not loaded'
    end if
    
    ! Convert input coordinates
    s = x(1)**2  ! s = r^2
    theta_v = x(2)
    zeta_v = x(3)
    
    print *, 'evaluate_booz_xform: s =', s
    
    ! Find radial grid point and interpolation weights
    js = minloc(abs(self%s_b - s), 1)
    if (js == self%ns_b) js = self%ns_b - 1
    if (js < 2) js = 2
    
    ! For simplicity, find the closest packed indices
    ! Since jlist goes from 2 to ns_b, and pack_rad = ns_b-1, we can map:
    js_pack = js - 1      ! Map surface index to packed index (2->1, 3->2, etc.)
    js_pack_m1 = js - 2   ! Previous surface
    
    ! Ensure we're within bounds
    if (js_pack > size(self%bmnc_b, 1)) js_pack = size(self%bmnc_b, 1)
    if (js_pack < 1) js_pack = 1
    if (js_pack_m1 < 1) js_pack_m1 = 1
    
    ! If we're at the same index, use neighboring indices
    if (js_pack == js_pack_m1) then
      if (js_pack > 1) then
        js_pack_m1 = js_pack - 1
      else
        js_pack = 2
        js_pack_m1 = 1
      end if
    end if
    
    whi = (s - self%s_b(js-1)) / max(1.0e-10_dp, self%s_b(js) - self%s_b(js-1))
    wlo = 1.0_dp - whi
    
    ! Interpolate radial profiles
    iota_local = wlo * self%iota_b(js-1) + whi * self%iota_b(js)
    I_local = wlo * self%buco_b(js-1) + whi * self%buco_b(js)
    G_local = wlo * self%bvco_b(js-1) + whi * self%bvco_b(js)
    phip_local = wlo * self%phip_b(js-1) + whi * self%phip_b(js)
    
    ! For now, use VMEC angles directly as Boozer angles (simplified)
    ! This is not the full transformation but should allow basic evaluation
    theta_b = theta_v
    zeta_b = zeta_v
    
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
      
      ! Check array bounds
      if (js_pack_m1 < 1 .or. js_pack_m1 > size(self%rmnc_b, 1) .or. &
          js_pack < 1 .or. js_pack > size(self%rmnc_b, 1)) cycle
      
      ! R, B, and Jacobian are cosine components
      R_booz = R_booz + cos(angle_arg) * &
               (wlo * self%rmnc_b(js_pack_m1, imn) + whi * self%rmnc_b(js_pack, imn))
      B_booz = B_booz + cos(angle_arg) * &
               (wlo * self%bmnc_b(js_pack_m1, imn) + whi * self%bmnc_b(js_pack, imn))
      Jac_booz = Jac_booz + cos(angle_arg) * &
                 (wlo * self%gmn_b(js_pack_m1, imn) + whi * self%gmn_b(js_pack, imn))
      
      ! Z is sine component
      Z_booz = Z_booz + sin(angle_arg) * &
               (wlo * self%zmns_b(js_pack_m1, imn) + whi * self%zmns_b(js_pack, imn))
      
      ! Derivatives for metric
      dR_dtheta = dR_dtheta - self%ixm_b(imn) * sin(angle_arg) * &
                  (wlo * self%rmnc_b(js_pack_m1, imn) + whi * self%rmnc_b(js_pack, imn))
      dR_dzeta = dR_dzeta + self%ixn_b(imn) * sin(angle_arg) * &
                 (wlo * self%rmnc_b(js_pack_m1, imn) + whi * self%rmnc_b(js_pack, imn))
      dZ_dtheta = dZ_dtheta + self%ixm_b(imn) * cos(angle_arg) * &
                  (wlo * self%zmns_b(js_pack_m1, imn) + whi * self%zmns_b(js_pack, imn))
      dZ_dzeta = dZ_dzeta - self%ixn_b(imn) * cos(angle_arg) * &
                 (wlo * self%zmns_b(js_pack_m1, imn) + whi * self%zmns_b(js_pack, imn))
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
    Acov(2) = -(wlo * self%phi_b(js-1) + whi * self%phi_b(js)) / twopi  ! A_theta = -Phi/(2*pi)
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