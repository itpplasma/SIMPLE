!> Module for reading and evaluating magnetic fields from BOOZXFORM output files
!> This allows direct input of pre-computed Boozer coordinate fields
module field_booz_xform
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_base, only: MagneticField
  implicit none

  private
  public :: BoozXformField

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
    use netcdf_utils, only: nc_open_r, nc_close, nc_get_dim, nc_get_var
    class(BoozXformField), intent(inout) :: this
    character(len=*), intent(in) :: filename
    
    integer :: ncid, lasym_int
    integer :: ns_dim, mn_dim, pack_dim
    
    ! Open NetCDF file
    call nc_open_r(filename, ncid)
    
    ! Read dimensions
    call nc_get_dim(ncid, 'radius', ns_dim)
    call nc_get_dim(ncid, 'mn_modes', mn_dim)
    call nc_get_dim(ncid, 'pack_rad', pack_dim)
    
    this%ns_b = ns_dim
    this%mnboz = mn_dim
    
    ! Read scalar variables
    call nc_get_var(ncid, 'ns_b', this%ns_b)
    call nc_get_var(ncid, 'nfp_b', this%nfp_b)
    call nc_get_var(ncid, 'mboz_b', this%mboz_b)
    call nc_get_var(ncid, 'nboz_b', this%nboz_b)
    call nc_get_var(ncid, 'mnboz_b', this%mnboz)
    call nc_get_var(ncid, 'aspect_b', this%aspect_b)
    call nc_get_var(ncid, 'rmax_b', this%rmax_b)
    call nc_get_var(ncid, 'rmin_b', this%rmin_b)
    call nc_get_var(ncid, 'betaxis_b', this%betaxis_b)
    
    ! Check stellarator symmetry
    call nc_get_var(ncid, 'lasym__logical__', lasym_int)
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
    call nc_get_var(ncid, 'iota_b', this%iota_b)
    call nc_get_var(ncid, 'buco_b', this%buco_b)
    call nc_get_var(ncid, 'bvco_b', this%bvco_b)
    call nc_get_var(ncid, 'beta_b', this%beta_b)
    call nc_get_var(ncid, 'phip_b', this%phip_b)
    call nc_get_var(ncid, 'chi_b', this%chi_b)
    call nc_get_var(ncid, 'pres_b', this%pres_b)
    call nc_get_var(ncid, 'phi_b', this%phi_b)
    
    ! Create s array (normalized toroidal flux)
    this%s_b = [(real(i-1, dp) / real(ns_dim-1, dp), i = 1, ns_dim)]
    
    ! Read mode arrays
    call nc_get_var(ncid, 'ixm_b', this%ixm_b)
    call nc_get_var(ncid, 'ixn_b', this%ixn_b)
    call nc_get_var(ncid, 'jlist', this%jlist)
    
    ! Read Fourier coefficients
    call nc_get_var(ncid, 'rmnc_b', this%rmnc_b)
    call nc_get_var(ncid, 'zmns_b', this%zmns_b)
    call nc_get_var(ncid, 'pmns_b', this%pmns_b)
    call nc_get_var(ncid, 'gmn_b', this%gmn_b)
    call nc_get_var(ncid, 'bmnc_b', this%bmnc_b)
    
    ! Close file
    call nc_close(ncid)
    
    print *, 'Loaded BOOZXFORM file:', trim(filename)
    print *, '  ns =', this%ns_b, ', nfp =', this%nfp_b
    print *, '  mboz =', this%mboz_b, ', nboz =', this%nboz_b
    print *, '  mnboz =', this%mnboz, ' Fourier modes'
    
  end subroutine load_booz_xform

  !> Evaluate magnetic field at given Boozer coordinates
  subroutine evaluate_booz_xform(this, s, theta_b, zeta_b, B_r, B_theta, B_zeta, &
                                  B_s, B_theta_s, B_zeta_s, dBds, dBdtheta, dBdzeta)
    class(BoozXformField), intent(in) :: this
    real(dp), intent(in) :: s          !< Normalized toroidal flux
    real(dp), intent(in) :: theta_b    !< Boozer poloidal angle
    real(dp), intent(in) :: zeta_b     !< Boozer toroidal angle
    real(dp), intent(out) :: B_r, B_theta, B_zeta
    real(dp), intent(out), optional :: B_s, B_theta_s, B_zeta_s
    real(dp), intent(out), optional :: dBds, dBdtheta, dBdzeta
    
    real(dp) :: B_mag, sqrtg
    real(dp) :: iota_local, buco_local, bvco_local
    integer :: js, j_idx, k
    real(dp) :: ds, w1, w2
    real(dp) :: cos_arg, sin_arg
    real(dp) :: angle_arg
    integer :: m, n
    
    ! Find radial index for interpolation
    js = min(max(1, int(s * (this%ns_b - 1)) + 1), this%ns_b - 1)
    ds = s - this%s_b(js)
    w1 = 1.0_dp - ds * (this%ns_b - 1)
    w2 = ds * (this%ns_b - 1)
    
    ! Interpolate radial quantities
    iota_local = w1 * this%iota_b(js) + w2 * this%iota_b(js+1)
    buco_local = w1 * this%buco_b(js) + w2 * this%buco_b(js+1)
    bvco_local = w1 * this%bvco_b(js) + w2 * this%bvco_b(js+1)
    
    ! Get packed index for Fourier arrays
    if (js <= size(this%jlist)) then
      j_idx = this%jlist(js)
    else
      j_idx = size(this%jlist)
    end if
    
    ! Evaluate |B| from Fourier series
    B_mag = 0.0_dp
    do k = 1, this%mnboz
      m = this%ixm_b(k)
      n = this%ixn_b(k) 
      angle_arg = m * theta_b - n * zeta_b
      cos_arg = cos(angle_arg)
      
      if (j_idx > 0 .and. j_idx <= size(this%bmnc_b, 1)) then
        B_mag = B_mag + this%bmnc_b(j_idx, k) * cos_arg
      end if
    end do
    
    ! For Boozer coordinates, the field components are:
    ! B^s = 0 (by definition of Boozer coordinates)
    ! B^theta = I / sqrt(g)
    ! B^zeta = G / sqrt(g)
    ! where sqrt(g) = (G + iota*I) / B
    
    sqrtg = (bvco_local + iota_local * buco_local) / B_mag
    
    ! Contravariant components
    B_r = 0.0_dp  ! B^s = 0 in Boozer coordinates
    B_theta = buco_local / sqrtg
    B_zeta = bvco_local / sqrtg
    
    ! Optional: covariant components
    if (present(B_s)) B_s = 0.0_dp
    if (present(B_theta_s)) B_theta_s = buco_local
    if (present(B_zeta_s)) B_zeta_s = bvco_local
    
    ! Optional: derivatives (simplified - would need more complete implementation)
    if (present(dBds)) dBds = 0.0_dp
    if (present(dBdtheta)) dBdtheta = 0.0_dp
    if (present(dBdzeta)) dBdzeta = 0.0_dp
    
  end subroutine evaluate_booz_xform

  !> Clean up allocated arrays
  subroutine cleanup_booz_xform(this)
    class(BoozXformField), intent(inout) :: this
    
    if (allocated(this%s_b)) deallocate(this%s_b)
    if (allocated(this%iota_b)) deallocate(this%iota_b)
    if (allocated(this%buco_b)) deallocate(this%buco_b)
    if (allocated(this%bvco_b)) deallocate(this%bvco_b)
    if (allocated(this%beta_b)) deallocate(this%beta_b)
    if (allocated(this%phip_b)) deallocate(this%phip_b)
    if (allocated(this%chi_b)) deallocate(this%chi_b)
    if (allocated(this%pres_b)) deallocate(this%pres_b)
    if (allocated(this%phi_b)) deallocate(this%phi_b)
    if (allocated(this%ixm_b)) deallocate(this%ixm_b)
    if (allocated(this%ixn_b)) deallocate(this%ixn_b)
    if (allocated(this%jlist)) deallocate(this%jlist)
    if (allocated(this%rmnc_b)) deallocate(this%rmnc_b)
    if (allocated(this%zmns_b)) deallocate(this%zmns_b)
    if (allocated(this%pmns_b)) deallocate(this%pmns_b)
    if (allocated(this%gmn_b)) deallocate(this%gmn_b)
    if (allocated(this%bmnc_b)) deallocate(this%bmnc_b)
    
  end subroutine cleanup_booz_xform

end module field_booz_xform