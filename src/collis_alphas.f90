  !
module collis_alp
  implicit none

  ! Define real(wp) kind parameter
  integer, parameter :: wp = kind(1.0d0)
  integer, parameter :: nsorts=3
  real(wp), dimension(nsorts) :: efcolf,velrat,enrat
  
  ! Physical constants
  real(wp), parameter :: pi = 3.14159265358979d0
  real(wp), parameter :: pmass = 1.6726d-24    ! Proton mass [g]
  real(wp), parameter :: emass = 9.1094d-28    ! Electron mass [g]
  real(wp), parameter :: e_charge = 4.8032d-10 ! Elementary charge [esu]
  real(wp), parameter :: ev_to_erg = 1.6022d-12 ! eV to erg conversion
  
  ! Collision constants
  real(wp), parameter :: sqrt_pi = 1.7724538d0
  real(wp), parameter :: cons_onsager = 0.75225278d0  ! 4/(3*sqrt(pi))
  real(wp), parameter :: pmin = 1.e-8  ! Minimum momentum module

  contains
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine coleff(p,dpp,dhh,fpeff)
  !
  !  Computes local values of dimensionless contravariant components
  !  of collisional diffusion tensor and friction force for nonrelativistic
  !  plasma. Backgound temperature is the same for all sorts.
  !
  !     Input variables:
  !        formal: p      - dimensionless momentum module (p/(sqrt(2)*p_T)
  !        common: efcolf - dmls collision frequencies
  !                velrat - ratio of test species thermal velocity to
  !                         background species thermal velocity
  !     Output variables:
  !        formal: dpp    - dimensionless momentum module diffusion
  !                         coefficient
  !                dhh    - dimensionless pitch angle diffusion coeff.
  !                fpeff  - effective dimensionless drag force (prop. to linear
  !                         deviation in Fokker-Planck eq.)
  !
  real(wp), intent(in) :: p
  real(wp), intent(out) :: dpp, dhh, fpeff
  
  real(wp) :: plim
  integer :: i
  
  plim = max(p, 1.d-8)
  
  ! Initialize outputs
  dpp = 0.0d0
  dhh = 0.0d0
  fpeff = 0.0d0
  
  ! Sum contributions from all species
  do i = 1, nsorts
    call add_species_contribution(p, plim, i, dpp, dhh, fpeff)
  end do
  
  ! Normalize pitch angle diffusion
  dhh = dhh / plim**2
  
  end subroutine coleff
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper subroutine to add contribution from a single species
  subroutine add_species_contribution(p, plim, species_idx, dpp, dhh, fpeff)
    real(wp), intent(in) :: p, plim
    integer, intent(in) :: species_idx
    real(wp), intent(inout) :: dpp, dhh, fpeff
    
    real(wp) :: xbeta, dp, dh, dpd
    
    xbeta = p * velrat(species_idx)
    call onseff(xbeta, dp, dh, dpd)
    
    dpp = dpp + dp * efcolf(species_idx)
    dhh = dhh + dh * efcolf(species_idx)
    fpeff = fpeff + (dpd/plim - 2.0*dp*p*enrat(species_idx)) * efcolf(species_idx)
  
  end subroutine add_species_contribution
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine onseff(v,dp,dh,dpd)
  !
  !  dp - dimensionless dpp
  !  dh - dhh*p^2     (p - dmls)
  !  dpd - (1/p)(d/dp)p^2*dp   (p - dmls)
  !
  real(wp), intent(in) :: v
  real(wp), intent(out) :: dp, dh, dpd
  
  real(wp) :: v2, v3
  
  v2 = v**2
  v3 = v2 * v
  
  if (v < 0.01d0) then
    ! Small velocity limit
    call onseff_small_v(v2, dp, dh, dpd)
  elseif (v > 6.d0) then
    ! Large velocity limit
    call onseff_large_v(v2, v3, v, dp, dh, dpd)
  else
    ! Intermediate velocity range
    call onseff_intermediate_v(v, v2, v3, dp, dh, dpd)
  endif
  
  end subroutine onseff
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper for small velocity limit
  subroutine onseff_small_v(v2, dp, dh, dpd)
    real(wp), intent(in) :: v2
    real(wp), intent(out) :: dp, dh, dpd
    
    dp = cons_onsager * (1.d0 - 0.6d0*v2)
    dh = cons_onsager * (1.d0 - 0.2d0*v2)
    dpd = 2.d0 * cons_onsager * (1.d0 - 1.2d0*v2)
  
  end subroutine onseff_small_v
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper for large velocity limit
  subroutine onseff_large_v(v2, v3, v, dp, dh, dpd)
    real(wp), intent(in) :: v2, v3, v
    real(wp), intent(out) :: dp, dh, dpd
    
    dp = 1.d0 / v3
    dh = (1.d0 - 0.5d0/v2) / v
    dpd = -1.d0 / v3
  
  end subroutine onseff_large_v
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper for intermediate velocity range
  subroutine onseff_intermediate_v(v, v2, v3, dp, dh, dpd)
    real(wp), intent(in) :: v, v2, v3
    real(wp), intent(out) :: dp, dh, dpd
    
    real(wp) :: ex, er
    
    ex = exp(-v2) / sqrt_pi
    er = erf(v)
    dp = er/v3 - 2.d0*ex/v2
    dh = er*(1.d0 - 0.5d0/v2)/v + ex/v2
    dpd = 4.d0*ex - dp
  
  end subroutine onseff_intermediate_v
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      FUNCTION ERF(X)
  !      PARAMETER  ( A1 = 0.07052 30784, A2 = 0.04228 20123,
  !     ,             A3 = 0.00927 05272, A4 = 0.00015 10143,
  !     ,             A5 = 0.00027 65672, A6 = 0.00004 30638 )
  !      F(T) = 1./((1.+T*(A1+T*(A2+T*(A3+T*(A4+T*(A5+T*A6))))))**4)**4
  !      W = 1. - F(ABS(X))
  !      ERF = SIGN(W,X)
  !      RETURN
  !      END
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha, &
                          v0,dchichi,slowrate,dchichi_norm,slowrate_norm)
  !
  !   Performs precomputation of the constants for Coulomb collision
  !   operator for alpha-particles colliding with 2 sorts of ions and electrons
  !
  real(wp), intent(in) :: am1, am2, Z1, Z2, densi1, densi2
  real(wp), intent(in) :: tempi1, tempi2, tempe, ealpha
  real(wp), intent(out) :: v0, dchichi, slowrate, dchichi_norm, slowrate_norm
  
  real(wp) :: dense, p, dpp, dhh, fpeff
  
  ! Compute energy ratios
  call compute_energy_ratios(ealpha, tempi1, tempi2, tempe)
  
  ! Compute velocities and velocity ratios
  call compute_velocities(ealpha, tempi1, tempi2, tempe, am1, am2, v0)
  
  ! Compute collision frequencies
  dense = densi1*Z1 + densi2*Z2
  call compute_collision_frequencies(am1, am2, Z1, Z2, densi1, densi2, &
                                     tempi1, tempi2, tempe, ealpha, dense, v0)
  
  ! Calculate collision effects at reference momentum
  p = 1.d0
  call coleff(p, dpp, dhh, fpeff)
  
  ! Set output quantities
  dchichi = dhh * v0
  slowrate = abs(fpeff) * v0
  dchichi_norm = dhh
  slowrate_norm = abs(fpeff)
  
  end subroutine loacol_alpha
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute energy ratios
  subroutine compute_energy_ratios(ealpha, tempi1, tempi2, tempe)
    real(wp), intent(in) :: ealpha, tempi1, tempi2, tempe
    
    enrat(1) = ealpha / tempi1
    enrat(2) = ealpha / tempi2
    enrat(3) = ealpha / tempe
  
  end subroutine compute_energy_ratios
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute velocities and velocity ratios
  subroutine compute_velocities(ealpha, tempi1, tempi2, tempe, am1, am2, v0)
    real(wp), intent(in) :: ealpha, tempi1, tempi2, tempe, am1, am2
    real(wp), intent(out) :: v0
    
    real(wp) :: vti1, vti2, vte
    
    ! Alpha particle birth velocity
    v0 = sqrt(2.d0 * ealpha * ev_to_erg / (4.d0 * pmass))
    
    ! Thermal velocities of background species
    vti1 = sqrt(2.d0 * tempi1 * ev_to_erg / (pmass * am1))
    vti2 = sqrt(2.d0 * tempi2 * ev_to_erg / (pmass * am2))
    vte = sqrt(2.d0 * tempe * ev_to_erg / emass)
    
    ! Velocity ratios
    velrat(1) = v0 / vti1
    velrat(2) = v0 / vti2
    velrat(3) = v0 / vte
  
  end subroutine compute_velocities
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute collision frequencies
  subroutine compute_collision_frequencies(am1, am2, Z1, Z2, densi1, densi2, &
                                           tempi1, tempi2, tempe, ealpha, dense, v0)
    real(wp), intent(in) :: am1, am2, Z1, Z2, densi1, densi2
    real(wp), intent(in) :: tempi1, tempi2, tempe, ealpha, dense, v0
    
    real(wp) :: alami1, alami2, alame, frecol_base
    
    ! Coulomb logarithms
    alami1 = compute_coulomb_log_ion(densi1, Z1, tempi1, am1, ealpha)
    alami2 = compute_coulomb_log_ion(densi2, Z2, tempi2, am2, ealpha)
    alame = 24.d0 - log(sqrt(dense) / tempe)
    
    ! Base collision frequency
    frecol_base = 2.d0 * pi * dense * e_charge**4 * 4.d0 / ((4.d0*pmass)**2 * v0**4)
    
    ! Species-specific collision frequencies
    efcolf(1) = frecol_base * Z1**2 * alami1 * densi1 / dense
    efcolf(2) = frecol_base * Z2**2 * alami2 * densi2 / dense
    efcolf(3) = frecol_base * alame
    
    ! Apply velocity ratio scaling
    efcolf = efcolf * velrat
  
  end subroutine compute_collision_frequencies
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute Coulomb logarithm for ions
  function compute_coulomb_log_ion(densi, Z, tempi, am, ealpha) result(alami)
    real(wp), intent(in) :: densi, Z, tempi, am, ealpha
    real(wp) :: alami
    
    real(wp) :: arg
    
    arg = sqrt(densi * Z**2 / tempi) * 2.d0 * Z * (4.d0 + am) / (4.d0*tempi + am*ealpha)
    alami = 23.d0 - log(max(epsilon(1.d0), arg))
  
  end function compute_coulomb_log_ion
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine stost(z,dtauc,iswmode,ierr)
  !
  !  Stochastic collision operator
  !
  real(wp), dimension(5) :: z                     ! Phase space coordinates (modified)
  real(wp) :: dtauc                               ! Normalized time step
  integer :: iswmode                              ! Collision mode switch
  integer :: ierr                                 ! Error code
  
  real(wp) :: p, dpp, dhh, fpeff
  
  ! Get collision coefficients
  p = z(4)
  call coleff(p, dpp, dhh, fpeff)
  
  ierr = 0
  
  ! Apply pitch angle scattering if requested
  if (iswmode == 1 .or. iswmode == 4) then
    call apply_pitch_scattering(z(5), dhh, dtauc, ierr)
    if (iswmode == 4) return
  endif
  
  ! Apply energy scattering and drag if requested
  if (iswmode < 3) then
    call apply_energy_scattering(z(4), dpp, fpeff, dtauc)
  else
    ! Drag only
    z(4) = z(4) + fpeff * dtauc
  endif
  
  ! Handle momentum boundary
  call enforce_momentum_boundary(z(4), ierr)
  
  end subroutine stost
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Apply pitch angle scattering
  subroutine apply_pitch_scattering(alam, dhh, dtauc, ierr)
    real(wp) :: alam
    real(wp) :: dhh, dtauc
    integer :: ierr
    
    real(wp) :: coala, dalam
    real :: ur
    
    coala = 1.d0 - alam**2
    
    ! Check for invalid pitch
    if (coala < 0.d0) then
      ierr = 1
      return
    endif
    
    ! Generate pitch angle change
    call getran(1, ur)
    dalam = sqrt(2.d0 * dhh * coala * dtauc) * dble(ur) - 2.d0 * alam * dhh * dtauc
    
    ! Handle large pitch changes
    if (abs(dalam) > 1.d0) then
      ierr = 2
      call random_number(ur)
      alam = 2.d0 * (dble(ur) - 0.5d0)
    else
      alam = alam + dalam
      ! Reflect from boundaries
      if (alam > 1.d0) then
        ierr = 3
        alam = 2.d0 - alam
      elseif (alam < -1.d0) then
        ierr = 3
        alam = -2.d0 - alam
      endif
    endif
  
  end subroutine apply_pitch_scattering
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Apply energy scattering
  subroutine apply_energy_scattering(p, dpp, fpeff, dtauc)
    real(wp) :: p
    real(wp) :: dpp, fpeff, dtauc
    
    real :: ur
    
    call getran(0, ur)
    p = p + sqrt(abs(2.d0 * dpp * dtauc)) * dble(ur) + fpeff * dtauc
  
  end subroutine apply_energy_scattering
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Enforce momentum boundary conditions
  subroutine enforce_momentum_boundary(p, ierr)
    real(wp) :: p
    integer :: ierr
    
    if (p < pmin) then
      ierr = ierr + 10
      p = pmin + abs(pmin - p)
    endif
  
  end subroutine enforce_momentum_boundary
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine getran(irand,ur)
  !
  !  Produces the random number with zero mean and unit variance
  !
  !  Input parameters: irand - 0 for continuous, 1 for discrete (+1/-1)
  !  Output parameters: ur   - random number
  !
    integer, intent(in) :: irand
    real, intent(out) :: ur
    
    real, parameter :: sqrt_12 = 3.464102  ! sqrt(12)
    
    call random_number(ur)
    
    if (irand == 0) then
      ! Continuous random number with zero mean and unit variance
      ur = sqrt_12 * (ur - 0.5)
    else
      ! Discrete random number (+1 or -1)
      ur = merge(1.0, -1.0, ur > 0.5)
    endif
    
  end subroutine getran
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module collis_alp
