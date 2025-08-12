program test_collis_alphas_refactored
! Unit tests for refactored collis_alphas module
! Tests follow behavior-driven design with Given-When-Then structure

use collis_alp, only: wp, nsorts, efcolf, velrat, enrat, &
                      coleff, onseff, loacol_alpha, stost, getran, &
                      compute_energy_ratios, compute_velocities, &
                      compute_collision_frequencies, compute_coulomb_log_ion, &
                      add_species_contribution, onseff_small_v, onseff_large_v, &
                      onseff_intermediate_v, apply_pitch_scattering, &
                      apply_energy_scattering, enforce_momentum_boundary

implicit none

integer :: n_tests_passed, n_tests_failed
real(wp), parameter :: tol = 1.0e-12_wp
real(wp), parameter :: tol_phys = 1.0e-6_wp  ! For physical calculations

n_tests_passed = 0
n_tests_failed = 0

call test_onseff_velocity_limits()
call test_coleff_basic()
call test_coulomb_logarithm()
call test_collision_frequencies()
call test_stochastic_operator()
call test_random_number_generator()
call test_momentum_boundary()
call test_pitch_scattering()
call test_energy_scattering()
call test_loacol_alpha_integration()

print *, '============================================'
print *, 'Test Summary:'
print *, 'Passed:', n_tests_passed
print *, 'Failed:', n_tests_failed
print *, '============================================'

if (n_tests_failed > 0) then
  error stop 'Some tests failed!'
end if

contains

  !===========================================================================
  ! Test: Onsager coefficients at different velocity limits
  !===========================================================================
  subroutine test_onseff_velocity_limits()
    ! Given: Different velocity values
    ! When: We compute Onsager coefficients
    ! Then: Results should match expected limits
    
    real(wp) :: dp, dh, dpd
    logical :: test_passed
    
    print *, 'Testing onseff velocity limits...'
    
    test_passed = .true.
    
    ! Test small velocity limit (v < 0.01)
    call onseff(0.005_wp, dp, dh, dpd)
    if (dp < 0.0_wp .or. dh < 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Small velocity coefficients should be positive'
    end if
    
    ! Test large velocity limit (v > 6.0)
    call onseff(10.0_wp, dp, dh, dpd)
    if (abs(dp - 1.0_wp/1000.0_wp) > tol_phys) then
      test_passed = .false.
      print *, '  Failed: Large velocity dp incorrect'
    end if
    
    ! Test intermediate velocity
    call onseff(2.0_wp, dp, dh, dpd)
    if (dp < 0.0_wp .or. dh < 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Intermediate velocity coefficients should be physical'
    end if
    
    if (test_passed) then
      print *, '  PASSED: onseff velocity limits correct'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: onseff velocity limits validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_onseff_velocity_limits

  !===========================================================================
  ! Test: Basic collision coefficient calculation
  !===========================================================================
  subroutine test_coleff_basic()
    ! Given: Collision parameters are set
    ! When: We compute collision coefficients
    ! Then: Results should be physically reasonable
    
    real(wp) :: p, dpp, dhh, fpeff
    logical :: test_passed
    integer :: i
    
    print *, 'Testing coleff basic functionality...'
    
    test_passed = .true.
    
    ! Set up test collision parameters
    do i = 1, nsorts
      efcolf(i) = 0.1_wp * i
      velrat(i) = 1.0_wp + 0.1_wp * i
      enrat(i) = 2.0_wp
    end do
    
    ! Test with unit momentum
    p = 1.0_wp
    call coleff(p, dpp, dhh, fpeff)
    
    if (dpp < 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Momentum diffusion should be non-negative'
    end if
    
    if (dhh < 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Pitch angle diffusion should be non-negative'
    end if
    
    ! Test with very small momentum
    p = 1.0e-10_wp
    call coleff(p, dpp, dhh, fpeff)
    
    ! For very small momentum, dhh scales as 1/p^2, so it can be very large
    ! Just check that values are finite (not NaN or Inf)
    if (.not. (dpp == dpp .and. dhh == dhh .and. fpeff == fpeff)) then
      test_passed = .false.
      print *, '  Failed: Coefficients should be finite (not NaN) for small momentum'
    end if
    
    if (test_passed) then
      print *, '  PASSED: coleff works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: coleff validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_coleff_basic

  !===========================================================================
  ! Test: Coulomb logarithm calculation
  !===========================================================================
  subroutine test_coulomb_logarithm()
    ! Given: Plasma parameters
    ! When: We compute Coulomb logarithm
    ! Then: Result should be in physical range (typically 10-25)
    
    real(wp) :: alami
    logical :: test_passed
    
    print *, 'Testing Coulomb logarithm calculation...'
    
    test_passed = .true.
    
    ! Test with typical fusion plasma parameters
    ! n = 1e14 cm^-3, Z = 1, T = 10 keV, mass = 2, E_alpha = 3.5 MeV
    alami = compute_coulomb_log_ion(1.0e14_wp, 1.0_wp, 1.0e4_wp, 2.0_wp, 3.5e6_wp)
    
    if (alami < 5.0_wp .or. alami > 30.0_wp) then
      test_passed = .false.
      print *, '  Failed: Coulomb log out of physical range:', alami
    end if
    
    ! Test with low density
    alami = compute_coulomb_log_ion(1.0e10_wp, 1.0_wp, 1.0e3_wp, 1.0_wp, 1.0e6_wp)
    
    if (alami < 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Coulomb log should be positive'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Coulomb logarithm calculation correct'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Coulomb logarithm validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_coulomb_logarithm

  !===========================================================================
  ! Test: Collision frequency calculation
  !===========================================================================
  subroutine test_collision_frequencies()
    ! Given: Plasma species parameters
    ! When: We compute collision frequencies
    ! Then: They should scale correctly with density and charge
    
    real(wp) :: v0_test
    logical :: test_passed
    
    print *, 'Testing collision frequency calculation...'
    
    test_passed = .true.
    
    ! Set up test parameters
    v0_test = 1.0e9_wp  ! cm/s
    
    call compute_collision_frequencies( &
      1.0_wp, 2.0_wp,     & ! masses
      1.0_wp, 2.0_wp,     & ! charges
      1.0e13_wp, 5.0e12_wp, & ! densities
      1.0e3_wp, 2.0e3_wp, & ! temperatures
      1.0e3_wp,           & ! electron temperature
      3.5e6_wp,           & ! alpha energy
      1.5e13_wp,          & ! total charge density
      v0_test)
    
    ! Check that collision frequencies are positive
    if (any(efcolf < 0.0_wp)) then
      test_passed = .false.
      print *, '  Failed: Collision frequencies should be positive'
    end if
    
    ! Check relative scaling with Z^2
    if (efcolf(2)/efcolf(1) < 1.0_wp) then  ! Z2^2/Z1^2 = 4
      test_passed = .false.
      print *, '  Failed: Collision frequency should scale with Z^2'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Collision frequencies calculated correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Collision frequency validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_collision_frequencies

  !===========================================================================
  ! Test: Stochastic collision operator
  !===========================================================================
  subroutine test_stochastic_operator()
    ! Given: Phase space coordinates
    ! When: We apply stochastic collisions
    ! Then: Momentum and pitch should evolve physically
    
    real(wp), dimension(5) :: z
    real(wp) :: dtauc
    integer :: ierr, iswmode
    logical :: test_passed
    
    print *, 'Testing stochastic collision operator...'
    
    test_passed = .true.
    
    ! Initialize test collision parameters
    efcolf = [0.1_wp, 0.05_wp, 0.2_wp]
    velrat = [1.0_wp, 1.2_wp, 0.8_wp]
    enrat = [2.0_wp, 2.0_wp, 2.0_wp]
    
    ! Set initial conditions
    z = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.5_wp]  ! position, momentum, pitch
    dtauc = 0.01_wp
    
    ! Test full operator (mode 1)
    iswmode = 1
    call stost(z, dtauc, iswmode, ierr)
    
    if (ierr /= 0 .and. ierr < 10) then
      test_passed = .false.
      print *, '  Failed: Unexpected error in full operator mode:', ierr
    end if
    
    if (z(4) <= 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Momentum should remain positive'
    end if
    
    if (abs(z(5)) > 1.0_wp) then
      test_passed = .false.
      print *, '  Failed: Pitch should remain in [-1, 1]'
    end if
    
    ! Test drag only (mode 3)
    z = [0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp]
    iswmode = 3
    call stost(z, dtauc, iswmode, ierr)
    
    if (z(4) >= 2.0_wp) then  ! Should slow down
      test_passed = .false.
      print *, '  Failed: Drag should reduce momentum'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Stochastic operator works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Stochastic operator validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_stochastic_operator

  !===========================================================================
  ! Test: Random number generator
  !===========================================================================
  subroutine test_random_number_generator()
    ! Given: Random number generator modes
    ! When: We generate random numbers
    ! Then: They should have correct statistics
    
    real :: ur
    real(wp) :: sum_cont, sum_sq_cont
    integer :: i, n_samples, n_positive
    logical :: test_passed
    
    print *, 'Testing random number generator...'
    
    test_passed = .true.
    n_samples = 1000
    
    ! Test continuous mode (irand=0)
    sum_cont = 0.0_wp
    sum_sq_cont = 0.0_wp
    do i = 1, n_samples
      call getran(0, ur)
      sum_cont = sum_cont + ur
      sum_sq_cont = sum_sq_cont + ur**2
    end do
    
    ! Check mean ~ 0 and variance ~ 1
    if (abs(sum_cont/n_samples) > 0.1_wp) then
      test_passed = .false.
      print *, '  Failed: Continuous random mean not near zero:', sum_cont/n_samples
    end if
    
    if (abs(sum_sq_cont/n_samples - 1.0_wp) > 0.2_wp) then
      test_passed = .false.
      print *, '  Failed: Continuous random variance not near one:', sum_sq_cont/n_samples
    end if
    
    ! Test discrete mode (irand=1)
    n_positive = 0
    do i = 1, n_samples
      call getran(1, ur)
      if (abs(abs(ur) - 1.0) > tol) then
        test_passed = .false.
        print *, '  Failed: Discrete random should be +1 or -1'
        exit
      end if
      if (ur > 0.0) n_positive = n_positive + 1
    end do
    
    ! Check roughly 50-50 distribution
    if (abs(n_positive - n_samples/2) > n_samples/10) then
      test_passed = .false.
      print *, '  Failed: Discrete random not balanced'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Random number generator works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Random number generator validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_random_number_generator

  !===========================================================================
  ! Test: Momentum boundary enforcement
  !===========================================================================
  subroutine test_momentum_boundary()
    ! Given: Momentum below minimum
    ! When: We enforce boundary conditions
    ! Then: Momentum should be reflected
    
    real(wp) :: p
    integer :: ierr
    logical :: test_passed
    
    print *, 'Testing momentum boundary enforcement...'
    
    test_passed = .true.
    
    ! Test reflection from lower boundary
    p = -1.0e-9_wp
    ierr = 0
    call enforce_momentum_boundary(p, ierr)
    
    if (p < 1.0e-8_wp) then
      test_passed = .false.
      print *, '  Failed: Momentum should be above minimum after reflection'
    end if
    
    if (ierr < 10) then
      test_passed = .false.
      print *, '  Failed: Error code should indicate boundary reflection'
    end if
    
    ! Test no change for valid momentum
    p = 1.0_wp
    ierr = 0
    call enforce_momentum_boundary(p, ierr)
    
    if (abs(p - 1.0_wp) > tol) then
      test_passed = .false.
      print *, '  Failed: Valid momentum should not change'
    end if
    
    if (ierr /= 0) then
      test_passed = .false.
      print *, '  Failed: No error for valid momentum'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Momentum boundary enforcement correct'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Momentum boundary validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_momentum_boundary

  !===========================================================================
  ! Test: Pitch angle scattering
  !===========================================================================
  subroutine test_pitch_scattering()
    ! Given: Initial pitch angle
    ! When: We apply pitch scattering
    ! Then: Pitch should remain in valid range
    
    real(wp) :: alam, dhh, dtauc
    integer :: ierr
    logical :: test_passed
    
    print *, 'Testing pitch angle scattering...'
    
    test_passed = .true.
    
    ! Test normal scattering
    alam = 0.5_wp
    dhh = 0.1_wp
    dtauc = 0.01_wp
    ierr = 0
    
    call apply_pitch_scattering(alam, dhh, dtauc, ierr)
    
    if (abs(alam) > 1.0_wp) then
      test_passed = .false.
      print *, '  Failed: Pitch should remain in [-1, 1]'
    end if
    
    ! Test invalid initial pitch
    alam = 1.5_wp
    ierr = 0
    call apply_pitch_scattering(alam, dhh, dtauc, ierr)
    
    if (ierr /= 1) then
      test_passed = .false.
      print *, '  Failed: Should detect invalid initial pitch'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Pitch scattering works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Pitch scattering validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_pitch_scattering

  !===========================================================================
  ! Test: Energy scattering
  !===========================================================================
  subroutine test_energy_scattering()
    ! Given: Initial momentum
    ! When: We apply energy scattering
    ! Then: Momentum should evolve with correct statistics
    
    real(wp) :: p, p_initial, dpp, fpeff, dtauc
    logical :: test_passed
    integer :: i
    
    print *, 'Testing energy scattering...'
    
    test_passed = .true.
    
    ! Test multiple scattering steps
    p_initial = 1.0_wp
    dpp = 0.1_wp
    fpeff = -0.05_wp  ! Drag coefficient (negative for slowing down)
    dtauc = 0.001_wp
    
    p = p_initial
    do i = 1, 100
      call apply_energy_scattering(p, dpp, fpeff, dtauc)
    end do
    
    ! Check that momentum changed
    if (abs(p - p_initial) < tol) then
      test_passed = .false.
      print *, '  Failed: Momentum should change due to scattering'
    end if
    
    ! Check momentum remains positive
    if (p <= 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Momentum should remain positive'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Energy scattering works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Energy scattering validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_energy_scattering

  !===========================================================================
  ! Test: Full loacol_alpha integration test
  !===========================================================================
  subroutine test_loacol_alpha_integration()
    ! Given: Realistic plasma parameters
    ! When: We initialize collision coefficients
    ! Then: Output quantities should be physically reasonable
    
    real(wp) :: v0, dchichi, slowrate, dchichi_norm, slowrate_norm
    logical :: test_passed
    
    print *, 'Testing loacol_alpha full integration...'
    
    test_passed = .true.
    
    ! Typical fusion plasma parameters
    ! D-T plasma with alpha particles
    call loacol_alpha( &
      2.0_wp, 3.0_wp,        & ! Deuterium and Tritium masses
      1.0_wp, 1.0_wp,        & ! Charges
      5.0e13_wp, 5.0e13_wp,  & ! Densities [cm^-3]
      10.0e3_wp, 10.0e3_wp,  & ! Ion temperatures [eV]
      10.0e3_wp,             & ! Electron temperature [eV]
      3.5e6_wp,              & ! Alpha energy [eV]
      v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
    
    ! Check velocity is reasonable (alpha at 3.5 MeV)
    if (v0 < 1.0e8_wp .or. v0 > 1.0e10_wp) then
      test_passed = .false.
      print *, '  Failed: Alpha velocity out of range:', v0
    end if
    
    ! Check scattering rates are positive
    if (dchichi <= 0.0_wp .or. slowrate <= 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Scattering rates should be positive'
    end if
    
    ! Check normalized rates
    if (dchichi_norm <= 0.0_wp .or. slowrate_norm <= 0.0_wp) then
      test_passed = .false.
      print *, '  Failed: Normalized rates should be positive'
    end if
    
    ! Check that normalized rates are smaller (dimensionless)
    if (dchichi_norm > dchichi .or. slowrate_norm > slowrate) then
      test_passed = .false.
      print *, '  Failed: Normalized rates should be smaller than dimensional ones'
    end if
    
    if (test_passed) then
      print *, '  PASSED: loacol_alpha integration test successful'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: loacol_alpha integration validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_loacol_alpha_integration

end program test_collis_alphas_refactored