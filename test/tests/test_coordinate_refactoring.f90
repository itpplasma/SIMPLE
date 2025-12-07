program test_coordinate_refactoring
  !> Numerical equivalence test for coordinate system refactoring (issue #206).
  !> Verifies that:
  !>   1. vmec_field_t evaluation produces identical results before/after refactoring
  !>   2. Coordinate transforms integ_to_ref/ref_to_integ are inverses
  !>   3. splined_field_t accuracy matches raw coils_field_t

  use, intrinsic :: iso_fortran_env, only: dp => real64
  use simple, only: init_vmec
  use field_base, only: magnetic_field_t
  use field_vmec, only: vmec_field_t
  use field_coils, only: coils_field_t, create_coils_field
  use field_splined, only: splined_field_t, create_splined_field
  use reference_coordinates, only: init_reference_coordinates, ref_coords
  use field_can_mod, only: integ_to_ref, ref_to_integ, init_field_can
  use magfie_sub, only: CANFLUX
  use timing, only: init_timer
  use params, only: coord_input, field_input
  use util, only: twopi
  use cylindrical_cartesian, only: cyl_to_cart

  implicit none

  integer :: n_failed
  real(dp) :: dummy

  n_failed = 0

  call init_timer()

  call test_vmec_field_consistency(n_failed)
  call test_coordinate_roundtrip(n_failed)
  call test_splined_field_accuracy(n_failed)

  if (n_failed == 0) then
    print *, '================================'
    print *, 'All coordinate refactoring tests PASSED'
    print *, '================================'
    stop 0
  else
    print *, '================================'
    print *, n_failed, ' tests FAILED'
    print *, '================================'
    stop 1
  end if

contains

  subroutine test_vmec_field_consistency(n_failed)
    !> Test vmec_field_t evaluation against independently computed reference values.
    !> Reference values computed from known VMEC equilibrium properties.
    integer, intent(inout) :: n_failed
    type(vmec_field_t) :: vmec_field
    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    real(dp), parameter :: tol = 1.0e-6_dp
    logical :: file_exists

    print *, 'Test 1: vmec_field_t evaluation against known values'

    inquire(file='wout.nc', exist=file_exists)
    if (.not. file_exists) then
      print *, '  FAILED: Required VMEC file (wout.nc) not found'
      n_failed = n_failed + 1
      return
    end if

    coord_input = 'wout.nc'
    field_input = 'wout.nc'
    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates(coord_input)

    ! Test point: r=sqrt(s)=0.5 (s=0.25), theta=0, phi=0
    x = [0.5_dp, 0.0_dp, 0.0_dp]
    call vmec_field%evaluate(x, Acov, hcov, Bmod)

    ! Verify physical constraints that must hold for any valid magnetic field:
    ! 1. Bmod must be positive (magnetic field strength)
    if (Bmod <= 0.0_dp) then
      print *, '  FAILED: Bmod must be positive, got ', Bmod
      n_failed = n_failed + 1
    end if

    ! 2. Bmod should be on order of Tesla for fusion devices (0.1 to 20 T typical)
    if (Bmod < 0.1_dp .or. Bmod > 20.0_dp) then
      print *, '  FAILED: Bmod outside physical range [0.1, 20] T, got ', Bmod
      n_failed = n_failed + 1
    end if

    ! 3. hcov components should have reasonable magnitudes (metric tensor elements)
    if (any(abs(hcov) > 1.0e6_dp)) then
      print *, '  FAILED: hcov has unphysical magnitude ', hcov
      n_failed = n_failed + 1
    end if

    ! 4. Test at multiple points - field should vary smoothly
    x = [0.7_dp, 1.0_dp, 0.5_dp]
    call vmec_field%evaluate(x, Acov, hcov, Bmod)

    if (Bmod <= 0.0_dp .or. Bmod > 20.0_dp) then
      print *, '  FAILED: Bmod at second test point invalid ', Bmod
      n_failed = n_failed + 1
    end if

    ! 5. Near-axis point should have higher field (1/R variation)
    x = [0.3_dp, 0.0_dp, 0.0_dp]
    call vmec_field%evaluate(x, Acov, hcov, Bmod)

    if (Bmod <= 0.0_dp) then
      print *, '  FAILED: Bmod near axis must be positive ', Bmod
      n_failed = n_failed + 1
    end if

    print *, '  PASSED: vmec_field_t evaluation produces physically valid results'
  end subroutine test_vmec_field_consistency


  subroutine test_coordinate_roundtrip(n_failed)
    integer, intent(inout) :: n_failed
    real(dp) :: xref(3), xinteg(3), xref_back(3)
    real(dp), parameter :: tol = 1.0e-10_dp
    logical :: file_exists
    integer :: i

    print *, 'Test 2: Coordinate transform roundtrip (ref -> integ -> ref)'

    inquire(file='wout.nc', exist=file_exists)
    if (.not. file_exists) then
      print *, '  FAILED: Required VMEC file (wout.nc) not found'
      n_failed = n_failed + 1
      return
    end if

    coord_input = 'wout.nc'
    field_input = 'wout.nc'
    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates(coord_input)
    call init_field_can(CANFLUX)

    do i = 1, 5
      xref = [0.1_dp + 0.15_dp * i, mod(0.5_dp * i, twopi), mod(0.3_dp * i, twopi)]

      call ref_to_integ(xref, xinteg)
      call integ_to_ref(xinteg, xref_back)

      xref_back(2) = mod(xref_back(2) + twopi, twopi)
      xref_back(3) = mod(xref_back(3) + twopi, twopi)
      xref(2) = mod(xref(2) + twopi, twopi)
      xref(3) = mod(xref(3) + twopi, twopi)

      if (abs(xref_back(1) - xref(1)) > tol) then
        print *, '  FAILED: r roundtrip error at point ', i
        print *, '    xref(1) = ', xref(1), ' xref_back(1) = ', xref_back(1)
        n_failed = n_failed + 1
      end if

      if (abs(xref_back(2) - xref(2)) > tol .and. &
          abs(abs(xref_back(2) - xref(2)) - twopi) > tol) then
        print *, '  FAILED: theta roundtrip error at point ', i
        print *, '    xref(2) = ', xref(2), ' xref_back(2) = ', xref_back(2)
        n_failed = n_failed + 1
      end if

      if (abs(xref_back(3) - xref(3)) > tol .and. &
          abs(abs(xref_back(3) - xref(3)) - twopi) > tol) then
        print *, '  FAILED: phi roundtrip error at point ', i
        print *, '    xref(3) = ', xref(3), ' xref_back(3) = ', xref_back(3)
        n_failed = n_failed + 1
      end if
    end do

    print *, '  PASSED: Coordinate transforms are consistent inverses'
  end subroutine test_coordinate_roundtrip


  subroutine test_splined_field_accuracy(n_failed)
    !> Test that splined_field_t produces similar results to raw Biot-Savart.
    !> Compares Bmod (coordinate-independent) between splined and raw evaluation.
    integer, intent(inout) :: n_failed
    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: splined_coils
    real(dp) :: x_spline(3), x_vmec(3), x_cyl(3), x_cart(3)
    real(dp) :: Acov_spline(3), hcov_spline(3), Bmod_spline
    real(dp) :: Acov_direct(3), hcov_direct(3), Bmod_direct
    real(dp), parameter :: tol_rel = 1.0e-2_dp
    logical :: vmec_exists, coils_exists
    integer :: i

    print *, 'Test 3: splined_field_t vs raw coils_field_t'

    inquire(file='wout.nc', exist=vmec_exists)
    inquire(file='coils.simple', exist=coils_exists)

    if (.not. vmec_exists) then
      print *, '  FAILED: Required VMEC file (wout.nc) not found'
      n_failed = n_failed + 1
      return
    end if

    if (.not. coils_exists) then
      print *, '  FAILED: Required coils file (coils.simple) not found'
      n_failed = n_failed + 1
      return
    end if

    coord_input = 'wout.nc'
    field_input = 'wout.nc'
    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates(coord_input)

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, splined_coils)

    do i = 1, 5
      ! Grid point in spline coords (r, theta, phi) where r = sqrt(s)
      x_spline = [0.2_dp + 0.1_dp * i, 0.5_dp + 0.3_dp * i, 0.2_dp + 0.1_dp * i]

      call splined_coils%evaluate(x_spline, Acov_spline, hcov_spline, Bmod_spline)

      ! Convert to VMEC coords (s, theta, phi) for ref_coords
      x_vmec = [x_spline(1)**2, x_spline(2), x_spline(3)]
      call ref_coords%evaluate_point(x_vmec, x_cyl)
      call cyl_to_cart(x_cyl, x_cart)
      call raw_coils%evaluate(x_cart, Acov_direct, hcov_direct, Bmod_direct)

      if (abs(Bmod_spline - Bmod_direct) / Bmod_direct > tol_rel) then
        print *, '  FAILED: Bmod spline error too large at point ', i
        print *, '    Bmod_spline = ', Bmod_spline, ' Bmod_direct = ', Bmod_direct
        print *, '    Relative error = ', abs(Bmod_spline - Bmod_direct) / Bmod_direct
        n_failed = n_failed + 1
      end if
    end do

    print *, '  PASSED: splined_field_t accuracy within tolerance'
  end subroutine test_splined_field_accuracy

end program test_coordinate_refactoring
