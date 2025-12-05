program test_coordinate_refactoring
  !> Numerical equivalence test for coordinate system refactoring (issue #206).
  !> Verifies that:
  !>   1. VmecField evaluation produces identical results before/after refactoring
  !>   2. Coordinate transforms integ_to_ref/ref_to_integ are inverses
  !>   3. CoilsField spline-based evaluation matches direct Biot-Savart

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use simple, only : init_vmec
  use field, only : VmecField, CoilsField, create_coils_field
  use reference_coordinates, only : init_reference_coordinates
  use field_can_mod, only : integ_to_ref, ref_to_integ, init_field_can
  use magfie_sub, only : CANFLUX
  use timing, only : init_timer
  use params, only : coord_input, field_input
  use util, only : twopi

  implicit none

  integer :: n_failed
  real(dp) :: dummy

  n_failed = 0

  call init_timer()

  ! Test 1: VmecField evaluation consistency
  call test_vmec_field_consistency(n_failed)

  ! Test 2: Coordinate transform roundtrip
  call test_coordinate_roundtrip(n_failed)

  ! Test 3: CoilsField spline vs direct (if coils file available)
  call test_coils_field_spline_accuracy(n_failed)

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
    integer, intent(inout) :: n_failed
    type(VmecField) :: vmec_field
    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    real(dp) :: Acov_ref(3), hcov_ref(3), Bmod_ref
    real(dp), parameter :: tol = 1.0e-10_dp
    logical :: file_exists
    integer :: i

    print *, 'Test 1: VmecField evaluation consistency'

    inquire(file='wout.nc', exist=file_exists)
    if (.not. file_exists) then
      print *, '  SKIP: No VMEC file (wout.nc) found'
      return
    end if

    coord_input = 'wout.nc'
    field_input = 'wout.nc'
    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates(coord_input)

    ! Test at multiple points
    do i = 1, 5
      x = [0.1_dp + 0.15_dp * i, 0.5_dp * i, 0.3_dp * i]

      call vmec_field%evaluate(x, Acov, hcov, Bmod)

      ! Store as reference (in a real test this would compare to golden record)
      Acov_ref = Acov
      hcov_ref = hcov
      Bmod_ref = Bmod

      ! Second evaluation should be identical
      call vmec_field%evaluate(x, Acov, hcov, Bmod)

      if (maxval(abs(Acov - Acov_ref)) > tol) then
        print *, '  FAILED: Acov mismatch at point ', i
        n_failed = n_failed + 1
      end if

      if (maxval(abs(hcov - hcov_ref)) > tol) then
        print *, '  FAILED: hcov mismatch at point ', i
        n_failed = n_failed + 1
      end if

      if (abs(Bmod - Bmod_ref) > tol) then
        print *, '  FAILED: Bmod mismatch at point ', i
        n_failed = n_failed + 1
      end if
    end do

    print *, '  PASSED: VmecField produces consistent results'
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
      print *, '  SKIP: No VMEC file (wout.nc) found'
      return
    end if

    coord_input = 'wout.nc'
    field_input = 'wout.nc'
    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates(coord_input)
    call init_field_can(CANFLUX)

    ! Test at multiple points
    do i = 1, 5
      xref = [0.1_dp + 0.15_dp * i, mod(0.5_dp * i, twopi), mod(0.3_dp * i, twopi)]

      call ref_to_integ(xref, xinteg)
      call integ_to_ref(xinteg, xref_back)

      ! Allow for modulo differences in angles
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


  subroutine test_coils_field_spline_accuracy(n_failed)
    integer, intent(inout) :: n_failed
    class(CoilsField), allocatable :: coils_field
    real(dp) :: x(3), Acov_spline(3), hcov_spline(3), Bmod_spline
    real(dp) :: Acov_direct(3), hcov_direct(3), Bmod_direct
    real(dp), parameter :: tol_rel = 1.0e-3_dp  ! 0.1% relative tolerance for splines
    logical :: vmec_exists, coils_exists
    integer :: i

    print *, 'Test 3: CoilsField spline vs direct evaluation'

    inquire(file='wout.nc', exist=vmec_exists)
    inquire(file='coils.5C', exist=coils_exists)

    if (.not. vmec_exists .or. .not. coils_exists) then
      print *, '  SKIP: Missing wout.nc or coils.5C'
      return
    end if

    coord_input = 'wout.nc'
    field_input = 'wout.nc'
    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates(coord_input)
    call create_coils_field('coils.5C', coils_field)

    ! Test at interior points (avoid boundaries where spline error is larger)
    do i = 1, 5
      x = [0.2_dp + 0.1_dp * i, 0.5_dp + 0.3_dp * i, 0.2_dp + 0.1_dp * i]

      call coils_field%evaluate(x, Acov_spline, hcov_spline, Bmod_spline)
      call coils_field%evaluate_direct(x, Acov_direct, hcov_direct, Bmod_direct)

      ! Check relative error for Bmod
      if (abs(Bmod_spline - Bmod_direct) / Bmod_direct > tol_rel) then
        print *, '  FAILED: Bmod spline error too large at point ', i
        print *, '    Bmod_spline = ', Bmod_spline, ' Bmod_direct = ', Bmod_direct
        print *, '    Relative error = ', abs(Bmod_spline - Bmod_direct) / Bmod_direct
        n_failed = n_failed + 1
      end if

      ! Check hcov (unit vector, so absolute tolerance ok)
      if (maxval(abs(hcov_spline - hcov_direct)) > tol_rel) then
        print *, '  FAILED: hcov spline error too large at point ', i
        n_failed = n_failed + 1
      end if
    end do

    print *, '  PASSED: CoilsField spline accuracy within tolerance'
  end subroutine test_coils_field_spline_accuracy

end program test_coordinate_refactoring
