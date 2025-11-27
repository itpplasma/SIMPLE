program test_sympl
  use, intrinsic :: iso_fortran_env, only : dp => real64
  use orbit_symplectic, only : SymplecticIntegrator, orbit_sympl_init, &
    orbit_timestep_sympl, EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT, GAUSS4
  use field_can_mod, only : FieldCan, FieldCan_init, eval_field => evaluate, &
    field_can_from_name
  use timing, only : get_wtime

  implicit none

  integer, parameter :: steps_per_bounce = 8
  integer, parameter :: nbounce = 30
  integer, parameter :: ntau = 1
  real(dp), parameter :: taub = 7800.0_dp

  type(FieldCan) :: f
  type(SymplecticIntegrator) :: integ
  real(dp) :: z0(4)
  real(dp) :: dt, h0, hmin, hmax, rel_drift, rel_tol, denom
  real(dp) :: starttime, endtime
  integer :: nt, k, ierr, unit, i
  integer, parameter :: nmode = 4
  integer, dimension(nmode) :: mode_ids
  character(len=12), dimension(nmode) :: labels
  character(len=40) :: outfile
  real(dp), allocatable :: orbit(:, :)

  dt = taub/real(steps_per_bounce, dp)
  nt = steps_per_bounce*nbounce
  mode_ids = (/EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT, GAUSS4/)
  labels = (/'euler1     ', 'euler2     ', 'midpoint   ', 'gauss4     '/)

  allocate(orbit(5, nt))

  do i = 1, nmode
    call field_can_from_name('test')
    call FieldCan_init(f, 1.0e-5_dp, 1.0_dp, 0.0_dp)

    z0 = (/0.1_dp, 1.5_dp, 0.0_dp, 0.0_dp/)
    call eval_field(f, z0(1), z0(2), z0(3), 0)
    z0(4) = f%vpar*f%hph + f%Aph/f%ro0

    rel_tol = select_tol(mode_ids(i))
    call orbit_sympl_init(integ, f, z0, dt, ntau, 1.0e-12_dp, mode_ids(i))

    orbit(:, 1) = (/z0, f%H/)
    h0 = f%H
    hmin = h0
    hmax = h0

    starttime = get_wtime()
    do k = 2, nt
      ierr = 0
      call orbit_timestep_sympl(integ, f, ierr)
      if (ierr /= 0) then
        write(*, '(a, i0, a)') 'orbit_timestep_sympl error ', ierr, &
          ' in '//trim(labels(i))
        error stop 1
      end if
      orbit(1:4, k) = integ%z
      orbit(5, k) = f%H
      hmax = max(hmax, f%H)
      hmin = min(hmin, f%H)
    end do
    endtime = get_wtime()

    denom = abs(h0)
    if (denom < 1.0e-14_dp) denom = 1.0_dp
    rel_drift = (hmax - hmin)/denom
    if (rel_drift > rel_tol) then
      write(*, '(a, a, a, 1pE10.3, a, 1pE10.3)') 'Hamiltonian spread for ', &
        trim(labels(i)), ' exceeds tolerance: ', rel_drift, ' > ', rel_tol
      error stop 1
    end if

    write(*, '(a, a, a, f9.5, a, 1pE10.3)') 'tokamak ', trim(labels(i)), ': ', &
      endtime - starttime, ' s, max dH/H=', rel_drift

    write(outfile, '(a, a)') 'tokamak_testfield_', trim(labels(i))//'.dat'
    open(newunit=unit, file=trim(outfile), action='write', status='replace', &
      recl=4096)
    do k = 1, nt
      write(unit, *) orbit(:, k)
    end do
    close(unit)
  end do

  deallocate(orbit)

contains

  function select_tol(mode) result(tol)
    integer, intent(in) :: mode
    real(dp) :: tol

    select case (mode)
    case (EXPL_IMPL_EULER)
      tol = 1.24e-1_dp
    case (IMPL_EXPL_EULER)
      tol = 1.20e-1_dp
    case (MIDPOINT)
      tol = 5.3e-3_dp
    case (GAUSS4)
      tol = 4.5e-5_dp
    case default
      tol = 1.0e-5_dp
    end select
  end function select_tol
end program test_sympl
