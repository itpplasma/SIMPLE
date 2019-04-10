module neo_orb
  use omp_lib
  use common, only: pi, c, e_charge, e_mass, p_mass, ev
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp

  use parmot_mod, only : rmu, ro0
  use velo_mod,   only : isw_field_type
  use orbit_symplectic, only : orbit_sympl_init, orbit_timestep_sympl
  use diag_mod, only : icounter

  implicit none

  double precision :: dtau, dtaumax, v0
  double precision, dimension(5) :: z
  integer          :: n_e, n_d

  integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet
  double precision :: relerr

  logical :: firstrun = .True.

contains

  subroutine init_field(ans_s, ans_tp, amultharm, aintegmode)
    ! initialize field geometry
    ! character*32, intent(in) :: vmec_file
    integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
    integer             :: ierr
    integer             :: L1i
    double precision    :: RT0, R0i, cbfi, bz0i, bf0

    netcdffile = 'wout.nc'  ! TODO: don't hardcode this
    ns_s = ans_s
    ns_tp = ans_tp
    multharm = amultharm
    integmode = aintegmode

    call spline_vmec_data ! initialize splines for VMEC field
    call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius
    print *, 'R0 = ', RT0, ' cm'
    isw_field_type = 1 ! evaluate fields in VMEC coords (0 = CAN, 1 = VMEC)
    if (integmode>=0) then
      call get_canonical_coordinates ! pre-compute transformation to canonical coords
      isw_field_type = 0 ! evaluate fields in canonical coords (0 = CAN, 1 = VMEC)
    end if

    ! initialize position and do first check if z is inside vacuum chamber
    z = 0.0d0
    call chamb_can(z(1:2), z(3), ierr)
    ! TODO: error handling
  end subroutine init_field

  subroutine init_params(Z_charge, m_mass, E_kin, adtau, adtaumax, arelerr)
    ! Initializes normalization for velocity and Larmor radius based on kinetic energy
    ! of plasma particles (= temperature for thermal particles).

    integer, intent(in) :: Z_charge, m_mass
    real(8), intent(in) :: E_kin, adtau, adtaumax
    real(8), intent(in) :: arelerr

    n_e = Z_charge
    n_d = m_mass
    relerr = arelerr

    ! Neglect relativistic effects by large inverse relativistic temperature
    rmu=1d8

    ! Reference velocity and normalized Larmor radius
    v0 = sqrt(2.d0*E_kin*ev/(n_d*p_mass))
    ro0 = v0*n_d*p_mass*c/(n_e*e_charge)

    dtau = adtau ! timestep where to get results
    dtaumax = adtaumax ! maximum timestep for adaptive integration

end subroutine init_params

  subroutine timestep(s, th, ph, lam, ierr)
    real(8), intent(inout) :: s, th, ph, lam
    integer, intent(out) :: ierr

    z(1) = s
    z(2) = th
    z(3) = ph
    z(4) = 1d0
    z(5) = lam

    if (integmode <= 0) then
      call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
    else
      if (firstrun) then
        call orbit_sympl_init(z, dtau, dtaumax, relerr, integmode)
        firstrun = .False.
      endif
      call orbit_timestep_sympl(z, ierr)
    endif

    s = z(1)
    th = z(2)
    ph = z(3)
    lam = z(5)
  end subroutine timestep

end module neo_orb
