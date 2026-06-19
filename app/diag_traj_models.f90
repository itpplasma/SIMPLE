program diag_traj_models
  ! Dump one trapped (banana) orbit traced by three orbit models for the
  ! cross-validation figures: GC (symplectic guiding center), CPP (flux-canonical
  ! Pauli particle), and the gyro-resolved full orbit (Boris). Two field cases:
  !
  !   analytic : circular-tokamak "test" canonical chart for GC/CPP
  !              (field_can_test, R0=1, a=0.5, B0=1, iota0=1) and the matching
  !              analytic cylindrical tokamak provider for the full orbit
  !              (orbit_full_tokamak, same equilibrium). No field file.
  !   vmec     : BOOZER chart of the QA wout for GC/CPP. The full-orbit VMEC
  !              curvilinear provider depends on libneo metric/Christoffel
  !              (PR #322) and is not wired here, so the vmec case dumps GC and
  !              CPP only.
  !
  ! Each model writes <outdir>/<case>_<model>.dat with columns
  !   t   R   Z   Hrel   pphi   mu
  ! where Hrel = H/H0 - 1, in R-Z [cm] for the figures (compare_gc_pauli.py).
  !
  ! Usage: ./diag_traj_models.x <analytic|vmec> <outdir> [nstep]
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: twopi, c, p_mass, e_charge
  use field_can_mod, only: field_can_t, field_can_from_name, evaluate, get_val, &
      integ_to_ref
  use orbit_symplectic, only: orbit_timestep_sympl, orbit_sympl_init
  use orbit_symplectic_base, only: symplectic_integrator_t, GAUSS1
  use orbit_cpp, only: orbit_timestep_cpp, orbit_cpp_init, cpp_stages_from_mode
  use orbit_full, only: FullOrbitState, init_full_orbit_state, &
      timestep_full_orbit, compute_energy, ORBIT_BORIS, COORD_CYL, FO_OK
  use orbit_full_tokamak, only: tokamak_provider_t
  use simple, only: tracer_t, init_sympl, init_params
  use simple_main, only: init_field
  use simple_coordinates, only: get_transform, transform_i
  use params, only: isw_field_type, field_input, coord_input, integmode
  use new_vmec_stuff_mod, only: rmajor
  use magfie_sub, only: BOOZER

  implicit none

  character(64)  :: case_arg, nstep_arg
  character(256) :: outdir
  integer :: nargs, nstep
  type(tracer_t) :: norb

  nargs = command_argument_count()
  if (nargs < 2) then
    print *, 'Usage: ./diag_traj_models.x <analytic|vmec> <outdir> [nstep]'
    error stop 1
  end if
  call get_command_argument(1, case_arg)
  call get_command_argument(2, outdir)
  nstep = 600
  if (nargs >= 3) then
    call get_command_argument(3, nstep_arg); read(nstep_arg, *) nstep
  end if

  select case (trim(case_arg))
  case ('analytic')
    call run_analytic(trim(outdir), nstep)
  case ('vmec')
    call run_vmec(norb, trim(outdir), nstep)
  case default
    print *, 'unknown case: ', trim(case_arg)
    error stop 1
  end select

contains

  ! Initialize GC on the analytic test chart with explicit Larmor scale ro0.
  subroutine init_test_chart(si, f, z0, ro0_in, dt, ntau, rtol, mode)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: z0(5), ro0_in, dt, rtol
    integer, intent(in) :: ntau, mode
    real(dp) :: z(4)

    call evaluate(f, z0(1), z0(2), z0(3), 0)
    f%mu  = 0.5_dp*z0(4)**2*(1.0_dp - z0(5)**2)/f%Bmod*2.0_dp
    f%ro0 = ro0_in
    f%vpar = z0(4)*z0(5)*sqrt(2.0_dp)
    z(1:3) = z0(1:3)
    z(4)   = f%vpar*f%hph + f%Aph/f%ro0
    call orbit_sympl_init(si, f, z, dt, ntau, rtol, mode)
  end subroutine init_test_chart

  subroutine run_analytic(outdir, nstep)
    character(*), intent(in) :: outdir
    integer, intent(in) :: nstep
    real(dp), parameter :: R0 = 1.0_dp, a = 0.5_dp, B0 = 1.0_dp, iota0 = 1.0_dp
    real(dp) :: z0(5), dt, rtol, ro0
    integer :: s

    call field_can_from_name('test')
    if (.not. associated(evaluate)) error stop 'evaluate not associated (test)'

    ! Trapped banana IC on the analytic chart: outer-mid-radius surface
    ! (r=0.45 of a=0.5), small pitch, deeply trapped so the banana is stable
    ! over the whole trace and wide enough to be visible.
    z0   = [0.45_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.20_dp]
    dt   = 2.0_dp
    rtol = 1.0e-13_dp
    ro0  = 1.0e-3_dp
    s = cpp_stages_from_mode(GAUSS1)

    call dump_gc_analytic(outdir, z0, ro0, dt, rtol, GAUSS1, nstep, R0)
    call dump_cpp_analytic(outdir, z0, ro0, dt, rtol, GAUSS1, s, nstep, R0)
    call dump_fullorbit_analytic(outdir, z0, ro0, dt, nstep, R0, a, B0, iota0)
  end subroutine run_analytic

  subroutine dump_gc_analytic(outdir, z0, ro0, dt, rtol, mode, nstep, R0)
    character(*), intent(in) :: outdir
    real(dp), intent(in) :: z0(5), ro0, dt, rtol, R0
    integer, intent(in) :: mode, nstep
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    integer :: it, ierr, u
    real(dp) :: H0, t

    call init_test_chart(si, f, z0, ro0, dt, 1, rtol, mode)
    call get_val(f, si%z(4)); H0 = f%H
    open(newunit=u, file=outdir//'/analytic_gc.dat', status='replace')
    write(u, '(A)') '# t  R  Z  Hrel  pphi  mu'
    call write_can_row(u, 0.0_dp, si%z, f%H, H0, si%z(4), f%mu, R0)
    do it = 1, nstep
      call orbit_timestep_sympl(si, f, ierr)
      if (ierr /= 0) exit
      call get_val(f, si%z(4))
      t = it*dt
      call write_can_row(u, t, si%z, f%H, H0, si%z(4), f%mu, R0)
    end do
    close(u)
    print '(A,I0,A)', 'analytic GC: ', it-1, ' steps -> analytic_gc.dat'
  end subroutine dump_gc_analytic

  subroutine dump_cpp_analytic(outdir, z0, ro0, dt, rtol, mode, s, nstep, R0)
    character(*), intent(in) :: outdir
    real(dp), intent(in) :: z0(5), ro0, dt, rtol, R0
    integer, intent(in) :: mode, s, nstep
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    integer :: it, ierr, u
    real(dp) :: H0, t

    call init_test_chart(si, f, z0, ro0, dt, 1, rtol, mode)
    call orbit_cpp_init(si, f)
    call get_val(f, si%z(4)); H0 = f%H
    open(newunit=u, file=outdir//'/analytic_cpp.dat', status='replace')
    write(u, '(A)') '# t  R  Z  Hrel  pphi  mu'
    call write_can_row(u, 0.0_dp, si%z, f%H, H0, si%z(4), f%mu, R0)
    do it = 1, nstep
      call orbit_timestep_cpp(si, f, s, ierr)
      if (ierr /= 0) exit
      call get_val(f, si%z(4))
      t = it*dt
      call write_can_row(u, t, si%z, f%H, H0, si%z(4), f%mu, R0)
    end do
    close(u)
    print '(A,I0,A)', 'analytic CPP: ', it-1, ' steps -> analytic_cpp.dat'
  end subroutine dump_cpp_analytic

  ! Full orbit (Boris) on the matching analytic cylindrical tokamak. The IC is
  ! placed at the banana's outer-midplane turning point with vpar from the same
  ! pitch lambda, so the gyro-resolved orbit gyrates around the GC banana.
  subroutine dump_fullorbit_analytic(outdir, z0, ro0, dt_gc, nstep, R0, a, B0, iota0)
    character(*), intent(in) :: outdir
    real(dp), intent(in) :: z0(5), ro0, dt_gc, R0, a, B0, iota0
    integer, intent(in) :: nstep
    type(tokamak_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, vtot, lambda, r_min, R, Z
    real(dp) :: Bvec(3), Bmod, hcov(3), vpar, vperp
    real(dp) :: u0(3), w0(3), b_hat(3), e1(3), e2(3)
    real(dp) :: Omega, period, dt, t, H0
    integer :: it, ierr, u, nsub, isub, ierrf, nout

    prov%R0 = R0; prov%a = a; prov%B0 = B0; prov%iota0 = iota0

    ! Physical scale: choose v0 so r_L/a = ro0 (matches GC ro0 = Larmor scale).
    ! r_L = m c vperp/(q B0). Pick mass/charge of an alpha, derive v from ro0.
    mass   = 4.0_dp * p_mass
    charge = 2.0_dp * e_charge

    r_min  = z0(1)              ! r is already the physical minor radius (<= a)
    lambda = z0(5)
    ! speed unit: r_L(at vperp=vtot) = ro0*a  =>  vtot = ro0*a*q*B0/(m c)
    vtot   = ro0 * a * charge * B0 / (mass * c)
    vpar   = lambda * vtot
    vperp  = sqrt(max(vtot*vtot - vpar*vpar, 0.0_dp))

    ! Start at outer midplane of the flux surface r_min: (R,Z) = (R0+r_min, 0).
    R = R0 + r_min
    Z = 0.0_dp
    u0 = [R, 0.0_dp, Z]

    call prov%eval_field(u0, Bvec, Bmod, hcov, ierrf)
    if (ierrf /= FO_OK) then
      print *, 'analytic full orbit: field eval failed at start'
      return
    end if
    ! orthonormal field unit vector and a perpendicular basis (e1 in R-Z plane).
    b_hat = Bvec / Bmod
    e1 = [0.0_dp, 0.0_dp, 1.0_dp]                 ! e_Z (in poloidal plane)
    e1 = e1 - dot_product(e1, b_hat)*b_hat
    e1 = e1 / sqrt(dot_product(e1, e1))
    e2 = cross(b_hat, e1)
    ! orthonormal velocity, then contravariant (v^phi = v_phi^orth / R).
    w0 = vpar*b_hat + vperp*e1
    w0 = [w0(1), w0(2)/R, w0(3)]                  ! contravariant in (R,phi,Z)

    Omega  = charge * Bmod / (mass * c)
    period = twopi / Omega
    ! resolve the gyration: many substeps per gyroperiod, dump per gyroperiod
    ! sample so the output stride ~ matches the GC step count.
    nsub   = 60
    dt     = period / nsub
    nout   = nstep                                 ! number of gyroperiods to dump

    call init_full_orbit_state(st, u0, w0, ORBIT_BORIS, COORD_CYL, &
                               mass, charge, dt, prov)
    H0 = compute_energy(st)

    open(newunit=u, file=outdir//'/analytic_fullorbit.dat', status='replace')
    write(u, '(A)') '# t  R  Z  Hrel  pphi  mu'
    t = 0.0_dp
    call write_fo_row(u, t, st, H0)
    do it = 1, nout
      do isub = 1, nsub
        call timestep_full_orbit(st, ierr)
        if (ierr /= FO_OK) then
          print '(A,I0)', 'analytic full orbit: stopped at gyroperiod ', it
          close(u)
          return
        end if
      end do
      t = it*period
      call write_fo_row(u, t, st, H0)
    end do
    close(u)
    print '(A,I0,A)', 'analytic full orbit: ', nout, &
        ' gyroperiods -> analytic_fullorbit.dat'
  end subroutine dump_fullorbit_analytic

  subroutine run_vmec(norb, outdir, nstep)
    type(tracer_t), intent(inout) :: norb
    character(*), intent(in) :: outdir
    integer, intent(in) :: nstep
    integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
    real(dp) :: z0(5), dtau, rtol, rbig
    integer :: s

    isw_field_type = BOOZER
    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, GAUSS1)
    call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0e-12_dp)
    if (.not. associated(evaluate)) error stop 'evaluate not associated (boozer)'

    rbig = rmajor*1.0e2_dp
    dtau = twopi*rbig/256.0_dp
    rtol = 1.0e-13_dp
    s = cpp_stages_from_mode(GAUSS1)
    integmode = GAUSS1
    norb%relerr = rtol

    ! Trapped IC (small pitch) on a mid-radius surface.
    z0 = [0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.1_dp]

    call dump_vmec_model(outdir, 'vmec_gc.dat',  z0, dtau, rtol, GAUSS1, s, nstep, .false.)
    call dump_vmec_model(outdir, 'vmec_cpp.dat', z0, dtau, rtol, GAUSS1, s, nstep, .true.)
  end subroutine run_vmec

  subroutine dump_vmec_model(outdir, fname, z0, dtau, rtol, mode, s, nstep, use_cpp)
    character(*), intent(in) :: outdir, fname
    real(dp), intent(in) :: z0(5), dtau, rtol
    integer, intent(in) :: mode, s, nstep
    logical, intent(in) :: use_cpp
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    procedure(transform_i), pointer :: integ2ref_cyl
    real(dp) :: H0, t
    integer :: it, ierr, u

    call init_sympl(si, f, z0, dtau, dtau, rtol, mode)
    if (use_cpp) call orbit_cpp_init(si, f)
    integ2ref_cyl => get_transform('vmec', 'cyl')
    call get_val(f, si%z(4)); H0 = f%H

    open(newunit=u, file=outdir//'/'//fname, status='replace')
    write(u, '(A)') '# t  R  Z  Hrel  pphi  mu'
    call write_vmec_row(u, 0.0_dp, si%z, f%H, H0, si%z(4), f%mu, integ2ref_cyl)
    do it = 1, nstep
      if (use_cpp) then
        call orbit_timestep_cpp(si, f, s, ierr)
      else
        call orbit_timestep_sympl(si, f, ierr)
      end if
      if (ierr /= 0) exit
      call get_val(f, si%z(4))
      t = it*dtau
      call write_vmec_row(u, t, si%z, f%H, H0, si%z(4), f%mu, integ2ref_cyl)
    end do
    close(u)
    print '(A,A,A,I0,A)', 'vmec ', trim(fname), ': ', it-1, ' steps'
  end subroutine dump_vmec_model

  ! --- row writers ---

  ! Analytic test chart: the radial coordinate r is the minor-radius-like
  ! variable used directly in field_can_test (0 <= r <= a = 0.5), so the
  ! circular cross-section maps as R = R0 + r cos th, Z = r sin th.
  subroutine write_can_row(u, t, z, H, H0, pphi, mu, R0)
    integer, intent(in) :: u
    real(dp), intent(in) :: t, z(4), H, H0, pphi, mu, R0
    real(dp) :: Rc, Zc
    Rc = R0 + z(1)*cos(z(2))
    Zc = z(1)*sin(z(2))
    write(u, '(6ES18.9)') t, Rc, Zc, H/H0 - 1.0_dp, pphi, mu
  end subroutine write_can_row

  ! VMEC BOOZER: integ -> ref (VMEC) -> cyl (R,phi,Z) [cm].
  subroutine write_vmec_row(u, t, z, H, H0, pphi, mu, to_cyl)
    integer, intent(in) :: u
    real(dp), intent(in) :: t, z(4), H, H0, pphi, mu
    procedure(transform_i), pointer, intent(in) :: to_cyl
    real(dp) :: ref(3), cyl(3)
    call integ_to_ref(z(1:3), ref)
    call to_cyl(ref, cyl)
    write(u, '(6ES18.9)') t, cyl(1), cyl(3), H/H0 - 1.0_dp, pphi, mu
  end subroutine write_vmec_row

  subroutine write_fo_row(u, t, st, H0)
    integer, intent(in) :: u
    real(dp), intent(in) :: t, H0
    type(FullOrbitState), intent(in) :: st
    real(dp) :: H
    H = compute_energy(st)
    ! pphi and mu not canonical here; report mu_full = m vperp^2/(2B), pphi=0.
    write(u, '(6ES18.9)') t, st%z(1), st%z(3), H/H0 - 1.0_dp, 0.0_dp, st%mu
  end subroutine write_fo_row

  pure function cross(a, b) result(cc)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: cc(3)
    cc(1) = a(2)*b(3) - a(3)*b(2)
    cc(2) = a(3)*b(1) - a(1)*b(3)
    cc(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

end program diag_traj_models
