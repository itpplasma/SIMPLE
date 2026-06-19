program diag_crossval_figures
  ! Dump the three-system cross-validation trajectories and invariant time series
  ! used by examples/plot_crossval_figures.py to regenerate the banana and
  ! conservation panels (analytic tokamak + VMEC). Real integrator output only.
  !
  ! Systems on the SHARED analytic tokamak (field_can_test, R0=1, a=0.5):
  !   GC      : SIMPLE's production symplectic guiding-centre integrator
  !   CPP_SYM : 6D Pauli canonical-midpoint (big dt=80, GC-scale macro-step)
  !   CP      : full charged particle, gyro-resolved (small dt=1)
  ! plus the CPP energy/mu/p_phi conservation sweep over dt for the CP/CPP models.
  !
  ! On VMEC flux coordinates (wout.nc, nfp=2): CP (small dt) and CPP_SYM (big dt),
  ! radial band + Hamiltonian energy time series.
  !
  ! Each file is whitespace columns with a one-line "# " header naming the columns.
  ! Output directory is argv(1) (default ./crossval_data).
  use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
  use field_can_mod, only: field_can_t, field_can_from_name, evaluate
  use orbit_symplectic_base, only: symplectic_integrator_t, GAUSS1
  use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl
  use orbit_cpp_canonical, only: cpp_canon_state_t, cpp_canon_init, cpp_canon_step, &
       cpp_canon_energy, MODEL_CP, MODEL_CPP_SYM, MODEL_CPP_VAR, &
       COORD_TOK, COORD_VMEC
  use orbit_cpp_vmec_metric, only: vmec_metric_init, vmec_metric_ready
  use util, only: c, twopi
  use field_pauli_cart, only: pauli_field_params_t
  use pauli_gc_drift, only: gc_drift_rhs
  use orbit_cpp_pauli, only: pauli6d_state_t, pauli6d_init, pauli6d_step, &
       pauli6d_energy, pauli6d_mu, pauli6d_to_gc
  implicit none

  real(dp), parameter :: mu = 1.0e-5_dp, mass = 1.0_dp, charge = 1.0_dp
  real(dp), parameter :: x0(3) = [0.1_dp, 1.5_dp, 0.0_dp]
  real(dp), parameter :: B0ic = 0.99749164651988482_dp
  character(len=256) :: outdir
  integer :: nargs

  nargs = command_argument_count()
  if (nargs >= 1) then
    call get_command_argument(1, outdir)
  else
    outdir = 'crossval_data'
  end if
  call execute_command_line('mkdir -p '//trim(outdir))

  call dump_analytic_cpp(MODEL_CP, 1.0_dp, 80000, sqrt(2.0_dp*mu*B0ic/mass), &
       trim(outdir)//'/analytic_cp.dat')
  call dump_analytic_cpp(MODEL_CPP_SYM, 80.0_dp, 1000, 0.0_dp, &
       trim(outdir)//'/analytic_cpp_sym.dat')
  call dump_analytic_cpp(MODEL_CPP_VAR, 800.0_dp, 400, 0.0_dp, &
       trim(outdir)//'/analytic_cpp_var.dat')
  call dump_analytic_gc(trim(outdir)//'/analytic_gc.dat')
  call dump_cpp_dt_sweep(trim(outdir)//'/analytic_cpp_sym_dtsweep.dat')

  ! Shared-field banana overlay: the genuine 6D Pauli particle (full gyration,
  ! gyro-fattened) vs the independent GC drift on the SAME analytic Cartesian
  ! tokamak and the SAME trapped launch (test_cpp_pauli_gc_banana parameters).
  ! The Pauli banana wraps the smooth GC banana to O(rho*).
  call dump_pauli_gc_banana(trim(outdir))

  call dump_vmec(trim(outdir))

  write(output_unit,'(A)') 'crossval data written to '//trim(outdir)
contains

  ! One 6D canonical model on the analytic tokamak. Columns:
  !   step  r  theta  phi  Rcyl-R0  Zcyl  E  mu  pphi
  ! Rcyl-R0 = r cos(theta), Zcyl = r sin(theta): the poloidal R-Z section.
  subroutine dump_analytic_cpp(model, dt, nsteps, vperp0, fname)
    integer, intent(in) :: model, nsteps
    real(dp), intent(in) :: dt, vperp0
    character(*), intent(in) :: fname
    type(cpp_canon_state_t) :: st
    real(dp) :: x_init(3), E0, E
    integer :: it, ierr, u

    x_init = x0
    if (model == MODEL_CPP_VAR) x_init(3) = 1.0_dp
    call cpp_canon_init(st, model, COORD_TOK, x_init, 0.0_dp, vperp0, mu, &
                        mass, charge, dt)
    E0 = cpp_canon_energy(st)
    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(A)') '# step r theta phi RmR0 Z E mu pphi'
    call write_row(u, 0, st, E0)
    do it = 1, nsteps
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) exit
      E = cpp_canon_energy(st)
      call write_row(u, it, st, E0)
    end do
    close(u)
    write(output_unit,'(A,I0,A)') '  wrote '//trim(fname)//' (', it-1, ' steps)'
  end subroutine dump_analytic_cpp

  subroutine write_row(u, it, st, E0)
    integer, intent(in) :: u, it
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: E0
    real(dp) :: r, th, E, dEr
    r = st%z(1); th = st%z(2)
    E = cpp_canon_energy(st)
    dEr = (E - E0)/abs(E0)
    write(u,'(I8,1X,8(ES22.14,1X))') it, r, th, st%z(3), &
         r*cos(th), r*sin(th), dEr, st%mu, st%z(6)
  end subroutine write_row

  ! Shared-field banana: genuine 6D Pauli particle (full gyration) and the
  ! independent GC drift on the same analytic Cartesian tokamak (R0=1, a=0.5,
  ! B0=1, iota0=1), same trapped launch as test_cpp_pauli_gc_banana
  ! (r0=0.3, lambda=0.3, rho*=0.04). Writes the GC-reduced minor-radius poloidal
  ! section for the Pauli, and the GC-drift section, plus the Pauli invariants.
  subroutine dump_pauli_gc_banana(outdir)
    character(*), intent(in) :: outdir
    real(dp), parameter :: R0 = 1.0_dp, B0 = 1.0_dp, iota0 = 1.0_dp, a = 0.5_dp
    real(dp), parameter :: r0p = 0.30_dp, lambda = 0.30_dp, v0 = 1.0_dp
    real(dp), parameter :: rhostar = 0.04_dp
    real(dp) :: Bmid, vpar, vperp, mu_b, rho, q

    Bmid  = B0*(1.0_dp - r0p/R0)
    vpar  = lambda*v0
    vperp = sqrt(max(v0*v0 - vpar*vpar, 0.0_dp))
    mu_b  = mass*vperp*vperp/(2.0_dp*Bmid)
    rho   = rhostar*a
    q     = mass*c*vperp/(Bmid*rho)   ! charge fixing rho_L/a = rho* at launch

    call dump_pauli6d(R0, B0, iota0, a, r0p, vpar, vperp, mass, q, &
         trim(outdir)//'/banana_pauli.dat')
    call dump_gc_drift(R0, B0, iota0, a, r0p, vpar, vperp, mu_b, mass, q, &
         trim(outdir)//'/banana_gc_drift.dat')
  end subroutine dump_pauli_gc_banana

  ! 6D Pauli particle, ~16 implicit-midpoint steps per gyration. Columns:
  !   step RmR0 Z dE dmu  (R-R0, Z from the GC-reduced minor radius/poloidal angle)
  subroutine dump_pauli6d(R0, B0, iota0, a, r0p, vpar, vperp, mass, charge, fname)
    real(dp), intent(in) :: R0, B0, iota0, a, r0p, vpar, vperp, mass, charge
    character(*), intent(in) :: fname
    type(pauli6d_state_t) :: st
    type(pauli_field_params_t) :: fp
    real(dp) :: xgc(3), dt, E0, E, mu0, mn, r, th, ph, vp, Omega, period
    integer :: it, ierr, nstep, u

    fp%R0 = R0; fp%B0 = B0; fp%iota0 = iota0; fp%a = a
    xgc = [R0 + r0p, 0.0_dp, 0.0_dp]
    Omega = charge*B0/(mass*c)
    period = twopi/Omega
    dt = period/16.0_dp
    nstep = 60000
    call pauli6d_init(st, fp, xgc, vpar, vperp, mass, charge, dt)
    E0 = pauli6d_energy(st); mu0 = pauli6d_mu(st)
    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(A)') '# step RmR0 Z dE dmu'
    do it = 0, nstep
      if (it > 0) then
        call pauli6d_step(st, ierr)
        if (ierr /= 0) exit
      end if
      call pauli6d_to_gc(st, r, th, ph, vp)
      E = pauli6d_energy(st); mn = pauli6d_mu(st)
      write(u,'(I8,1X,4(ES22.14,1X))') it, r*cos(th), r*sin(th), &
           (E-E0)/abs(E0), (mn-mu0)/abs(mu0)
    end do
    close(u)
    write(output_unit,'(A,I0,A)') '  wrote '//trim(fname)//' (', it-1, ' steps)'
  end subroutine dump_pauli6d

  ! Independent GC drift RK4 on the same Cartesian field. Columns: step RmR0 Z
  subroutine dump_gc_drift(R0, B0, iota0, a, r0p, vpar, vperp, mu_b, mass, charge, fname)
    real(dp), intent(in) :: R0, B0, iota0, a, r0p, vpar, vperp, mu_b, mass, charge
    character(*), intent(in) :: fname
    type(pauli_field_params_t) :: fp
    real(dp) :: Y(4), k1(4), k2(4), k3(4), k4(4), dt, Rcyl, dR, r, th
    integer :: it, nstep, u

    fp%R0 = R0; fp%B0 = B0; fp%iota0 = iota0; fp%a = a
    Y(1:3) = [R0 + r0p, 0.0_dp, 0.0_dp]; Y(4) = vpar
    dt = 2.0e-4_dp; nstep = 200000
    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(A)') '# step RmR0 Z'
    do it = 0, nstep
      if (it > 0) then
        call gc_drift_rhs(fp, mass, charge, vperp, mu_b, Y, k1)
        call gc_drift_rhs(fp, mass, charge, vperp, mu_b, Y + 0.5_dp*dt*k1, k2)
        call gc_drift_rhs(fp, mass, charge, vperp, mu_b, Y + 0.5_dp*dt*k2, k3)
        call gc_drift_rhs(fp, mass, charge, vperp, mu_b, Y + dt*k3, k4)
        Y = Y + dt/6.0_dp*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
      end if
      Rcyl = sqrt(Y(1)*Y(1) + Y(2)*Y(2)); dR = Rcyl - R0
      r = sqrt(dR*dR + Y(3)*Y(3)); th = atan2(Y(3), dR)
      if (mod(it, 20) == 0) write(u,'(I8,1X,2(ES22.14,1X))') it, r*cos(th), r*sin(th)
    end do
    close(u)
    write(output_unit,'(A)') '  wrote '//trim(fname)
  end subroutine dump_gc_drift

  ! SIMPLE production symplectic GC on the same analytic field (field_can_test
  ! chart). Same trapped launch as test_cpp_pauli_gc_banana so the GC banana is
  ! the smooth slow-manifold reference the CPP/CP bananas wrap. Columns:
  !   step r theta phi RmR0 Z dE pphi
  subroutine dump_analytic_gc(fname)
    character(*), intent(in) :: fname
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    real(dp), parameter :: r0p = 0.30_dp, lambda = 0.30_dp, v0 = 1.0_dp
    real(dp), parameter :: a = 0.5_dp, rhostar = 0.04_dp
    real(dp) :: z(4), dt, ro0_gc, H0, H, r, th, pph
    integer :: it, ierr, nstep, u

    call field_can_from_name('test')
    ro0_gc = rhostar*a
    call evaluate(f, r0p, 0.0_dp, 0.0_dp, 0)
    f%mu  = 0.5_dp*v0*v0*(1.0_dp - lambda*lambda)/f%Bmod
    f%ro0 = ro0_gc
    f%vpar = v0*lambda
    z(1) = r0p; z(2) = 0.0_dp; z(3) = 0.0_dp
    z(4) = f%vpar*f%hph + f%Aph/f%ro0
    dt = 1.0e-3_dp
    nstep = 20000
    call orbit_sympl_init(si, f, z, dt, 1, 1.0e-13_dp, GAUSS1)
    call gc_invariants(f, si%z, ro0_gc, H0, pph)
    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(A)') '# step r theta phi RmR0 Z dE pphi'
    do it = 0, nstep
      if (it > 0) then
        call orbit_timestep_sympl(si, f, ierr)
        if (ierr /= 0) exit
      end if
      r = si%z(1); th = si%z(2)
      call gc_invariants(f, si%z, ro0_gc, H, pph)
      write(u,'(I8,1X,7(ES22.14,1X))') it, r, th, si%z(3), &
           r*cos(th), r*sin(th), (H-H0)/abs(H0), pph
    end do
    close(u)
    write(output_unit,'(A,I0,A)') '  wrote '//trim(fname)//' (', it-1, ' steps)'
  end subroutine dump_analytic_gc

  ! GC Hamiltonian H = vpar^2/2 + mu B and canonical p_phi = pph (z(4)*ro0 form).
  subroutine gc_invariants(f, z, ro0, H, pph)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: z(4), ro0
    real(dp), intent(out) :: H, pph
    real(dp) :: vpar
    call evaluate(f, z(1), z(2), z(3), 0)
    vpar = (z(4) - f%Aph/ro0)/f%hph
    H = 0.5_dp*vpar*vpar + f%mu*f%Bmod
    pph = z(4)
  end subroutine gc_invariants

  ! CP/CPP-sym energy + p_phi conservation as a function of dt: max relative
  ! energy band and max |dp_phi| over a fixed physical time. Columns:
  !   dt  maxdE  maxdpphi
  subroutine dump_cpp_dt_sweep(fname)
    character(*), intent(in) :: fname
    real(dp) :: dts(5), dt, E0, E, Emin, Emax, pph0, pphdev
    type(cpp_canon_state_t) :: st
    integer :: i, it, ierr, n, u
    dts = [10.0_dp, 20.0_dp, 40.0_dp, 80.0_dp, 160.0_dp]
    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(A)') '# dt maxdE maxdpphi'
    do i = 1, 5
      dt = dts(i)
      call cpp_canon_init(st, MODEL_CPP_SYM, COORD_TOK, x0, 0.0_dp, 0.0_dp, mu, &
                          mass, charge, dt)
      E0 = cpp_canon_energy(st); Emin = E0; Emax = E0
      pph0 = st%z(6); pphdev = 0.0_dp
      n = nint(80000.0_dp/dt)
      do it = 1, n
        call cpp_canon_step(st, ierr)
        if (ierr /= 0) exit
        E = cpp_canon_energy(st); Emin = min(Emin, E); Emax = max(Emax, E)
        pphdev = max(pphdev, abs(st%z(6) - pph0))
      end do
      write(u,'(3(ES22.14,1X))') dt, (Emax-Emin)/abs(E0), pphdev
    end do
    close(u)
    write(output_unit,'(A)') '  wrote '//trim(fname)
  end subroutine dump_cpp_dt_sweep

  ! VMEC flux-coordinate CP (small dt) and CPP-sym (big dt) on wout.nc. Columns:
  !   step s theta phi sqrts_cosT sqrts_sinT dE
  ! sqrt(s) cos/sin theta gives a flux-surface poloidal proxy section.
  subroutine dump_vmec(outdir)
    character(*), intent(in) :: outdir
    type(cpp_canon_state_t) :: st
    real(dp), parameter :: xv0(3) = [0.3_dp, 0.6_dp, 0.2_dp]
    real(dp) :: E0, E, s, th
    integer :: it, ierr, u
    logical :: have_wout

    inquire(file='wout.nc', exist=have_wout)
    if (.not. have_wout) then
      write(output_unit,'(A)') '  wout.nc not in cwd; skipping VMEC dump'
      return
    end if
    call vmec_metric_init('wout.nc')
    if (.not. vmec_metric_ready()) then
      write(output_unit,'(A)') '  VMEC metric not ready; skip VMEC'
      return
    end if

    ! CP gyro-resolved (small dt).
    call cpp_canon_init(st, MODEL_CP, COORD_VMEC, xv0, 0.0_dp, 3.0e5_dp, 0.0_dp, &
                        mass, charge, 2.0e-8_dp)
    E0 = cpp_canon_energy(st)
    open(newunit=u, file=trim(outdir)//'/vmec_cp.dat', status='replace', action='write')
    write(u,'(A)') '# step s theta phi sqrts_cosT sqrts_sinT dE'
    do it = 0, 2000
      if (it > 0) then
        call cpp_canon_step(st, ierr)
        if (ierr /= 0) exit
      end if
      s = st%z(1); th = st%z(2); E = cpp_canon_energy(st)
      write(u,'(I8,1X,6(ES22.14,1X))') it, s, th, st%z(3), &
           sqrt(max(s,0.0_dp))*cos(th), sqrt(max(s,0.0_dp))*sin(th), (E-E0)/abs(E0)
    end do
    close(u)
    write(output_unit,'(A,I0,A)') '  wrote '//trim(outdir)//'/vmec_cp.dat (', it-1, ' steps)'

    ! CPP-sym Pauli (big GC-scale dt).
    call cpp_canon_init(st, MODEL_CPP_SYM, COORD_VMEC, xv0, 1.0e5_dp, 0.0_dp, &
                        1.0e-3_dp, mass, charge, 5.0e-7_dp)
    E0 = cpp_canon_energy(st)
    open(newunit=u, file=trim(outdir)//'/vmec_cpp_sym.dat', status='replace', action='write')
    write(u,'(A)') '# step s theta phi sqrts_cosT sqrts_sinT dE'
    do it = 0, 1000
      if (it > 0) then
        call cpp_canon_step(st, ierr)
        if (ierr /= 0) exit
      end if
      s = st%z(1); th = st%z(2); E = cpp_canon_energy(st)
      write(u,'(I8,1X,6(ES22.14,1X))') it, s, th, st%z(3), &
           sqrt(max(s,0.0_dp))*cos(th), sqrt(max(s,0.0_dp))*sin(th), (E-E0)/abs(E0)
    end do
    close(u)
    write(output_unit,'(A,I0,A)') '  wrote '//trim(outdir)//'/vmec_cpp_sym.dat (', it-1, ' steps)'
  end subroutine dump_vmec

end program diag_crossval_figures
