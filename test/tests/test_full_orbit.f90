program test_full_orbit
  ! Behavioral tests for the full-orbit Boris pusher using only the analytic
  ! Cartesian and cylindrical mock providers (no libneo field). Oracles:
  !   1. uniform-B exact gyration: |v|, vpar, energy, mu conserved; closed
  !      circle of Larmor radius r_L = m c vperp/(q B), period T = 2pi m c/(qB).
  !   2. Cartesian linear grad-B drift vs analytic v_gradB.
  !   3. cylindrical 1/R curvature + grad-B drift vs analytic v_d, separated
  !      into pure curvature (vperp=0) and pure grad-B (vpar=0).
  !   4. mu adiabatic invariance over many gyroperiods.
  !   5. cylindrical Christoffel mock vs closed form.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: c, twopi, pi, p_mass, e_charge
  use orbit_full, only: FullOrbitState, init_full_orbit_state, &
      timestep_full_orbit, convert_full_to_gc, compute_energy, &
      ORBIT_BORIS, COORD_CART, COORD_CYL, FO_OK
  use orbit_full_mock_cart, only: cartesian_provider_t, FIELD_UNIFORM, FIELD_LINGRAD
  use orbit_full_mock_cyl, only: cylindrical_provider_t
  implicit none

  integer :: nfail
  nfail = 0

  call test_uniform_gyration(nfail)
  call test_cart_gradb_drift(nfail)
  call test_cyl_curvature_drift(nfail)
  call test_cyl_gradb_drift(nfail)
  call test_mu_invariance(nfail)
  call test_cyl_christoffel(nfail)

  if (nfail == 0) then
    print *, 'ALL FULL-ORBIT TESTS PASSED'
  else
    print *, 'FULL-ORBIT TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine check(name, ok, nfail)
    character(*), intent(in) :: name
    logical, intent(in) :: ok
    integer, intent(inout) :: nfail
    if (ok) then
      print '(A,A)', 'PASS  ', name
    else
      print '(A,A)', 'FAIL  ', name
      nfail = nfail + 1
    end if
  end subroutine check

  subroutine test_uniform_gyration(nfail)
    integer, intent(inout) :: nfail
    type(cartesian_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, B0, vperp, vpar, speed0
    real(dp) :: Omega, period, dt, rL
    real(dp) :: x0(3), v0(3), xstart(3)
    real(dp) :: e0, e1, mu0, vz0, errpos, errv
    integer :: i, nstep, ierr

    mass   = 4.0_dp * p_mass
    charge = 2.0_dp * e_charge
    B0     = 1.0d4
    prov%field_kind = FIELD_UNIFORM
    prov%B0 = [0.0_dp, 0.0_dp, B0]

    vperp = 1.0d7
    vpar  = 3.0d6
    speed0 = sqrt(vperp**2 + vpar**2)
    Omega = charge * B0 / (mass * c)
    period = twopi / Omega
    rL = mass * c * vperp / (charge * B0)
    nstep = 400
    dt = period / nstep

    x0 = [0.0_dp, 0.0_dp, 0.0_dp]
    v0 = [vperp, 0.0_dp, vpar]
    call init_full_orbit_state(st, x0, v0, ORBIT_BORIS, COORD_CART, &
                               mass, charge, dt, prov)
    xstart = st%z(1:3)
    e0  = compute_energy(st)
    mu0 = st%mu
    vz0 = st%z(6)

    do i = 1, nstep
      call timestep_full_orbit(st, ierr)
      if (ierr /= FO_OK) then
        call check('uniform: timestep ierr', .false., nfail)
        return
      end if
    end do

    e1 = compute_energy(st)
    errv = abs(sqrt(dot_product(st%z(4:6), st%z(4:6))) - speed0) / speed0
    ! transverse return: z advances by vpar*period, so only (x,y) close.
    errpos = sqrt((st%z(1)-xstart(1))**2 + (st%z(2)-xstart(2))**2)

    print '(A,ES12.4,A,ES12.4)', '  uniform: r_L=', rL, ' period=', period
    print '(A,ES12.4,A,ES12.4)', '  uniform: |v| relerr=', errv, &
        ' return-pos err=', errpos
    print '(A,ES12.4)', '  uniform: dE/E=', abs(e1-e0)/e0
    print '(A,ES12.4)', '  uniform: dmu/mu=', abs(st%mu-mu0)/mu0

    call check('uniform: |v| constant', errv < 1d-10, nfail)
    call check('uniform: energy constant', abs(e1-e0)/e0 < 1d-9, nfail)
    call check('uniform: vpar=vz constant', abs(st%z(6)-vz0) < 1d-6*speed0, nfail)
    ! closed circle: error scales with dt; Boris is 2nd order -> bound on rL.
    call check('uniform: return to start', errpos < 1d-3*rL, nfail)
    call check('uniform: mu constant', abs(st%mu-mu0)/mu0 < 1d-9, nfail)
  end subroutine test_uniform_gyration

  subroutine test_cart_gradb_drift(nfail)
    integer, intent(inout) :: nfail
    type(cartesian_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, B0, g, vperp, vpar
    real(dp) :: Omega, period, dt, vd_exact, vd_meas
    real(dp) :: x0(3), v0(3), ygc0, ygc1, t_total
    integer :: i, nper_run, nstep_per, nstep, ierr

    mass   = 4.0_dp * p_mass
    charge = 2.0_dp * e_charge
    B0     = 1.0d4
    g      = 1.0d2           ! G/cm gradient of B_z along x
    prov%field_kind = FIELD_LINGRAD
    prov%B0 = [0.0_dp, 0.0_dp, B0]
    prov%gradB = 0.0_dp
    prov%gradB(3,1) = g       ! B_z = B0 + g*x

    vperp = 1.0d7
    vpar  = 0.0_dp
    Omega = charge * B0 / (mass * c)
    period = twopi / Omega
    nstep_per = 200
    nper_run = 2000
    nstep = nstep_per * nper_run
    dt = period / nstep_per

    ! analytic grad-B drift: v = (m c vperp^2)/(2 q B0) * (g/B0), along +e_y for
    ! q>0, B along +z, grad|B| along +x.
    vd_exact = mass * c * vperp**2 / (2.0_dp * charge * B0) * (g / B0)

    x0 = [0.0_dp, 0.0_dp, 0.0_dp]
    v0 = [vperp, 0.0_dp, vpar]
    call init_full_orbit_state(st, x0, v0, ORBIT_BORIS, COORD_CART, &
                               mass, charge, dt, prov)
    ygc0 = guiding_center_y_cart(st)
    do i = 1, nstep
      call timestep_full_orbit(st, ierr)
      if (ierr /= FO_OK) then
        call check('cart gradB: timestep ierr', .false., nfail)
        return
      end if
    end do
    ygc1 = guiding_center_y_cart(st)
    t_total = nstep * dt
    ! guiding-center y-drift removes the gyration, leaving the secular drift.
    vd_meas = (ygc1 - ygc0) / t_total

    print '(A,ES12.4,A,ES12.4)', '  cart gradB: vd_exact=', vd_exact, &
        ' vd_meas=', vd_meas
    call check('cart gradB: drift sign/magnitude', &
        abs(vd_meas - vd_exact) < 0.05_dp*abs(vd_exact), nfail)
  end subroutine test_cart_gradb_drift

  ! Guiding-center y from the instantaneous Cartesian state:
  ! x_gc = x - (1/Omega) (b x v), Omega = qB/(mc), b = B/|B|.
  function guiding_center_y_cart(st) result(ygc)
    type(FullOrbitState), intent(in) :: st
    real(dp) :: ygc
    real(dp) :: Bvec(3), Bmod, hcov(3), Omega, rho(3)
    integer :: ierr
    call st%prov%eval_field(st%z(1:3), Bvec, Bmod, hcov, ierr)
    Omega = st%qm * Bmod / c
    rho = cross_local(hcov, st%z(4:6)) / Omega
    ygc = st%z(2) - rho(2)
  end function guiding_center_y_cart

  pure function cross_local(a, b) result(cc)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: cc(3)
    cc(1) = a(2)*b(3) - a(3)*b(2)
    cc(2) = a(3)*b(1) - a(1)*b(3)
    cc(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_local

  subroutine test_cyl_curvature_drift(nfail)
    ! Pure curvature: vperp -> 0 (small), vpar finite. v_d = (mc/qBR) vpar^2.
    integer, intent(inout) :: nfail
    call run_cyl_drift(nfail, 'cyl curvature', vpar_in=1.0d7, vperp_in=1.0d5)
  end subroutine test_cyl_curvature_drift

  subroutine test_cyl_gradb_drift(nfail)
    ! Pure grad-B: vpar -> 0 (small), vperp finite. v_d = (mc/qBR)(vperp^2/2).
    integer, intent(inout) :: nfail
    call run_cyl_drift(nfail, 'cyl gradB', vpar_in=1.0d5, vperp_in=1.0d7)
  end subroutine test_cyl_gradb_drift

  subroutine run_cyl_drift(nfail, tag, vpar_in, vperp_in)
    integer, intent(inout) :: nfail
    character(*), intent(in) :: tag
    real(dp), intent(in) :: vpar_in, vperp_in
    type(cylindrical_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, B0, R0, R, Bloc
    real(dp) :: Omega, period, dt, vd_exact, vd_meas
    real(dp) :: u0(3), w0(3), z0, t_total
    integer :: i, nstep_per, nper_run, nstep, ierr

    mass   = 4.0_dp * p_mass
    charge = 2.0_dp * e_charge
    B0     = 1.0d4
    R0     = 200.0_dp
    R      = 200.0_dp
    prov%B0 = B0
    prov%R0 = R0

    Bloc = B0 * R0 / R       ! = B0 here since R=R0
    Omega = charge * Bloc / (mass * c)
    period = twopi / Omega
    nstep_per = 300
    nper_run = 40
    nstep = nstep_per * nper_run
    dt = period / nstep_per

    ! v_d = (m c)/(q B R) * (vpar^2 + vperp^2/2), along +Z for q>0, B toroidal.
    vd_exact = (mass * c) / (charge * Bloc * R) * &
               (vpar_in**2 + 0.5_dp * vperp_in**2)

    ! contravariant velocity in (R,phi,Z): orthonormal vphi_phys=vpar (toroidal
    ! is field direction); vperp put into v^R (radial). v^phi = vpar/R.
    u0 = [R, 0.0_dp, 0.0_dp]
    w0 = [vperp_in, vpar_in / R, 0.0_dp]
    call init_full_orbit_state(st, u0, w0, ORBIT_BORIS, COORD_CYL, &
                               mass, charge, dt, prov)
    z0 = st%z(3)
    do i = 1, nstep
      call timestep_full_orbit(st, ierr)
      if (ierr /= FO_OK) then
        call check(tag//': timestep ierr', .false., nfail)
        return
      end if
    end do
    t_total = nstep * dt
    vd_meas = (st%z(3) - z0) / t_total

    print '(A,A,A,ES12.4,A,ES12.4)', '  ', tag, ': vd_exact=', vd_exact, &
        ' vd_meas=', vd_meas
    call check(tag//': vertical drift', &
        abs(vd_meas - vd_exact) < 0.10_dp*abs(vd_exact), nfail)
  end subroutine run_cyl_drift

  subroutine test_mu_invariance(nfail)
    integer, intent(inout) :: nfail
    type(cartesian_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, B0, g, vperp, vpar
    real(dp) :: Omega, period, dt, mu0, mumax_dev, mu_now
    real(dp) :: Bvec(3), Bmod, hcov(3), vperp2, vpar_now
    real(dp) :: x0(3), v0(3)
    integer :: i, nstep, ierr

    mass   = 4.0_dp * p_mass
    charge = 2.0_dp * e_charge
    B0     = 1.0d4
    g      = 1.0d0
    prov%field_kind = FIELD_LINGRAD
    prov%B0 = [0.0_dp, 0.0_dp, B0]
    prov%gradB = 0.0_dp
    prov%gradB(3,1) = g

    vperp = 1.0d7
    vpar  = 2.0d6
    Omega = charge * B0 / (mass * c)
    period = twopi / Omega
    nstep = 200 * 60
    dt = period / 200

    x0 = [0.0_dp, 0.0_dp, 0.0_dp]
    v0 = [vperp, 0.0_dp, vpar]
    call init_full_orbit_state(st, x0, v0, ORBIT_BORIS, COORD_CART, &
                               mass, charge, dt, prov)
    mu0 = st%mu
    mumax_dev = 0.0_dp
    do i = 1, nstep
      call timestep_full_orbit(st, ierr)
      if (ierr /= FO_OK) then
        call check('mu inv: timestep ierr', .false., nfail)
        return
      end if
      call prov%eval_field(st%z(1:3), Bvec, Bmod, hcov, ierr)
      vpar_now = dot_product(st%z(4:6), hcov)
      vperp2 = dot_product(st%z(4:6), st%z(4:6)) - vpar_now**2
      mu_now = mass * vperp2 / (2.0_dp * Bmod)
      mumax_dev = max(mumax_dev, abs(mu_now - mu0) / mu0)
    end do

    print '(A,ES12.4)', '  mu inv: max rel deviation=', mumax_dev
    call check('mu adiabatic invariance', mumax_dev < 1d-3, nfail)
  end subroutine test_mu_invariance

  subroutine test_cyl_christoffel(nfail)
    integer, intent(inout) :: nfail
    type(cylindrical_provider_t) :: prov
    real(dp) :: Gamma(3,3,3), x(3), R
    real(dp) :: gfd(3,3,3), err
    logical :: ok

    prov%B0 = 1.0_dp
    prov%R0 = 1.0_dp
    R = 1.7_dp
    x = [R, 0.3_dp, -0.5_dp]
    call prov%christoffel(x, Gamma)

    ! closed-form check
    ok = abs(Gamma(1,2,2) - (-R)) < 1d-12 .and. &
         abs(Gamma(2,1,2) - 1.0_dp/R) < 1d-12 .and. &
         abs(Gamma(2,2,1) - 1.0_dp/R) < 1d-12
    call check('cyl christoffel: closed form entries', ok, nfail)

    ! symmetry Gamma^i_{mn} = Gamma^i_{nm}
    ok = christoffel_symmetric(Gamma)
    call check('cyl christoffel: symmetry', ok, nfail)

    ! all other entries zero
    call check('cyl christoffel: only known nonzeros', &
        only_known_nonzero(Gamma, R), nfail)

    ! finite-difference metric -> Gamma, compare to closed form
    call christoffel_fd(prov, x, gfd)
    err = maxval(abs(gfd - Gamma))
    print '(A,ES12.4)', '  cyl christoffel: max FD-vs-analytic err=', err
    call check('cyl christoffel: FD agrees', err < 1d-5, nfail)
  end subroutine test_cyl_christoffel

  logical function christoffel_symmetric(Gamma) result(ok)
    real(dp), intent(in) :: Gamma(3,3,3)
    integer :: i, m, n
    ok = .true.
    do i = 1, 3
      do m = 1, 3
        do n = 1, 3
          if (abs(Gamma(i,m,n) - Gamma(i,n,m)) > 1d-14) ok = .false.
        end do
      end do
    end do
  end function christoffel_symmetric

  logical function only_known_nonzero(Gamma, R) result(ok)
    real(dp), intent(in) :: Gamma(3,3,3), R
    real(dp) :: g(3,3,3)
    g = Gamma
    g(1,2,2) = 0.0_dp
    g(2,1,2) = 0.0_dp
    g(2,2,1) = 0.0_dp
    ok = maxval(abs(g)) < 1d-14
  end function only_known_nonzero

  subroutine christoffel_fd(prov, x, gfd)
    type(cylindrical_provider_t), intent(in) :: prov
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: gfd(3,3,3)
    real(dp) :: h, gp(3,3), gm(3,3), ginv(3,3), sqrtg
    real(dp) :: dg(3,3,3), gnow(3,3), xp(3)
    integer :: i, j, k, l
    h = 1d-4

    call prov%metric(x, gnow, ginv, sqrtg)
    ! dg(k,i,j) = d g_ij / d x_k
    do k = 1, 3
      xp = x; xp(k) = x(k) + h
      call prov%metric(xp, gp, ginv, sqrtg)
      xp = x; xp(k) = x(k) - h
      call prov%metric(xp, gm, ginv, sqrtg)
      dg(k,:,:) = (gp - gm) / (2.0_dp*h)
    end do
    call prov%metric(x, gnow, ginv, sqrtg)

    gfd = 0.0_dp
    do i = 1, 3
      do j = 1, 3       ! j = m
        do k = 1, 3     ! k = n
          do l = 1, 3
            gfd(i,j,k) = gfd(i,j,k) + 0.5_dp*ginv(i,l) * &
                (dg(j,l,k) + dg(k,l,j) - dg(l,j,k))
          end do
        end do
      end do
    end do
  end subroutine christoffel_fd

end program test_full_orbit
