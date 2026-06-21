module orbit_cpp_boris
  ! Cartesian explicit pusher for the classical particle (CP) and classical Pauli
  ! particle (CPP), the production 6D full-orbit / large-step path (issue #420).
  ! The particle advances in Cartesian (x, v): the magnetic rotation is exact for
  ! constant B over a step, the kinetic metric is the identity, the geodesic terms
  ! vanish, and the magnetic axis is an ordinary point. CP (pauli=.false.) is the
  ! plain Boris drift-rotate-drift; CPP (pauli=.true.) adds the half mirror kicks
  ! v -= 0.5*dt*(mu/m)*grad|B| around the rotation on the regular Pauli Lagrangian
  ! H = |v|^2/2 + mu|B|, with frozen mu and an optional rotation-angle filter.
  !
  ! Field and geometry come from the chartmap (the Cartesian-side representation,
  ! issue #420): at the Cartesian point we invert to the logical chart
  ! u=(rho, theta_B, phi_B) with the chartmap forward map (ref_coords%from_cart,
  ! rho=sqrt(s)), evaluate the production Boozer field there (chartmap_eval_field:
  ! |B|, the covariant field direction h_i, d|B|/du_i), and push the physical
  ! vectors to Cartesian with the chartmap Jacobian Jc = d(x)/du (covariant_basis)
  ! and inverse metric g^{ij} (metric_tensor):
  !   B_cart = Jc (|B| g^{ij} h_j),   grad|B|_cart = Jc^{-T} d|B|/du.
  ! The chartmap also owns the loss boundary: from_cart flags rho>=1 (out of the
  ! s<1 plasma) -- the ONLY confinement loss. A field-locate non-convergence is a
  ! numerical fault, retried/reported, never a loss (#419, #420).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use reference_coordinates, only: ref_coords
  use orbit_cpp_chartmap_metric, only: chartmap_eval_field
  use libneo_coordinates, only: chartmap_coordinate_system_t
  implicit none
  private

  real(dp), parameter :: c = 1.0_dp
  real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

  ! cart_field / locate status: regular interior point, physical edge loss, or a
  ! numerical locate fault (NOT a loss).
  integer, parameter, public :: CPB_OK = 0, CPB_LOSS = 1, CPB_LOCATE_FAIL = 2

  public :: cpp_boris_state_t, cpp_boris_init, cpp_boris_step, cpp_boris_energy, &
            cpp_boris_mu, cpp_boris_to_gc

  type :: cpp_boris_state_t
    real(dp) :: x(3)   = 0.0_dp    ! Cartesian position (scaled cm)
    real(dp) :: v(3)   = 0.0_dp    ! Cartesian velocity (normalized)
    real(dp) :: u(3)   = 0.0_dp    ! last logical (rho, theta_B, phi_B)
    real(dp) :: mu     = 0.0_dp    ! magnetic moment parameter
    real(dp) :: dt     = 0.0_dp
    real(dp) :: mass   = 1.0_dp
    real(dp) :: charge = 1.0_dp
    real(dp) :: ro0    = 1.0_dp
    real(dp) :: pabs   = 0.0_dp    ! normalized speed (carried for z(4) write-back)
    logical  :: pauli  = .true.    ! .true. CPP (+mu|B|); .false. CP (full orbit)
    logical  :: filtered = .false. ! HLW large-step rotation filter
  end type cpp_boris_state_t

contains

  ! Cartesian -> logical chart (rho, theta_B, phi_B) by a warm-started Newton on the
  ! chartmap forward map x(u) = evaluate_cart(u), Jc = covariant_basis(u). The carried
  ! u_guess (previous substep) is a Larmor-step away, so Newton converges in 1-2 cheap
  ! iterations -- far more robust and faster than the cold multi-seed from_cart, and
  ! thread-safe (read-only spline evaluation). status: CPB_OK interior, CPB_LOSS when
  ! the point maps to rho>=1 (out of the s<1 plasma), CPB_LOCATE_FAIL on no convergence.
  subroutine invert_cart_warm(x, u_guess, u, status)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: u(3)
    integer, intent(out) :: status
    integer, parameter :: maxit = 30, maxls = 30
    ! The forward map is a deterministic spline, so a damped Newton on the wedge
    ! point converges to ~machine precision; tol targets that. accept_tol only
    ! classifies a Newton that has stalled at the spline floor.
    real(dp), parameter :: tol = 1.0e-9_dp, accept_tol = 1.0e-6_dp, rho_edge = 1.0_dp
    real(dp) :: xc(3), Jc(3,3), Jinv(3,3), res(3), du(3), ut(3), rn, rnew, alpha
    integer :: it, ls, i

    u = u_guess
    call ref_coords%evaluate_cart(u, xc)
    res = xc - x
    rn = sqrt(res(1)**2 + res(2)**2 + res(3)**2)
    do it = 1, maxit
      if (rn < tol) then
        status = accept_or_fail(u(1), rn, accept_tol, rho_edge)
        return
      end if
      call ref_coords%covariant_basis(u, Jc)
      if (.not. jacobian_ok(Jc)) then   ! near-axis singular chart: bail, not a loss
        status = CPB_LOCATE_FAIL; return
      end if
      call inv3(Jc, Jinv)
      do i = 1, 3
        du(i) = -(Jinv(i,1)*res(1) + Jinv(i,2)*res(2) + Jinv(i,3)*res(3))
      end do
      ! backtracking line search: Newton is not monotonic for a finite offset.
      alpha = 1.0_dp
      do ls = 1, maxls
        ut = u + alpha*du
        if (ut(1) < 0.0_dp) ut(1) = -ut(1)        ! reflect through the axis
        ! A trial overshoot past rho=1 is NOT a loss: evaluate_cart clamps rho to
        ! the grid edge so an interior target yields a large residual and the step
        ! backtracks. Loss is decided only on the converged rho (accept_or_fail);
        ! a genuinely-outside particle converges to rho~1 and is flagged there.
        call ref_coords%evaluate_cart(ut, xc)
        res = xc - x
        rnew = sqrt(res(1)**2 + res(2)**2 + res(3)**2)
        if (rnew < rn) exit
        alpha = 0.5_dp*alpha
      end do
      if (rnew >= rn) then   ! line search could not improve -> stalled at the floor
        status = accept_or_fail(u(1), rn, accept_tol, rho_edge)
        return
      end if
      u = ut
      rn = rnew
    end do
    status = accept_or_fail(u(1), rn, accept_tol, rho_edge)
  end subroutine invert_cart_warm

  ! Classify a finished Newton. The inverse itself never declares a confinement
  ! loss: a converged locate (interior or just past the clamped edge) is CPB_OK and
  ! reports the radius through u, so the caller decides loss on the guiding-centre
  ! radius. A point that stalls at the edge is still OK (its rho>=1 is reported); an
  ! interior point that cannot converge is a numerical fault (counted confined).
  pure integer function accept_or_fail(rho, rn, accept_tol, rho_edge) result(status)
    real(dp), intent(in) :: rho, rn, accept_tol, rho_edge
    if (rn < accept_tol .or. rho >= rho_edge - 1.0e-3_dp) then
      status = CPB_OK
    else
      status = CPB_LOCATE_FAIL
    end if
  end function accept_or_fail

  ! Geometric field period 2*pi/nfp; the device is exactly nfp-fold symmetric about
  ! Z, so a rotation by this angle maps one field period onto the next.
  real(dp) function field_period()
    integer :: nfp
    nfp = 1
    select type (cs => ref_coords)
    class is (chartmap_coordinate_system_t)
      nfp = cs%num_field_periods
    end select
    field_period = twopi/real(max(nfp, 1), dp)
  end function field_period

  pure function rotz(v, ca, sa) result(w)
    real(dp), intent(in) :: v(3), ca, sa
    real(dp) :: w(3)
    w(1) = ca*v(1) - sa*v(2)
    w(2) = sa*v(1) + ca*v(2)
    w(3) = v(3)
  end function rotz

  ! Map a global Cartesian point into the fundamental field-period wedge by an
  ! integer rotation about Z. The chartmap stores geometry only on one period, so
  ! the inversion and field evaluation run in the wedge; (ca, sa) rotate the wedge
  ! field vectors back to the global frame. This is what lets the Cartesian inverse
  ! converge to machine precision on a multi-period (nfp>1) device.
  subroutine to_wedge(x, xw, ca, sa)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: xw(3), ca, sa
    real(dp) :: phi, period, alpha
    period = field_period()
    phi = atan2(x(2), x(1))
    alpha = period*floor(phi/period)
    ca = cos(alpha); sa = sin(alpha)
    xw = rotz(x, ca, -sa)     ! rotate by -alpha into the wedge
  end subroutine to_wedge

  ! Cartesian B, |B|, grad|B| at logical u from the chartmap field (field_can) and
  ! geometry (ref_coords): B^i = |B| g^{ij} h_j, B_cart = Jc B^i; grad|B| covariant
  ! d|B|/du -> Cartesian by Jc^{-T}. Jc returned for downstream Larmor offsets.
  subroutine field_at_logical(u, Bvec, Bmod, gradB, Jc, status)
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: Bvec(3), Bmod, gradB(3), Jc(3,3)
    integer, intent(out) :: status
    real(dp) :: ue(3), Acov(3), dA(3,3), dBmod(3), hcov(3)
    real(dp) :: g(3,3), ginv(3,3), sqrtg, Bctr(3), Jinv(3,3)
    integer :: i

    ! A particle may gyro-excurse a Larmor radius past s=1; evaluate the field at
    ! the clamped edge there (field_can is undefined past the last closed surface).
    ue = u
    ue(1) = min(ue(1), 1.0_dp - 1.0e-9_dp)
    call chartmap_eval_field(ue, Acov, dA, Bmod, dBmod, hcov)
    call ref_coords%metric_tensor(ue, g, ginv, sqrtg)
    call ref_coords%covariant_basis(ue, Jc)
    if (.not. jacobian_ok(Jc)) then   ! near-axis singular chart
      status = CPB_LOCATE_FAIL; return
    end if
    ! B = curl A with the Boozer flux-function potential A = (0, A_theta(s),
    ! A_phi(s)). Then B^s = d_theta A_phi - d_phi A_theta = 0 EXACTLY (B tangent to
    ! the flux surface), B^theta = -dA_phi/drho / sqrtg, B^phi = dA_theta/drho /
    ! sqrtg, carrying the exact equilibrium pitch. This is divergence-free by
    ! construction; raising the unit field h with the metric instead left a spurious
    ! B^s (radial streaming) that drove the CP over-loss. Renormalize to |B|.
    Bctr(1) = 0.0_dp
    Bctr(2) = -dA(3,1)/sqrtg
    Bctr(3) = dA(2,1)/sqrtg
    do i = 1, 3
      Bvec(i) = Jc(i,1)*Bctr(1) + Jc(i,2)*Bctr(2) + Jc(i,3)*Bctr(3)
    end do
    Bvec = Bvec*(Bmod/max(sqrt(Bvec(1)**2 + Bvec(2)**2 + Bvec(3)**2), 1.0e-30_dp))
    call inv3(Jc, Jinv)
    do i = 1, 3
      gradB(i) = Jinv(1,i)*dBmod(1) + Jinv(2,i)*dBmod(2) + Jinv(3,i)*dBmod(3)
    end do
    status = CPB_OK
  end subroutine field_at_logical

  ! Cartesian B vector, |B|, and grad|B| at Cartesian x, from the chartmap field.
  ! status: CPB_OK (interior, u_out valid), CPB_LOSS (rho>=1 edge loss),
  ! CPB_LOCATE_FAIL (numerical inversion fault). On loss/fault Bvec etc. are
  ! undefined and the caller must not push.
  subroutine cart_field(x, u_guess, Bvec, Bmod, gradB, u_out, status)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: Bvec(3), Bmod, gradB(3), u_out(3)
    integer, intent(out) :: status
    real(dp) :: xw(3), ca, sa, u(3), Jc(3,3), Bw(3), gw(3)

    call to_wedge(x, xw, ca, sa)
    call invert_cart_warm(xw, u_guess, u, status)
    if (status /= CPB_OK) return
    u_out = u
    call field_at_logical(u, Bw, Bmod, gw, Jc, status)
    if (status /= CPB_OK) return
    Bvec = rotz(Bw, ca, sa)       ! wedge field vector -> global frame
    gradB = rotz(gw, ca, sa)
  end subroutine cart_field

  ! Logical chart of a Cartesian point and the local covariant frame, for seeding
  ! and Larmor offsets. status as in cart_field. u_guess warm-starts the inversion.
  subroutine locate(x, u_guess, u_out, bhat, eperp, Bmod, status)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: u_out(3), bhat(3), eperp(3), Bmod
    integer, intent(out) :: status
    real(dp) :: xw(3), ca, sa, u(3), Jc(3,3), g(3,3), ginv(3,3), sqrtg
    real(dp) :: Acov(3), dA(3,3), dBmod(3), hcov(3), hctr(3), eperp_u(3)
    real(dp) :: bw(3), ew(3)
    integer :: i

    call to_wedge(x, xw, ca, sa)
    call invert_cart_warm(xw, u_guess, u, status)
    if (status /= CPB_OK) return
    u_out = u

    call chartmap_eval_field(u, Acov, dA, Bmod, dBmod, hcov)
    call ref_coords%covariant_basis(u, Jc)
    call ref_coords%metric_tensor(u, g, ginv, sqrtg)
    if (.not. jacobian_ok(Jc)) then   ! near-axis singular chart
      status = CPB_LOCATE_FAIL; return
    end if
    ! bhat (wedge Cartesian) from B = curl A (B^s = 0 exactly, exact pitch); same
    ! construction as field_at_logical so seed/readout match the push.
    hctr = [0.0_dp, -dA(3,1)/sqrtg, dA(2,1)/sqrtg]
    do i = 1, 3
      bw(i) = Jc(i,1)*hctr(1) + Jc(i,2)*hctr(2) + Jc(i,3)*hctr(3)
    end do
    bw = bw/max(sqrt(bw(1)**2 + bw(2)**2 + bw(3)**2), 1.0e-30_dp)
    ! perpendicular gyrophase reference: contravariant flux direction -> Cartesian.
    call perp_unit_dir_flux(g, ginv, hcov, eperp_u)
    do i = 1, 3
      ew(i) = Jc(i,1)*eperp_u(1) + Jc(i,2)*eperp_u(2) + Jc(i,3)*eperp_u(3)
    end do
    ew = ew/max(sqrt(ew(1)**2 + ew(2)**2 + ew(3)**2), 1.0e-30_dp)
    bhat = rotz(bw, ca, sa)       ! wedge -> global frame
    eperp = rotz(ew, ca, sa)
    status = CPB_OK
  end subroutine locate

  ! Seed from a guiding-centre start record u0=(s, theta_B, phi_B) with parallel
  ! speed vpar0 and perpendicular speed vperp0. pauli=.true. (CPP) keeps v_perp=0
  ! on the Pauli slow manifold and carries mu|B|. pauli=.false. (full-orbit CP)
  ! places the particle a Larmor vector off the guiding centre in Cartesian (regular
  ! through the axis) and seeds v = vpar0 bhat + vperp0 e_perp with e_perp the same
  ! gyrophase reference the position offset uses.
  subroutine cpp_boris_init(st, pauli, x0_boozer, vpar0, vperp0, mu_in, mass, &
                            charge, dt, ro0_in, pabs, filtered)
    type(cpp_boris_state_t), intent(out) :: st
    logical, intent(in) :: pauli
    real(dp), intent(in) :: x0_boozer(3), vpar0, vperp0, mu_in, mass, charge, &
                            dt, ro0_in, pabs
    logical, intent(in), optional :: filtered
    real(dp) :: u_gc(3), xyz_gc(3), u_p(3), x_p(3), qc
    real(dp) :: bhat(3), eperp(3), Bmod
    integer :: status

    st%pauli = pauli
    st%mass = mass; st%charge = charge; st%dt = dt; st%ro0 = ro0_in
    st%mu = mu_in; st%pabs = pabs
    if (present(filtered)) st%filtered = filtered

    ! GC logical coords: chartmap radial label is rho = sqrt(s).
    u_gc = [sqrt(max(x0_boozer(1), 0.0_dp)), x0_boozer(2), x0_boozer(3)]
    call ref_coords%evaluate_cart(u_gc, xyz_gc)
    qc = charge/ro0_in

    x_p = xyz_gc
    u_p = u_gc
    if (.not. pauli .and. vperp0 > 0.0_dp) then
      ! Larmor offset off the guiding centre. If the offset point falls outside the
      ! chart (a near-edge marker whose gyro-circle pokes past s=1) or fails to
      ! locate, fall back to seeding at the guiding centre: the offset is O(rho_L),
      ! and a genuine edge orbit is then lost during integration, not at init. Never
      ! abort -- this runs per particle inside the OpenMP loop.
      call gc_to_particle(xyz_gc, u_gc, vperp0, mass, qc, x_p, u_p, status)
      if (status /= CPB_OK) then
        x_p = xyz_gc
        u_p = u_gc
      end if
    end if

    st%x = x_p
    st%u = u_p
    call locate(x_p, u_p, u_p, bhat, eperp, Bmod, status)
    if (status /= CPB_OK) then
      ! Cannot even seed the frame at the guiding centre: leave v=0 so the first
      ! orbit step reports a locate fault (counted confined), never a crash.
      st%v = 0.0_dp
      return
    end if
    st%u = u_p
    st%v = vpar0*bhat
    if (.not. pauli .and. vperp0 > 0.0_dp) st%v = st%v + vperp0*eperp
  end subroutine cpp_boris_init

  ! Cartesian guiding centre x_gc -> particle position a Larmor vector off it,
  ! solved by the fixed point x_p with cart(x_p) - rho(x_p) = x_gc, rho the Larmor
  ! vector built from the perpendicular speed at x_p (same gyrophase reference as
  ! the velocity seed), so the seed offset and the GC reconstruction are inverses.
  subroutine gc_to_particle(xyz_gc, u_gc, vperp0, mass, qc, x_p, u_p, status)
    real(dp), intent(in) :: xyz_gc(3), u_gc(3), vperp0, mass, qc
    real(dp), intent(out) :: x_p(3), u_p(3)
    integer, intent(out) :: status
    integer, parameter :: maxfp = 50
    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: bhat(3), eperp(3), Bmod
    real(dp) :: rho_l(3), xnew(3)
    integer :: it

    x_p = xyz_gc
    do it = 1, maxfp
      call locate(x_p, u_gc, u_p, bhat, eperp, Bmod, status)
      if (status /= CPB_OK) return
      ! rho = (m/(qc|B|)) bhat x v_perp, v_perp = vperp0 eperp (Cartesian).
      rho_l = (mass/(qc*Bmod))*cross(bhat, vperp0*eperp)
      xnew = xyz_gc + rho_l
      if (maxval(abs(xnew - x_p)) < tol) then
        x_p = xnew
        call locate(x_p, u_gc, u_p, bhat, eperp, Bmod, status)
        return
      end if
      x_p = xnew
    end do
    status = CPB_OK
  end subroutine gc_to_particle

  subroutine cpp_boris_step(st, status)
    type(cpp_boris_state_t), intent(inout) :: st
    integer, intent(out) :: status
    real(dp) :: x(3), v(3), Bvec(3), Bmod, gradB(3), u(3)
    real(dp) :: tvec(3), svec(3), vp(3), tmag2, qcm, fac

    x = st%x
    v = st%v
    qcm = st%charge/(c*st%ro0*st%mass)   ! rotation: dv/dt = qcm v x B

    x = x + 0.5_dp*st%dt*v
    call cart_field(x, st%u, Bvec, Bmod, gradB, u, status)
    if (status /= CPB_OK) return
    st%u = u

    ! half mirror kick (Pauli only): m dv = -mu grad|B|.
    if (st%pauli) v = v - 0.5_dp*st%dt*(st%mu/st%mass)*gradB

    ! exact magnetic rotation; optional HLW large-step filter on the angle.
    tvec = qcm*Bvec*0.5_dp*st%dt
    if (st%filtered) then
      tmag2 = sqrt(tvec(1)**2 + tvec(2)**2 + tvec(3)**2)
      if (tmag2 > 1.0e-30_dp) then
        fac = tan(tmag2)/tmag2          ! HLW: replace t by tan(theta/2)/(theta/2) t
        tvec = fac*tvec
      end if
    end if
    tmag2 = tvec(1)**2 + tvec(2)**2 + tvec(3)**2
    svec = 2.0_dp*tvec/(1.0_dp + tmag2)
    vp = v + cross(v, tvec)
    v = v + cross(vp, svec)

    if (st%pauli) v = v - 0.5_dp*st%dt*(st%mu/st%mass)*gradB

    x = x + 0.5_dp*st%dt*v

    st%x = x
    st%v = v
    status = CPB_OK
  end subroutine cpp_boris_step

  function cpp_boris_energy(st) result(energy)
    type(cpp_boris_state_t), intent(in) :: st
    real(dp) :: energy, Bvec(3), Bmod, gradB(3), u(3)
    integer :: status
    call cart_field(st%x, st%u, Bvec, Bmod, gradB, u, status)
    energy = 0.5_dp*st%mass*(st%v(1)**2 + st%v(2)**2 + st%v(3)**2)
    if (st%pauli .and. status == CPB_OK) energy = energy + st%mu*Bmod
  end function cpp_boris_energy

  ! Guiding-centre magnetic moment mu = m v_perp^2/(2|B_gc|) (#421): evaluate at
  ! the Larmor-corrected guiding centre, not the raw particle point, so the gyro
  ! ripple O(rho/L) is removed and mu is conserved to O((rho/L)^2). Diagnostic only.
  function cpp_boris_mu(st) result(mu)
    type(cpp_boris_state_t), intent(in) :: st
    real(dp) :: mu, s, th, ph, vpar, Bgc
    integer :: status
    call cpp_boris_to_gc(st, s, th, ph, vpar, status, Bmod_gc=Bgc)
    if (status /= CPB_OK) then
      mu = 0.0_dp; return
    end if
    mu = 0.5_dp*st%mass*max(st%v(1)**2 + st%v(2)**2 + st%v(3)**2 - vpar**2, &
                            0.0_dp)/max(Bgc, 1.0e-30_dp)
  end function cpp_boris_mu

  ! Guiding-centre reduction for output (#421): remove the Larmor vector
  ! (particle_to_gc, Cartesian) and report the centre in (s, theta_B, phi_B) with
  ! the parallel speed at the centre. status: CPB_OK / CPB_LOSS / CPB_LOCATE_FAIL.
  subroutine cpp_boris_to_gc(st, s, th, ph, vpar, status, Bmod_gc)
    type(cpp_boris_state_t), intent(in) :: st
    real(dp), intent(out) :: s, th, ph, vpar
    integer, intent(out) :: status
    real(dp), intent(out), optional :: Bmod_gc
    real(dp) :: u_p(3), x_gc(3), u_gc(3), qc
    real(dp) :: bhat(3), eperp(3), Bmod
    real(dp) :: vpar_p, vperp_cart(3), rho_l(3)

    s = 0.0_dp; th = 0.0_dp; ph = 0.0_dp; vpar = 0.0_dp
    if (present(Bmod_gc)) Bmod_gc = 0.0_dp

    call locate(st%x, st%u, u_p, bhat, eperp, Bmod, status)
    if (status /= CPB_OK) return
    qc = st%charge/st%ro0

    ! Larmor vector from the particle's perpendicular velocity at x (Cartesian):
    ! rho = (m/(qc|B|)) bhat x v_perp; x_gc = x - rho.
    vpar_p = st%v(1)*bhat(1) + st%v(2)*bhat(2) + st%v(3)*bhat(3)
    vperp_cart = st%v - vpar_p*bhat
    rho_l = (st%mass/(qc*Bmod))*cross(bhat, vperp_cart)
    x_gc = st%x - rho_l

    call locate(x_gc, u_p, u_gc, bhat, eperp, Bmod, status)
    if (status /= CPB_OK) return
    s = u_gc(1)**2               ! chart rho -> s
    th = u_gc(2); ph = u_gc(3)
    vpar = st%v(1)*bhat(1) + st%v(2)*bhat(2) + st%v(3)*bhat(3)
    if (present(Bmod_gc)) Bmod_gc = Bmod
    ! Confinement loss is decided here, on the Larmor-corrected guiding centre
    ! (#421): the particle may gyro-excurse past s=1 and return, as in ASCOT5; only
    ! the guiding centre crossing the last closed surface is a loss.
    if (u_gc(1) >= 1.0_dp) status = CPB_LOSS
  end subroutine cpp_boris_to_gc

  ! Unit perpendicular direction in contravariant flux components: raised radial
  ! covector projected off the field-parallel part, normalized in the metric.
  subroutine perp_unit_dir_flux(g, ginv, hcov, eperp)
    real(dp), intent(in) :: g(3,3), ginv(3,3), hcov(3)
    real(dp), intent(out) :: eperp(3)
    real(dp) :: er(3), hcon(3), hpar, nrm
    integer :: i, j
    er = [ginv(1,1), ginv(2,1), ginv(3,1)]
    do i = 1, 3
      hcon(i) = ginv(i,1)*hcov(1) + ginv(i,2)*hcov(2) + ginv(i,3)*hcov(3)
    end do
    hpar = hcov(1)*er(1) + hcov(2)*er(2) + hcov(3)*er(3)
    eperp = er - hpar*hcon
    nrm = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        nrm = nrm + g(i,j)*eperp(i)*eperp(j)
      end do
    end do
    if (nrm <= 0.0_dp) error stop 'perp_unit_dir_flux: degenerate direction'
    eperp = eperp/sqrt(nrm)
  end subroutine perp_unit_dir_flux

  ! A chart Jacobian is usable when its determinant is well above the round-off
  ! floor relative to its size (the chartmap is singular at the magnetic axis,
  ! rho->0). Rejecting near-singular Jc keeps the field push and the inversion off
  ! the axis singularity instead of producing Inf/NaN.
  pure logical function jacobian_ok(Jc)
    real(dp), intent(in) :: Jc(3,3)
    real(dp) :: det, scale
    det = Jc(1,1)*(Jc(2,2)*Jc(3,3) - Jc(2,3)*Jc(3,2)) &
        - Jc(1,2)*(Jc(2,1)*Jc(3,3) - Jc(2,3)*Jc(3,1)) &
        + Jc(1,3)*(Jc(2,1)*Jc(3,2) - Jc(2,2)*Jc(3,1))
    scale = sqrt(sum(Jc**2))
    jacobian_ok = (det == det) .and. abs(det) > 1.0e-8_dp*max(scale, 1.0e-30_dp)**3
  end function jacobian_ok

  pure function cross(a, b) result(cr)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: cr(3)
    cr(1) = a(2)*b(3) - a(3)*b(2)
    cr(2) = a(3)*b(1) - a(1)*b(3)
    cr(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  pure subroutine inv3(A, Ainv)
    real(dp), intent(in) :: A(3,3)
    real(dp), intent(out) :: Ainv(3,3)
    real(dp) :: det
    det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
        - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
        + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
    Ainv(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2))/det
    Ainv(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3))/det
    Ainv(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2))/det
    Ainv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))/det
    Ainv(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1))/det
    Ainv(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3))/det
    Ainv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))/det
    Ainv(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2))/det
    Ainv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/det
  end subroutine inv3

end module orbit_cpp_boris
