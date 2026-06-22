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
  ! u=(rho, theta_B, phi_B) with the chartmap forward map (rho=sqrt(s)), evaluate
  ! the production Boozer flux potential there (chartmap_eval_field: A_theta(s),
  ! A_phi(s), |B|, d|B|/du), and build the Cartesian field as B = curl A so B^s = 0
  ! exactly (B tangent to the flux surface) with grad|B|_cart = Jc^{-T} d|B|/du.
  ! The magnetic axis (rho->0) is healed by a pseudo-Cartesian chart
  ! w=(X,Y,phi)=(rho cos th, rho sin th, phi): the polar basis dx/dtheta ~ rho
  ! makes det(Jc)->0, but the (X,Y) chart Jacobian Jw is regular through rho=0, so
  ! the inverse Newton and the field assembly both stay well-conditioned there.
  ! The chartmap also owns the loss boundary: rho>=1 (out of the s<1 plasma) is
  ! the ONLY confinement loss. A field-locate non-convergence is a numerical fault,
  ! retried/reported, never a loss (#419, #420).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use reference_coordinates, only: ref_coords
  use orbit_cpp_chartmap_metric, only: chartmap_eval_field
  use libneo_coordinates, only: chartmap_coordinate_system_t, chartmap_from_cyl_ok
  implicit none
  private

  real(dp), parameter :: c = 1.0_dp
  real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

  ! Inverse-Newton classification: a converged locate has residual below
  ! NEWTON_ACCEPT_TOL; a loosely converged point counts as the edge only when its
  ! residual is below EDGE_FRAC of a radial cell (a genuine gyro-overshoot loss sits
  ! within a Larmor radius of RHO_EDGE), otherwise it is a stalled interior fault.
  real(dp), parameter :: NEWTON_ACCEPT_TOL = 1.0e-6_dp, RHO_EDGE = 1.0_dp, &
                         EDGE_FRAC = 0.05_dp

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

  ! Cartesian (wedge) -> logical chart (rho, theta_B, phi_B). Warm damped Newton on
  ! the chartmap forward map x(u)=evaluate_cart(u) from the carried guess (a Larmor
  ! step away, 1-2 iters, thread-safe read-only spline eval); on stall -- e.g. the
  ! guess went stale across a field-period seam -- fall back to the chartmap's robust
  ! multi-seed from_cart, which seeds zeta across [0, 2pi/nfp). status: CPB_OK (located,
  ! rho reported through u for the caller's guiding-centre loss test) or CPB_LOCATE_FAIL
  ! (no converged root -> numerical fault, never itself a loss).
  subroutine invert_cart_warm(x, u_guess, u, status)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: u(3)
    integer, intent(out) :: status
    real(dp) :: xc(3), Jc(3,3), rn
    integer :: ierr

    call invert_warm_newton(x, u_guess, u, status)
    if (status /= CPB_LOCATE_FAIL) return
    select type (cs => ref_coords)               ! robust multi-seed fallback
    class is (chartmap_coordinate_system_t)
      call cs%from_cart(x, u, ierr)
    class default
      return
    end select
    if (ierr /= chartmap_from_cyl_ok) return     ! genuine no-root: keep LOCATE_FAIL
    ! Re-verify the fallback root with the same residual-vs-radial-cell criterion as
    ! the warm path: from_cart clamps rho to [0,1], so a seam point it cannot solve
    ! comes back pinned at rho=1 with a large residual. Accepting that as the edge
    ! fakes a loss from mid-radius, so classify by the actual residual, not by rho.
    call ref_coords%evaluate_cart(u, xc)
    call ref_coords%covariant_basis(u, Jc)
    rn = sqrt((xc(1) - x(1))**2 + (xc(2) - x(2))**2 + (xc(3) - x(3))**2)
    status = accept_or_fail(u(1), rn, radial_scale(Jc), NEWTON_ACCEPT_TOL, RHO_EDGE)
  end subroutine invert_cart_warm

  subroutine invert_warm_newton(x, u_guess, u, status)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: u(3)
    integer, intent(out) :: status
    integer, parameter :: maxit = 30, maxls = 30
    ! The forward map is a deterministic spline, so a damped Newton on the wedge
    ! point converges to ~machine precision; tol targets that. accept_tol only
    ! classifies a Newton that has stalled at the spline floor.
    real(dp), parameter :: tol = 1.0e-9_dp
    real(dp) :: xc(3), Jc(3,3), Jw(3,3), Jinv(3,3), res(3)
    real(dp) :: w(3), wt(3), ut(3), dw(3), cth, sth, rho, rn, rnew, alpha
    integer :: it, ls, i

    ! Iterate in the pseudo-Cartesian chart w=(X,Y,phi)=(rho cos th, rho sin th,
    ! phi): the polar (rho,theta) Newton is singular at the axis (dx/dtheta ~ rho,
    ! det(Jc)->0), whereas the (X,Y) step stays regular and crosses the axis
    ! without the reflect hack. w_to_u recovers rho>=0, theta=atan2 automatically.
    ! Warm the radial/poloidal guess (rho, theta do not jump between substeps), but
    ! seed the toroidal coordinate from the in-wedge geometric angle atan2(y,x). The
    ! carried u_guess(3) is the Boozer phi on the global multi-period sheet; across a
    ! field-period seam to_wedge rotates x into the next wedge while u_guess(3) still
    ! sits a full period 2*pi/nfp away, stalling the Newton and tripping a spurious
    ! loss. The wedge geometric angle differs from logical phi only by the Boozer
    ! shift O(0.1 rad), so Newton converges in 1-2 iters at and away from the seam.
    u = u_guess
    u(3) = atan2(x(2), x(1))
    w(1) = u(1)*cos(u(2)); w(2) = u(1)*sin(u(2)); w(3) = u(3)
    call ref_coords%evaluate_cart(u, xc)
    res = xc - x
    rn = sqrt(res(1)**2 + res(2)**2 + res(3)**2)
    do it = 1, maxit
      if (rn < tol) then
        status = accept_or_fail(u(1), rn, 0.0_dp, NEWTON_ACCEPT_TOL, RHO_EDGE)
        return
      end if
      call ref_coords%covariant_basis(u, Jc)
      call pseudocart_basis(u, Jc, Jw, cth, sth, rho)
      if (.not. jacobian_ok(Jw)) then   ! genuinely degenerate (off the chart)
        status = CPB_LOCATE_FAIL; return
      end if
      call inv3(Jw, Jinv)
      do i = 1, 3
        dw(i) = -(Jinv(i,1)*res(1) + Jinv(i,2)*res(2) + Jinv(i,3)*res(3))
      end do
      ! backtracking line search: Newton is not monotonic for a finite offset.
      ! A trial overshoot past rho=1 is NOT a loss: evaluate_cart clamps rho to the
      ! grid edge so an interior target yields a large residual and the step
      ! backtracks. Loss is decided only on the converged rho (accept_or_fail).
      alpha = 1.0_dp
      do ls = 1, maxls
        wt = w + alpha*dw
        call w_to_u(wt, ut)
        call ref_coords%evaluate_cart(ut, xc)
        res = xc - x
        rnew = sqrt(res(1)**2 + res(2)**2 + res(3)**2)
        if (rnew < rn) exit
        alpha = 0.5_dp*alpha
      end do
      if (rnew >= rn) then   ! line search could not improve -> stalled at the floor
        status = accept_or_fail(u(1), rn, radial_scale(Jc), NEWTON_ACCEPT_TOL, RHO_EDGE)
        return
      end if
      w = wt
      u = ut
      rn = rnew
    end do
    status = accept_or_fail(u(1), rn, radial_scale(Jc), NEWTON_ACCEPT_TOL, RHO_EDGE)
  end subroutine invert_warm_newton

  ! Length of one unit-rho radial step |dx/drho| = |Jc(:,1)|, the chart scale used to
  ! judge a stalled Newton: a residual that is a small fraction of a radial cell means
  ! the target is essentially at the edge, a large fraction means an interior stall.
  pure real(dp) function radial_scale(Jc) result(s)
    real(dp), intent(in) :: Jc(3,3)
    s = sqrt(Jc(1,1)**2 + Jc(2,1)**2 + Jc(3,1)**2)
  end function radial_scale

  ! Classify a finished Newton. A converged locate (rn below accept_tol) is CPB_OK and
  ! reports the radius through u, so the caller decides loss on the guiding-centre
  ! radius. A loosely converged point AT the clamped edge is CPB_OK only when the
  ! residual is a small fraction of a radial cell (a genuine gyro-overshoot loss sits
  ! within a Larmor radius of rho=1). A Newton that stalls at the clamped edge while
  ! its true radius is well inside the plasma has a residual of order a radial cell:
  ! that is a numerical fault (counted confined), NOT a loss -- otherwise an inversion
  ! that clamps to rho=1 at a field-period seam fakes an edge loss from mid-radius.
  pure integer function accept_or_fail(rho, rn, scale, accept_tol, rho_edge) result(status)
    real(dp), intent(in) :: rho, rn, scale, accept_tol, rho_edge
    if (rn < accept_tol) then
      status = CPB_OK
    else if (rho >= rho_edge - 1.0e-3_dp .and. rn < EDGE_FRAC*scale) then
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
    real(dp) :: Jw(3,3), Jwinv(3,3), cth, sth, rho, detJw, Bw(3), dBw(3), Bn
    integer :: i

    ! A particle may gyro-excurse a Larmor radius past s=1; evaluate the field at
    ! the clamped edge there (field_can is undefined past the last closed surface).
    ue = u
    ue(1) = min(ue(1), 1.0_dp - 1.0e-9_dp)
    call chartmap_eval_field(ue, Acov, dA, Bmod, dBmod, hcov)
    call ref_coords%covariant_basis(ue, Jc)
    ! Heal the axis: work in the pseudo-Cartesian chart w=(X,Y,phi). The polar
    ! basis column Jc(:,2)=dx/dtheta ~ rho makes det(Jc)->0 at the axis; the (X,Y)
    ! chart Jacobian Jw is regular there (det(Jw)=det(Jc)/rho, same sign).
    call pseudocart_basis(ue, Jc, Jw, cth, sth, rho)
    if (.not. jacobian_ok(Jw)) then   ! genuinely degenerate (off the chart)
      status = CPB_LOCATE_FAIL; return
    end if
    detJw = Jw(1,1)*(Jw(2,2)*Jw(3,3) - Jw(2,3)*Jw(3,2)) &
          - Jw(1,2)*(Jw(2,1)*Jw(3,3) - Jw(2,3)*Jw(3,1)) &
          + Jw(1,3)*(Jw(2,1)*Jw(3,2) - Jw(2,2)*Jw(3,1))
    ! B = curl A with the Boozer flux-function potential A = (0, A_theta(s),
    ! A_phi(s)): B^s = d_theta A_phi - d_phi A_theta = 0 EXACTLY (B tangent to the
    ! flux surface). In the regular (X,Y,phi) chart the contravariant components
    ! stay finite through the axis (det(Jw) is bounded, and dA_theta/drho ~ rho
    ! cancels the surviving 1/rho). The signed det(Jw) keeps the left-handed chart
    ! orientation: the unsigned sqrt(|det|) flips B and sends trapped bananas
    ! outward. B^X = sin th dA_phi/drho /detJw, B^Y = -cos th dA_phi/drho /detJw,
    ! B^phi = (dA_theta/drho)/rho /detJw. Renormalize to |B|.
    Bw(1) =  sth*dA(3,1)/detJw
    Bw(2) = -cth*dA(3,1)/detJw
    Bw(3) = (dA(2,1)/rho)/detJw
    do i = 1, 3
      Bvec(i) = Jw(i,1)*Bw(1) + Jw(i,2)*Bw(2) + Jw(i,3)*Bw(3)
    end do
    Bn = sqrt(Bvec(1)**2 + Bvec(2)**2 + Bvec(3)**2)
    Bvec = Bvec*(Bmod/max(Bn, 1.0e-30_dp))
    ! grad|B| in Cartesian: Jw^{-T} d|B|/dw, with d|B|/dw mapped from d|B|/du by the
    ! pseudo-Cartesian chain rule (d|B|/dtheta ~ rho cancels the 1/rho).
    dBw(1) = cth*dBmod(1) - (sth/rho)*dBmod(2)
    dBw(2) = sth*dBmod(1) + (cth/rho)*dBmod(2)
    dBw(3) = dBmod(3)
    call inv3(Jw, Jwinv)
    do i = 1, 3
      gradB(i) = Jwinv(1,i)*dBw(1) + Jwinv(2,i)*dBw(2) + Jwinv(3,i)*dBw(3)
    end do
    status = CPB_OK
  end subroutine field_at_logical

  ! Cartesian B vector, |B|, and grad|B| at Cartesian x, from the chartmap field.
  ! status: CPB_OK (located, u_out valid) or CPB_LOCATE_FAIL (numerical inversion
  ! fault). On fault Bvec etc. are undefined and the caller must not push. Loss is
  ! not decided here -- the field is defined through the clamped edge, and only the
  ! guiding-centre crossing rho>=1 in cpp_boris_to_gc is a confinement loss.
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
    real(dp) :: xw(3), ca, sa, u(3), Jc(3,3), Bw(3), gw(3), bw_hat(3), ew(3)

    call to_wedge(x, xw, ca, sa)
    call invert_cart_warm(xw, u_guess, u, status)
    if (status /= CPB_OK) return
    u_out = u
    ! Same axis-healed field assembly as the push (field_at_logical), so the seed
    ! frame and the orbit-step field are identical.
    call field_at_logical(u, Bw, Bmod, gw, Jc, status)
    if (status /= CPB_OK) return
    bw_hat = Bw/max(sqrt(Bw(1)**2 + Bw(2)**2 + Bw(3)**2), 1.0e-30_dp)
    call perp_ref(bw_hat, ew)     ! arbitrary unit vector perpendicular to b
    bhat = rotz(bw_hat, ca, sa)   ! wedge -> global frame
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

  ! Pseudo-Cartesian near-axis chart w=(X,Y,phi)=(rho cos th, rho sin th, phi).
  ! The chartmap polar chart (rho,theta) is singular at the magnetic axis: the
  ! covariant basis column Jc(:,2)=dx/dtheta ~ rho vanishes, so det(Jc)->0 and both
  ! the inverse Newton (ill-conditioned in theta) and the field assembly degrade.
  ! The (X,Y) basis stays regular through rho=0 (Pfefferle et al.,
  ! arXiv:1412.5464; libneo flux_pseudocartesian). Returns the regular chart
  ! Jacobian Jw(a,i)=dx_a/dw_i and the trig used to map field components.
  subroutine pseudocart_basis(u, Jc, Jw, cth, sth, rho)
    real(dp), intent(in) :: u(3), Jc(3,3)
    real(dp), intent(out) :: Jw(3,3), cth, sth, rho
    integer :: a
    rho = max(u(1), 1.0e-30_dp)
    cth = cos(u(2)); sth = sin(u(2))
    do a = 1, 3
      Jw(a,1) = Jc(a,1)*cth - Jc(a,2)*(sth/rho)   ! e_X = dx/dX
      Jw(a,2) = Jc(a,1)*sth + Jc(a,2)*(cth/rho)   ! e_Y = dx/dY
      Jw(a,3) = Jc(a,3)                           ! e_phi
    end do
  end subroutine pseudocart_basis

  ! Pseudo-Cartesian w=(X,Y,phi) -> logical u=(rho,theta,phi). rho>=0 and the
  ! atan2 branch make the axis an ordinary point (no reflect hack on the inverse).
  pure subroutine w_to_u(w, u)
    real(dp), intent(in) :: w(3)
    real(dp), intent(out) :: u(3)
    u(1) = sqrt(w(1)**2 + w(2)**2)
    u(2) = atan2(w(2), w(1))
    u(3) = w(3)
  end subroutine w_to_u

  ! An arbitrary unit vector perpendicular to b, regular everywhere. The gyrophase
  ! reference is gauge: only b and |v_perp| are physical, and the guiding-centre
  ! reduction recovers v_perp from v directly, not from this choice. Gram-Schmidt
  ! off the least-aligned axis so the subtraction never cancels.
  pure subroutine perp_ref(b, e)
    real(dp), intent(in) :: b(3)
    real(dp), intent(out) :: e(3)
    real(dp) :: a(3), d, n
    if (abs(b(3)) < 0.9_dp) then
      a = [0.0_dp, 0.0_dp, 1.0_dp]
    else
      a = [1.0_dp, 0.0_dp, 0.0_dp]
    end if
    d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
    e = a - d*b
    n = sqrt(e(1)**2 + e(2)**2 + e(3)**2)
    e = e/max(n, 1.0e-30_dp)
  end subroutine perp_ref

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
