module orbit_fo_boris
  ! Full-orbit (FO) Boris pusher: a gyro-resolved classical charged particle in
  ! Cartesian (x, v), the full-orbit counterpart to SIMPLE's guiding-center (GC)
  ! model and the ASCOT-style reference for it. The particle advances by the
  ! explicit Boris drift-rotate-drift: the magnetic rotation is exact for constant
  ! B over a step, the kinetic metric is the identity, the geodesic terms vanish,
  ! and the magnetic axis is an ordinary point.
  !
  ! Field and geometry come from the chartmap (the Cartesian-side representation):
  ! at the Cartesian point we invert to the logical chart u=(rho, theta_B, phi_B)
  ! with the chartmap forward map (rho=sqrt(s)), evaluate the production Boozer
  ! flux potential there (fo_eval_field: A_theta(s), A_phi(s), |B|, d|B|/du), and
  ! build the Cartesian field as B = curl A so B^s = 0 exactly (B tangent to the
  ! flux surface) with grad|B|_cart = Jc^{-T} d|B|/du. The magnetic axis (rho->0)
  ! is healed by a pseudo-Cartesian chart w=(X,Y,phi)=(rho cos th, rho sin th, phi):
  ! the polar basis dx/dtheta ~ rho makes det(Jc)->0, but the (X,Y) chart Jacobian
  ! Jw is regular through rho=0, so the inverse Newton and the field assembly both
  ! stay well-conditioned there. The chartmap also owns the loss boundary: the
  ! guiding-centre crossing rho>=1 (out of the s<1 plasma) is the ONLY confinement
  ! loss. A field-locate non-convergence is a numerical fault, retried/reported,
  ! never a loss.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use reference_coordinates, only: ref_coords
  use orbit_fo_field, only: fo_eval_field
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

  ! A guiding-centre loss (u_gc >= 1) is confirmed only when the robustly-located
  ! particle radius u_p is within this gap of the edge: the GC and particle differ by
  ! one Larmor radius (rho*/a ~ 0.005-0.01 here), so 0.05 is several Larmor radii of
  ! margin and rejects field-period-seam reconstruction glitches that put a mid-radius
  ! particle's reconstructed GC spuriously at rho >= 1.
  real(dp), parameter :: GC_PARTICLE_GAP = 0.05_dp

  ! cart_field / locate status: regular interior point, physical edge loss, or a
  ! numerical locate fault (NOT a loss).
  integer, parameter, public :: FO_OK = 0, FO_LOSS = 1, FO_LOCATE_FAIL = 2

  public :: fo_state_t, fo_init, fo_step, fo_energy, fo_mu, fo_to_gc, accept_or_fail

  type :: fo_state_t
    real(dp) :: x(3)   = 0.0_dp    ! Cartesian position (scaled cm)
    real(dp) :: v(3)   = 0.0_dp    ! Cartesian velocity (normalized)
    real(dp) :: u(3)   = 0.0_dp    ! last logical (rho, theta_B, phi_B)
    real(dp) :: mu     = 0.0_dp    ! guiding-centre magnetic moment (diagnostic)
    real(dp) :: dt     = 0.0_dp
    real(dp) :: mass   = 1.0_dp
    real(dp) :: charge = 1.0_dp
    real(dp) :: ro0    = 1.0_dp
    real(dp) :: pabs   = 0.0_dp    ! normalized speed (carried for z(4) write-back)
  end type fo_state_t

contains

  ! Cartesian (wedge) -> logical chart (rho, theta_B, phi_B). Warm damped Newton on
  ! the chartmap forward map x(u)=evaluate_cart(u) from the carried guess (a Larmor
  ! step away, 1-2 iters, thread-safe read-only spline eval); on stall -- e.g. the
  ! guess went stale across a field-period seam -- fall back to the chartmap's robust
  ! multi-seed from_cart, which seeds zeta across [0, 2pi/nfp). status: FO_OK (located,
  ! rho reported through u for the caller's guiding-centre loss test) or FO_LOCATE_FAIL
  ! (no converged root -> numerical fault, never itself a loss).
  subroutine invert_cart_warm(x, u_guess, u, status)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: u(3)
    integer, intent(out) :: status
    real(dp) :: xc(3), Jc(3,3), rn
    integer :: ierr

    call invert_warm_newton(x, u_guess, u, status)
    if (status /= FO_LOCATE_FAIL) return
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
    status = accept_or_fail(u(1), rn, radial_scale(Jc), NEWTON_ACCEPT_TOL, RHO_EDGE, &
                            u_guess(1))
  end subroutine invert_cart_warm

  ! Cartesian (wedge) -> logical Newton. Seed rho and theta from the carried guess
  ! (they do not jump between substeps) and the toroidal coordinate from the in-wedge
  ! geometric angle atan2(y,x). The carried Boozer phi lives on the global
  ! multi-period sheet and goes a full period 2*pi/nfp stale across a field-period
  ! seam; the geometric angle is always in-wedge and differs from logical phi only by
  ! the Boozer shift O(0.1 rad). One robust seed, no seam special case. (A warm phi
  ! seed is faster away from seams but cannot be trusted at them: evaluate_cart wraps
  ! phi mod 2*pi/nfp, so a stale guess can converge to a clamped-edge root and fake a
  ! loss.) On stall the caller (invert_cart_warm) runs the multi-seed from_cart.
  subroutine invert_warm_newton(x, u_guess, u, status)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: u(3)
    integer, intent(out) :: status

    call newton_from(x, [u_guess(1), u_guess(2), atan2(x(2), x(1))], u, status)
  end subroutine invert_warm_newton

  ! Damped Newton on the chartmap forward map x(u)=evaluate_cart(u) from an explicit
  ! seed. Iterate in the pseudo-Cartesian chart w=(X,Y,phi)=(rho cos th, rho sin th,
  ! phi): the polar (rho,theta) Newton is singular at the axis (dx/dtheta ~ rho,
  ! det(Jc)->0), whereas the (X,Y) step stays regular and crosses the axis without
  ! the reflect hack. w_to_u recovers rho>=0, theta=atan2 automatically.
  subroutine newton_from(x, u_seed, u, status)
    real(dp), intent(in) :: x(3), u_seed(3)
    real(dp), intent(out) :: u(3)
    integer, intent(out) :: status
    integer, parameter :: maxit = 30, maxls = 30
    ! The forward map is a deterministic spline, so a damped Newton converges to
    ! ~machine precision; tol targets that. accept_or_fail classifies a stall.
    real(dp), parameter :: tol = 1.0e-9_dp
    real(dp) :: xc(3), Jc(3,3), Jw(3,3), Jinv(3,3), res(3)
    real(dp) :: w(3), wt(3), ut(3), dw(3), cth, sth, rho, rn, rnew, alpha
    integer :: it, ls, i

    u = u_seed
    w(1) = u(1)*cos(u(2)); w(2) = u(1)*sin(u(2)); w(3) = u(3)
    call ref_coords%evaluate_cart(u, xc)
    res = xc - x
    rn = sqrt(res(1)**2 + res(2)**2 + res(3)**2)
    do it = 1, maxit
      if (rn < tol) then
        status = accept_or_fail(u(1), rn, 0.0_dp, NEWTON_ACCEPT_TOL, RHO_EDGE, u_seed(1))
        return
      end if
      call ref_coords%covariant_basis(u, Jc)
      call pseudocart_basis(u, Jc, Jw, cth, sth, rho)
      if (.not. jacobian_ok(Jw)) then   ! genuinely degenerate (off the chart)
        status = FO_LOCATE_FAIL; return
      end if
      call inv3(Jw, Jinv)
      do i = 1, 3
        dw(i) = -(Jinv(i,1)*res(1) + Jinv(i,2)*res(2) + Jinv(i,3)*res(3))
      end do
      ! Backtracking line search: Newton is not monotonic for a finite offset.
      ! A trial that overshoots past rho=1 must be rejected, not evaluated: the
      ! forward map clamps rho to the grid edge, so a past-edge trial returns a
      ! point ON the edge whose residual can be smaller than the current interior
      ! point -- the line search would accept it and the next step stalls on the
      ! flat clamped region (the failure seen for mid-radius orbits crossing a
      ! field-period seam, where the seam-corrupted step points outward). w_to_u
      ! gives the trial rho before the clamp, so reject rho > 1 and keep backtracking
      ! toward the true interior root. Loss is decided only on the guiding centre.
      alpha = 1.0_dp
      do ls = 1, maxls
        wt = w + alpha*dw
        call w_to_u(wt, ut)
        if (ut(1) > 1.0_dp) then
          alpha = 0.5_dp*alpha
          cycle
        end if
        call ref_coords%evaluate_cart(ut, xc)
        res = xc - x
        rnew = sqrt(res(1)**2 + res(2)**2 + res(3)**2)
        if (rnew < rn) exit
        alpha = 0.5_dp*alpha
      end do
      if (rnew >= rn) then   ! line search could not improve -> stalled at the floor
        status = accept_or_fail(u(1), rn, radial_scale(Jc), NEWTON_ACCEPT_TOL, RHO_EDGE, &
                                u_seed(1))
        return
      end if
      w = wt
      u = ut
      rn = rnew
    end do
    status = accept_or_fail(u(1), rn, radial_scale(Jc), NEWTON_ACCEPT_TOL, RHO_EDGE, &
                            u_seed(1))
  end subroutine newton_from

  ! Length of one unit-rho radial step |dx/drho| = |Jc(:,1)|, the chart scale used to
  ! judge a stalled Newton: a residual that is a small fraction of a radial cell means
  ! the target is essentially at the edge, a large fraction means an interior stall.
  pure real(dp) function radial_scale(Jc) result(s)
    real(dp), intent(in) :: Jc(3,3)
    s = sqrt(Jc(1,1)**2 + Jc(2,1)**2 + Jc(3,1)**2)
  end function radial_scale

  ! Classify a finished Newton by the Cartesian residual rn, judged against the local
  ! radial cell |dx/drho| (scale), not an absolute length. A point is located when rn
  ! is at machine tolerance OR a small fraction (EDGE_FRAC) of a radial cell. The
  ! relative test is what makes this scale-correct: on a reactor-size chartmap
  ! (positions ~1e3 cm) and near the magnetic axis, where the chart is barely resolved
  ! (innermost rho grid point ~1e-3) and |dx/drho| is large, the residual floors well
  ! above any fixed accept_tol while the point still sits a tiny fraction of a cell
  ! from its target. A genuine stall -- a residual that is a sizable fraction of a
  ! radial cell -- is a fault, EXCEPT when it clamps to the edge (rho ~ rho_edge) for a
  ! marker whose warm guess rho_guess was already within GC_PARTICLE_GAP of the edge:
  ! that is a marker leaving the plasma, whose true position is past rho=1 where the
  ! forward map cannot represent it, so the inverse stalls on the clamped edge. Accept
  ! it as located at the edge so the push continues on the clamped-edge field and
  ! fo_to_gc decides the loss on the guiding centre -- never a blanket confined fault
  ! that would silently keep an exiting marker. A mid-radius seam glitch that clamps to
  ! rho=1 has rho_guess well inside, fails the guard, and stays a fault, so it is never
  ! turned into a spurious loss.
  pure integer function accept_or_fail(rho, rn, scale, accept_tol, rho_edge, rho_guess) &
                                       result(status)
    real(dp), intent(in) :: rho, rn, scale, accept_tol, rho_edge, rho_guess
    if (rn < accept_tol .or. (scale > 0.0_dp .and. rn < EDGE_FRAC*scale)) then
      status = FO_OK
    else if (rho >= rho_edge - GC_PARTICLE_GAP .and. &
             rho_guess >= rho_edge - GC_PARTICLE_GAP) then
      status = FO_OK
    else
      status = FO_LOCATE_FAIL
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
    call fo_eval_field(ue, Acov, dA, Bmod, dBmod, hcov)
    call ref_coords%covariant_basis(ue, Jc)
    ! Heal the axis: work in the pseudo-Cartesian chart w=(X,Y,phi). The polar
    ! basis column Jc(:,2)=dx/dtheta ~ rho makes det(Jc)->0 at the axis; the (X,Y)
    ! chart Jacobian Jw is regular there (det(Jw)=det(Jc)/rho, same sign).
    call pseudocart_basis(ue, Jc, Jw, cth, sth, rho)
    if (.not. jacobian_ok(Jw)) then   ! genuinely degenerate (off the chart)
      status = FO_LOCATE_FAIL; return
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
    status = FO_OK
  end subroutine field_at_logical

  ! Cartesian B vector, |B|, and grad|B| at Cartesian x, from the chartmap field.
  ! status: FO_OK (located, u_out valid) or FO_LOCATE_FAIL (numerical inversion
  ! fault). On fault Bvec etc. are undefined and the caller must not push. Loss is
  ! not decided here -- the field is defined through the clamped edge, and only the
  ! guiding-centre crossing rho>=1 in fo_to_gc is a confinement loss.
  subroutine cart_field(x, u_guess, Bvec, Bmod, gradB, u_out, status)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: Bvec(3), Bmod, gradB(3), u_out(3)
    integer, intent(out) :: status
    real(dp) :: xw(3), ca, sa, u(3), Jc(3,3), Bw(3), gw(3)

    call to_wedge(x, xw, ca, sa)
    call invert_cart_warm(xw, u_guess, u, status)
    if (status /= FO_OK) return
    u_out = u
    call field_at_logical(u, Bw, Bmod, gw, Jc, status)
    if (status /= FO_OK) return
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
    if (status /= FO_OK) return
    u_out = u
    ! Same axis-healed field assembly as the push (field_at_logical), so the seed
    ! frame and the orbit-step field are identical.
    call field_at_logical(u, Bw, Bmod, gw, Jc, status)
    if (status /= FO_OK) return
    bw_hat = Bw/max(sqrt(Bw(1)**2 + Bw(2)**2 + Bw(3)**2), 1.0e-30_dp)
    call perp_ref(bw_hat, ew)     ! arbitrary unit vector perpendicular to b
    bhat = rotz(bw_hat, ca, sa)   ! wedge -> global frame
    eperp = rotz(ew, ca, sa)
    status = FO_OK
  end subroutine locate

  ! Seed from a guiding-centre start record u0=(s, theta_B, phi_B) with parallel
  ! speed vpar0 and perpendicular speed vperp0. Place the particle a Larmor vector
  ! off the guiding centre in Cartesian (regular through the axis) and seed
  ! v = vpar0 bhat + vperp0 e_perp with e_perp the same gyrophase reference the
  ! position offset uses.
  subroutine fo_init(st, x0_boozer, vpar0, vperp0, mu_in, mass, &
                            charge, dt, ro0_in, pabs)
    type(fo_state_t), intent(out) :: st
    real(dp), intent(in) :: x0_boozer(3), vpar0, vperp0, mu_in, mass, charge, &
                            dt, ro0_in, pabs
    real(dp) :: u_gc(3), xyz_gc(3), u_p(3), x_p(3), qc
    real(dp) :: bhat(3), eperp(3), Bmod
    integer :: status

    st%mass = mass; st%charge = charge; st%dt = dt; st%ro0 = ro0_in
    st%mu = mu_in; st%pabs = pabs

    ! GC logical coords: chartmap radial label is rho = sqrt(s).
    u_gc = [sqrt(max(x0_boozer(1), 0.0_dp)), x0_boozer(2), x0_boozer(3)]
    call ref_coords%evaluate_cart(u_gc, xyz_gc)
    qc = charge/ro0_in

    x_p = xyz_gc
    u_p = u_gc
    if (vperp0 > 0.0_dp) then
      ! Larmor offset off the guiding centre. If the offset point falls outside the
      ! chart (a near-edge marker whose gyro-circle pokes past s=1) or fails to
      ! locate, fall back to seeding at the guiding centre: the offset is O(rho_L),
      ! and a genuine edge orbit is then lost during integration, not at init. Never
      ! abort -- this runs per particle inside the OpenMP loop.
      call gc_to_particle(xyz_gc, u_gc, vperp0, mass, qc, x_p, u_p, status)
      if (status /= FO_OK) then
        x_p = xyz_gc
        u_p = u_gc
      end if
    end if

    st%x = x_p
    st%u = u_p
    call locate(x_p, u_p, u_p, bhat, eperp, Bmod, status)
    if (status /= FO_OK) then
      ! Cannot even seed the frame at the guiding centre: leave v=0 so the first
      ! orbit step reports a locate fault (counted confined), never a crash.
      st%v = 0.0_dp
      return
    end if
    st%u = u_p
    st%v = vpar0*bhat
    if (vperp0 > 0.0_dp) st%v = st%v + vperp0*eperp
  end subroutine fo_init

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
      if (status /= FO_OK) return
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
    status = FO_OK
  end subroutine gc_to_particle

  subroutine fo_step(st, status)
    type(fo_state_t), intent(inout) :: st
    integer, intent(out) :: status
    real(dp) :: x(3), v(3), Bvec(3), Bmod, gradB(3), u(3)
    real(dp) :: tvec(3), svec(3), vp(3), tmag2, qcm

    x = st%x
    v = st%v
    qcm = st%charge/(c*st%ro0*st%mass)   ! rotation: dv/dt = qcm v x B

    x = x + 0.5_dp*st%dt*v
    call cart_field(x, st%u, Bvec, Bmod, gradB, u, status)
    if (status /= FO_OK) return    ! unresolved: leave st at the last resolved state
    st%u = u

    ! exact magnetic rotation (constant B over the step).
    tvec = qcm*Bvec*0.5_dp*st%dt
    tmag2 = tvec(1)**2 + tvec(2)**2 + tvec(3)**2
    svec = 2.0_dp*tvec/(1.0_dp + tmag2)
    vp = v + cross(v, tvec)
    v = v + cross(vp, svec)

    x = x + 0.5_dp*st%dt*v

    st%x = x
    st%v = v
    status = FO_OK
  end subroutine fo_step

  function fo_energy(st) result(energy)
    type(fo_state_t), intent(in) :: st
    real(dp) :: energy
    energy = 0.5_dp*st%mass*(st%v(1)**2 + st%v(2)**2 + st%v(3)**2)
  end function fo_energy

  ! Guiding-centre magnetic moment mu = m v_perp^2/(2|B_gc|): evaluate at the
  ! Larmor-corrected guiding centre, not the raw particle point, so the gyro
  ! ripple O(rho/L) is removed and mu is conserved to O((rho/L)^2). Diagnostic only.
  function fo_mu(st) result(mu)
    type(fo_state_t), intent(in) :: st
    real(dp) :: mu, s, th, ph, vpar, Bgc
    integer :: status
    call fo_to_gc(st, s, th, ph, vpar, status, Bmod_gc=Bgc)
    if (status /= FO_OK) then
      mu = 0.0_dp; return
    end if
    mu = 0.5_dp*st%mass*max(st%v(1)**2 + st%v(2)**2 + st%v(3)**2 - vpar**2, &
                            0.0_dp)/max(Bgc, 1.0e-30_dp)
  end function fo_mu

  ! Guiding-centre reduction for output: remove the Larmor vector in Cartesian and
  ! report the centre in (s, theta_B, phi_B) with the parallel speed at the centre.
  ! status: FO_OK / FO_LOSS / FO_LOCATE_FAIL.
  subroutine fo_to_gc(st, s, th, ph, vpar, status, Bmod_gc)
    type(fo_state_t), intent(in) :: st
    real(dp), intent(out) :: s, th, ph, vpar
    integer, intent(out) :: status
    real(dp), intent(out), optional :: Bmod_gc
    real(dp) :: u_p(3), x_gc(3), u_gc(3), qc
    real(dp) :: bhat(3), eperp(3), Bmod
    real(dp) :: vpar_p, vperp_cart(3), rho_l(3)

    s = 0.0_dp; th = 0.0_dp; ph = 0.0_dp; vpar = 0.0_dp
    if (present(Bmod_gc)) Bmod_gc = 0.0_dp

    call locate(st%x, st%u, u_p, bhat, eperp, Bmod, status)
    if (status /= FO_OK) return
    qc = st%charge/st%ro0

    ! Larmor vector from the particle's perpendicular velocity at x (Cartesian):
    ! rho = (m/(qc|B|)) bhat x v_perp; x_gc = x - rho.
    vpar_p = st%v(1)*bhat(1) + st%v(2)*bhat(2) + st%v(3)*bhat(3)
    vperp_cart = st%v - vpar_p*bhat
    rho_l = (st%mass/(qc*Bmod))*cross(bhat, vperp_cart)
    x_gc = st%x - rho_l

    call locate(x_gc, u_p, u_gc, bhat, eperp, Bmod, status)
    if (status /= FO_OK) return
    s = u_gc(1)**2               ! chart rho -> s
    th = u_gc(2); ph = u_gc(3)
    vpar = st%v(1)*bhat(1) + st%v(2)*bhat(2) + st%v(3)*bhat(3)
    if (present(Bmod_gc)) Bmod_gc = Bmod
    ! Confinement loss: the Larmor-corrected guiding centre crosses the last closed
    ! surface (u_gc >= 1). The GC must be locatable, so the loss keys on u_gc rather
    ! than the particle (a particle gyro-excursed past s=1 is off-chart and cannot be
    ! inverted). But u_gc is a second, cold-guess locate of x_gc and carries the
    ! residual field-period-seam noise, which occasionally returns rho >= 1 for a
    ! particle that is in fact at mid-radius. Reject that with the robust warm-started
    ! particle locate u_p: a real loss has the particle within ~a Larmor radius of the
    ! edge (|x_gc - x| = |rho_l| is a Larmor radius), so u_gc >= 1 while u_p is well
    ! inside is a reconstruction glitch, not a loss. The field, integrator and energy
    ! match the ASCOT full-orbit reference, so the loss detector must be this clean.
    if (u_gc(1) >= 1.0_dp .and. u_p(1) >= 1.0_dp - GC_PARTICLE_GAP) status = FO_LOSS
  end subroutine fo_to_gc

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

end module orbit_fo_boris
