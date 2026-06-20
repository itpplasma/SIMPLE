module orbit_cp_explicit
  ! Classical charged particle (CP) pusher in NATIVE VMEC flux coordinates
  ! u = (s, vartheta, varphi), gyro-resolved, for the production CP6D loss path
  ! (FACTS design item 2). It integrates the curvilinear Lorentz ODE in canonical
  ! Hamiltonian form WITHOUT any Newton iteration or Jacobian: the earlier implicit
  ! canonical-midpoint CP6D path (orbit_cpp_canonical MODEL_CP + COORD_VMEC) used a
  ! finite-difference Newton Jacobian that goes noisy at banana turning points
  ! (v_par -> 0) and spuriously ejects all trapped particles. Here there is no
  ! Jacobian, so a turning point is a smooth point of the iteration.
  !
  ! SYMPLECTIC IMPLICIT MIDPOINT, fixed-point (Picard) iterated -- not RK4 and not
  ! a Newton solve. A non-symplectic RK4 heats the orbit (the magnetic force does
  ! no work, but RK4 leaks energy secularly) and over a 10 ms gyro-resolved trace
  ! -- O(1e6) gyrations -- the accumulated heating widens the banana until the
  ! orbit drifts to the axis or edge and is spuriously lost. The implicit midpoint
  ! is symplectic, so energy stays bounded with no secular drift over the whole
  ! trace; and because the gyro-resolved step makes dt*Omega < 1, the midpoint
  ! fixed-point map is a contraction and Picard converges in a few cheap iterations
  ! -- no Jacobian, no LU, robust through v_par -> 0.
  !
  ! State and Hamiltonian. z = (x^k, p_k), p the canonical covariant momentum:
  !   p_k     = m g_kj v^j + qc A_k
  !   v^k     = g^kj (p_j - qc A_j) / m            (contravariant velocity)
  !   H       = (1/2m)(p - qc A) g^ij (p - qc A)
  ! Hamilton's equations (the curvilinear Lorentz equation in canonical form):
  !   dx^k/dt = v^k
  !   dp_k/dt = (m/2) g_ij,k v^i v^j + qc A_i,k v^i
  ! The geodesic (Christoffel) part is folded into g_ij,k v^i v^j and the metric in
  ! v^k; the magnetic part into qc A_i,k v^i. No magnetic moment mu enters the EOM
  ! (the full particle resolves the gyration); mu seeds vperp from the GC pitch.
  !
  ! Normalization: the SIMPLE GC sqrt(2) convention, identical to init_sympl /
  ! init_cpp. mass = 1, qc = charge/(c ro0) = 1/ro0_bar with ro0_bar = ro0/sqrt(2),
  ! step dt = dtaumin/sqrt(2). mass = 1 keeps v ~ O(vpar_bar) ~ O(1).
  !
  ! GPU portability: vmec_field_metric_eval is !$acc routine seq (single libneo
  ! spline evaluation, no class() dispatch), and the midpoint/Picard arithmetic is
  ! pure fixed-size, so cp_explicit_step is device-callable -- one particle per
  ! thread.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use vmec_field_metric, only: vmec_field_metric_eval
  implicit none
  private

  ! Thesis normalization speed of light c = 1 (the same c orbit_cpp_canonical
  ! uses, NOT the physical CGS c that would make the magnetic coupling vanish);
  ! qc = charge/(c ro0) reads as 1/ro0 with charge = 1.
  real(dp), parameter :: c = 1.0_dp

  public :: cp_explicit_state_t, cp_explicit_init, cp_explicit_step, &
            cp_explicit_to_gc, cp_explicit_energy

  type :: cp_explicit_state_t
    real(dp) :: x(3)   = 0.0_dp   ! (s, vartheta, varphi)
    real(dp) :: p(3)   = 0.0_dp   ! covariant canonical momentum p_i
    real(dp) :: mass   = 1.0_dp
    real(dp) :: ro0    = 1.0_dp   ! ro0_bar = ro0/sqrt(2); qc = 1/(c ro0)
    real(dp) :: dt     = 0.0_dp
    real(dp) :: mu     = 0.0_dp   ! GC magnetic moment (seed + energy diag only)
    real(dp) :: pabs   = 0.0_dp   ! normalized speed (GC z(4)), carried for write-back
  end type cp_explicit_state_t

contains

  ! Seed the explicit CP state from the SAME GC start as init_sympl: position x0,
  ! parallel speed vpar0 (vpar_bar), magnetic moment mu_in. The perpendicular seed
  ! direction is the metric-unit radial direction projected perpendicular to the
  ! field at a fixed gyrophase; vperp = sqrt(2 mu |B|) from the GC pitch. The
  ! resulting orbit gyrates about a center within O(rho*) of the GC start (the FLR
  ! offset is the physics). p_i = m g_ij v^j + qc A_i.
  !$acc routine seq
  subroutine cp_explicit_init(st, x0, vpar0, mu_in, mass, ro0_in, dt)
    type(cp_explicit_state_t), intent(out) :: st
    real(dp), intent(in) :: x0(3), vpar0, mu_in, mass, ro0_in, dt

    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: hcon(3), eperp(3), vcon(3), qc, vperp, hpar, nrm
    integer :: i, j

    st%mass = mass
    st%ro0  = ro0_in
    st%dt   = dt
    st%mu   = mu_in
    st%x    = x0

    call vmec_field_metric_eval(x0, g, ginv, sqrtg, dg, Acov, dA, &
                                Bctr, Bcov, Bmod, dBmod, hcov)
    qc = 1.0_dp/(c*st%ro0)

    ! Parallel velocity v_par^i = vpar0 g^ij h_j.
    call raise(ginv, hcov, hcon)
    do i = 1, 3
      vcon(i) = vpar0*hcon(i)
    end do

    ! Perpendicular seed: radial covariant direction e_r = (1,0,0) raised to
    ! e_r^i = g^i1, projected perpendicular to h (|h|^2 = 1), normalized in g.
    eperp = [ginv(1,1), ginv(2,1), ginv(3,1)]
    hpar = hcov(1)*eperp(1) + hcov(2)*eperp(2) + hcov(3)*eperp(3)
    do i = 1, 3
      eperp(i) = eperp(i) - hpar*hcon(i)
    end do
    nrm = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        nrm = nrm + g(i,j)*eperp(i)*eperp(j)
      end do
    end do
    if (nrm > 0.0_dp) then
      eperp = eperp/sqrt(nrm)
    else
      eperp = [1.0_dp, 0.0_dp, 0.0_dp]
    end if

    vperp = sqrt(2.0_dp*st%mu*Bmod/mass)
    do i = 1, 3
      vcon(i) = vcon(i) + vperp*eperp(i)
    end do

    ! p_i = m g_ij v^j + qc A_i.
    do i = 1, 3
      st%p(i) = qc*Acov(i)
      do j = 1, 3
        st%p(i) = st%p(i) + mass*g(i,j)*vcon(j)
      end do
    end do
  end subroutine cp_explicit_init

  ! One symplectic implicit-midpoint step, Picard (fixed-point) iterated. The
  ! midpoint map z_{n+1} = z_n + dt f((z_n + z_{n+1})/2) is solved by iterating
  ! z^{m+1} = z_n + dt f((z_n + z^m)/2); for the gyro-resolved step (dt*Omega < 1)
  ! this is a contraction and converges in a handful of iterations. No Jacobian, so
  ! v_par -> 0 is a smooth point. Boundary guard keeps s in (0,1): ierr = 2 on
  ! s >= 1 (loss), ierr = 0 otherwise.
  !$acc routine seq
  subroutine cp_explicit_step(st, ierr)
    type(cp_explicit_state_t), intent(inout) :: st
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 200
    real(dp), parameter :: tol = 1.0e-12_dp
    real(dp) :: zold(6), z(6), znew(6), zmid(6), f(6), smid
    real(dp) :: scale(6), relchg, pscale
    integer :: kit, i

    ierr = 0
    zold(1:3) = st%x
    zold(4:6) = st%p
    z = zold

    ! Per-component convergence scale. The state mixes angles (s ~ 1, theta/phi can
    ! be hundreds of radians) with covariant momenta of a different magnitude, so a
    ! single mixed norm would declare convergence while the small radial momentum
    ! p_1 -- the component that drives the s-drift -- is still moving. A loosely
    ! converged midpoint is no longer symplectic and leaks a spurious radial drift.
    ! Scale s and the two angles by 1, the momenta by their own magnitude.
    pscale = max(abs(zold(4)) + abs(zold(5)) + abs(zold(6)), 1.0e-30_dp)
    scale(1:3) = 1.0_dp
    scale(4:6) = pscale

    do kit = 1, maxit
      zmid = 0.5_dp*(zold + z)
      smid = zmid(1)
      ! Keep the midpoint evaluation inside the domain; an s>=1 midpoint means the
      ! orbit has reached the edge -- flag the loss.
      if (smid >= 1.0_dp) then
        ierr = 2; return
      end if
      if (smid <= 0.0_dp) smid = 1.0e-8_dp
      zmid(1) = smid
      call cp_rhs(st, zmid, f)
      znew = zold + st%dt*f
      relchg = 0.0_dp
      do i = 1, 6
        relchg = max(relchg, abs(znew(i) - z(i))/scale(i))
      end do
      z = znew
      if (relchg <= tol) exit
    end do

    if (z(1) >= 1.0_dp) then
      ierr = 2; return
    end if
    if (z(1) <= 0.0_dp) z(1) = 1.0e-8_dp

    st%x = z(1:3)
    st%p = z(4:6)
  end subroutine cp_explicit_step

  ! RHS of the (x, p) Hamiltonian ODE: dx^k/dt = v^k, dp_k/dt = (m/2) g_ij,k v^i
  ! v^j + qc A_i,k v^i, with v^k = g^kj (p_j - qc A_j)/m.
  !$acc routine seq
  subroutine cp_rhs(st, z, dzdt)
    type(cp_explicit_state_t), intent(in) :: st
    real(dp), intent(in) :: z(6)
    real(dp), intent(out) :: dzdt(6)
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: xx(3), pcov(3), vcov(3), vcon(3), qc, geo, em
    integer :: i, j, k

    xx = z(1:3)
    pcov = z(4:6)
    call vmec_field_metric_eval(xx, g, ginv, sqrtg, dg, Acov, dA, &
                                Bctr, Bcov, Bmod, dBmod, hcov)
    qc = 1.0_dp/(c*st%ro0)

    do k = 1, 3
      vcov(k) = pcov(k) - qc*Acov(k)
    end do
    call raise(ginv, vcov, vcon)
    do k = 1, 3
      vcon(k) = vcon(k)/st%mass
    end do

    do k = 1, 3
      dzdt(k) = vcon(k)
    end do
    do k = 1, 3
      geo = 0.0_dp
      do j = 1, 3
        do i = 1, 3
          geo = geo + dg(i,j,k)*vcon(i)*vcon(j)
        end do
      end do
      em = 0.0_dp
      do i = 1, 3
        em = em + dA(i,k)*vcon(i)
      end do
      dzdt(3+k) = 0.5_dp*st%mass*geo + qc*em
    end do
  end subroutine cp_rhs

  ! Guiding-center reduction at the current state: position is x; vpar = h_i v^i.
  !$acc routine seq
  subroutine cp_explicit_to_gc(st, s, th, ph, vpar)
    type(cp_explicit_state_t), intent(in) :: st
    real(dp), intent(out) :: s, th, ph, vpar
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: vcov(3), vcon(3), qc
    integer :: k

    call vmec_field_metric_eval(st%x, g, ginv, sqrtg, dg, Acov, dA, &
                                Bctr, Bcov, Bmod, dBmod, hcov)
    qc = 1.0_dp/(c*st%ro0)
    do k = 1, 3
      vcov(k) = (st%p(k) - qc*Acov(k))/st%mass
    end do
    call raise(ginv, vcov, vcon)
    s = st%x(1); th = st%x(2); ph = st%x(3)
    vpar = hcov(1)*vcon(1) + hcov(2)*vcon(2) + hcov(3)*vcon(3)
  end subroutine cp_explicit_to_gc

  ! Kinetic energy H = (1/2m)(p - qc A) g^ij (p - qc A). The full charged particle
  ! resolves the gyration, so its kinetic term already carries the perpendicular
  ! energy: no separate mu|B| term (matches cpp_canon_energy for MODEL_CP). The
  ! symplectic midpoint keeps it bounded with no secular drift.
  !$acc routine seq
  function cp_explicit_energy(st) result(energy)
    type(cp_explicit_state_t), intent(in) :: st
    real(dp) :: energy
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: vcov(3), vcon(3), qc
    integer :: k

    call vmec_field_metric_eval(st%x, g, ginv, sqrtg, dg, Acov, dA, &
                                Bctr, Bcov, Bmod, dBmod, hcov)
    qc = 1.0_dp/(c*st%ro0)
    do k = 1, 3
      vcov(k) = st%p(k) - qc*Acov(k)
    end do
    call raise(ginv, vcov, vcon)
    energy = 0.0_dp
    do k = 1, 3
      energy = energy + 0.5_dp/st%mass*vcov(k)*vcon(k)
    end do
  end function cp_explicit_energy

  ! Raise a covariant vector: v^i = g^ij v_j.
  !$acc routine seq
  pure subroutine raise(ginv, vcov, vcon)
    real(dp), intent(in) :: ginv(3,3), vcov(3)
    real(dp), intent(out) :: vcon(3)
    integer :: i
    do i = 1, 3
      vcon(i) = ginv(i,1)*vcov(1) + ginv(i,2)*vcov(2) + ginv(i,3)*vcov(3)
    end do
  end subroutine raise

end module orbit_cp_explicit
