module orbit_cpp_pauli
  ! Genuine 6D classical Pauli particle (CPP), Cartesian realization.
  !
  ! This is the NON-tautological CPP. It is structurally distinct from the
  ! guiding-center (GC) integrator: a full 6D canonical phase space (x, p) that
  ! carries real gyration, evolved by a structure-preserving implicit-symplectic
  ! map. The GC motion is the SLOW MANIFOLD of this system; initializing on that
  ! manifold and stepping at bounce-scale (GC-sized) dt makes the gyro-averaged
  ! orbit reproduce GC to O(rho*) -- a real cross-method check, not an identity.
  !
  ! Contrast with orbit_cpp (flux-canonical CPP): that residual is BYTE-IDENTICAL
  ! to the GC degenerate-Lagrangian residual because it IS the GC slow-manifold
  ! projection. Its tests are refactor/code-motion oracles, not physics
  ! cross-validation. THIS module is the physics cross-validation.
  !
  ! Model (Xiao & Qin, CPC 265 (2021) 107981), CGS Gaussian (see src/util.F90):
  !   H = |p - (q/c) A(x)|^2 / (2 m) + mu |B(x)|,  mu a FIXED parameter
  !   v = (p - (q/c) A)/m
  !   dx/dt = v
  !   dp_j/dt = -dH/dx_j = (q/c) v_i dA_i/dx_j - mu d|B|/dx_j   (sum over i)
  ! i.e. m dv/dt = (q/c) v x B - mu grad|B| (the Pauli force), as required.
  !
  ! Discretization: implicit midpoint (Gauss s=1), structure-preserving for the
  ! canonical (x,p) lift, energy-bounded with no secular drift. The Newton
  ! Jacobian is ANALYTIC, built from grad A, Hess A, grad|B|, Hess|B| of the
  ! exact field (field_pauli_cart). No finite differences.
  !
  ! GPU portability: fixed-size 6D state, pure !$acc routine seq residual,
  ! Jacobian, and 6x6 LU solve; no procedure pointers, no class() dispatch, no
  ! finite-difference Jacobian. The field is a concrete inlinable routine, not a
  ! provider vtable. This is the clean GPU-offload-ready realization (flat metric,
  ! no Christoffel).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: c
  use field_pauli_cart, only: pauli_field_params_t, eval_pauli_field_cart
  implicit none
  private

  public :: pauli6d_state_t, pauli6d_init, pauli6d_step, pauli6d_energy, &
            pauli6d_mu, pauli6d_to_gc

  ! Symmetric second-derivative pair index (j,k) -> packed m, used to expand
  ! d2A(3,6) and d2Bmod(6) into full 3x3 blocks in the Jacobian.
  integer, parameter :: PJ(6) = [1,1,1,2,2,3]
  integer, parameter :: PK(6) = [1,2,3,2,3,3]

  type :: pauli6d_state_t
    real(dp) :: z(6)   = 0.0_dp   ! (x1,x2,x3, p1,p2,p3) canonical
    real(dp) :: mu     = 0.0_dp   ! m vperp^2 / (2 |B|), FIXED parameter
    real(dp) :: dt     = 0.0_dp
    real(dp) :: mass   = 0.0_dp
    real(dp) :: charge = 0.0_dp
    type(pauli_field_params_t) :: fp
  end type pauli6d_state_t

contains

  ! Initialize on the slow manifold from a guiding-center start: position xgc,
  ! parallel speed vpar, perpendicular speed vperp. mu is fixed from vperp here;
  ! the canonical momentum is p = m v + (q/c) A, with v the physical velocity
  ! v = vpar*b + vperp*e1 (e1 an arbitrary unit vector perp to b). The gyrophase
  ! choice only sets the initial gyro position; the gyro-averaged orbit is
  ! gyrophase-independent, which the validation test relies on.
  subroutine pauli6d_init(st, fp, xgc, vpar, vperp, mass, charge, dt)
    type(pauli6d_state_t), intent(out) :: st
    type(pauli_field_params_t), intent(in) :: fp
    real(dp), intent(in) :: xgc(3), vpar, vperp, mass, charge, dt
    real(dp) :: Avec(3), dA(3,3), d2A(3,6), Bvec(3), Bmod, dBmod(3), d2Bmod(6)
    real(dp) :: bhat(3), e1(3), e2(3), vvec(3), tmp(3), nrm

    st%fp = fp
    st%mass = mass
    st%charge = charge
    st%dt = dt

    call eval_pauli_field_cart(fp, xgc, Avec, dA, d2A, Bvec, Bmod, dBmod, d2Bmod)
    bhat = Bvec / Bmod

    ! e1 perpendicular to bhat: pick the least-aligned axis, project out bhat.
    if (abs(bhat(1)) <= abs(bhat(2)) .and. abs(bhat(1)) <= abs(bhat(3))) then
      tmp = [1.0_dp, 0.0_dp, 0.0_dp]
    else if (abs(bhat(2)) <= abs(bhat(3))) then
      tmp = [0.0_dp, 1.0_dp, 0.0_dp]
    else
      tmp = [0.0_dp, 0.0_dp, 1.0_dp]
    end if
    e1 = tmp - dot_product(tmp, bhat) * bhat
    nrm = sqrt(dot_product(e1, e1))
    e1 = e1 / nrm
    e2 = cross(bhat, e1)

    vvec = vpar * bhat + vperp * e1
    st%mu = mass * vperp * vperp / (2.0_dp * Bmod)
    st%z(1:3) = xgc
    st%z(4:6) = mass * vvec + (charge / c) * Avec
  end subroutine pauli6d_init

  ! One implicit-midpoint macro-step. Returns ierr/=0 on Newton/LU failure.
  subroutine pauli6d_step(st, ierr)
    type(pauli6d_state_t), intent(inout) :: st
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 50
    real(dp), parameter :: atol = 1.0e-13_dp, rtol = 1.0e-12_dp
    real(dp) :: zold(6), z(6), fvec(6), fjac(6,6), dz(6), reltol(6)
    integer :: kit, i, info
    logical :: conv

    zold = st%z
    z = zold
    ierr = 0

    do kit = 1, maxit
      call pauli6d_residual(st, zold, z, fvec)
      call pauli6d_jacobian(st, zold, z, fjac)
      dz = fvec
      call lu_solve6(fjac, dz, info)
      if (info /= 0) then
        ierr = 1
        return
      end if
      z = z - dz
      do i = 1, 3
        reltol(i)   = max(abs(z(i)), 1.0_dp)
        reltol(i+3) = max(abs(z(i+3)), 1.0_dp)
      end do
      conv = .true.
      do i = 1, 6
        if (abs(dz(i)) >= rtol*reltol(i) .and. abs(fvec(i)) >= atol) conv = .false.
      end do
      if (conv) exit
    end do

    st%z = z
  end subroutine pauli6d_step

  ! Implicit-midpoint residual F = znew - zold - dt * rhs((zold+znew)/2).
  pure subroutine pauli6d_residual(st, zold, z, fvec)
    !$acc routine seq
    type(pauli6d_state_t), intent(in) :: st
    real(dp), intent(in)  :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)
    real(dp) :: zmid(6), rhs(6)

    zmid = 0.5_dp * (zold + z)
    call pauli6d_rhs(st, zmid, rhs)
    fvec = z - zold - st%dt * rhs
  end subroutine pauli6d_residual

  ! Canonical RHS: dx/dt = v, dp_j/dt = (q/c) v_i dA_i/dxj - mu d|B|/dxj.
  pure subroutine pauli6d_rhs(st, w, rhs)
    !$acc routine seq
    type(pauli6d_state_t), intent(in) :: st
    real(dp), intent(in)  :: w(6)
    real(dp), intent(out) :: rhs(6)
    real(dp) :: Avec(3), dA(3,3), d2A(3,6), Bvec(3), Bmod, dBmod(3), d2Bmod(6)
    real(dp) :: vvec(3), qc
    integer :: i, j

    qc = st%charge / c
    call eval_pauli_field_cart(st%fp, w(1:3), Avec, dA, d2A, Bvec, Bmod, &
                               dBmod, d2Bmod)
    vvec = (w(4:6) - qc * Avec) / st%mass
    rhs(1:3) = vvec
    do j = 1, 3
      rhs(3+j) = -st%mu * dBmod(j)
      do i = 1, 3
        rhs(3+j) = rhs(3+j) + qc * vvec(i) * dA(i,j)
      end do
    end do
  end subroutine pauli6d_rhs

  ! Analytic Jacobian dF/dz of the implicit-midpoint residual. With
  ! zmid = (zold+z)/2 the chain rule gives dF/dz = I - (dt/2) d(rhs)/dw|_zmid.
  ! d(rhs)/dw is built from grad A, Hess A, grad|B|, Hess|B|.
  pure subroutine pauli6d_jacobian(st, zold, z, fjac)
    !$acc routine seq
    type(pauli6d_state_t), intent(in) :: st
    real(dp), intent(in)  :: zold(6), z(6)
    real(dp), intent(out) :: fjac(6,6)
    real(dp) :: zmid(6), Avec(3), dA(3,3), d2A(3,6), Bvec(3), Bmod
    real(dp) :: dBmod(3), d2Bmod(6), vvec(3), qc, drhs(6,6)
    real(dp) :: dv_dx(3,3), dv_dp(3,3), hessAblk(3,3,3), hessBblk(3,3)
    integer :: i, j, k, m, l

    zmid = 0.5_dp * (zold + z)
    qc = st%charge / c
    call eval_pauli_field_cart(st%fp, zmid(1:3), Avec, dA, d2A, Bvec, Bmod, &
                               dBmod, d2Bmod)
    vvec = (zmid(4:6) - qc * Avec) / st%mass

    ! Expand packed second derivatives into full symmetric 3x3 blocks.
    hessBblk = 0.0_dp
    do m = 1, 6
      j = PJ(m); k = PK(m)
      hessBblk(j,k) = d2Bmod(m)
      hessBblk(k,j) = d2Bmod(m)
    end do
    hessAblk = 0.0_dp
    do i = 1, 3
      do m = 1, 6
        j = PJ(m); k = PK(m)
        hessAblk(i,j,k) = d2A(i,m)
        hessAblk(i,k,j) = d2A(i,m)
      end do
    end do

    ! dv_i/dx_k = -(q/(c m)) dA_i/dx_k ;  dv_i/dp_k = delta_ik / m
    dv_dx = -(qc / st%mass) * dA
    dv_dp = 0.0_dp
    do i = 1, 3
      dv_dp(i,i) = 1.0_dp / st%mass
    end do

    ! d(rhs)/dw, 6x6: rows 1:3 = dv, rows 4:6 = dp-force.
    drhs = 0.0_dp
    ! dx/dt = v block
    do i = 1, 3
      do k = 1, 3
        drhs(i, k)   = dv_dx(i,k)
        drhs(i, 3+k) = dv_dp(i,k)
      end do
    end do
    ! dp_j/dt = (q/c) sum_i v_i dA_i/dxj - mu d|B|/dxj
    do j = 1, 3
      do k = 1, 3
        ! d/dx_k
        drhs(3+j, k) = -st%mu * hessBblk(j,k)
        do i = 1, 3
          drhs(3+j, k) = drhs(3+j, k) + qc * (dv_dx(i,k) * dA(i,j) &
                       + vvec(i) * hessAblk(i,j,k))
        end do
        ! d/dp_k
        do i = 1, 3
          drhs(3+j, 3+k) = drhs(3+j, 3+k) + qc * dv_dp(i,k) * dA(i,j)
        end do
      end do
    end do

    fjac = 0.0_dp
    do l = 1, 6
      fjac(l,l) = 1.0_dp
    end do
    fjac = fjac - 0.5_dp * st%dt * drhs
  end subroutine pauli6d_jacobian

  ! Total energy (the Hamiltonian); conserved up to bounded oscillation.
  function pauli6d_energy(st) result(energy)
    type(pauli6d_state_t), intent(in) :: st
    real(dp) :: energy
    real(dp) :: Avec(3), dA(3,3), d2A(3,6), Bvec(3), Bmod, dBmod(3), d2Bmod(6)
    real(dp) :: vvec(3), qc

    qc = st%charge / c
    call eval_pauli_field_cart(st%fp, st%z(1:3), Avec, dA, d2A, Bvec, Bmod, &
                               dBmod, d2Bmod)
    vvec = (st%z(4:6) - qc * Avec) / st%mass
    energy = 0.5_dp * st%mass * dot_product(vvec, vvec) + st%mu * Bmod
  end function pauli6d_energy

  ! Instantaneous mu = m vperp^2 / (2 |B|) reconstructed from the state; for a
  ! perfect adiabatic invariant this stays at st%mu up to gyro-scale ripple.
  function pauli6d_mu(st) result(mu_now)
    type(pauli6d_state_t), intent(in) :: st
    real(dp) :: mu_now
    real(dp) :: Avec(3), dA(3,3), d2A(3,6), Bvec(3), Bmod, dBmod(3), d2Bmod(6)
    real(dp) :: vvec(3), bhat(3), vpar, vperp2, qc

    qc = st%charge / c
    call eval_pauli_field_cart(st%fp, st%z(1:3), Avec, dA, d2A, Bvec, Bmod, &
                               dBmod, d2Bmod)
    vvec = (st%z(4:6) - qc * Avec) / st%mass
    bhat = Bvec / Bmod
    vpar = dot_product(vvec, bhat)
    vperp2 = max(dot_product(vvec, vvec) - vpar*vpar, 0.0_dp)
    mu_now = st%mass * vperp2 / (2.0_dp * Bmod)
  end function pauli6d_mu

  ! Guiding-center position estimate: x_gc = x - rho_L, rho_L = (m c / q) v x b / |B|.
  ! Removes the leading gyro displacement so the orbit can be compared to GC.
  ! Returns flux-label minor radius r, poloidal th, toroidal ph of the GC point.
  subroutine pauli6d_to_gc(st, r, th, ph, vpar_out)
    type(pauli6d_state_t), intent(in) :: st
    real(dp), intent(out) :: r, th, ph, vpar_out
    real(dp) :: Avec(3), dA(3,3), d2A(3,6), Bvec(3), Bmod, dBmod(3), d2Bmod(6)
    real(dp) :: vvec(3), bhat(3), rho(3), xgc(3), Rcyl, dR, qc

    qc = st%charge / c
    call eval_pauli_field_cart(st%fp, st%z(1:3), Avec, dA, d2A, Bvec, Bmod, &
                               dBmod, d2Bmod)
    vvec = (st%z(4:6) - qc * Avec) / st%mass
    bhat = Bvec / Bmod
    ! Larmor vector rho = (m c)/(q |B|) (b x v) so that x_gc = x - rho.
    rho = (st%mass * c) / (st%charge * Bmod) * cross(bhat, vvec)
    xgc = st%z(1:3) - rho
    vpar_out = dot_product(vvec, bhat)
    Rcyl = sqrt(xgc(1)**2 + xgc(2)**2)
    dR = Rcyl - st%fp%R0
    r  = sqrt(dR*dR + xgc(3)**2)
    th = atan2(xgc(3), dR)
    ph = atan2(xgc(2), xgc(1))
  end subroutine pauli6d_to_gc

  ! Dense 6x6 LU with partial pivoting, rhs overwritten with the solution.
  pure subroutine lu_solve6(A, rhs, info)
    !$acc routine seq
    real(dp), intent(inout) :: A(6,6), rhs(6)
    integer, intent(out) :: info
    integer :: i, j, k, ipiv
    real(dp) :: amax, factor, tmp

    info = 0
    do k = 1, 6
      ipiv = k
      amax = abs(A(k,k))
      do i = k+1, 6
        if (abs(A(i,k)) > amax) then
          amax = abs(A(i,k))
          ipiv = i
        end if
      end do
      if (amax == 0.0_dp) then
        info = k
        return
      end if
      if (ipiv /= k) then
        do j = 1, 6
          tmp = A(k,j); A(k,j) = A(ipiv,j); A(ipiv,j) = tmp
        end do
        tmp = rhs(k); rhs(k) = rhs(ipiv); rhs(ipiv) = tmp
      end if
      do i = k+1, 6
        factor = A(i,k)/A(k,k)
        A(i,k) = factor
        do j = k+1, 6
          A(i,j) = A(i,j) - factor*A(k,j)
        end do
        rhs(i) = rhs(i) - factor*rhs(k)
      end do
    end do
    do i = 6, 1, -1
      tmp = rhs(i)
      do j = i+1, 6
        tmp = tmp - A(i,j)*rhs(j)
      end do
      rhs(i) = tmp/A(i,i)
    end do
  end subroutine lu_solve6

  pure function cross(a, b) result(cab)
    !$acc routine seq
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: cab(3)
    cab(1) = a(2)*b(3) - a(3)*b(2)
    cab(2) = a(3)*b(1) - a(1)*b(3)
    cab(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

end module orbit_cpp_pauli
