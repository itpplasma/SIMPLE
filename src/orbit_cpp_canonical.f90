module orbit_cpp_canonical
  ! Curvilinear canonical-midpoint 6D port of the Egger-Feiel thesis integrators
  ! (DVI_python: cp_sym_midpoint.py, cpp_sym_midpoint.py, cpp_var_midpoint.py),
  ! generalized to ARBITRARY curvilinear coordinates with a full (non-diagonal)
  ! metric g_ij/g^ij and its direction derivatives dg_ij,k.
  !
  ! Hamiltonian  H = (1/2m)(p_i - qc A_i) g^ij (p_j - qc A_j) [+ mu|B|].
  !   q_dot^k = (1/m) g^kj (p_j - qc A_j) = v^k/m
  !   p_dot_k = qc A_j,k v^j + (m/2) g_ij,k v^i v^j [- mu |B|,k]
  ! discretized with the implicit midpoint (SIMPLE's existing scheme).
  !
  ! Three models, integer-dispatched:
  !   MODEL_CP     full classical charged particle, gyro-resolved (dt=1)
  !   MODEL_CPP_SYM Pauli symplectic midpoint, H_full + mu|B|         (dt=80)
  !   MODEL_CPP_VAR Pauli variational midpoint, discrete Euler-Lagrange (dt=800)
  ! Two coordinate blocks, integer-dispatched:
  !   COORD_TOK  analytic toroidal metric + exact-curl tokamak field, fully
  !              inline, !$acc routine seq, class-free, analytic Jacobian.
  !   COORD_VMEC real VMEC flux coordinates: SINGLE-SOURCE full metric g_ij/g^ij,
  !              analytic dg_ij,k, covariant A_i and dA, and |B|=sqrt(g_ij B^i B^j)
  !              with analytic dBmod -- all from one vmec_field_metric_eval (libneo
  !              splint_vmec_data_d2), so h_i g^ij h_j = 1 identically and dg is the
  !              true derivative of g. Host-side (splines); Jacobian is a SIMPLIFIED
  !              first-derivative analytic Jacobian of the same residual (g, ginv,
  !              dg, dA, dBmod; d2 terms dropped) -- self-consistent (needs dg=d g),
  !              smooth through v_par -> 0 where the old finite-difference Jacobian
  !              went noisy and spuriously ejected trapped particles.
  ! The diagonal toroidal metric is the special case of the general full-metric
  ! arithmetic (off-diagonals zero), so COORD_TOK reproduces the validated python
  ! oracle bit-for-bit while the same residual runs on a stellarator metric.
  !
  ! 6D state z = (q1,q2,q3, p1,p2,p3). q canonical, p canonical covariant. The
  ! position rows (1:3) solve the thesis midpoint; the momentum rows (4:6) carry p
  ! as explicit residual rows p_state - p_new(x), giving a square 6x6 Newton
  ! system solved with the device LU rk_solve from linalg_lu_device.
  !
  ! GPU portability: cpp_canon_step_tok is the device entry. The whole COORD_TOK
  ! chain (cpp_canon_step_tok -> residual_tok -> eval_block_tok / dLdq / raise /
  ! residual_blk, jacobian_analytic -> grad_jacobian_tok, rk_solve) is
  ! !$acc routine seq with fixed-size 6 state, integer model dispatch, analytic
  ! Jacobian, no class()/proc-ptr -- one particle per GPU thread. The host
  ! cpp_canon_step keeps the coord dispatcher so COORD_VMEC (libneo class dispatch
  ! + spline reads, host-only) shares the same residual math and Newton LU.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: twopi
  use field_can_base, only: field_can_t
  use field_can_test, only: eval_field_correct_test
  use linalg_lu_device, only: rk_solve
  implicit none
  private

  integer, parameter, public :: MODEL_CP = 0, MODEL_CPP_SYM = 1, MODEL_CPP_VAR = 2
  integer, parameter, public :: COORD_TOK = 0, COORD_VMEC = 1, COORD_CHARTMAP = 2
  ! COORD_BOOZER: single-source Boozer field+metric (boozer_field_metric), the
  ! straight-field-line chart the production GC runs in, so the 6D state shares the
  ! GC's Boozer angles and field. Same first-derivative analytic Jacobian as VMEC.
  integer, parameter, public :: COORD_BOOZER = 3

  ! Thesis normalization: e = m = c = 1. qe/c uses this c, not the physical
  ! CGS speed of light in util (which would make the magnetic coupling vanish).
  real(dp), parameter :: c = 1.0_dp

  public :: cpp_canon_state_t, cpp_canon_init, cpp_canon_step, cpp_canon_step_tok, &
            cpp_canon_energy, cpp_canon_to_gc, cpp_canon_boozer_guiding_center
  public :: residual, jacobian   ! exposed for the Jacobian FD self-check in tests

  type :: cpp_canon_state_t
    real(dp) :: z(6)      = 0.0_dp   ! (q1,q2,q3, p1,p2,p3)
    real(dp) :: pold(3)   = 0.0_dp   ! carried covariant p_i of the previous step
    real(dp) :: dpdtold(3) = 0.0_dp  ! variational carry: dL/dq_i of previous step
    real(dp) :: mu        = 0.0_dp
    real(dp) :: dt        = 0.0_dp
    real(dp) :: mass      = 1.0_dp
    real(dp) :: charge    = 1.0_dp
    ! Magnetic-coupling length normalization: the canonical momentum couples to A
    ! through qc = charge/(c*ro0). ro0=1 (default) reproduces the thesis e=m=c=1
    ! coupling qc=charge/c for COORD_TOK/COORD_VMEC. The SIMPLE-normalized GC wire
    ! sets ro0 = ro0_bar = ro0/sqrt(2) so that p_i = vpar*h_i + A_i/ro0_bar.
    real(dp) :: ro0       = 1.0_dp
    real(dp) :: pabs      = 0.0_dp   ! normalized particle speed (GC z(4)), carried
    integer  :: model     = MODEL_CP
    integer  :: coord     = COORD_TOK
  end type cpp_canon_state_t

  ! Full metric + field block at a point. The toroidal block fills the diagonal
  ! entries and leaves off-diagonals zero, so the general arithmetic reduces to
  ! the validated diagonal case.
  type :: block_t
    real(dp) :: g(3,3)    = 0.0_dp   ! covariant metric g_ij
    real(dp) :: ginv(3,3) = 0.0_dp   ! contravariant metric g^ij
    real(dp) :: dg(3,3,3) = 0.0_dp   ! dg(i,j,k) = d g_ij / d q_k
    real(dp) :: Acov(3)   = 0.0_dp   ! covariant vector potential A_i (A_1 = 0)
    real(dp) :: dA(3,3)   = 0.0_dp   ! dA(i,k) = d A_i / d q_k
    real(dp) :: Bmod      = 0.0_dp   ! field modulus |B|
    real(dp) :: dBmod(3)  = 0.0_dp   ! d|B|/dq_k
    real(dp) :: d2Bmod(6) = 0.0_dp   ! packed Hessian of |B| (1=rr,2=rth,3=rph,4=thth,5=thph,6=phph)
    real(dp) :: hcov(3)   = 0.0_dp   ! covariant unit field h_i
  end type block_t

contains

  ! Evaluate the full metric + field block at q. COORD_TOK fills the analytic
  ! diagonal block (its analytic Jacobian also uses d2g/d2A/d2Bmod); COORD_VMEC
  ! fills the single-source block whose first derivatives (dg, dA, dBmod) drive
  ! the simplified first-derivative analytic Jacobian.
  subroutine eval_block(coord, q, blk)
    integer, intent(in) :: coord
    real(dp), intent(in) :: q(3)
    type(block_t), intent(out) :: blk

    select case (coord)
    case (COORD_VMEC)
      call eval_block_vmec(q, blk)
    case (COORD_BOOZER)
      call eval_block_boozer(q, blk)
    case (COORD_CHARTMAP)
      call eval_block_chartmap(q, blk)
    case default
      call eval_block_tok(q, blk)
    end select
  end subroutine eval_block

  ! Single-source Boozer block (host-side): full metric g_ij/g^ij and its analytic
  ! derivative dg pulled back from the VMEC R,Z geometry through the Boozer angle
  ! map, with the field (A_i, |B|, dBmod, h_i) taken directly from the production
  ! Boozer splines so the 6D field equals the GC field. dg is the genuine
  ! derivative of g (test_boozer_field_metric: dg vs FD ~1e-11), so the
  ! first-derivative analytic Jacobian is self-consistent like the VMEC path.
  subroutine eval_block_boozer(q, blk)
    use boozer_field_metric, only: boozer_field_metric_eval
    real(dp), intent(in) :: q(3)
    type(block_t), intent(out) :: blk
    real(dp) :: sqrtg, Bctr(3), Bcov(3)

    call boozer_field_metric_eval(q, blk%g, blk%ginv, sqrtg, blk%dg, blk%Acov, &
         blk%dA, Bctr, Bcov, blk%Bmod, blk%dBmod, blk%hcov)
  end subroutine eval_block_boozer

  ! Analytic toroidal metric (R0=1) + exact-curl tokamak field. Diagonal metric;
  ! the only nonzero metric derivatives are dg22/dr, dg33/dr, dg33/dth (the latter
  ! carries the factor r: dg33/dth = -2 r (R0+r cos th) sin th). !$acc routine seq,
  ! class-free: the GPU-portable block.
  subroutine eval_block_tok(q, blk)
    !$acc routine seq
    real(dp), intent(in) :: q(3)
    type(block_t), intent(out) :: blk
    type(field_can_t) :: fc
    real(dp) :: r, cth, sth, Rr

    r = q(1); cth = cos(q(2)); sth = sin(q(2)); Rr = 1.0_dp + r*cth
    blk%g = 0.0_dp; blk%ginv = 0.0_dp; blk%dg = 0.0_dp
    blk%g(1,1) = 1.0_dp;  blk%g(2,2) = r*r;        blk%g(3,3) = Rr*Rr
    blk%ginv(1,1) = 1.0_dp; blk%ginv(2,2) = 1.0_dp/(r*r); blk%ginv(3,3) = 1.0_dp/(Rr*Rr)
    blk%dg(2,2,1) = 2.0_dp*r              ! dg22/dr
    blk%dg(3,3,1) = 2.0_dp*Rr*cth         ! dg33/dr
    blk%dg(3,3,2) = -2.0_dp*r*Rr*sth      ! dg33/dth  (CORRECT: factor r)

    call eval_field_correct_test(fc, q(1), q(2), q(3), 1)
    blk%Acov = [0.0_dp, fc%Ath, fc%Aph]
    blk%dA = 0.0_dp
    blk%dA(2,:) = fc%dAth
    blk%dA(3,:) = fc%dAph
    blk%Bmod = fc%Bmod
    blk%dBmod = fc%dBmod
    blk%d2Bmod = fc%d2Bmod
    blk%hcov = [0.0_dp, fc%hth, fc%hph]
  end subroutine eval_block_tok

  ! Real VMEC flux block (host-side). SINGLE-SOURCE: the full non-diagonal metric
  ! g_ij/g^ij, its analytic derivatives dg_ij,k, the covariant A_i and its gradient
  ! dA, and |B| = sqrt(g_ij B^i B^j) with its analytic gradient dBmod ALL come from
  ! one vmec_field_metric_eval call (libneo splint_vmec_data_d2). Two reasons this
  ! must be single-source: (1) |B| from the SAME g gives h_i g^ij h_j = 1 to
  ! round-off (the dual-source path gave 1.009); (2) the analytic Jacobian below
  ! requires dg to be the genuine derivative of g -- with the dual-source split
  ! (g from libneo metric_tensor, dg from a separate Christoffel call) dg is NOT
  ! the derivative of g, so an analytic Jacobian is inconsistent and Newton fails.
  ! Here dg is analytic from the same R,Z derivatives (test_vmec_field_metric: dg
  ! vs FD ~1e-8), so the first-derivative analytic Jacobian is self-consistent and
  ! Newton converges smoothly through v_par -> 0.
  subroutine eval_block_vmec(q, blk)
    use vmec_field_metric, only: vmec_field_metric_eval
    real(dp), intent(in) :: q(3)
    type(block_t), intent(out) :: blk
    real(dp) :: sqrtg, Bctr(3), Bcov(3)

    call vmec_field_metric_eval(q, blk%g, blk%ginv, sqrtg, blk%dg, blk%Acov, &
         blk%dA, Bctr, Bcov, blk%Bmod, blk%dBmod, blk%hcov)
  end subroutine eval_block_vmec

  ! Production Boozer/chartmap block (host-side). The 6D state runs in the chartmap
  ! coordinates u=(rho,theta_B,phi_B): the libneo chartmap metric/Christoffel from
  ! ref_coords is native there, and the production field_can (Boozer) field is
  ! reparametrized from s=rho^2 with dF/drho=2 rho dF/ds. This is THE chart whose
  ! metric matches the production field_can chart (libneo #322), so it backs the
  ! production macrostep. NOT GPU-portable (class-dispatched metric + spline field).
  subroutine eval_block_chartmap(q, blk)
    use orbit_cpp_chartmap_metric, only: chartmap_eval_metric, chartmap_eval_field
    real(dp), intent(in) :: q(3)
    type(block_t), intent(out) :: blk

    call chartmap_eval_metric(q, blk%g, blk%ginv, blk%dg)
    call chartmap_eval_field(q, blk%Acov, blk%dA, blk%Bmod, blk%dBmod, blk%hcov)
  end subroutine eval_block_chartmap

  ! Raise a covariant vector: v^i = g^ij v_j.
  pure subroutine raise(ginv, vcov, vcon)
    !$acc routine seq
    real(dp), intent(in) :: ginv(3,3), vcov(3)
    real(dp), intent(out) :: vcon(3)
    integer :: i
    do i = 1, 3
      vcon(i) = ginv(i,1)*vcov(1) + ginv(i,2)*vcov(2) + ginv(i,3)*vcov(3)
    end do
  end subroutine raise

  ! Metric-unit perpendicular direction for the CP gyration seed: the radial
  ! covariant direction e_r = (1,0,0) raised to e_r^i = g^i1, projected
  ! perpendicular to the field (subtract its h-parallel part using |h|^2=1), then
  ! normalized in the metric so g_ij eperp^i eperp^j = 1. On the diagonal tokamak
  ! (h_1 = 0, g^11 = 1) this reduces to eperp = (1,0,0). A fixed gyrophase: the
  ! O(rho*) FLR offset of the seeded gyro-center is the physics, not an error.
  subroutine perp_unit_dir(blk, eperp)
    type(block_t), intent(in) :: blk
    real(dp), intent(out) :: eperp(3)
    real(dp) :: er(3), hcon(3), hpar, nrm
    integer :: i, j

    er = [blk%ginv(1,1), blk%ginv(2,1), blk%ginv(3,1)]   ! e_r^i = g^i1
    call raise(blk%ginv, blk%hcov, hcon)                 ! h^i

    ! Parallel component along h: (h_i e_r^i) with |h|^2 = h_i h^i = 1.
    hpar = blk%hcov(1)*er(1) + blk%hcov(2)*er(2) + blk%hcov(3)*er(3)
    do i = 1, 3
      eperp(i) = er(i) - hpar*hcon(i)
    end do

    ! Normalize in the metric: |eperp|_g^2 = g_ij eperp^i eperp^j.
    nrm = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        nrm = nrm + blk%g(i,j)*eperp(i)*eperp(j)
      end do
    end do
    if (nrm > 0.0_dp) then
      eperp = eperp/sqrt(nrm)
    else
      eperp = [1.0_dp, 0.0_dp, 0.0_dp]
    end if
  end subroutine perp_unit_dir

  subroutine boozer_larmor_offset(g, sqrtg, hcov, Bmod, vperp_con, mass, qc, rho)
    real(dp), intent(in) :: g(3,3), sqrtg, hcov(3), Bmod, vperp_con(3)
    real(dp), intent(in) :: mass, qc
    real(dp), intent(out) :: rho(3)
    real(dp) :: vcov(3), factor
    integer :: i

    do i = 1, 3
      vcov(i) = g(i,1)*vperp_con(1) + g(i,2)*vperp_con(2) &
        + g(i,3)*vperp_con(3)
    end do

    factor = mass/(qc*Bmod*sqrtg)
    rho(1) = factor*(hcov(2)*vcov(3) - hcov(3)*vcov(2))
    rho(2) = factor*(hcov(3)*vcov(1) - hcov(1)*vcov(3))
    rho(3) = factor*(hcov(1)*vcov(2) - hcov(2)*vcov(1))
  end subroutine boozer_larmor_offset

  ! Lagrangian gradient dL/dq_k at (vmid, midpoint block), general full metric:
  !   dL/dq_k = (m/2) g_ij,k vmid^i vmid^j + qc A_i,k vmid^i [- mu |B|,k].
  ! mu_active gates the Pauli +mu|B| term so MODEL_CP folds it out.
  pure subroutine dLdq(mass, charge, ro0, mu, mu_active, vmid, blk, out)
    !$acc routine seq
    real(dp), intent(in) :: mass, charge, ro0, mu, vmid(3)
    logical, intent(in) :: mu_active
    type(block_t), intent(in) :: blk
    real(dp), intent(out) :: out(3)
    real(dp) :: qc, geo, em
    integer :: k, i, j

    qc = charge/(c*ro0)
    do k = 1, 3
      geo = 0.0_dp
      do j = 1, 3
        do i = 1, 3
          geo = geo + blk%dg(i,j,k)*vmid(i)*vmid(j)
        end do
      end do
      em = 0.0_dp
      do i = 1, 3
        em = em + blk%dA(i,k)*vmid(i)
      end do
      out(k) = 0.5_dp*mass*geo + qc*em
      if (mu_active) out(k) = out(k) - mu*blk%dBmod(k)
    end do
  end subroutine dLdq

  ! Symplectic-midpoint residual math on a pre-evaluated block, shared by
  ! MODEL_CP (mu_active=.false.) and MODEL_CPP_SYM (.true.). q rows:
  ! q-qold - dt/m g^kj (pmid_j - qc Amid_j). p rows: p_state - (pold + dt dLdq).
  ! Block-as-argument so the same math runs host (dispatcher) and device (TOK).
  pure subroutine sym_residual_blk(st, mu_active, zold, z, blk, fvec)
    !$acc routine seq
    type(cpp_canon_state_t), intent(in) :: st
    logical, intent(in) :: mu_active
    real(dp), intent(in) :: zold(6), z(6)
    type(block_t), intent(in) :: blk
    real(dp), intent(out) :: fvec(6)
    real(dp) :: vmid(3), grad(3), pmid(3), vcov(3), vcon(3), qc
    integer :: k

    vmid = (z(1:3) - zold(1:3))/st%dt
    call dLdq(st%mass, st%charge, st%ro0, st%mu, mu_active, vmid, blk, grad)

    qc = st%charge/(c*st%ro0)
    pmid = st%pold + 0.5_dp*st%dt*grad
    do k = 1, 3
      vcov(k) = pmid(k) - qc*blk%Acov(k)
    end do
    call raise(blk%ginv, vcov, vcon)
    do k = 1, 3
      fvec(k) = z(k) - zold(k) - st%dt/st%mass*vcon(k)
      fvec(3+k) = z(3+k) - (st%pold(k) + st%dt*grad(k))
    end do
  end subroutine sym_residual_blk

  ! Variational-midpoint residual math on a pre-evaluated block (MODEL_CPP_VAR):
  ! discrete Euler-Lagrange. p rows carry p = m g_ij vmid^j + qc Amid; q rows:
  ! (dpdt + dLdxold) dt/2 - (p - dLdxdotold).
  pure subroutine var_residual_blk(st, zold, z, blk, fvec)
    !$acc routine seq
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    type(block_t), intent(in) :: blk
    real(dp), intent(out) :: fvec(6)
    real(dp) :: vmid(3), dpdt(3), pnew(3), qc
    integer :: k, j

    vmid = (z(1:3) - zold(1:3))/st%dt
    call dLdq(st%mass, st%charge, st%ro0, st%mu, .true., vmid, blk, dpdt)

    qc = st%charge/(c*st%ro0)
    do k = 1, 3
      pnew(k) = qc*blk%Acov(k)
      do j = 1, 3
        pnew(k) = pnew(k) + st%mass*blk%g(k,j)*vmid(j)
      end do
      fvec(k) = (dpdt(k) + st%dpdtold(k))*0.5_dp*st%dt - (pnew(k) - st%pold(k))
      fvec(3+k) = z(3+k) - pnew(k)
    end do
  end subroutine var_residual_blk

  ! Model-dispatched residual math, block-as-argument. Integer dispatch only, no
  ! class()/proc-ptr: !$acc routine seq, the device residual core.
  pure subroutine residual_blk(st, zold, z, blk, fvec)
    !$acc routine seq
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    type(block_t), intent(in) :: blk
    real(dp), intent(out) :: fvec(6)

    select case (st%model)
    case (MODEL_CP)
      call sym_residual_blk(st, .false., zold, z, blk, fvec)
    case (MODEL_CPP_SYM)
      call sym_residual_blk(st, .true., zold, z, blk, fvec)
    case (MODEL_CPP_VAR)
      call var_residual_blk(st, zold, z, blk, fvec)
    case default
      fvec = 0.0_dp
    end select
  end subroutine residual_blk

  ! Host residual dispatcher: evaluate the block (TOK or VMEC) then the math.
  subroutine residual(st, zold, z, fvec)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)
    type(block_t) :: blk

    call eval_block(st%coord, 0.5_dp*(zold(1:3) + z(1:3)), blk)
    call residual_blk(st, zold, z, blk, fvec)
  end subroutine residual

  ! Device residual (COORD_TOK only): inline analytic block, no VMEC dispatch.
  subroutine residual_tok(st, zold, z, fvec)
    !$acc routine seq
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)
    type(block_t) :: blk

    call eval_block_tok(0.5_dp*(zold(1:3) + z(1:3)), blk)
    call residual_blk(st, zold, z, blk, fvec)
  end subroutine residual_tok

  ! Jacobian dF/dz. COORD_TOK uses the analytic diagonal-metric Jacobian (with
  ! d2g/d2A/d2Bmod, validated by the analytic-vs-FD self-check). COORD_VMEC uses a
  ! simplified FIRST-derivative analytic Jacobian built from the same block (g,
  ! ginv, dg, dA, dBmod) the residual uses, dropping the d2g/d2A/d2Bmod
  ! force-gradient terms. The dropped terms make it an APPROXIMATE Jacobian, but it
  ! is SMOOTH (the finite-difference Jacobian it replaces went noisy at banana
  ! turning points v_par -> 0 and spuriously ejected all trapped particles); Newton
  ! converges to the residual root with a smooth approximate Jacobian. Both feed
  ! the same portable Newton LU.
  subroutine jacobian(st, zold, z, jac)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: jac(6,6)

    if (st%coord == COORD_VMEC .or. st%coord == COORD_BOOZER) then
      call jacobian_vmec_analytic(st, zold, z, jac)   ! metric-based, dg-consistent
    else
      call jacobian_analytic(st, zold, z, jac)
    end if
  end subroutine jacobian

  ! Simplified first-derivative analytic Jacobian for the full-metric sym residual
  ! (COORD_VMEC, MODEL_CP / MODEL_CPP_SYM). It uses the SAME single-source block
  ! (g, ginv, dg, dA, Acov, dBmod) at qmid = (zold+z)/2 that the residual
  ! evaluates, where dg is the genuine derivative of g, so every term below uses
  ! ONLY first derivatives -- the second derivatives of g, A and |B| (d2g, d2A,
  ! d2Bmod) are dropped, the agreed simplification. The dropped terms make it an
  ! APPROXIMATE Jacobian, but it is self-consistent and SMOOTH (the FD Jacobian it
  ! replaces went noisy at v_par -> 0); Newton converges to the residual root.
  !
  ! sym residual:
  !   grad_k  = (m/2) dg_ij,k v^i v^j + qc dA_i,k v^i [- mu dBmod_k],  v=(z-zold)/dt
  !   pmid_l  = pold_l + (dt/2) grad_l
  !   vcov_l  = pmid_l - qc Acov_l,   vcon_k = ginv_kl vcov_l
  !   Fq_k    = z_k - zold_k - (dt/m) vcon_k
  !   Fp_k    = z_(3+k) - (pold_k + dt grad_k)
  ! With block first derivatives w.r.t. z_j = (1/2) d/dq_j (qmid carries the 1/2)
  ! and the explicit v dependence dv^i/dz_j = delta_ij/dt:
  !   dgrad_dz(k,j) = (m sum_l dg_jl,k v^l + qc dA_j,k)/dt          (d2 terms dropped)
  !   dginv_dz(k,l,j) = -(1/2) ginv_ka dg_ab,j ginv_bl              (from dg only)
  ! giving
  !   dFq_k/dz_j   = delta_kj - (dt/m)[ dginv_dz(k,l,j) vcov_l
  !                  + ginv_kl ( (dt/2) dgrad_dz(l,j) - (qc/2) dA_l,j ) ]
  !   dFq_k/dp_m   = 0      (pmid uses pold, not z(4:6))
  !   dFp_k/dz_j   = -dt dgrad_dz(k,j)
  !   dFp_k/dp_m   = delta_km
  subroutine jacobian_vmec_analytic(st, zold, z, jac)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: jac(6,6)
    type(block_t) :: blk
    real(dp) :: qmid(3), vmid(3), grad(3), vcov(3), qc, mu_use
    real(dp) :: dgrad_dz(3,3), dginv_dz(3,3,3)
    integer :: k, j, l, a, b
    logical :: mu_active

    qmid = 0.5_dp*(zold(1:3) + z(1:3))
    vmid = (z(1:3) - zold(1:3))/st%dt
    call eval_block(st%coord, qmid, blk)   ! VMEC or BOOZER single-source block
    qc = st%charge/(c*st%ro0)
    mu_active = (st%model /= MODEL_CP)
    mu_use = merge(st%mu, 0.0_dp, mu_active)

    ! dgrad_dz(k,j): explicit v dependence only (block d2 terms dropped). dLdq is
    ! symmetric in dg's first two indices, so the v-derivative collapses to one sum.
    do k = 1, 3
      do j = 1, 3
        dgrad_dz(k,j) = 0.0_dp
        do l = 1, 3
          dgrad_dz(k,j) = dgrad_dz(k,j) + blk%dg(j,l,k)*vmid(l)
        end do
        dgrad_dz(k,j) = (st%mass*dgrad_dz(k,j) + qc*blk%dA(j,k))/st%dt
      end do
    end do

    ! dginv_dz(k,l,j) = d g^kl / d z_j = -(1/2) g^ka (dg_ab,j) g^bl.
    do j = 1, 3
      do l = 1, 3
        do k = 1, 3
          dginv_dz(k,l,j) = 0.0_dp
          do a = 1, 3
            do b = 1, 3
              dginv_dz(k,l,j) = dginv_dz(k,l,j) &
                  - 0.5_dp*blk%ginv(k,a)*blk%dg(a,b,j)*blk%ginv(b,l)
            end do
          end do
        end do
      end do
    end do

    ! grad and vcov at the current iterate (vcov = pmid - qc Acov).
    call dLdq(st%mass, st%charge, st%ro0, mu_use, mu_active, vmid, blk, grad)
    do l = 1, 3
      vcov(l) = st%pold(l) + 0.5_dp*st%dt*grad(l) - qc*blk%Acov(l)
    end do

    jac = 0.0_dp
    do k = 1, 3
      do j = 1, 3
        jac(k,j) = 0.0_dp
        do l = 1, 3
          jac(k,j) = jac(k,j) + dginv_dz(k,l,j)*vcov(l) &
              + blk%ginv(k,l)*(0.5_dp*st%dt*dgrad_dz(l,j) - 0.5_dp*qc*blk%dA(l,j))
        end do
        jac(k,j) = -st%dt/st%mass*jac(k,j)
        jac(3+k,j) = -st%dt*dgrad_dz(k,j)
      end do
      jac(k,k) = jac(k,k) + 1.0_dp
      jac(3+k,3+k) = 1.0_dp
    end do
  end subroutine jacobian_vmec_analytic

  ! Analytic 6x6 Jacobian for the diagonal toroidal block (COORD_TOK). The
  ! position rows depend on z(1:3) only, so the p rows are linear: [Jqq 0; Jpq I].
  ! Metric/field first derivatives are analytic (in block_t); d2g and d2A come from
  ! central differences of the block's own dg/dA, while the mu|B| force gradient
  ! uses the block's analytic Hessian d2Bmod -- a true Hessian of the corrected
  ! |B|, validated by the analytic-vs-FD self-check. The diagonal metric keeps
  ! g^kj = ginv_kk delta_kj, so the q-row k couples to z(1:3) only through qmid.
  subroutine jacobian_analytic(st, zold, z, jac)
    !$acc routine seq
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: jac(6,6)
    type(block_t) :: blk
    real(dp) :: qmid(3), vmid(3), grad(3), qc, mu_use
    real(dp) :: dgrad_dx(3,3), dginv_dx(3,3), dgii_dx(3,3)
    integer :: k, j
    logical :: is_var

    is_var = (st%model == MODEL_CPP_VAR)
    mu_use = st%mu
    if (st%model == MODEL_CP) mu_use = 0.0_dp

    qmid = 0.5_dp*(zold(1:3) + z(1:3))
    vmid = (z(1:3) - zold(1:3))/st%dt
    call eval_block_tok(qmid, blk)   ! analytic Jacobian path is COORD_TOK only
    qc = st%charge/(c*st%ro0)

    ! Diagonal-metric derivative blocks: d(g_kk)/dx_j and d(g^kk)/dx_j.
    do k = 1, 3
      do j = 1, 3
        dgii_dx(k,j) = blk%dg(k,k,j)
        dginv_dx(k,j) = -blk%dg(k,k,j)/(blk%g(k,k)*blk%g(k,k))
      end do
    end do

    call grad_jacobian_tok(qmid, st%mass, qc, mu_use, vmid, blk, st%dt, dgrad_dx)
    call dLdq(st%mass, st%charge, st%ro0, mu_use, st%model /= MODEL_CP, vmid, blk, grad)

    jac = 0.0_dp
    if (.not. is_var) then
      do k = 1, 3
        do j = 1, 3
          jac(k,j) = -st%dt/st%mass*( &
              0.5_dp*dginv_dx(k,j)*(st%pold(k) + 0.5_dp*st%dt*grad(k) - qc*blk%Acov(k)) &
              + blk%ginv(k,k)*(0.5_dp*st%dt*dgrad_dx(k,j) - qc*0.5_dp*blk%dA(k,j)) )
        end do
        jac(k,k) = jac(k,k) + 1.0_dp
        do j = 1, 3
          jac(3+k,j) = -st%dt*dgrad_dx(k,j)
        end do
        jac(3+k,3+k) = 1.0_dp
      end do
    else
      do k = 1, 3
        do j = 1, 3
          jac(k,j) = 0.5_dp*st%dt*dgrad_dx(k,j) &
                     - ( 0.5_dp*st%mass*dgii_dx(k,j)*vmid(k) + 0.5_dp*qc*blk%dA(k,j) )
        end do
        jac(k,k) = jac(k,k) - st%mass*blk%g(k,k)/st%dt
        do j = 1, 3
          jac(3+k,j) = -( 0.5_dp*st%mass*dgii_dx(k,j)*vmid(k) + 0.5_dp*qc*blk%dA(k,j) )
        end do
        jac(3+k,k) = jac(3+k,k) - st%mass*blk%g(k,k)/st%dt
        jac(3+k,3+k) = 1.0_dp
      end do
    end if
  end subroutine jacobian_analytic

  ! d(dLdq_k)/dx_j for the diagonal toroidal block. vmid=(z-zold)/dt scales 1/dt;
  ! qmid=(z+zold)/2 scales 1/2. d2g and d2A are central differences of the block's
  ! own dg/dA at qmid; the mu|B| force gradient uses the block's TRUE analytic
  ! Hessian d2Bmod (closed form of |B|=sqrt(W)). All GPU-portable (block evals).
  subroutine grad_jacobian_tok(qmid, mass, qc, mu, vmid, blk, dt, dgrad_dx)
    !$acc routine seq
    real(dp), intent(in) :: qmid(3), mass, qc, mu, vmid(3), dt
    type(block_t), intent(in) :: blk
    real(dp), intent(out) :: dgrad_dx(3,3)
    type(block_t) :: bp, bm
    real(dp) :: qp(3), qm(3), d2g(3,3,3), d2A(3,3,3), dBgrad(3,3)
    real(dp), parameter :: h = 1.0e-7_dp
    integer :: k, j, i

    ! Central differences of dg, dA give the metric/A second derivatives.
    d2g = 0.0_dp; d2A = 0.0_dp
    do j = 1, 3
      qp = qmid; qm = qmid; qp(j) = qp(j) + h; qm(j) = qm(j) - h
      call eval_block_tok(qp, bp)
      call eval_block_tok(qm, bm)
      do k = 1, 3
        do i = 1, 3
          d2g(i,k,j) = (bp%dg(i,i,k) - bm%dg(i,i,k))/(2.0_dp*h)
        end do
        d2A(2,k,j) = (bp%dA(2,k) - bm%dA(2,k))/(2.0_dp*h)
        d2A(3,k,j) = (bp%dA(3,k) - bm%dA(3,k))/(2.0_dp*h)
      end do
    end do

    ! dBgrad(k,j) = d(d|B|/dq_k)/dq_j = analytic Hessian of |B| (packed -> dense).
    dBgrad(1,1) = blk%d2Bmod(1); dBgrad(1,2) = blk%d2Bmod(2); dBgrad(1,3) = blk%d2Bmod(3)
    dBgrad(2,1) = blk%d2Bmod(2); dBgrad(2,2) = blk%d2Bmod(4); dBgrad(2,3) = blk%d2Bmod(5)
    dBgrad(3,1) = blk%d2Bmod(3); dBgrad(3,2) = blk%d2Bmod(5); dBgrad(3,3) = blk%d2Bmod(6)

    do k = 1, 3
      do j = 1, 3
        dgrad_dx(k,j) = mass*blk%dg(j,j,k)*vmid(j)/dt
        do i = 1, 3
          dgrad_dx(k,j) = dgrad_dx(k,j) + 0.25_dp*mass*d2g(i,k,j)*vmid(i)*vmid(i)
        end do
        if (j == 2) dgrad_dx(k,j) = dgrad_dx(k,j) + qc*blk%dA(2,k)/dt
        if (j == 3) dgrad_dx(k,j) = dgrad_dx(k,j) + qc*blk%dA(3,k)/dt
        dgrad_dx(k,j) = dgrad_dx(k,j) &
            + 0.5_dp*qc*(d2A(2,k,j)*vmid(2) + d2A(3,k,j)*vmid(3))
        if (mu > 0.0_dp) dgrad_dx(k,j) = dgrad_dx(k,j) - 0.5_dp*mu*dBgrad(k,j)
      end do
    end do
  end subroutine grad_jacobian_tok

  ! Initialize the 6D state. CP: FULL velocity v^i = vpar0 h^i + vperp e_perp^i,
  ! vperp = sqrt(2 mu B/m), e_perp the metric-unit radial direction projected
  ! perpendicular to h (a fixed gyrophase); p = m g_ij v^j + qc A places the
  ! gyro-center within O(rho*) of the GC start. CPP-sym: vel along h; CPP-var:
  ! vel=0, p=qc A, dpdt0=-mu dB.
  subroutine cpp_canon_init(st, model, coord, x0, vpar0, vperp0, mu_in, &
                            mass, charge, dt, ro0_in)
    type(cpp_canon_state_t), intent(out) :: st
    integer, intent(in) :: model, coord
    real(dp), intent(in) :: x0(3), vpar0, vperp0, mu_in, mass, charge, dt
    real(dp), intent(in), optional :: ro0_in
    type(block_t) :: blk
    real(dp) :: vcon(3), eperp(3), qc, vperp
    integer :: i, j

    vcon = 0.0_dp
    st%model = model
    st%coord = coord
    st%mass = mass
    st%charge = charge
    st%dt = dt
    st%z(1:3) = x0
    if (present(ro0_in)) st%ro0 = ro0_in
    qc = charge/(c*st%ro0)

    call eval_block(coord, x0, blk)

    select case (model)
    case (MODEL_CP)
      ! Full classical particle: resolve the gyration, so seed the full velocity.
      ! mu from vperp0 (vperp = sqrt(2 mu B/m)); on the diagonal tokamak with
      ! h_1 = 0 the perpendicular direction reduces to the bare radial seed
      ! [vperp,0,0], so the COORD_TOK oracle is reproduced bit-for-bit.
      st%mu = mass*vperp0*vperp0/(2.0_dp*blk%Bmod)
      vperp = sqrt(2.0_dp*st%mu*blk%Bmod/mass)
      call perp_unit_dir(blk, eperp)
      call raise(blk%ginv, vpar0*blk%hcov, vcon)   ! parallel piece v_par^i
      do i = 1, 3
        vcon(i) = vcon(i) + vperp*eperp(i)         ! + perpendicular gyration
      end do
    case (MODEL_CPP_SYM)
      st%mu = mu_in
      ! Parallel start: v^i = vpar0 g^ij h_j (raise the covariant field direction).
      call raise(blk%ginv, vpar0*blk%hcov, vcon)
    case (MODEL_CPP_VAR)
      st%mu = mu_in
    end select

    ! p_i = m g_ij v^j + qc A_i.
    do i = 1, 3
      st%pold(i) = qc*blk%Acov(i)
      do j = 1, 3
        st%pold(i) = st%pold(i) + mass*blk%g(i,j)*vcon(j)
      end do
    end do
    if (model == MODEL_CPP_VAR) st%dpdtold = -st%mu*blk%dBmod
    st%z(4:6) = st%pold
  end subroutine cpp_canon_init

  ! One canonical-midpoint macro-step. Newton with the 6x6 Jacobian and the
  ! device LU (rk_solve). Boundary guard keeps q(1) in (0,1). Returns ierr/=0 on
  ! LU failure or non-convergence. Updates carried pold/dpdtold.
  subroutine cpp_canon_step(st, ierr)
    type(cpp_canon_state_t), intent(inout) :: st
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 50
    real(dp), parameter :: atol = 1.0e-13_dp, rtol = 1.0e-12_dp
    real(dp) :: zold(6), z(6), fvec(6), fjac(6,6), dz(6), reltol(6)
    type(block_t) :: blk
    real(dp) :: vmid(3), qc
    integer :: kit, i, info, j
    logical :: res_conv, step_conv

    zold = st%z
    z = zold
    ierr = 0

    do kit = 1, maxit
      if (z(1) <= 0.0_dp) z(1) = 1.0e-3_dp
      if (z(1) >= 1.0_dp) then
        ierr = 2
        return
      end if
      call residual(st, zold, z, fvec)
      call jacobian(st, zold, z, fjac)
      dz = fvec
      call rk_solve(6, fjac, dz, info)
      if (info /= 0) then
        ierr = 1
        return
      end if
      z = z - dz
      reltol(1) = 1.0_dp; reltol(2) = twopi; reltol(3) = twopi
      do i = 1, 3
        reltol(3+i) = max(abs(z(3+i)), 1.0_dp)
      end do
      res_conv = .true.; step_conv = .true.
      do i = 1, 6
        if (abs(fvec(i)) >= atol) res_conv = .false.
        if (abs(dz(i)) >= rtol*reltol(i)) step_conv = .false.
      end do
      if (res_conv .or. step_conv) exit
    end do

    if (kit > maxit) ierr = 3

    if (st%model == MODEL_CPP_VAR) then
      vmid = (z(1:3) - zold(1:3))/st%dt
      call eval_block(st%coord, 0.5_dp*(zold(1:3)+z(1:3)), blk)
      qc = st%charge/(c*st%ro0)
      call dLdq(st%mass, st%charge, st%ro0, st%mu, .true., vmid, blk, st%dpdtold)
      do i = 1, 3
        st%pold(i) = qc*blk%Acov(i)
        do j = 1, 3
          st%pold(i) = st%pold(i) + st%mass*blk%g(i,j)*vmid(j)
        end do
      end do
    else
      st%pold = z(4:6)
    end if
    st%z = z
  end subroutine cpp_canon_step

  ! Device COORD_TOK macro-step (!$acc routine seq): identical Newton iteration to
  ! cpp_canon_step, but hardwired to the analytic toroidal block so the whole
  ! kernel chain (residual_tok -> eval_block_tok/dLdq/raise, jacobian_analytic ->
  ! grad_jacobian_tok, rk_solve) is device-callable. Integer model dispatch only;
  ! no class()/proc-ptr; no VMEC branch. Runs one particle per GPU thread.
  subroutine cpp_canon_step_tok(st, ierr)
    !$acc routine seq
    type(cpp_canon_state_t), intent(inout) :: st
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 50
    real(dp), parameter :: atol = 1.0e-13_dp, rtol = 1.0e-12_dp
    real(dp) :: zold(6), z(6), fvec(6), fjac(6,6), dz(6), reltol(6)
    type(block_t) :: blk
    real(dp) :: vmid(3), qc
    integer :: kit, i, info, j
    logical :: res_conv, step_conv

    zold = st%z
    z = zold
    ierr = 0

    do kit = 1, maxit
      if (z(1) <= 0.0_dp) z(1) = 1.0e-3_dp
      if (z(1) >= 1.0_dp) then
        ierr = 2
        return
      end if
      call residual_tok(st, zold, z, fvec)
      call jacobian_analytic(st, zold, z, fjac)
      dz = fvec
      call rk_solve(6, fjac, dz, info)
      if (info /= 0) then
        ierr = 1
        return
      end if
      z = z - dz
      reltol(1) = 1.0_dp; reltol(2) = twopi; reltol(3) = twopi
      do i = 1, 3
        reltol(3+i) = max(abs(z(3+i)), 1.0_dp)
      end do
      res_conv = .true.; step_conv = .true.
      do i = 1, 6
        if (abs(fvec(i)) >= atol) res_conv = .false.
        if (abs(dz(i)) >= rtol*reltol(i)) step_conv = .false.
      end do
      if (res_conv .or. step_conv) exit
    end do

    if (kit > maxit) ierr = 3

    if (st%model == MODEL_CPP_VAR) then
      vmid = (z(1:3) - zold(1:3))/st%dt
      call eval_block_tok(0.5_dp*(zold(1:3)+z(1:3)), blk)
      qc = st%charge/(c*st%ro0)
      call dLdq(st%mass, st%charge, st%ro0, st%mu, .true., vmid, blk, st%dpdtold)
      do i = 1, 3
        st%pold(i) = qc*blk%Acov(i)
        do j = 1, 3
          st%pold(i) = st%pold(i) + st%mass*blk%g(i,j)*vmid(j)
        end do
      end do
    else
      st%pold = z(4:6)
    end if
    st%z = z
  end subroutine cpp_canon_step_tok

  ! Hamiltonian H = (1/2m)(p-qc A) g^ij (p-qc A) [+ mu|B|]. MODEL_CP omits the
  ! mu|B| term because the full charged particle resolves the perpendicular
  ! gyromotion directly: its kinetic energy (1/2m)(p-qcA)g(p-qcA) already contains
  ! the perpendicular kinetic energy. The Pauli models (CPP_SYM/CPP_VAR) drop the
  ! resolved gyromotion and reinstate it as the guiding-center mu|B| (the magnetic
  ! moment is the Pauli kinetic piece), so they add it. This is a model difference,
  ! not a discretization detail of midpoint vs stored p.
  function cpp_canon_energy(st) result(energy)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp) :: energy
    type(block_t) :: blk
    real(dp) :: vcov(3), vcon(3), qc
    integer :: k

    call eval_block(st%coord, st%z(1:3), blk)
    qc = st%charge/(c*st%ro0)
    do k = 1, 3
      vcov(k) = st%z(3+k) - qc*blk%Acov(k)
    end do
    call raise(blk%ginv, vcov, vcon)
    energy = 0.0_dp
    do k = 1, 3
      energy = energy + 0.5_dp/st%mass*vcov(k)*vcon(k)
    end do
    if (st%model /= MODEL_CP) energy = energy + st%mu*blk%Bmod
  end function cpp_canon_energy

  ! Guiding-center reduction: position is q itself; vpar = h_i v^i with
  ! v^i = g^ij (p_j - qc A_j)/m.
  subroutine cpp_canon_to_gc(st, r, th, ph, vpar)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(out) :: r, th, ph, vpar
    type(block_t) :: blk
    real(dp) :: vcov(3), vcon(3), qc
    integer :: k

    call eval_block(st%coord, st%z(1:3), blk)
    qc = st%charge/(c*st%ro0)
    do k = 1, 3
      vcov(k) = (st%z(3+k) - qc*blk%Acov(k))/st%mass
    end do
    call raise(blk%ginv, vcov, vcon)
    r = st%z(1); th = st%z(2); ph = st%z(3)
    vpar = blk%hcov(1)*vcon(1) + blk%hcov(2)*vcon(2) + blk%hcov(3)*vcon(3)
  end subroutine cpp_canon_to_gc

  subroutine cpp_canon_boozer_guiding_center(st, xgc)
    use boozer_field_metric, only: boozer_field_metric_eval
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(out) :: xgc(3)

    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: qc, vcov(3), vcon(3), hcon(3), vpar, vperp_con(3), rho(3)
    integer :: i

    if (st%coord /= COORD_BOOZER) error stop &
      'CP guiding-center reconstruction requires COORD_BOOZER'

    call boozer_field_metric_eval(st%z(1:3), g, ginv, sqrtg, dg, Acov, dA, &
      Bctr, Bcov, Bmod, dBmod, hcov)

    qc = st%charge/(c*st%ro0)
    do i = 1, 3
      vcov(i) = (st%z(3+i) - qc*Acov(i))/st%mass
      vcon(i) = ginv(i,1)*vcov(1) + ginv(i,2)*vcov(2) + ginv(i,3)*vcov(3)
      hcon(i) = ginv(i,1)*hcov(1) + ginv(i,2)*hcov(2) + ginv(i,3)*hcov(3)
    end do

    vpar = hcov(1)*vcon(1) + hcov(2)*vcon(2) + hcov(3)*vcon(3)
    vperp_con = vcon - vpar*hcon
    call boozer_larmor_offset(g, sqrtg, hcov, Bmod, vperp_con, st%mass, qc, rho)
    xgc = st%z(1:3) - rho
  end subroutine cpp_canon_boozer_guiding_center

end module orbit_cpp_canonical
