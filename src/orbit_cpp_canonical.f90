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
  !   COORD_VMEC real VMEC flux coordinates: full metric g_ij/g^ij + Christoffel
  !              from libneo (#322) via orbit_cpp_vmec_metric, covariant A_i and
  !              |B| from SIMPLE's native VMEC field. Host-side (libneo class +
  !              splines); Jacobian by finite difference of the same residual.
  ! The diagonal toroidal metric is the special case of the general full-metric
  ! arithmetic (off-diagonals zero), so COORD_TOK reproduces the validated python
  ! oracle bit-for-bit while the same residual runs on a stellarator metric.
  !
  ! 6D state z = (q1,q2,q3, p1,p2,p3). q canonical, p canonical covariant. The
  ! position rows (1:3) solve the thesis midpoint; the momentum rows (4:6) carry p
  ! as explicit residual rows p_state - p_new(x), giving a square 6x6 Newton
  ! system solved with the device LU rk_solve from orbit_rk_core.
  !
  ! GPU portability: COORD_TOK keeps fixed-size 6 state, integer dispatch, !$acc
  ! routine seq, analytic Jacobian, no class()/proc-ptr. COORD_VMEC is host-side
  ! by necessity (libneo class dispatch + spline reads); the Newton LU is the same
  ! portable kernel.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: twopi
  use field_can_base, only: field_can_t
  use field_can_test, only: eval_field_correct_test
  use orbit_rk_core, only: rk_solve
  implicit none
  private

  integer, parameter, public :: MODEL_CP = 0, MODEL_CPP_SYM = 1, MODEL_CPP_VAR = 2
  integer, parameter, public :: COORD_TOK = 0, COORD_VMEC = 1

  ! Thesis normalization: e = m = c = 1. qe/c uses this c, not the physical
  ! CGS speed of light in util (which would make the magnetic coupling vanish).
  real(dp), parameter :: c = 1.0_dp

  public :: cpp_canon_state_t, cpp_canon_init, cpp_canon_step, &
            cpp_canon_energy, cpp_canon_to_gc
  public :: residual, jacobian   ! exposed for the Jacobian FD self-check in tests

  type :: cpp_canon_state_t
    real(dp) :: z(6)      = 0.0_dp   ! (q1,q2,q3, p1,p2,p3)
    real(dp) :: pold(3)   = 0.0_dp   ! carried covariant p_i of the previous step
    real(dp) :: dpdtold(3) = 0.0_dp  ! variational carry: dL/dq_i of previous step
    real(dp) :: mu        = 0.0_dp
    real(dp) :: dt        = 0.0_dp
    real(dp) :: mass      = 1.0_dp
    real(dp) :: charge    = 1.0_dp
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
    real(dp) :: hcov(3)   = 0.0_dp   ! covariant unit field h_i
  end type block_t

contains

  ! Evaluate the full metric + field block at q. mode_secders unused here (the
  ! Jacobian uses analytic dg/dA for COORD_TOK and finite differences for the
  ! mu|B| force and the whole COORD_VMEC path).
  subroutine eval_block(coord, q, blk)
    integer, intent(in) :: coord
    real(dp), intent(in) :: q(3)
    type(block_t), intent(out) :: blk

    select case (coord)
    case (COORD_VMEC)
      call eval_block_vmec(q, blk)
    case default
      call eval_block_tok(q, blk)
    end select
  end subroutine eval_block

  ! Analytic toroidal metric (R0=1) + exact-curl tokamak field. Diagonal metric;
  ! the only nonzero metric derivatives are dg22/dr, dg33/dr, dg33/dth (the latter
  ! with the CORRECT factor r the python listing drops). !$acc routine seq,
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

    call eval_field_correct_test(fc, q(1), q(2), q(3), 0)
    blk%Acov = [0.0_dp, fc%Ath, fc%Aph]
    blk%dA = 0.0_dp
    blk%dA(2,:) = fc%dAth
    blk%dA(3,:) = fc%dAph
    blk%Bmod = fc%Bmod
    blk%dBmod = fc%dBmod
    blk%hcov = [0.0_dp, fc%hth, fc%hph]
  end subroutine eval_block_tok

  ! Real VMEC flux block (host-side). Full non-diagonal metric + Christoffel from
  ! libneo; covariant A_i and |B| from the native VMEC field. dA is taken by a
  ! central difference of A_i (the native evaluator returns analytic dA only in s).
  subroutine eval_block_vmec(q, blk)
    use orbit_cpp_vmec_metric, only: vmec_eval_metric, vmec_eval_field
    real(dp), intent(in) :: q(3)
    type(block_t), intent(out) :: blk
    real(dp) :: Ap(3), Am(3), Bmp, dBmp(3), hp(3), qp(3), qm(3)
    real(dp), parameter :: h = 1.0e-6_dp
    integer :: k

    call vmec_eval_metric(q, blk%g, blk%ginv, blk%dg)
    call vmec_eval_field(q, blk%Acov, blk%Bmod, blk%dBmod, blk%hcov)
    blk%dA = 0.0_dp
    do k = 1, 3
      qp = q; qm = q; qp(k) = qp(k) + h; qm(k) = qm(k) - h
      call vmec_eval_field(qp, Ap, Bmp, dBmp, hp)
      call vmec_eval_field(qm, Am, Bmp, dBmp, hp)
      blk%dA(:,k) = (Ap - Am)/(2.0_dp*h)
    end do
  end subroutine eval_block_vmec

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

  ! Lagrangian gradient dL/dq_k at (vmid, midpoint block), general full metric:
  !   dL/dq_k = (m/2) g_ij,k vmid^i vmid^j + qc A_i,k vmid^i [- mu |B|,k].
  ! mu_active gates the Pauli +mu|B| term so MODEL_CP folds it out.
  pure subroutine dLdq(mass, charge, mu, mu_active, vmid, blk, out)
    !$acc routine seq
    real(dp), intent(in) :: mass, charge, mu, vmid(3)
    logical, intent(in) :: mu_active
    type(block_t), intent(in) :: blk
    real(dp), intent(out) :: out(3)
    real(dp) :: qc, geo, em
    integer :: k, i, j

    qc = charge/c
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

  ! Symplectic-midpoint residual shared by MODEL_CP (mu_active=.false.) and
  ! MODEL_CPP_SYM (.true.). q rows: q-qold - dt/m g^kj (pmid_j - qc Amid_j).
  ! p rows: p_state - p_new with p_new = pold + dt dLdq(vmid).
  subroutine sym_residual(st, mu_active, zold, z, fvec)
    type(cpp_canon_state_t), intent(in) :: st
    logical, intent(in) :: mu_active
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)
    type(block_t) :: blk
    real(dp) :: qmid(3), vmid(3), grad(3), pmid(3), vcov(3), vcon(3), qc
    integer :: k

    qmid = 0.5_dp*(zold(1:3) + z(1:3))
    vmid = (z(1:3) - zold(1:3))/st%dt
    call eval_block(st%coord, qmid, blk)
    call dLdq(st%mass, st%charge, st%mu, mu_active, vmid, blk, grad)

    qc = st%charge/c
    pmid = st%pold + 0.5_dp*st%dt*grad
    do k = 1, 3
      vcov(k) = pmid(k) - qc*blk%Acov(k)
    end do
    call raise(blk%ginv, vcov, vcon)
    do k = 1, 3
      fvec(k) = z(k) - zold(k) - st%dt/st%mass*vcon(k)
      fvec(3+k) = z(3+k) - (st%pold(k) + st%dt*grad(k))
    end do
  end subroutine sym_residual

  ! Variational-midpoint residual (MODEL_CPP_VAR): discrete Euler-Lagrange.
  ! p rows carry p = m g_ij vmid^j + qc Amid; q rows:
  ! (dpdt + dLdxold) dt/2 - (p - dLdxdotold). Carries dpdt->dpdtold, p->pold.
  subroutine var_residual(st, zold, z, fvec)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)
    type(block_t) :: blk
    real(dp) :: qmid(3), vmid(3), dpdt(3), pnew(3), qc
    integer :: k, j

    qmid = 0.5_dp*(zold(1:3) + z(1:3))
    vmid = (z(1:3) - zold(1:3))/st%dt
    call eval_block(st%coord, qmid, blk)
    call dLdq(st%mass, st%charge, st%mu, .true., vmid, blk, dpdt)

    qc = st%charge/c
    do k = 1, 3
      pnew(k) = qc*blk%Acov(k)
      do j = 1, 3
        pnew(k) = pnew(k) + st%mass*blk%g(k,j)*vmid(j)
      end do
      fvec(k) = (dpdt(k) + st%dpdtold(k))*0.5_dp*st%dt - (pnew(k) - st%pold(k))
      fvec(3+k) = z(3+k) - pnew(k)
    end do
  end subroutine var_residual

  ! Model-dispatched residual.
  subroutine residual(st, zold, z, fvec)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)

    select case (st%model)
    case (MODEL_CP)
      call sym_residual(st, .false., zold, z, fvec)
    case (MODEL_CPP_SYM)
      call sym_residual(st, .true., zold, z, fvec)
    case (MODEL_CPP_VAR)
      call var_residual(st, zold, z, fvec)
    case default
      fvec = 0.0_dp
    end select
  end subroutine residual

  ! Jacobian dF/dz. COORD_TOK uses the analytic full-metric Jacobian (validated by
  ! the analytic-vs-FD self-check); COORD_VMEC uses a central-difference Jacobian
  ! of the same residual (the host metric/field are spline+FD based, so a closed
  ! Hessian would be inconsistent). Both feed the same portable Newton LU.
  subroutine jacobian(st, zold, z, jac)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: jac(6,6)

    if (st%coord == COORD_VMEC) then
      call jacobian_fd(st, zold, z, jac)
    else
      call jacobian_analytic(st, zold, z, jac)
    end if
  end subroutine jacobian

  ! Finite-difference Jacobian of the residual (host path).
  subroutine jacobian_fd(st, zold, z, jac)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: jac(6,6)
    real(dp) :: zp(6), zm(6), rp(6), rm(6), h
    integer :: j

    do j = 1, 6
      h = 1.0e-7_dp*max(abs(z(j)), 1.0_dp)
      zp = z; zm = z; zp(j) = zp(j) + h; zm(j) = zm(j) - h
      call residual(st, zold, zp, rp)
      call residual(st, zold, zm, rm)
      jac(:,j) = (rp - rm)/(2.0_dp*h)
    end do
  end subroutine jacobian_fd

  ! Analytic 6x6 Jacobian for the diagonal toroidal block (COORD_TOK). The
  ! position rows depend on z(1:3) only, so the p rows are linear: [Jqq 0; Jpq I].
  ! Metric/field first derivatives are analytic (in block_t); the second
  ! derivatives d2g, d2A and the mu|B| force gradient come from central
  ! differences of the block's own dg/dA/dBmod -- exact-consistent with the
  ! residual, GPU-portable (just block evals). The diagonal metric keeps
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
    qc = st%charge/c

    ! Diagonal-metric derivative blocks: d(g_kk)/dx_j and d(g^kk)/dx_j.
    do k = 1, 3
      do j = 1, 3
        dgii_dx(k,j) = blk%dg(k,k,j)
        dginv_dx(k,j) = -blk%dg(k,k,j)/(blk%g(k,k)*blk%g(k,k))
      end do
    end do

    call grad_jacobian_tok(qmid, st%mass, qc, mu_use, vmid, blk, st%dt, dgrad_dx)
    call dLdq(st%mass, st%charge, mu_use, st%model /= MODEL_CP, vmid, blk, grad)

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
  ! qmid=(z+zold)/2 scales 1/2. d2g, d2A and the mu|B| gradient are central
  ! differences of the block's own dg/dA/dBmod at qmid -- consistent with the
  ! residual whichever (oracle-faithful) form it uses, GPU-portable.
  subroutine grad_jacobian_tok(qmid, mass, qc, mu, vmid, blk, dt, dgrad_dx)
    !$acc routine seq
    real(dp), intent(in) :: qmid(3), mass, qc, mu, vmid(3), dt
    type(block_t), intent(in) :: blk
    real(dp), intent(out) :: dgrad_dx(3,3)
    type(block_t) :: bp, bm
    real(dp) :: qp(3), qm(3), d2g(3,3,3), d2A(3,3,3), dBgrad(3,3)
    real(dp), parameter :: h = 1.0e-7_dp
    integer :: k, j, i

    ! Central differences of dg, dA, dBmod give the diagonal second derivatives.
    d2g = 0.0_dp; d2A = 0.0_dp; dBgrad = 0.0_dp
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
        dBgrad(k,j) = (bp%dBmod(k) - bm%dBmod(k))/(2.0_dp*h)
      end do
    end do

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

  ! Initialize the 6D state. CP: vel=(v^r=sqrt(2 mu B/ (m g_rr)),0,0) so the
  ! radial gyration energy is mu B; p=g_ij v^j + qc A. CPP-sym: vel along h;
  ! CPP-var: vel=0, p=qc A, dpdt0=-mu dB.
  subroutine cpp_canon_init(st, model, coord, x0, vpar0, vperp0, mu_in, &
                            mass, charge, dt)
    type(cpp_canon_state_t), intent(out) :: st
    integer, intent(in) :: model, coord
    real(dp), intent(in) :: x0(3), vpar0, vperp0, mu_in, mass, charge, dt
    type(block_t) :: blk
    real(dp) :: vcon(3), qc
    integer :: i, j

    vcon = 0.0_dp
    st%model = model
    st%coord = coord
    st%mass = mass
    st%charge = charge
    st%dt = dt
    st%z(1:3) = x0
    qc = charge/c

    call eval_block(coord, x0, blk)

    select case (model)
    case (MODEL_CP)
      st%mu = mass*vperp0*vperp0/(2.0_dp*blk%Bmod)
      vcon = [sqrt(blk%ginv(1,1)*2.0_dp*st%mu*blk%Bmod), 0.0_dp, 0.0_dp]
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
      qc = st%charge/c
      call dLdq(st%mass, st%charge, st%mu, .true., vmid, blk, st%dpdtold)
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

  ! Hamiltonian H = (1/2m)(p-qc A) g^ij (p-qc A) [+ mu|B|]. CP has no mu term.
  function cpp_canon_energy(st) result(energy)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp) :: energy
    type(block_t) :: blk
    real(dp) :: vcov(3), vcon(3), qc
    integer :: k

    call eval_block(st%coord, st%z(1:3), blk)
    qc = st%charge/c
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
    qc = st%charge/c
    do k = 1, 3
      vcov(k) = (st%z(3+k) - qc*blk%Acov(k))/st%mass
    end do
    call raise(blk%ginv, vcov, vcon)
    r = st%z(1); th = st%z(2); ph = st%z(3)
    vpar = blk%hcov(1)*vcon(1) + blk%hcov(2)*vcon(2) + blk%hcov(3)*vcon(3)
  end subroutine cpp_canon_to_gc

end module orbit_cpp_canonical
