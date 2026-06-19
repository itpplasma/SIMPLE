module orbit_cpp_canonical
  ! Curvilinear canonical-midpoint 6D port of the Egger-Feiel thesis integrators
  ! (DVI_python: cp_sym_midpoint.py, cpp_sym_midpoint.py, cpp_var_midpoint.py).
  !
  ! This SUPERSEDES the Cartesian orbit_cpp_pauli discretization. The thesis
  ! scheme works in curvilinear (r,theta,phi) with the contravariant metric in
  ! the position equation and the geodesic metric-derivative force in the
  ! momentum equation. Three models, integer-dispatched:
  !   MODEL_CP     full classical charged particle, gyro-resolved (dt=1)
  !   MODEL_CPP_SYM Pauli symplectic midpoint, H_full + mu|B|         (dt=80)
  !   MODEL_CPP_VAR Pauli variational midpoint, discrete Euler-Lagrange (dt=800)
  ! Coordinate block: COORD_TOK = inline analytic toroidal metric (only one wired
  ! here); COORD_VMEC reserved for the libneo metric_tensor/christoffel path.
  !
  ! 6D state z = (q1,q2,q3, p1,p2,p3) = (r,theta,phi, p_r,p_th,p_ph). q canonical,
  ! p canonical covariant. The position rows (1:3) solve the thesis midpoint; the
  ! momentum rows (4:6) carry p as explicit residual rows p_state - p_new(x), so
  ! the Jacobian is square 6x6 and the carried p (the python global side effect of
  ! F) becomes part of the root. p_new is linear in p_state, so the (4:6) rows
  ! decouple and Newton converges in the same iterations as the python 3D root.
  !
  ! GPU portability: fixed-size 6 state, integer (model,coord) dispatch, !$acc
  ! routine seq residual/Jacobian/LU, no class()/proc-ptr in the hot loop. The
  ! Jacobian is analytic except the tiny O(mu) |B| force, whose gradient is a
  ! 2-eval central difference of the field's own dBmod (the oracle-faithful dBmod
  ! is not a true gradient, so a closed Hessian would be inconsistent). Reuses
  ! rk_solve (device LU) from orbit_rk_core.
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
    real(dp) :: z(6)      = 0.0_dp   ! (r,th,ph, p_r,p_th,p_ph)
    real(dp) :: pold(3)   = 0.0_dp   ! carried covariant p_i of the previous step
    real(dp) :: dpdtold(3) = 0.0_dp  ! variational carry: dL/dq_i of previous step
    real(dp) :: mu        = 0.0_dp
    real(dp) :: dt        = 0.0_dp
    real(dp) :: mass      = 1.0_dp
    real(dp) :: charge    = 1.0_dp
    integer  :: model     = MODEL_CP
    integer  :: coord     = COORD_TOK
  end type cpp_canon_state_t

contains

  ! Metric + field block at q=(r,th,ph). Returns the contravariant/covariant
  ! diagonal metric, the metric direction-derivatives d_g(i,k)=dg_ii/dq_k and
  ! d2_g(i,k,l)=d2 g_ii/dq_k dq_l, and the field_can_t carrying A,dA,d2A,
  ! Bmod,dBmod,d2Bmod. mode_secders>0 fills the d2 blocks for the Jacobian.
  subroutine eval_block(coord, q, mode_secders, fc, gii, ginv, d_g, d2_g)
    !$acc routine seq
    integer, intent(in) :: coord
    real(dp), intent(in) :: q(3)
    integer, intent(in) :: mode_secders
    type(field_can_t), intent(inout) :: fc
    real(dp), intent(out) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3)
    real(dp) :: r, cth, sth, Rr

    select case (coord)
    case default  ! COORD_TOK: analytic toroidal metric, R0=1.
      r = q(1); cth = cos(q(2)); sth = sin(q(2)); Rr = 1.0_dp + r*cth
      gii = [1.0_dp, r*r, Rr*Rr]
      ginv = [1.0_dp, 1.0_dp/(r*r), 1.0_dp/(Rr*Rr)]
      d_g = 0.0_dp
      d_g(2,1) = 2.0_dp*r
      d_g(3,1) = 2.0_dp*Rr*cth                 ! dg33/dr
      d_g(3,2) = -2.0_dp*r*Rr*sth              ! dg33/dth  (CORRECT: factor r)
      call eval_field_correct_test(fc, q(1), q(2), q(3), mode_secders)
      if (mode_secders <= 0) return
      d2_g = 0.0_dp
      d2_g(2,1,1) = 2.0_dp                    ! d2 g22 / dr dr
      d2_g(3,1,1) = 2.0_dp*cth*cth            ! d2 g33 / dr dr
      d2_g(3,1,2) = -2.0_dp*(2.0_dp*r*cth + 1.0_dp)*sth  ! d2 g33 / dr dth
      d2_g(3,2,1) = d2_g(3,1,2)
      d2_g(3,2,2) = 2.0_dp*r*(r*sth*sth - Rr*cth)         ! d2 g33 / dth dth
    end select
  end subroutine eval_block

  ! Lagrangian gradient dL/dq_k at (vmid, midpoint field/metric block). The
  ! +mu|B| Pauli term is gated by mu_active so MODEL_CP folds it out at compile
  ! time of the branch (the only branch inside the residual arithmetic).
  pure subroutine dLdq(mass, charge, mu, mu_active, vmid, fc, d_g, out)
    !$acc routine seq
    real(dp), intent(in) :: mass, charge, mu, vmid(3)
    logical, intent(in) :: mu_active
    type(field_can_t), intent(in) :: fc
    real(dp), intent(in) :: d_g(3,3)
    real(dp), intent(out) :: out(3)
    real(dp) :: qc
    integer :: k

    qc = charge/c
    do k = 1, 3
      out(k) = 0.5_dp*mass*(d_g(1,k)*vmid(1)**2 + d_g(2,k)*vmid(2)**2 &
                            + d_g(3,k)*vmid(3)**2) &
               + qc*(fc%dAth(k)*vmid(2) + fc%dAph(k)*vmid(3))
      if (mu_active) out(k) = out(k) - mu*fc%dBmod(k)
    end do
  end subroutine dLdq

  ! Symplectic-midpoint residual shared by MODEL_CP (mu_active=.false.) and
  ! MODEL_CPP_SYM (.true.). q rows: q-qold - dt/m g^ii (pmid - qe/c Amid).
  ! p rows: p_state - p_new with p_new = pold + dt dLdq(vmid) (linear, decoupled).
  subroutine sym_residual(st, mu_active, zold, z, fvec)
    !$acc routine seq
    type(cpp_canon_state_t), intent(in) :: st
    logical, intent(in) :: mu_active
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)
    type(field_can_t) :: fc
    real(dp) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3)
    real(dp) :: qmid(3), vmid(3), grad(3), pmid(3), Amid(3), qc
    integer :: k

    qmid = 0.5_dp*(zold(1:3) + z(1:3))
    vmid = (z(1:3) - zold(1:3))/st%dt
    call eval_block(st%coord, qmid, 0, fc, gii, ginv, d_g, d2_g)
    call dLdq(st%mass, st%charge, st%mu, mu_active, vmid, fc, d_g, grad)

    qc = st%charge/c
    Amid = [0.0_dp, fc%Ath, fc%Aph]
    pmid = st%pold + 0.5_dp*st%dt*grad
    do k = 1, 3
      fvec(k) = z(k) - zold(k) - st%dt/st%mass*ginv(k)*(pmid(k) - qc*Amid(k))
      fvec(3+k) = z(3+k) - (st%pold(k) + st%dt*grad(k))
    end do
  end subroutine sym_residual

  ! Variational-midpoint residual (MODEL_CPP_VAR): discrete Euler-Lagrange.
  ! p rows carry p = m g_ii vmid + qe/c Amid; residual q rows:
  ! (dpdt + dLdxold) dt/2 - (p - dLdxdotold). Carries dpdt->dpdtold, p->pold.
  subroutine var_residual(st, zold, z, fvec)
    !$acc routine seq
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)
    type(field_can_t) :: fc
    real(dp) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3)
    real(dp) :: qmid(3), vmid(3), dpdt(3), pnew(3), Amid(3), qc
    integer :: k

    qmid = 0.5_dp*(zold(1:3) + z(1:3))
    vmid = (z(1:3) - zold(1:3))/st%dt
    call eval_block(st%coord, qmid, 0, fc, gii, ginv, d_g, d2_g)
    call dLdq(st%mass, st%charge, st%mu, .true., vmid, fc, d_g, dpdt)

    qc = st%charge/c
    Amid = [0.0_dp, fc%Ath, fc%Aph]
    do k = 1, 3
      pnew(k) = st%mass*gii(k)*vmid(k) + qc*Amid(k)
      fvec(k) = (dpdt(k) + st%dpdtold(k))*0.5_dp*st%dt - (pnew(k) - st%pold(k))
      fvec(3+k) = z(3+k) - pnew(k)
    end do
  end subroutine var_residual

  ! Model-dispatched residual (integer select, GPU-portable).
  subroutine residual(st, zold, z, fvec)
    !$acc routine seq
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

  ! Analytic 6x6 Jacobian dF/dz. The position rows depend on z(1:3) only (pold is
  ! fixed during the step), so dF_q/dp_state = 0 and the p rows are linear, giving
  ! the block structure [Jqq 0; Jpq I]. Jqq and Jpq come from the field/metric
  ! second derivatives; assembled here from eval_block(mode=2).
  subroutine jacobian(st, zold, z, jac)
    !$acc routine seq
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(in) :: zold(6), z(6)
    real(dp), intent(out) :: jac(6,6)
    type(field_can_t) :: fc
    real(dp) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3)
    real(dp) :: qmid(3), vmid(3), grad(3), Amid(3), dA(3,3), qc, mu_use
    real(dp) :: dgrad_dx(3,3), dginv_dx(3,3), dgii_dx(3,3)
    integer :: k, j
    logical :: is_var

    is_var = (st%model == MODEL_CPP_VAR)
    mu_use = st%mu
    if (st%model == MODEL_CP) mu_use = 0.0_dp

    qmid = 0.5_dp*(zold(1:3) + z(1:3))
    vmid = (z(1:3) - zold(1:3))/st%dt
    call eval_block(st%coord, qmid, 2, fc, gii, ginv, d_g, d2_g)
    qc = st%charge/c
    Amid = [0.0_dp, fc%Ath, fc%Aph]
    ! dA(i,k) = dA_i / dq_k (A_r = 0). dAth/dAph are gradients in fc.
    dA = 0.0_dp
    do k = 1, 3
      dA(2,k) = fc%dAth(k)
      dA(3,k) = fc%dAph(k)
    end do

    ! dGii/dx_k and dgii/dx_k along each q (only k=1,2 nonzero here).
    ! Diagonal metric is positive-definite, so gii > 0; d(g^ii)/dx = -d(g_ii)/x/g_ii^2.
    do k = 1, 3
      do j = 1, 3
        dgii_dx(j,k) = d_g(j,k)
        dginv_dx(j,k) = -d_g(j,k)/(gii(j)*gii(j))
      end do
    end do

    ! d(dLdq_k)/dx_j: vmid depends on z via 1/dt, qmid (field/metric) via 1/2.
    call grad_jacobian(st%coord, qmid, st%mass, qc, mu_use, vmid, fc, d_g, d2_g, &
                       st%dt, dgrad_dx)
    call dLdq(st%mass, st%charge, mu_use, st%model /= MODEL_CP, vmid, fc, d_g, grad)

    jac = 0.0_dp
    if (.not. is_var) then
      ! Symplectic: q-row k = z(k) - zold(k) - dt/m Gii_k (pmid_k - qc Amid_k),
      ! pmid_k = pold_k + dt/2 dLdq_k.  d/dx_j:
      do k = 1, 3
        do j = 1, 3
          jac(k,j) = -st%dt/st%mass*( &
              0.5_dp*dginv_dx(k,j)*(st%pold(k) + 0.5_dp*st%dt*grad(k) - qc*Amid(k)) &
              + ginv(k)*(0.5_dp*st%dt*dgrad_dx(k,j) - qc*0.5_dp*dA(k,j)) )
        end do
        jac(k,k) = jac(k,k) + 1.0_dp
        ! p-row: p_state_k - (pold_k + dt dLdq_k); d/dx_j = -dt dgrad_dx(k,j).
        do j = 1, 3
          jac(3+k,j) = -st%dt*dgrad_dx(k,j)
        end do
        jac(3+k,3+k) = 1.0_dp
      end do
    else
      ! Variational: q-row k = (dpdt_k + dpdtold_k) dt/2 - (pnew_k - pold_k),
      ! pnew_k = m gii_k vmid_k + qc Amid_k.
      do k = 1, 3
        do j = 1, 3
          ! d pnew_k / dx_j = m (dgii_k/dx_j) vmid_k (*1/2 for qmid)
          !                 + m gii_k (d vmid_k/dx_j) + qc dA_k/dx_j (*1/2)
          jac(k,j) = 0.5_dp*st%dt*dgrad_dx(k,j) &
                     - ( 0.5_dp*st%mass*dgii_dx(k,j)*vmid(k) &
                         + 0.5_dp*qc*dA(k,j) )
        end do
        ! d vmid_k/dz(k) = 1/dt enters pnew: -(m gii_k / dt)
        jac(k,k) = jac(k,k) - st%mass*gii(k)/st%dt
        ! p-row: p_state_k - pnew_k.
        do j = 1, 3
          jac(3+k,j) = -( 0.5_dp*st%mass*dgii_dx(k,j)*vmid(k) + 0.5_dp*qc*dA(k,j) )
        end do
        jac(3+k,k) = jac(3+k,k) - st%mass*gii(k)/st%dt
        jac(3+k,3+k) = 1.0_dp
      end do
    end if
  end subroutine jacobian

  ! d(dLdq_k)/dx_j with vmid=(z-zold)/dt and qmid=(z+zold)/2 both depending on z.
  ! Explicit (vmid) part scales 1/dt; midpoint (qmid) part scales 1/2. The mu|B|
  ! term uses d(dBmod_k)/dq_j; field_correct_test's dBmod is the (intentionally
  ! oracle-faithful) python form whose Hessian is not symmetric, so its gradient
  ! block is taken by a tight central difference of the field's own dBmod at qmid
  ! -- fully consistent with the residual and GPU-portable (just field evals, no
  ! class/proc-ptr), not a Jacobian approximation of a different field.
  subroutine grad_jacobian(coord, qmid, mass, qc, mu, vmid, fc, d_g, d2_g, dt, dgrad_dx)
    !$acc routine seq
    integer, intent(in) :: coord
    real(dp), intent(in) :: qmid(3), mass, qc, mu, vmid(3), dt
    type(field_can_t), intent(in) :: fc
    real(dp), intent(in) :: d_g(3,3), d2_g(3,3,3)
    real(dp), intent(out) :: dgrad_dx(3,3)
    real(dp) :: d2A(3,3,3), dBgrad(3,3)
    integer :: k, j, i

    ! Expand packed d2A (order drdr,drdth,drdph,dthdth,dthdph,dphdph) into 3x3.
    call expand_sym(fc%d2Ath, d2A(2,:,:))
    call expand_sym(fc%d2Aph, d2A(3,:,:))
    d2A(1,:,:) = 0.0_dp
    dBgrad = 0.0_dp
    if (mu > 0.0_dp) call dBmod_grad(coord, qmid, dBgrad)

    do k = 1, 3
      do j = 1, 3
        ! Geodesic term m/2 sum_i d_g(i,k) vmid_i^2:
        !   explicit (vmid) part: m/2 d_g(j,k) * 2 vmid_j * (1/dt) [i=j only]
        !   midpoint  (qmid) part: m/2 sum_i d2_g(i,k,j) vmid_i^2 * (1/2)
        dgrad_dx(k,j) = mass*d_g(j,k)*vmid(j)/dt
        do i = 1, 3
          dgrad_dx(k,j) = dgrad_dx(k,j) + 0.25_dp*mass*d2_g(i,k,j)*vmid(i)*vmid(i)
        end do
        ! EM term qc (dAth_k vmid2 + dAph_k vmid3):
        !   explicit: qc (dAth_k [j==2] + dAph_k [j==3]) / dt
        !   midpoint: qc (d2Ath_kj vmid2 + d2Aph_kj vmid3) * 1/2
        if (j == 2) dgrad_dx(k,j) = dgrad_dx(k,j) + qc*fc%dAth(k)/dt
        if (j == 3) dgrad_dx(k,j) = dgrad_dx(k,j) + qc*fc%dAph(k)/dt
        dgrad_dx(k,j) = dgrad_dx(k,j) &
            + 0.5_dp*qc*(d2A(2,k,j)*vmid(2) + d2A(3,k,j)*vmid(3))
        ! mu term -mu d(dBmod_k)/dq_j * 1/2 (midpoint only).
        dgrad_dx(k,j) = dgrad_dx(k,j) - 0.5_dp*mu*dBgrad(k,j)
      end do
    end do
  end subroutine grad_jacobian

  ! Gradient block dBgrad(k,j) = d(dBmod_k)/dq_j of the field's own dBmod, by a
  ! tight central difference. Consistent with whatever dBmod the residual uses.
  subroutine dBmod_grad(coord, q, dBgrad)
    !$acc routine seq
    integer, intent(in) :: coord
    real(dp), intent(in) :: q(3)
    real(dp), intent(out) :: dBgrad(3,3)
    type(field_can_t) :: fp, fm
    real(dp) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3), qp(3), qm(3)
    real(dp), parameter :: h = 1.0e-7_dp
    integer :: j

    dBgrad = 0.0_dp
    do j = 1, 2   ! phi-derivatives vanish (axisymmetric field)
      qp = q; qm = q; qp(j) = qp(j) + h; qm(j) = qm(j) - h
      call eval_block(coord, qp, 0, fp, gii, ginv, d_g, d2_g)
      call eval_block(coord, qm, 0, fm, gii, ginv, d_g, d2_g)
      dBgrad(1,j) = (fp%dBmod(1) - fm%dBmod(1))/(2.0_dp*h)
      dBgrad(2,j) = (fp%dBmod(2) - fm%dBmod(2))/(2.0_dp*h)
      dBgrad(3,j) = (fp%dBmod(3) - fm%dBmod(3))/(2.0_dp*h)
    end do
  end subroutine dBmod_grad

  ! Expand packed symmetric second-derivative (6) into a full 3x3 block.
  pure subroutine expand_sym(packed, blk)
    !$acc routine seq
    real(dp), intent(in) :: packed(6)
    real(dp), intent(out) :: blk(3,3)
    blk(1,1) = packed(1)
    blk(1,2) = packed(2); blk(2,1) = packed(2)
    blk(1,3) = packed(3); blk(3,1) = packed(3)
    blk(2,2) = packed(4)
    blk(2,3) = packed(5); blk(3,2) = packed(5)
    blk(3,3) = packed(6)
  end subroutine expand_sym

  ! Initialize the 6D state. CP: vel=(sqrt(g^rr 2 mu B),0,0), p=(g_rr vr, qe/c Ath,
  ! qe/c Aph). CPP-sym/var: vel=0, p=(0, qe/c Ath, qe/c Aph); var also sets
  ! dpdt0 = -mu dB. mu is taken from vperp0 (CP) or passed for the CPP models.
  subroutine cpp_canon_init(st, model, coord, x0, vpar0, vperp0, mu_in, &
                            mass, charge, dt)
    type(cpp_canon_state_t), intent(out) :: st
    integer, intent(in) :: model, coord
    real(dp), intent(in) :: x0(3), vpar0, vperp0, mu_in, mass, charge, dt
    type(field_can_t) :: fc
    real(dp) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3), vel(3), qc

    st%model = model
    st%coord = coord
    st%mass = mass
    st%charge = charge
    st%dt = dt
    st%z(1:3) = x0
    qc = charge/c

    call eval_block(coord, x0, 0, fc, gii, ginv, d_g, d2_g)

    select case (model)
    case (MODEL_CP)
      ! mu from the perpendicular (radial) gyration energy: mu = m vperp^2/(2|B|).
      st%mu = mass*vperp0*vperp0/(2.0_dp*fc%Bmod)
      vel = [sqrt(ginv(1)*2.0_dp*st%mu*fc%Bmod), 0.0_dp, 0.0_dp]
      st%pold = [gii(1)*vel(1), gii(2)*vel(2) + qc*fc%Ath, gii(3)*vel(3) + qc*fc%Aph]
    case (MODEL_CPP_SYM)
      st%mu = mu_in
      vel = vpar0*[0.0_dp, fc%hth, fc%hph]
      st%pold = [vel(1), vel(2) + qc*fc%Ath, vel(3) + qc*fc%Aph]
    case (MODEL_CPP_VAR)
      st%mu = mu_in
      st%pold = [0.0_dp, qc*fc%Ath, qc*fc%Aph]
      st%dpdtold = -st%mu*fc%dBmod
    end select
    st%z(4:6) = st%pold
  end subroutine cpp_canon_init

  ! One canonical-midpoint macro-step. Newton with the analytic 6x6 Jacobian and
  ! the device LU (rk_solve). Boundary guard keeps r in (0,1). Returns ierr/=0 on
  ! LU failure or non-convergence. Updates carried pold/dpdtold for the next step.
  subroutine cpp_canon_step(st, ierr)
    type(cpp_canon_state_t), intent(inout) :: st
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 50
    real(dp), parameter :: atol = 1.0e-13_dp, rtol = 1.0e-12_dp
    real(dp) :: zold(6), z(6), fvec(6), fjac(6,6), dz(6), reltol(6)
    real(dp) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3), vmid(3), Amid(3), qc
    type(field_can_t) :: fc
    integer :: kit, i, info
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
      ! Converged when ALL residuals are below atol (the Newton root), or when the
      ! full step is below the relative floor (no further progress possible). The
      ! two are independent criteria; a per-component mix would accept a stalled
      ! component with a large residual.
      res_conv = .true.; step_conv = .true.
      do i = 1, 6
        if (abs(fvec(i)) >= atol) res_conv = .false.
        if (abs(dz(i)) >= rtol*reltol(i)) step_conv = .false.
      end do
      if (res_conv .or. step_conv) exit
    end do

    if (kit > maxit) ierr = 3

    ! Carry forward. For the variational model dpdtold/pold are recomputed at the
    ! converged midpoint, mirroring the python globals updated inside F.
    if (st%model == MODEL_CPP_VAR) then
      vmid = (z(1:3) - zold(1:3))/st%dt
      call eval_block(st%coord, 0.5_dp*(zold(1:3)+z(1:3)), 0, fc, gii, ginv, d_g, d2_g)
      qc = st%charge/c
      Amid = [0.0_dp, fc%Ath, fc%Aph]
      call dLdq(st%mass, st%charge, st%mu, .true., vmid, fc, d_g, st%dpdtold)
      do i = 1, 3
        st%pold(i) = st%mass*gii(i)*vmid(i) + qc*Amid(i)
      end do
    else
      st%pold = z(4:6)
    end if
    st%z = z
  end subroutine cpp_canon_step

  ! Hamiltonian H = (1/2m)(p-qe/c A) g^ii (p-qe/c A) [+ mu|B|]. CP has no mu term.
  function cpp_canon_energy(st) result(energy)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp) :: energy
    type(field_can_t) :: fc
    real(dp) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3), vcov(3), qc
    integer :: k

    call eval_block(st%coord, st%z(1:3), 0, fc, gii, ginv, d_g, d2_g)
    qc = st%charge/c
    vcov = [st%z(4) - 0.0_dp, st%z(5) - qc*fc%Ath, st%z(6) - qc*fc%Aph]
    energy = 0.0_dp
    do k = 1, 3
      energy = energy + 0.5_dp/st%mass*ginv(k)*vcov(k)*vcov(k)
    end do
    if (st%model /= MODEL_CP) energy = energy + st%mu*fc%Bmod
  end function cpp_canon_energy

  ! Guiding-center reduction: position is q itself (canonical curvilinear chart);
  ! vpar = h_i v^i with v^i = g^ii (p_i - qe/c A_i)/m.
  subroutine cpp_canon_to_gc(st, r, th, ph, vpar)
    type(cpp_canon_state_t), intent(in) :: st
    real(dp), intent(out) :: r, th, ph, vpar
    type(field_can_t) :: fc
    real(dp) :: gii(3), ginv(3), d_g(3,3), d2_g(3,3,3), vcon(3), qc

    call eval_block(st%coord, st%z(1:3), 0, fc, gii, ginv, d_g, d2_g)
    qc = st%charge/c
    vcon(1) = ginv(1)*(st%z(4) - 0.0_dp)/st%mass
    vcon(2) = ginv(2)*(st%z(5) - qc*fc%Ath)/st%mass
    vcon(3) = ginv(3)*(st%z(6) - qc*fc%Aph)/st%mass
    r = st%z(1); th = st%z(2); ph = st%z(3)
    vpar = fc%hth*vcon(2) + fc%hph*vcon(3)
  end subroutine cpp_canon_to_gc

end module orbit_cpp_canonical
