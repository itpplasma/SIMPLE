module orbit_full_device
  ! GPU-offload-ready full-orbit pushers (B2). This is the device path: NO
  ! class() vtable dispatch and NO finite-difference Jacobian in the per-particle
  ! hot loop. Concrete field evaluation is selected by an integer field code via
  ! select case to inlinable !$acc routine seq helpers; the state is fixed-size;
  ! the symplectic (implicit-midpoint Lorentz) Newton uses an ANALYTIC Jacobian
  ! built from the analytic field gradient.
  !
  ! Contrast with orbit_full (the CPU path): that module keeps the abstract
  ! field_metric_provider_t seam for mock-based unit tests and curvilinear
  ! geometry with Christoffel symbols. It is NOT device-offloadable because it
  ! dispatches through a class() pointer and differentiates the residual by
  ! finite differences. This module is the clean Cartesian (flat-metric)
  ! realization that needs no Christoffel symbols and inlines onto the device.
  !
  ! GPU-offload-ready models here (Cartesian, flat metric):
  !   FOFIELD_UNIFORM  - constant B
  !   FOFIELD_LINGRAD  - linear B_i = B0_i + gradB(i,j) x_j
  !   FOFIELD_TOKAMAK  - analytic divergence-free circular tokamak (field_pauli_cart)
  ! All three carry an exact analytic grad B, so both the Boris rotation and the
  ! implicit-midpoint analytic Jacobian are device-pure. Curvilinear / provider
  ! mock models remain CPU-only in orbit_full.
  !
  ! Units: CGS Gaussian (see src/util.F90), m dv/dt = (q/c) v x B.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: c
  use field_pauli_cart, only: pauli_field_params_t, eval_pauli_field_cart
  implicit none
  private

  integer, parameter, public :: FOFIELD_UNIFORM = 1
  integer, parameter, public :: FOFIELD_LINGRAD = 2
  integer, parameter, public :: FOFIELD_TOKAMAK = 3

  ! Integrator codes (mirror orbit_full's ORBIT_* for the device subset).
  integer, parameter, public :: FODEV_BORIS   = 2
  integer, parameter, public :: FODEV_FOSYMPL = 3

  integer, parameter, public :: FODEV_OK = 0
  integer, parameter, public :: FODEV_ERR_NO_CONVERGE = 2

  ! Fixed-size device state. No pointers, no polymorphism: the whole struct is
  ! copyable to the device and every field eval is a pure select-case call.
  type, public :: fo_device_state_t
    real(dp) :: z(6)       = 0.0_dp     ! (x1,x2,x3, v1,v2,v3) Cartesian phys.
    real(dp) :: dt         = 0.0_dp
    real(dp) :: mass       = 0.0_dp
    real(dp) :: charge     = 0.0_dp
    integer  :: field_code = FOFIELD_UNIFORM
    integer  :: integrator = FODEV_BORIS
    ! Field parameters (only the ones the selected field_code uses are read).
    real(dp) :: B0(3)      = [0.0_dp, 0.0_dp, 1.0_dp]   ! uniform / lingrad bias
    real(dp) :: gradB(3,3) = 0.0_dp                     ! lingrad: dB_i/dx_j
    type(pauli_field_params_t) :: tok                   ! tokamak field params
  end type fo_device_state_t

  public :: fo_device_init, fo_device_step, fo_device_eval_field, &
            fo_device_energy

contains

  ! Initialize the device state and seed mu from the launch (for diagnostics).
  subroutine fo_device_init(st, x0, v0, field_code, integrator, mass, charge, dt)
    type(fo_device_state_t), intent(out) :: st
    real(dp), intent(in) :: x0(3), v0(3)
    integer,  intent(in) :: field_code, integrator
    real(dp), intent(in) :: mass, charge, dt

    st%z(1:3)    = x0
    st%z(4:6)    = v0
    st%field_code = field_code
    st%integrator = integrator
    st%mass      = mass
    st%charge    = charge
    st%dt        = dt
  end subroutine fo_device_init

  ! Concrete field evaluation by integer code. Returns B and analytic grad B
  ! (gB(i,j) = dB_i/dx_j). Pure and device-inlinable: the select case resolves
  ! to one of three straight-line bodies, no indirection.
  pure subroutine fo_device_eval_field(st, x, Bvec, Bmod, gB)
    !$acc routine seq
    type(fo_device_state_t), intent(in) :: st
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Bvec(3), Bmod, gB(3,3)
    real(dp) :: Av(3), dA(3,3), d2A(3,6), dBm(3), d2Bm(6)
    integer :: i

    select case (st%field_code)
    case (FOFIELD_UNIFORM)
      Bvec = st%B0
      gB = 0.0_dp
    case (FOFIELD_LINGRAD)
      do i = 1, 3
        Bvec(i) = st%B0(i) + st%gradB(i,1)*x(1) + st%gradB(i,2)*x(2) &
                + st%gradB(i,3)*x(3)
      end do
      gB = st%gradB
    case (FOFIELD_TOKAMAK)
      ! Reuse the exact analytic A; B = curl A and grad B come from A's second
      ! derivatives (pair index packing: 1:(x,x)2:(x,y)3:(x,z)4:(y,y)5:(y,z)6:(z,z)).
      call eval_pauli_field_cart(st%tok, x, Av, dA, d2A, Bvec, Bmod, dBm, d2Bm)
      ! dB_x/dxj = d2A_z/(dy dxj) - d2A_y/(dz dxj), etc.
      gB(1,1) = d2A(3,2) - d2A(2,3)
      gB(1,2) = d2A(3,4) - d2A(2,5)
      gB(1,3) = d2A(3,5) - d2A(2,6)
      gB(2,1) = d2A(1,3) - d2A(3,1)
      gB(2,2) = d2A(1,5) - d2A(3,2)
      gB(2,3) = d2A(1,6) - d2A(3,3)
      gB(3,1) = d2A(2,1) - d2A(1,2)
      gB(3,2) = d2A(2,2) - d2A(1,4)
      gB(3,3) = d2A(2,3) - d2A(1,5)
    case default
      Bvec = 0.0_dp
      gB = 0.0_dp
    end select

    Bmod = sqrt(Bvec(1)**2 + Bvec(2)**2 + Bvec(3)**2)
  end subroutine fo_device_eval_field

  ! One macro-step. Integer select case to a concrete pusher, no procedure
  ! pointers. Both pushers are device-pure.
  subroutine fo_device_step(st, ierr)
    type(fo_device_state_t), intent(inout) :: st
    integer, intent(out) :: ierr

    select case (st%integrator)
    case (FODEV_BORIS)
      call boris_step_dev(st, ierr)
    case (FODEV_FOSYMPL)
      call fosympl_step_dev(st, ierr)
    case default
      ierr = FODEV_ERR_NO_CONVERGE
    end select
  end subroutine fo_device_step

  ! Cartesian Boris, drift-kick-drift, CGS. Pure device kernel.
  subroutine boris_step_dev(st, ierr)
    !$acc routine seq
    type(fo_device_state_t), intent(inout) :: st
    integer, intent(out) :: ierr
    real(dp) :: x(3), v(3), Bvec(3), Bmod, gB(3,3)
    real(dp) :: tvec(3), svec(3), vprime(3), tmag2, dt, qmc

    dt  = st%dt
    qmc = st%charge / (st%mass * c)
    x = st%z(1:3)
    v = st%z(4:6)

    x = x + 0.5_dp * dt * v
    call fo_device_eval_field(st, x, Bvec, Bmod, gB)

    tvec = qmc * Bvec * 0.5_dp * dt
    tmag2 = dot_product(tvec, tvec)
    svec = 2.0_dp * tvec / (1.0_dp + tmag2)
    vprime = v + cross(v, tvec)
    v = v + cross(vprime, svec)

    x = x + 0.5_dp * dt * v

    st%z(1:3) = x
    st%z(4:6) = v
    ierr = FODEV_OK
  end subroutine boris_step_dev

  ! Implicit-midpoint Lorentz full orbit with an ANALYTIC Jacobian. Residual
  !   F(1:3) = xn - xo - dt vmid
  !   F(4:6) = vn - vo - dt (q/(m c)) vmid x B(xmid)
  ! with xmid=(xo+xn)/2, vmid=(vo+vn)/2. Newton on the 6x6 system; the Jacobian
  ! uses grad B(xmid), no finite differences. Structure-preserving (symplectic
  ! for the canonical lift), energy-bounded.
  subroutine fosympl_step_dev(st, ierr)
    !$acc routine seq
    type(fo_device_state_t), intent(inout) :: st
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 50
    real(dp), parameter :: atol = 1.0e-13_dp, rtol = 1.0e-12_dp
    real(dp) :: zold(6), z(6), fvec(6), fjac(6,6), dz(6), reltol(6)
    integer :: kit, i, info
    logical :: conv

    zold = st%z
    z = zold
    ierr = FODEV_OK

    do kit = 1, maxit
      call fosympl_residual_dev(st, zold, z, fvec)
      call fosympl_jacobian_dev(st, zold, z, fjac)
      dz = fvec
      call lu_solve6(fjac, dz, info)
      if (info /= 0) then
        ierr = FODEV_ERR_NO_CONVERGE
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
  end subroutine fosympl_step_dev

  pure subroutine fosympl_residual_dev(st, zold, z, fvec)
    !$acc routine seq
    type(fo_device_state_t), intent(in) :: st
    real(dp), intent(in)  :: zold(6), z(6)
    real(dp), intent(out) :: fvec(6)
    real(dp) :: zmid(6), Bvec(3), Bmod, gB(3,3), qmc, dt

    dt = st%dt
    qmc = st%charge / (st%mass * c)
    zmid = 0.5_dp * (zold + z)
    call fo_device_eval_field(st, zmid(1:3), Bvec, Bmod, gB)
    fvec(1:3) = z(1:3) - zold(1:3) - dt * zmid(4:6)
    fvec(4:6) = z(4:6) - zold(4:6) - dt * qmc * cross(zmid(4:6), Bvec)
  end subroutine fosympl_residual_dev

  ! Analytic Jacobian dF/dz. With zmid=(zold+z)/2 the chain rule gives
  !   dF/dz = I - (dt/2) d(rhs)/dw|_zmid,  rhs = (v, (q/mc) v x B(x)).
  ! d(v x B)/dx_k = v x (dB/dx_k); d(v x B)/dv_k = e_k x B.
  pure subroutine fosympl_jacobian_dev(st, zold, z, fjac)
    !$acc routine seq
    type(fo_device_state_t), intent(in) :: st
    real(dp), intent(in)  :: zold(6), z(6)
    real(dp), intent(out) :: fjac(6,6)
    real(dp) :: zmid(6), Bvec(3), Bmod, gB(3,3), vmid(3), qmc, dt
    real(dp) :: dBk(3), term(3), ek(3), drhs(6,6)
    integer :: i, k, l

    dt = st%dt
    qmc = st%charge / (st%mass * c)
    zmid = 0.5_dp * (zold + z)
    vmid = zmid(4:6)
    call fo_device_eval_field(st, zmid(1:3), Bvec, Bmod, gB)

    drhs = 0.0_dp
    ! dx/dt = v block: d(rhs_x)/dv = I
    do i = 1, 3
      drhs(i, 3+i) = 1.0_dp
    end do
    ! dv/dt = (q/mc) v x B: derivatives wrt x_k and v_k.
    do k = 1, 3
      ! d(v x B)/dx_k = v x (dB/dx_k); dB/dx_k = gB(:,k)
      dBk = gB(:, k)
      term = qmc * cross(vmid, dBk)
      do i = 1, 3
        drhs(3+i, k) = term(i)
      end do
      ! d(v x B)/dv_k = e_k x B
      ek = 0.0_dp; ek(k) = 1.0_dp
      term = qmc * cross(ek, Bvec)
      do i = 1, 3
        drhs(3+i, 3+k) = term(i)
      end do
    end do

    fjac = 0.0_dp
    do l = 1, 6
      fjac(l,l) = 1.0_dp
    end do
    fjac = fjac - 0.5_dp * dt * drhs
  end subroutine fosympl_jacobian_dev

  ! 0.5 m |v|^2 (no electrostatic potential on the device path).
  function fo_device_energy(st) result(energy)
    type(fo_device_state_t), intent(in) :: st
    real(dp) :: energy
    energy = 0.5_dp * st%mass * dot_product(st%z(4:6), st%z(4:6))
  end function fo_device_energy

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

end module orbit_full_device
