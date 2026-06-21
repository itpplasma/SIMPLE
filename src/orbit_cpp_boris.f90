module orbit_cpp_boris
  ! Experimental large-step 6D classical Pauli particle by a Boris-type pusher
  ! (Xiao-Qin BAP2), an EXPLICIT, structure-preserving alternative to the implicit
  ! midpoint (orbit_cpp_canonical), which has no root at trapped bounces at large
  ! dt (issue #417). Boris has no nonlinear solve, hence no convergence floor: the
  ! step is exact rotation + force arithmetic. It is the same physics as CPP --
  ! H = |p - qc A|^2/2m + mu|B|, seeded with v_perp = 0 on the slow manifold.
  !
  ! The particle is advanced in Cartesian (x, v), where the magnetic rotation is
  ! exact for constant B over a step. The field is the production Boozer field:
  ! at the Cartesian point we invert to Boozer (cart_to_boozer), evaluate
  ! boozer_field_metric (contravariant B^i, |B|, d|B|/du_i), and push the physical
  ! vectors to Cartesian with the chart Jacobian Jc = d(xyz)/du:
  !   B_cart = Jc B^ctr,   grad|B|_cart = Jc^{-T} d|B|/du.
  ! The Pauli mirror force enters as the "electric" half-kick -mu grad|B|/m; the
  ! full charged particle (MODEL_CP) drops it. Energy H and the parameter mu are
  ! the validation invariants. set filtered=.true. for the Hairer-Lubich-Wang
  ! large-step filter on the rotation angle (keeps the modified mu bounded when
  ! dt*Omega_c >> 1).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use boozer_cartesian, only: boozer_to_cart, cart_to_boozer
  use boozer_field_metric, only: boozer_field_metric_eval
  implicit none
  private

  real(dp), parameter :: c = 1.0_dp

  public :: cpp_boris_state_t, cpp_boris_init, cpp_boris_step, cpp_boris_energy, &
            cpp_boris_to_gc

  type :: cpp_boris_state_t
    real(dp) :: x(3)   = 0.0_dp    ! Cartesian position (cm)
    real(dp) :: v(3)   = 0.0_dp    ! Cartesian velocity (normalized)
    real(dp) :: u(3)   = 0.0_dp    ! carried Boozer (s, vth, vph) = cart_to_boozer guess
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

  ! Cartesian B vector, |B|, and grad|B| at Cartesian x, from the Boozer field.
  ! u_guess seeds the cart->Boozer inversion and is updated to the found u.
  subroutine cart_field(x, u_guess, Bvec, Bmod, gradB, u_out, ierr)
    real(dp), intent(in) :: x(3), u_guess(3)
    real(dp), intent(out) :: Bvec(3), Bmod, gradB(3), u_out(3)
    integer, intent(out) :: ierr
    real(dp) :: u(3), xyz(3), Jc(3,3), Jinv(3,3)
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3), Acov(3), dA(3,3)
    real(dp) :: Bctr(3), Bcov(3), dBmod(3), hcov(3)
    integer :: i

    call cart_to_boozer(x, u_guess, u, ierr)
    if (ierr /= 0) return
    u_out = u
    call boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
         Bctr, Bcov, Bmod, dBmod, hcov)
    call boozer_to_cart(u, xyz, Jc)
    ! B_cart = Jc B^ctr (contravariant field pushed to Cartesian).
    do i = 1, 3
      Bvec(i) = Jc(i,1)*Bctr(1) + Jc(i,2)*Bctr(2) + Jc(i,3)*Bctr(3)
    end do
    ! grad|B|_cart = Jc^{-T} d|B|/du.
    call inv3(Jc, Jinv)
    do i = 1, 3
      gradB(i) = Jinv(1,i)*dBmod(1) + Jinv(2,i)*dBmod(2) + Jinv(3,i)*dBmod(3)
    end do
  end subroutine cart_field

  ! One Boris-Pauli macro-step in Cartesian: half drift, half mirror kick, exact
  ! magnetic rotation, half mirror kick, half drift. ierr/=0 on field-inversion
  ! failure (treated as a lost/aborted orbit by the caller).
  subroutine cpp_boris_step(st, ierr)
    type(cpp_boris_state_t), intent(inout) :: st
    integer, intent(out) :: ierr
    real(dp) :: x(3), v(3), Bvec(3), Bmod, gradB(3), u(3)
    real(dp) :: tvec(3), svec(3), vp(3), tmag2, qcm, fac
    integer :: i

    x = st%x
    v = st%v
    qcm = st%charge/(c*st%ro0*st%mass)   ! rotation: dv/dt = qcm v x B

    x = x + 0.5_dp*st%dt*v
    call cart_field(x, st%u, Bvec, Bmod, gradB, u, ierr)
    if (ierr /= 0) return
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
    ierr = 0
  end subroutine cpp_boris_step

  ! Seed from a guiding-centre start record (s,th,ph,vpar) with v_perp = 0 (CPP
  ! slow manifold) and the fixed mu. Cartesian position from boozer_to_cart; the
  ! parallel velocity along the Cartesian field direction.
  subroutine cpp_boris_init(st, pauli, x0_boozer, vpar0, mu_in, mass, charge, &
                            dt, ro0_in, pabs, filtered)
    type(cpp_boris_state_t), intent(out) :: st
    logical, intent(in) :: pauli
    real(dp), intent(in) :: x0_boozer(3), vpar0, mu_in, mass, charge, dt, ro0_in, pabs
    logical, intent(in), optional :: filtered
    real(dp) :: xyz(3), Jc(3,3), Bvec(3), Bmod, gradB(3), u(3)
    integer :: ierr

    st%pauli = pauli
    st%mass = mass; st%charge = charge; st%dt = dt; st%ro0 = ro0_in
    st%mu = mu_in; st%pabs = pabs
    if (present(filtered)) st%filtered = filtered
    st%u = x0_boozer
    call boozer_to_cart(x0_boozer, xyz, Jc)
    st%x = xyz
    call cart_field(xyz, x0_boozer, Bvec, Bmod, gradB, u, ierr)
    st%u = u
    ! v = vpar0 * b_hat (parallel only; v_perp = 0 on the slow manifold).
    st%v = vpar0*Bvec/max(Bmod_of(Bvec), 1.0e-30_dp)
  end subroutine cpp_boris_init

  function cpp_boris_energy(st) result(energy)
    type(cpp_boris_state_t), intent(in) :: st
    real(dp) :: energy, Bvec(3), Bmod, gradB(3), u(3)
    integer :: ierr
    call cart_field(st%x, st%u, Bvec, Bmod, gradB, u, ierr)
    energy = 0.5_dp*st%mass*(st%v(1)**2 + st%v(2)**2 + st%v(3)**2)
    if (st%pauli) energy = energy + st%mu*Bmod
  end function cpp_boris_energy

  ! Guiding-centre reduction for output: Boozer (s,th,ph) of the current point and
  ! the parallel speed lambda = vpar/|v|.
  subroutine cpp_boris_to_gc(st, s, th, ph, vpar, ierr)
    type(cpp_boris_state_t), intent(inout) :: st
    real(dp), intent(out) :: s, th, ph, vpar
    integer, intent(out) :: ierr
    real(dp) :: u(3), Bvec(3), Bmod, gradB(3), uf(3), vmag
    call cart_to_boozer(st%x, st%u, u, ierr)
    if (ierr /= 0) return
    st%u = u
    s = u(1); th = u(2); ph = u(3)
    call cart_field(st%x, u, Bvec, Bmod, gradB, uf, ierr)
    vmag = max(Bmod_of(Bvec), 1.0e-30_dp)
    vpar = (st%v(1)*Bvec(1) + st%v(2)*Bvec(2) + st%v(3)*Bvec(3))/vmag
  end subroutine cpp_boris_to_gc

  pure function cross(a, b) result(cr)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: cr(3)
    cr(1) = a(2)*b(3) - a(3)*b(2)
    cr(2) = a(3)*b(1) - a(1)*b(3)
    cr(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  pure function Bmod_of(B) result(m)
    real(dp), intent(in) :: B(3)
    real(dp) :: m
    m = sqrt(B(1)**2 + B(2)**2 + B(3)**2)
  end function Bmod_of

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
