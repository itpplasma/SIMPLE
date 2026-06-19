module pauli_gc_drift
  ! Independent guiding-centre drift model in the SAME Cartesian normalized units
  ! on the SAME analytic field (field_pauli_cart) as the 6D Pauli particle. This
  ! is the asymptotic GC slow manifold:
  !   dX/dt   = vpar bhat + (m c)/(q |B|) [ vpar^2 (bhat x kappa)
  !                                       + (vperp^2/2)(bhat x grad|B|)/|B| ]
  !   dvpar/dt = -(mu/m) bhat . grad|B|,   mu = m vperp^2/(2|B|) fixed,
  ! with curvature kappa = (bhat . grad) bhat. It is a DIFFERENT method from the
  ! 6D Pauli (drift-reduced vs full gyration), so their banana orbits agree only
  ! to O(rho*) -- a genuine cross-method check, not a tautology. Field gradients
  ! and grad|B| come from the analytic A and its second derivatives, so this GC
  ! oracle and the Pauli particle see literally the same field.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: c
  use field_pauli_cart, only: pauli_field_params_t, eval_pauli_field_cart
  implicit none
  private
  public :: gc_drift_rhs
  ! Symmetric pair index (j,k) -> packed second-derivative slot.
  integer, parameter :: MIDX(3,3) = reshape([1,2,3, 2,4,5, 3,5,6], [3,3])

contains

  subroutine gc_drift_rhs(fp, mass, charge, vperp, mu, Y, dY)
    type(pauli_field_params_t), intent(in) :: fp
    real(dp), intent(in)  :: mass, charge, vperp, mu, Y(4)
    real(dp), intent(out) :: dY(4)
    real(dp) :: Av(3), dA(3,3), d2A(3,6), Bv(3), Bm, dBm(3), d2Bm(6)
    real(dp) :: bhat(3), gB(3,3), gradb(3,3), kappa(3), drift(3), vp, coef
    integer :: i, j

    call eval_pauli_field_cart(fp, Y(1:3), Av, dA, d2A, Bv, Bm, dBm, d2Bm)
    vp = Y(4)
    bhat = Bv / Bm
    ! grad B_i/dx_j from B = curl A, i.e. second derivatives of A.
    do j = 1, 3
      gB(1,j) = d2A(3,MIDX(2,j)) - d2A(2,MIDX(3,j))
      gB(2,j) = d2A(1,MIDX(3,j)) - d2A(3,MIDX(1,j))
      gB(3,j) = d2A(2,MIDX(1,j)) - d2A(1,MIDX(2,j))
    end do
    do i = 1, 3
      do j = 1, 3
        gradb(i,j) = gB(i,j)/Bm - Bv(i)*dBm(j)/Bm**2
      end do
    end do
    do i = 1, 3
      kappa(i) = 0.0_dp
      do j = 1, 3
        kappa(i) = kappa(i) + bhat(j)*gradb(i,j)
      end do
    end do
    coef = mass*c/(charge*Bm)
    drift = coef*( vp*vp*cross(bhat,kappa) + (vperp*vperp/2.0_dp)*cross(bhat,dBm)/Bm )
    dY(1:3) = vp*bhat + drift
    dY(4)   = -(mu/mass)*dot_product(bhat, dBm)
  end subroutine gc_drift_rhs

  pure function cross(u, w) result(z)
    real(dp), intent(in) :: u(3), w(3)
    real(dp) :: z(3)
    z(1) = u(2)*w(3) - u(3)*w(2)
    z(2) = u(3)*w(1) - u(1)*w(3)
    z(3) = u(1)*w(2) - u(2)*w(1)
  end function cross

end module pauli_gc_drift
