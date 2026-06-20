module vmec_field_metric
  ! Single-source, device-callable VMEC metric + field evaluator in NATIVE VMEC
  ! flux coordinates u = (s, theta, varphi). Everything is assembled from ONE
  ! libneo evaluation (splint_vmec_data_d2, issue #322: R,Z map plus 1st and 2nd
  ! derivatives) so the metric and the field share the same g_ij. There is no
  ! separate field source and no class() dispatch, so the core routine is a plain
  ! subroutine marked !$acc routine seq.
  !
  ! Why a single source: the dual-source path (libneo coordinate_system_t metric
  ! plus the separate native VMEC field) gives h_i g^ij h_j = 1.009 because the
  ! two metrics differ. Here h_i = g_ij B^j / |B| with |B| = sqrt(g_ij B^i B^j)
  ! from the SAME g, so h_i g^ij h_j = 1 identically (to round-off).
  !
  ! Field in native VMEC flux coordinates (s, theta, varphi):
  !   A_i      = (0, A_theta(s), A_phi(s))           (flux functions of s)
  !   B^i      = physical field in the VMEC poloidal angle theta, carrying the
  !              lambda stream function (lmns) that converts the symmetry-flux
  !              (straight-field-line) curl-of-flux-A field to the VMEC angle:
  !                B^theta = (-dA_phi_ds - lam_p dA_theta_ds)/sqrtg
  !                B^phi   = (1 + lam_t) dA_theta_ds/sqrtg     (B^s = 0)
  !              Dropping lambda (the earlier curl-of-flux-A-only form) leaves
  !              h_i g^ij h_j = 1 but the WRONG field direction, so the grad-B
  !              and curvature drift are wrong and trapped orbits drift out.
  !   |B|      = sqrt(g_ij B^i B^j)
  ! Metric and metric derivatives:
  !   g_ij     = metric_tensor_vmec(R, dR, dZ)
  !   dg_ij,k  = analytic, from the same R,Z first and second derivatives
  !              (NOT finite difference, NOT the symflux C^T gV C path)
  ! All gradients (dg, dsqrtg, dB, d|B|) are analytic; only the radial second
  ! derivative of A_phi comes from splint_iota (the sA_phi spline second
  ! derivative), the same spline the field uses for dA_phi_ds.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: vmec_field_metric_eval

contains

  ! Single evaluation at u = (s, theta, varphi). Returns the full metric block
  ! (g, ginv, sqrtg, dg), the covariant vector potential and its gradient
  ! (Acov, dA), the contravariant and covariant field (Bctr, Bcov), the field
  ! modulus and its gradient (Bmod, dBmod) and the covariant unit field hcov.
  ! dg(i,j,k) = d g_ij / du_k, dA(i,k) = d A_i / du_k, dBmod(k) = d|B| / du_k.
  !$acc routine seq
  subroutine vmec_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
                                    Bctr, Bcov, Bmod, dBmod, hcov)
    use spline_vmec_sub, only: splint_vmec_data_d2, splint_iota
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp), intent(out) :: Acov(3), dA(3,3)
    real(dp), intent(out) :: Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)

    real(dp) :: s, theta, varphi
    real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, daiota_ds
    real(dp) :: R, Z, alam
    real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
    real(dp) :: dl_ds, dl_dt, dl_dp
    real(dp) :: d2R(6), d2Z(6), d2l(6)
    real(dp) :: dR(3), dZ(3), hR(3,3), hZ(3,3), hl(3,3)
    real(dp) :: dsqrtg(3), d2A_phi_ds2, dBctr(3,3)
    real(dp) :: det, B2, dB2(3), num2, num3, dnum2(3), dnum3(3)
    integer :: i, j, k, idx6(3,3)

    idx6 = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

    s = u(1); theta = u(2); varphi = u(3)

    call splint_vmec_data_d2(s, theta, varphi, A_phi, A_theta, &
         dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, &
         dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
         dl_ds, dl_dt, dl_dp, d2R, d2Z, d2l)

    dR = [dR_ds, dR_dt, dR_dp]
    dZ = [dZ_ds, dZ_dt, dZ_dp]
    do i = 1, 3
      do k = 1, 3
        hR(i,k) = d2R(idx6(i,k))
        hZ(i,k) = d2Z(idx6(i,k))
        hl(i,k) = d2l(idx6(i,k))
      end do
    end do

    ! Native VMEC metric g_ij = e_i . e_j with the (R, phi, Z) embedding.
    g(1,1) = dR(1)**2 + dZ(1)**2
    g(1,2) = dR(1)*dR(2) + dZ(1)*dZ(2)
    g(1,3) = dR(1)*dR(3) + dZ(1)*dZ(3)
    g(2,2) = dR(2)**2 + dZ(2)**2
    g(2,3) = dR(2)*dR(3) + dZ(2)*dZ(3)
    g(3,3) = R**2 + dR(3)**2 + dZ(3)**2
    g(2,1) = g(1,2); g(3,1) = g(1,3); g(3,2) = g(2,3)

    ! Analytic metric derivatives dg_ij,k from the same R,Z first/second
    ! derivatives (no finite difference). g_33 carries the extra R^2 term.
    do k = 1, 3
      dg(1,1,k) = 2.0_dp*(dR(1)*hR(1,k) + dZ(1)*hZ(1,k))
      dg(1,2,k) = hR(1,k)*dR(2) + dR(1)*hR(2,k) + hZ(1,k)*dZ(2) + dZ(1)*hZ(2,k)
      dg(1,3,k) = hR(1,k)*dR(3) + dR(1)*hR(3,k) + hZ(1,k)*dZ(3) + dZ(1)*hZ(3,k)
      dg(2,2,k) = 2.0_dp*(dR(2)*hR(2,k) + dZ(2)*hZ(2,k))
      dg(2,3,k) = hR(2,k)*dR(3) + dR(2)*hR(3,k) + hZ(2,k)*dZ(3) + dZ(2)*hZ(3,k)
      dg(3,3,k) = 2.0_dp*(dR(3)*hR(3,k) + dZ(3)*hZ(3,k)) + 2.0_dp*R*dR(k)
      dg(2,1,k) = dg(1,2,k); dg(3,1,k) = dg(1,3,k); dg(3,2,k) = dg(2,3,k)
    end do

    ! Inverse metric by cofactors.
    det = g(1,1)*(g(2,2)*g(3,3) - g(2,3)*g(3,2)) &
        - g(1,2)*(g(2,1)*g(3,3) - g(2,3)*g(3,1)) &
        + g(1,3)*(g(2,1)*g(3,2) - g(2,2)*g(3,1))
    ginv(1,1) = (g(2,2)*g(3,3) - g(2,3)*g(3,2))/det
    ginv(1,2) = (g(1,3)*g(3,2) - g(1,2)*g(3,3))/det
    ginv(1,3) = (g(1,2)*g(2,3) - g(1,3)*g(2,2))/det
    ginv(2,1) = (g(2,3)*g(3,1) - g(2,1)*g(3,3))/det
    ginv(2,2) = (g(1,1)*g(3,3) - g(1,3)*g(3,1))/det
    ginv(2,3) = (g(1,3)*g(2,1) - g(1,1)*g(2,3))/det
    ginv(3,1) = (g(2,1)*g(3,2) - g(2,2)*g(3,1))/det
    ginv(3,2) = (g(1,2)*g(3,1) - g(1,1)*g(3,2))/det
    ginv(3,3) = (g(1,1)*g(2,2) - g(1,2)*g(2,1))/det

    ! Native VMEC Jacobian sqrtg = R (dR_dt dZ_ds - dR_ds dZ_dt) and its gradient.
    sqrtg = R*(dR_dt*dZ_ds - dR_ds*dZ_dt)
    do k = 1, 3
      dsqrtg(k) = dR(k)*(dR_dt*dZ_ds - dR_ds*dZ_dt) &
                + R*(hR(2,k)*dZ_ds + dR_dt*hZ(1,k) - hR(1,k)*dZ_dt - dR_ds*hZ(2,k))
    end do

    ! Covariant vector potential (flux functions of s) and its gradient.
    Acov = [0.0_dp, A_theta, A_phi]
    dA = 0.0_dp
    dA(2,1) = dA_theta_ds
    dA(3,1) = dA_phi_ds

    ! Contravariant field in NATIVE VMEC angles. The physical field is the simple
    ! curl-of-flux-A form only in the symmetry-flux (straight-field-line) angle
    ! vartheta = theta + lambda(s,theta,phi); transformed to the VMEC poloidal angle
    ! theta it carries the lambda stream function (lmns), without which the field
    ! direction (and hence the grad-B/curvature drift) is wrong even though
    ! |B| stays a unit-h field. With sqrtg the VMEC Jacobian (B^s = 0):
    !   B^theta = (-dA_phi_ds - lam_p dA_theta_ds)/sqrtg
    !   B^phi   = (1 + lam_t) dA_theta_ds/sqrtg
    ! (lam_t = dl_dt, lam_p = dl_dp). dA_theta_ds = torflux is constant in s.
    num2 = -dA_phi_ds - dl_dp*dA_theta_ds
    num3 = (1.0_dp + dl_dt)*dA_theta_ds
    Bctr(1) = 0.0_dp
    Bctr(2) = num2/sqrtg
    Bctr(3) = num3/sqrtg

    ! Covariant field B_i = g_ij B^j (same metric).
    do i = 1, 3
      Bcov(i) = g(i,2)*Bctr(2) + g(i,3)*Bctr(3)
    end do

    ! |B| = sqrt(g_ij B^i B^j) = sqrt(B^i B_i) from the SAME g.
    B2 = Bctr(2)*Bcov(2) + Bctr(3)*Bcov(3)
    Bmod = sqrt(B2)

    ! Gradient of B^i = num/sqrtg. d(num)/du_k uses d2A_phi_ds2 (radial, via
    ! splint_iota) and the lambda second derivatives hl(.,.) (dA_theta_ds is
    ! constant so d(dA_theta_ds)=0).
    call splint_iota(s, aiota, daiota_ds)
    d2A_phi_ds2 = -daiota_ds*dA_theta_ds   ! aiota = -dA_phi_ds/torflux
    do k = 1, 3
      dnum2(k) = -hl(3,k)*dA_theta_ds      ! d(-dA_phi_ds - lam_p dA_theta_ds)
      dnum3(k) = hl(2,k)*dA_theta_ds       ! d((1+lam_t) dA_theta_ds)
    end do
    dnum2(1) = dnum2(1) - d2A_phi_ds2      ! radial part of -dA_phi_ds
    dBctr = 0.0_dp
    do k = 1, 3
      dBctr(2,k) = dnum2(k)/sqrtg - num2*dsqrtg(k)/sqrtg**2
      dBctr(3,k) = dnum3(k)/sqrtg - num3*dsqrtg(k)/sqrtg**2
    end do

    ! d(|B|^2)/du_k = dg_ij,k B^i B^j + 2 g_ij B^i dB^j/du_k, then chain to |B|.
    do k = 1, 3
      dB2(k) = 0.0_dp
      do i = 2, 3
        do j = 2, 3
          dB2(k) = dB2(k) + dg(i,j,k)*Bctr(i)*Bctr(j)
        end do
      end do
      dB2(k) = dB2(k) + 2.0_dp*(Bcov(2)*dBctr(2,k) + Bcov(3)*dBctr(3,k))
      dBmod(k) = 0.5_dp*dB2(k)/Bmod
    end do

    ! Covariant unit field h_i = B_i/|B|; h_i g^ij h_j = 1 by construction.
    do i = 1, 3
      hcov(i) = Bcov(i)/Bmod
    end do
  end subroutine vmec_field_metric_eval

end module vmec_field_metric
