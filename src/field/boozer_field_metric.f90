module boozer_field_metric
  ! Single-source, device-callable Boozer metric + field evaluator in NATIVE
  ! Boozer flux coordinates u = (s, vartheta_B, varphi_B). This is the Boozer
  ! analogue of vmec_field_metric: it returns the same block (g, ginv, sqrtg,
  ! dg, Acov, dA, Bctr, Bcov, Bmod, dBmod, hcov) with an IDENTICAL signature so
  ! the 6D CP/CPP integrators can switch coordinate charts without changing the
  ! call site. No class() dispatch, fixed-size arrays, !$acc routine seq.
  !
  ! Geometry: the VMEC angles are obtained from the Boozer angle map deltas,
  !   theta_V   = vartheta_B - deltheta_BV(s, vartheta_B, varphi_B)
  !   varphi_V  = varphi_B   - delphi_BV (s, vartheta_B, varphi_B)
  ! taken from the SAME B-side spline that supplies the Jacobian (delthe_delphi_BV_d2),
  ! so the geometry point and the Jacobian belong to one consistent angle map (a
  ! separate boozer_to_vmec Newton on the V-side spline would disagree at the
  ! grid-error level and spoil dg). The metric g_ij is built from
  ! the R,Z map and its first/second derivatives at the VMEC angles (the same
  ! splint_vmec_data_d2 that vmec_field_metric uses), giving the metric in the
  ! VMEC-ANGLE chart g_V and its analytic gradient dg_V w.r.t. (s, theta_V,
  ! varphi_V). It is then pulled back to the Boozer chart with the angle
  ! Jacobian J = d(s, theta_V, varphi_V)/d(s, vartheta_B, varphi_B):
  !   g_B   = J^T g_V J
  !   dg_B  = analytic, using dJm (the SECOND derivatives of the angle map,
  !           from delthe_delphi_BV_d2) and dg_V chained through J (NOT finite
  !           difference).
  !
  ! Field: taken DIRECTLY from the production Boozer splines (splint_boozer_coord)
  ! so the 6D field equals the guiding-centre Boozer field bit-for-bit:
  !   Acov = (0, A_theta = torflux*s, A_phi(s))
  !   Bcov = (B_r, B_vartheta_B = I(s), B_varphi_B = g(s))
  !   Bmod = splined |B|,   dBmod = splined gradient w.r.t. (s, vartheta_B, varphi_B)
  ! Bctr is raised with the Boozer ginv, hcov = Bcov/Bmod. With field and metric
  ! both in the Boozer chart, h_i g^ij h_j = 1 is the consistency gate between
  ! the transformed metric and the splined |B|.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: boozer_field_metric_eval

contains

  !$acc routine seq
  subroutine boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
                                      Bctr, Bcov, Bmod, dBmod, hcov)
    use spline_vmec_sub, only: splint_vmec_data_d2
    use boozer_sub, only: delthe_delphi_BV_d2, splint_boozer_coord
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp), intent(out) :: Acov(3), dA(3,3)
    real(dp), intent(out) :: Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)

    real(dp) :: s, vartheta_B, varphi_B, theta_V, varphi_V
    integer :: idx6(3,3), i, j, k, l, m

    ! Angle map deltas, first and second derivatives w.r.t. (s, vartheta_B, varphi_B)
    real(dp) :: del_t, del_p, ddel_t(3), ddel_p(3), d2del_t(6), d2del_p(6)

    ! VMEC-chart geometry from splint_vmec_data_d2
    real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
    real(dp) :: R, Zc, alam
    real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
    real(dp) :: dl_ds, dl_dt, dl_dp
    real(dp) :: d2R(6), d2Z(6), d2l(6)
    real(dp) :: dR(3), dZ(3), hR(3,3), hZ(3,3)

    ! Metric and its gradient in the VMEC-angle chart
    real(dp) :: gV(3,3), dgV(3,3,3)

    ! Angle Jacobian J = d(s,theta_V,varphi_V)/d(s,vartheta_B,varphi_B) and its
    ! gradient dJm(i,j,k) = d Jm(i,j) / d u_k (u the Boozer coordinate).
    real(dp) :: Jm(3,3), dJm(3,3,3)

    ! Boozer-chart field from production splines
    real(dp) :: A_theta_B, A_phi_B, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3
    real(dp) :: B_vartheta_B, dB_vartheta_B, d2B_vartheta_B
    real(dp) :: B_varphi_B, dB_varphi_B, d2B_varphi_B
    real(dp) :: Bmod_B, dBmod_B(3), d2Bmod_B(6)
    real(dp) :: B_r, dB_r(3), d2B_r(6)

    real(dp) :: det, tmp, dgVtot(3,3,3)

    idx6 = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

    s = u(1); vartheta_B = u(2); varphi_B = u(3)

    ! Angle map deltas and their first/second derivatives in Boozer coordinates,
    ! all from the SAME B-side spline. theta_V, varphi_V are derived from these
    ! deltas (NOT from a separate boozer_to_vmec Newton on the V-side spline) so
    ! the geometry point and the Jacobian belong to one consistent angle map.
    call delthe_delphi_BV_d2(s, vartheta_B, varphi_B, del_t, del_p, &
                             ddel_t, ddel_p, d2del_t, d2del_p)
    theta_V = vartheta_B - del_t
    varphi_V = varphi_B - del_p

    ! VMEC-chart geometry (R,Z map plus 1st and 2nd derivatives) at the VMEC angles.
    call splint_vmec_data_d2(s, theta_V, varphi_V, A_phi, A_theta, &
         dA_phi_ds, dA_theta_ds, aiota, R, Zc, alam, &
         dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
         dl_ds, dl_dt, dl_dp, d2R, d2Z, d2l)

    dR = [dR_ds, dR_dt, dR_dp]
    dZ = [dZ_ds, dZ_dt, dZ_dp]
    do i = 1, 3
      do k = 1, 3
        hR(i,k) = d2R(idx6(i,k))
        hZ(i,k) = d2Z(idx6(i,k))
      end do
    end do

    ! Metric in the VMEC-angle chart (R, phi, Z embedding), g_33 carries R^2.
    gV(1,1) = dR(1)**2 + dZ(1)**2
    gV(1,2) = dR(1)*dR(2) + dZ(1)*dZ(2)
    gV(1,3) = dR(1)*dR(3) + dZ(1)*dZ(3)
    gV(2,2) = dR(2)**2 + dZ(2)**2
    gV(2,3) = dR(2)*dR(3) + dZ(2)*dZ(3)
    gV(3,3) = R**2 + dR(3)**2 + dZ(3)**2
    gV(2,1) = gV(1,2); gV(3,1) = gV(1,3); gV(3,2) = gV(2,3)

    ! Analytic gradient of g_V w.r.t. the VMEC coordinates (s, theta_V, varphi_V).
    do k = 1, 3
      dgV(1,1,k) = 2.0_dp*(dR(1)*hR(1,k) + dZ(1)*hZ(1,k))
      dgV(1,2,k) = hR(1,k)*dR(2) + dR(1)*hR(2,k) + hZ(1,k)*dZ(2) + dZ(1)*hZ(2,k)
      dgV(1,3,k) = hR(1,k)*dR(3) + dR(1)*hR(3,k) + hZ(1,k)*dZ(3) + dZ(1)*hZ(3,k)
      dgV(2,2,k) = 2.0_dp*(dR(2)*hR(2,k) + dZ(2)*hZ(2,k))
      dgV(2,3,k) = hR(2,k)*dR(3) + dR(2)*hR(3,k) + hZ(2,k)*dZ(3) + dZ(2)*hZ(3,k)
      dgV(3,3,k) = 2.0_dp*(dR(3)*hR(3,k) + dZ(3)*hZ(3,k)) + 2.0_dp*R*dR(k)
      dgV(2,1,k) = dgV(1,2,k); dgV(3,1,k) = dgV(1,3,k); dgV(3,2,k) = dgV(2,3,k)
    end do

    ! Angle Jacobian Jm(a,b) = d x_V^a / d u_B^b, with
    !   x_V = (s, theta_V, varphi_V),  u_B = (s, vartheta_B, varphi_B).
    ! theta_V  = vartheta_B - del_t,  varphi_V = varphi_B - del_p, so
    !   Jm(1,:) = (1, 0, 0)
    !   Jm(2,:) = (-ddel_t(1), 1 - ddel_t(2),  -ddel_t(3))
    !   Jm(3,:) = (-ddel_p(1),  -ddel_p(2),  1 - ddel_p(3))
    Jm = 0.0_dp
    Jm(1,1) = 1.0_dp
    Jm(2,1) = -ddel_t(1); Jm(2,2) = 1.0_dp - ddel_t(2); Jm(2,3) = -ddel_t(3)
    Jm(3,1) = -ddel_p(1); Jm(3,2) = -ddel_p(2);          Jm(3,3) = 1.0_dp - ddel_p(3)

    ! dJm(a,b,k) = d Jm(a,b) / d u_B^k. Row 1 is constant -> 0. Rows 2,3 are minus
    ! the second derivatives of the angle map (packed idx6 over (s,t,p)).
    dJm = 0.0_dp
    do k = 1, 3
      dJm(2,1,k) = -d2del_t(idx6(1,k))
      dJm(2,2,k) = -d2del_t(idx6(2,k))
      dJm(2,3,k) = -d2del_t(idx6(3,k))
      dJm(3,1,k) = -d2del_p(idx6(1,k))
      dJm(3,2,k) = -d2del_p(idx6(2,k))
      dJm(3,3,k) = -d2del_p(idx6(3,k))
    end do

    ! Pull back the metric: g_B(i,j) = sum_{a,b} Jm(a,i) g_V(a,b) Jm(b,j).
    do i = 1, 3
      do j = 1, 3
        tmp = 0.0_dp
        do l = 1, 3
          do m = 1, 3
            tmp = tmp + Jm(l,i)*gV(l,m)*Jm(m,j)
          end do
        end do
        g(i,j) = tmp
      end do
    end do

    ! Total derivative of g_V along the Boozer coordinates: g_V depends on u_B
    ! both explicitly (none) and through x_V(u_B), so
    !   d g_V(a,b)/d u_B^k = sum_c dgV(a,b,c) * Jm(c,k).
    do k = 1, 3
      do i = 1, 3
        do j = 1, 3
          tmp = 0.0_dp
          do l = 1, 3
            tmp = tmp + dgV(i,j,l)*Jm(l,k)
          end do
          dgVtot(i,j,k) = tmp
        end do
      end do
    end do

    ! Gradient of the pulled-back metric (product rule on J^T g_V J):
    !   dg_B(i,j,k) = sum_{a,b} [ dJm(a,i,k) g_V(a,b) Jm(b,j)
    !                            + Jm(a,i) dgVtot(a,b,k) Jm(b,j)
    !                            + Jm(a,i) g_V(a,b) dJm(b,j,k) ].
    do k = 1, 3
      do i = 1, 3
        do j = 1, 3
          tmp = 0.0_dp
          do l = 1, 3
            do m = 1, 3
              tmp = tmp + dJm(l,i,k)*gV(l,m)*Jm(m,j) &
                        + Jm(l,i)*dgVtot(l,m,k)*Jm(m,j) &
                        + Jm(l,i)*gV(l,m)*dJm(m,j,k)
            end do
          end do
          dg(i,j,k) = tmp
        end do
      end do
    end do

    ! Inverse Boozer metric by cofactors.
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

    ! Boozer Jacobian sqrt(g) = sqrt(det g_B).
    sqrtg = sqrt(det)

    ! Field directly from the production Boozer splines (mode_secders=0 is enough;
    ! all needed first derivatives are returned). Abscissa r = s.
    call splint_boozer_coord(s, vartheta_B, varphi_B, 0, &
                             A_theta_B, A_phi_B, dA_theta_dr, dA_phi_dr, &
                             d2A_phi_dr2, d3A_phi_dr3, &
                             B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                             B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                             Bmod_B, dBmod_B, d2Bmod_B, &
                             B_r, dB_r, d2B_r)

    ! Covariant vector potential (flux functions of s) and its gradient.
    Acov = [0.0_dp, A_theta_B, A_phi_B]
    dA = 0.0_dp
    dA(2,1) = dA_theta_dr
    dA(3,1) = dA_phi_dr

    ! Covariant Boozer field B_i = (B_s, I(s), g(s)). B_s = B_r (radial covariant).
    Bcov(1) = B_r
    Bcov(2) = B_vartheta_B
    Bcov(3) = B_varphi_B

    ! |B| and its gradient straight from the spline.
    Bmod = Bmod_B
    dBmod = dBmod_B

    ! Contravariant field by raising with the Boozer metric.
    do i = 1, 3
      Bctr(i) = ginv(i,1)*Bcov(1) + ginv(i,2)*Bcov(2) + ginv(i,3)*Bcov(3)
    end do

    ! Covariant unit field h_i = B_i/|B|.
    do i = 1, 3
      hcov(i) = Bcov(i)/Bmod
    end do
  end subroutine boozer_field_metric_eval

end module boozer_field_metric
