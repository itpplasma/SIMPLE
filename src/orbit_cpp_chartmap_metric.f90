module orbit_cpp_chartmap_metric
  ! Host-side metric + field provider for the genuine 6D canonical-midpoint
  ! integrator on the PRODUCTION Boozer/chartmap chart. This is the chart whose
  ! metric matches the production field_can chart (libneo #322): field_can_boozer
  ! integrates in (s, theta_B, phi_B) with the chartmap reference coordinate
  ! rho = sqrt(s), SAME angles. The 6D state runs in u = (rho, theta_B, phi_B) so
  ! the libneo chartmap metric/Christoffel (from reference_coordinates%ref_coords)
  ! is native; the production field_can (Boozer) field is reparametrized from
  ! s = rho^2 with the radial chain rule dF/drho = 2 rho dF/ds.
  !
  ! ref_coords is the scaled chartmap coordinate system built by
  ! init_reference_coordinates; both chartmap_coordinate_system_t and its scaled
  ! extension provide metric_tensor/christoffel, so the base coordinate_system_t
  ! interface dispatches to the active chart.
  !
  ! NOT GPU-portable: libneo metric_tensor/christoffel are class()-dispatched and
  ! read 3D splines; field_can%evaluate is a procedure pointer over spline reads.
  ! Metric derivatives come from Christoffel via metric compatibility:
  !   dg_ij/du_k = g_il Gamma^l_jk + g_jl Gamma^l_ik.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: chartmap_metric_active, chartmap_eval_metric, chartmap_eval_field

contains

  ! True only when the production reference chart is a chartmap coordinate system
  ! (chartmap_coordinate_system_t or its scaled extension). The 6D Boozer/chartmap
  ! CPP path requires this; the generic-BOOZER-on-VMEC chart has no matching metric.
  logical function chartmap_metric_active()
    use reference_coordinates, only: ref_coords
    use libneo_coordinates, only: chartmap_coordinate_system_t

    chartmap_metric_active = .false.
    if (.not. allocated(ref_coords)) return
    select type (ref_coords)
    class is (chartmap_coordinate_system_t)
      chartmap_metric_active = .true.
    end select
  end function chartmap_metric_active

  ! Full metric block at u=(rho,theta_B,phi_B): g_ij, g^ij, and
  ! dg(i,j,k)=dg_ij/du_k from Christoffel + metric compatibility. Native chartmap
  ! coordinates, no s=rho^2 reparametrization (the metric IS in rho).
  subroutine chartmap_eval_metric(u, g, ginv, dg)
    use reference_coordinates, only: ref_coords
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: g(3,3), ginv(3,3), dg(3,3,3)
    real(dp) :: sqrtg, Gamma(3,3,3)
    integer :: i, j, k, l

    call ref_coords%metric_tensor(u, g, ginv, sqrtg)
    call ref_coords%christoffel(u, Gamma)
    do k = 1, 3
      do j = 1, 3
        do i = 1, 3
          dg(i,j,k) = 0.0_dp
          do l = 1, 3
            dg(i,j,k) = dg(i,j,k) + g(i,l)*Gamma(l,j,k) + g(j,l)*Gamma(l,i,k)
          end do
        end do
      end do
    end do
  end subroutine chartmap_eval_metric

  ! Production Boozer field at u=(rho,theta_B,phi_B), reparametrized from
  ! field_can(s=rho^2). field_can returns covariant A_theta,A_phi (A_s=0),
  ! covariant h_theta,h_phi (h_s=0), |B| and the s-derivatives. In the chartmap
  ! chart only the radial coordinate differs: F(rho)=F(s(rho)), dF/drho=2 rho dF/ds.
  ! Angular derivatives are unchanged. dA(i,k)=dA_i/du_k carries the same chain rule
  ! on its radial column.
  subroutine chartmap_eval_field(u, Acov, dA, Bmod, dBmod, hcov)
    use field_can_mod, only: eval_field => evaluate, field_can_t
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: Acov(3), dA(3,3), Bmod, dBmod(3), hcov(3)
    type(field_can_t) :: f
    real(dp) :: rho, drho_ds

    rho = u(1)
    drho_ds = 2.0_dp*rho   ! ds/drho = 2 rho (s = rho^2)

    call eval_field(f, rho*rho, u(2), u(3), 0)

    Acov = [0.0_dp, f%Ath, f%Aph]
    hcov = [0.0_dp, f%hth, f%hph]
    Bmod = f%Bmod

    ! dA(i,k) = dA_i/du_k. A_s = 0 (row 1 all zero). Rows 2,3 (A_theta, A_phi):
    ! radial column k=1 scales by ds/drho; angular columns unchanged.
    dA = 0.0_dp
    dA(2,1) = f%dAth(1)*drho_ds
    dA(2,2) = f%dAth(2)
    dA(2,3) = f%dAth(3)
    dA(3,1) = f%dAph(1)*drho_ds
    dA(3,2) = f%dAph(2)
    dA(3,3) = f%dAph(3)

    ! d|B|/du_k: radial column scales, angular columns unchanged.
    dBmod(1) = f%dBmod(1)*drho_ds
    dBmod(2) = f%dBmod(2)
    dBmod(3) = f%dBmod(3)
  end subroutine chartmap_eval_field

end module orbit_cpp_chartmap_metric
