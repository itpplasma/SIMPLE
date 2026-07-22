module orbit_fo_field
  ! Field provider for the full-orbit (FO) Boris pusher on the production
  ! Boozer/chartmap chart. The chartmap reference coordinate is rho = sqrt(s)
  ! with the same angles as the production field_can_boozer chart, so the libneo
  ! chartmap geometry (reference_coordinates%ref_coords) is native. The Boozer
  ! flux potential A_theta(s), A_phi(s), |B| and the field |B|-derivatives come
  ! from field_can, reparametrized from s = rho^2 by the radial chain rule
  ! dF/drho = 2 rho dF/ds; the angular derivatives are unchanged.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: fo_eval_field, fo_eval_reference_field, &
    fo_reference_field_available

contains

  logical function fo_reference_field_available()
    use magfie_sub, only: has_magfie_refcoords_field

    fo_reference_field_available = has_magfie_refcoords_field()
  end function fo_reference_field_available

  subroutine fo_eval_reference_field(u, hcov, Bmod, status)
    use magfie_sub, only: evaluate_magfie_refcoords_field
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: hcov(3), Bmod
    integer, intent(out) :: status
    real(dp) :: Acov(3)

    call evaluate_magfie_refcoords_field(u, Acov, hcov, Bmod, status)
  end subroutine fo_eval_reference_field

  ! Production Boozer field at u=(rho,theta_B,phi_B), reparametrized from
  ! field_can(s=rho^2). field_can returns covariant A_theta,A_phi (A_s=0),
  ! covariant h_theta,h_phi (h_s=0), |B| and the s-derivatives. Only the radial
  ! coordinate differs: F(rho)=F(s(rho)), dF/drho = 2 rho dF/ds; angular
  ! derivatives are unchanged. dA(i,k)=dA_i/du_k carries the same chain rule on
  ! its radial column.
  subroutine fo_eval_field(u, Acov, dA, Bmod, dBmod, hcov)
    use field_can_mod, only: eval_field => evaluate, field_can_t
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: Acov(3), dA(3,3), Bmod, dBmod(3), hcov(3)
    type(field_can_t) :: f
    real(dp) :: rho, ds_drho

    rho = u(1)
    ds_drho = 2.0_dp*rho   ! ds/drho = 2 rho (s = rho^2); dF/drho = ds/drho dF/ds

    call eval_field(f, rho*rho, u(2), u(3), 0)

    Acov = [0.0_dp, f%Ath, f%Aph]
    hcov = [0.0_dp, f%hth, f%hph]
    Bmod = f%Bmod

    ! dA(i,k) = dA_i/du_k. A_s = 0 (row 1 all zero). Rows 2,3 (A_theta, A_phi):
    ! radial column k=1 scales by ds/drho; angular columns unchanged.
    dA = 0.0_dp
    dA(2,1) = f%dAth(1)*ds_drho
    dA(2,2) = f%dAth(2)
    dA(2,3) = f%dAth(3)
    dA(3,1) = f%dAph(1)*ds_drho
    dA(3,2) = f%dAph(2)
    dA(3,3) = f%dAph(3)

    ! d|B|/du_k: radial column scales, angular columns unchanged.
    dBmod(1) = f%dBmod(1)*ds_drho
    dBmod(2) = f%dBmod(2)
    dBmod(3) = f%dBmod(3)
  end subroutine fo_eval_field

end module orbit_fo_field
