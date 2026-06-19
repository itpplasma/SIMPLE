module orbit_cpp_vmec_metric
  ! Host-side metric + field provider for the 6D canonical-midpoint integrator on
  ! REAL VMEC equilibria, in NATIVE VMEC flux coordinates u = (s, vartheta,
  ! varphi). Pairs SIMPLE's native VMEC field (covariant A_i, contravariant B^i,
  ! |B|) with libneo's coordinate_system metric tensor and Christoffel symbols
  ! (issue #322, libneo feature/metric-christoffel).
  !
  ! This is the general curvilinear path the analytic-tokamak block reduces to:
  ! the full (non-diagonal) metric g_ij/g^ij and its direction derivatives dg_ij,k
  ! drive the same geodesic momentum equation that orbit_cpp_canonical solves.
  !
  ! NOT GPU-portable: libneo's metric_tensor/christoffel are class()-dispatched
  ! and read 3D splines. The analytic-tokamak COORD_TOK block in
  ! orbit_cpp_canonical stays !$acc routine seq and class-free; only this VMEC
  ! block is host-side. The Newton LU solve (rk_solve) is portable for both.
  !
  ! Metric derivatives come from Christoffel via metric compatibility:
  !   dg_ij/du_k = g_il Gamma^l_jk + g_jl Gamma^l_ik.
  ! Field gradients of |B| use a central difference of the native |B| (the native
  ! evaluator returns covariant B_i and B^i but not analytic d|B|, the same
  ! central-difference convention orbit_cpp_canonical uses for the tokamak block).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use libneo_coordinates_base, only: coordinate_system_t
  use field_can_base, only: n_field_evaluations
  implicit none
  private

  public :: vmec_metric_init, vmec_metric_attach, vmec_metric_ready
  public :: vmec_eval_metric, vmec_eval_field, vmec_bmod

  class(coordinate_system_t), allocatable :: cs
  logical :: ready = .false.

contains

  ! Load VMEC splines from a wout file and build the libneo VMEC coordinate
  ! system. Idempotent guard via vmec_metric_ready. Stand-alone entry for tests
  ! that have not already splined a VMEC equilibrium.
  subroutine vmec_metric_init(wout_file)
    use new_vmec_stuff_mod, only: netcdffile, multharm, ns_s, ns_tp
    use spline_vmec_sub, only: spline_vmec_data
    character(*), intent(in) :: wout_file

    netcdffile = wout_file
    ns_s = 5
    ns_tp = 5
    multharm = 3
    call spline_vmec_data
    call vmec_metric_attach
  end subroutine vmec_metric_init

  ! Build the libneo VMEC coordinate system from VMEC splines that the caller has
  ! already loaded (production init_vmec/init_field). No re-splining, so the
  ! production ns_s/ns_tp/multharm and the equilibrium scaling are preserved.
  subroutine vmec_metric_attach
    use libneo_coordinates_vmec, only: make_vmec_coordinate_system

    if (allocated(cs)) deallocate(cs)
    call make_vmec_coordinate_system(cs)
    ready = .true.
  end subroutine vmec_metric_attach

  logical function vmec_metric_ready()
    vmec_metric_ready = ready
  end function vmec_metric_ready

  ! Full metric block at u=(s,th,ph): g_ij, g^ij, and dg(i,j,k)=dg_ij/du_k from
  ! Christoffel + metric compatibility.
  subroutine vmec_eval_metric(u, g, ginv, dg)
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: g(3,3), ginv(3,3), dg(3,3,3)
    real(dp) :: sqrtg, Gamma(3,3,3)
    integer :: i, j, k, l

    call cs%metric_tensor(u, g, ginv, sqrtg)
    call cs%christoffel(u, Gamma)
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
  end subroutine vmec_eval_metric

  ! Native VMEC field block at u=(s,th,ph): covariant A (A_s=0), |B|, and the
  ! covariant gradient d|B|/du via central difference. h_i = g_il B^l / |B|.
  subroutine vmec_eval_field(u, Acov, Bmod, dBmod, hcov)
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: Acov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: g(3,3), ginv(3,3), sqrtg
    real(dp) :: Bctr(3), Bcov(3)
    real(dp), parameter :: h = 1.0e-6_dp
    real(dp) :: up(3), um(3), bp, bm
    integer :: k, i

    ! Count one 6D field evaluation per primary call (mirrors field_can: one count
    ! per evaluate). The FD-perturbation vmec_bmod calls below are not counted, the
    ! same convention the field_can evaluators use for their derivative stencils.
    n_field_evaluations = n_field_evaluations + 1
    call native_field(u, Acov, Bctr, Bcov, Bmod)
    call cs%metric_tensor(u, g, ginv, sqrtg)
    ! h_i = B_i/|B| (covariant unit field; B_i already covariant from VMEC).
    do i = 1, 3
      hcov(i) = Bcov(i)/Bmod
    end do
    do k = 1, 3
      up = u; um = u; up(k) = up(k) + h; um(k) = um(k) - h
      bp = vmec_bmod(up); bm = vmec_bmod(um)
      dBmod(k) = (bp - bm)/(2.0_dp*h)
    end do
  end subroutine vmec_eval_field

  ! |B| only, for the central-difference gradient.
  real(dp) function vmec_bmod(u)
    real(dp), intent(in) :: u(3)
    real(dp) :: Acov(3), Bctr(3), Bcov(3)
    call native_field(u, Acov, Bctr, Bcov, vmec_bmod)
  end function vmec_bmod

  ! Native VMEC field: covariant A_i (A_s=0), contravariant B^i, covariant B_i,
  ! |B|. Uses SIMPLE's vmec_field_evaluate (libneo splint_vmec_data underneath).
  subroutine native_field(u, Acov, Bctr, Bcov, Bmod)
    use vmec_field_eval, only: vmec_field_evaluate
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: Acov(3), Bctr(3), Bcov(3), Bmod
    real(dp) :: s, theta, varphi
    real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, sqg, alam
    real(dp) :: dl_ds, dl_dt, dl_dp
    real(dp) :: Bctrvr_vartheta, Bctrvr_varphi
    real(dp) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi

    s = u(1); theta = u(2); varphi = u(3)
    call vmec_field_evaluate(s, theta, varphi, A_theta, A_phi, &
         dA_theta_ds, dA_phi_ds, aiota, sqg, alam, dl_ds, dl_dt, dl_dp, &
         Bctrvr_vartheta, Bctrvr_varphi, Bcovar_r, Bcovar_vartheta, Bcovar_varphi)

    Acov = [0.0_dp, A_theta, A_phi]
    Bctr = [0.0_dp, Bctrvr_vartheta, Bctrvr_varphi]
    Bcov = [Bcovar_r, Bcovar_vartheta, Bcovar_varphi]
    Bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
  end subroutine native_field

end module orbit_cpp_vmec_metric
