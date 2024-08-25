module field_can_meiss

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_base, only: FieldCan, n_field_evaluations
use simple_magfie, only: MagneticField
use interpolate, only: SplineData3D, construct_splines_3d, &
    evaluate_splines_3d, evaluate_splines_3d_der2

implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

class(MagneticField), allocatable :: magfie
integer :: n_r=63, n_phi=64, n_th=65
real(dp) :: xmin(3) = [1d-12, 0d0, 0d0]
real(dp) :: xmax(3) = [1d0, twopi, twopi]

real(dp) :: h_r, h_phi, h_th

! For splining covariant vector potential, h=B/Bmod and Bmod over canonical coordinates
type(SplineData3D) :: spl_Bmod, spl_A2, spl_A3, spl_h2, spl_h3

! For splining lambda (difference between canonical and toroidal cylinder angle)
! and chi (gauge transformation)
real(dp), dimension(:,:,:), allocatable :: lam_phi, chi_gauge
type(SplineData3D) :: spl_lam_phi, spl_chi_gauge

integer, parameter :: order(3) = [5, 5, 5]
logical, parameter :: periodic(3) = [.False., .True., .True.]

contains

subroutine init(magfie_, n_r_, n_phi_, n_th_, rmin, rmax, thmin, thmax)

  class(MagneticField), intent(in) :: magfie_
  integer, intent(in), optional :: n_r_, n_phi_, n_th_
  real(dp), intent(in), optional :: rmin, rmax, thmin, thmax

  if (allocated(magfie)) deallocate(magfie)
  allocate(magfie, source=magfie_)

  if (present(n_r_)) n_r = n_r_
  if (present(n_phi_)) n_phi = n_phi_
  if (present(n_th_)) n_th = n_th_

  if (present(rmin)) xmin(1) = rmin
  if (present(rmax)) xmax(1) = rmax
  if (present(thmin)) xmax(3) = thmin
  if (present(thmax)) xmax(3) = thmax

  h_r = (xmax(1)-xmin(1))/(n_r-1)
  h_phi = (xmax(2)-xmin(2))/(n_phi-1)
  h_th = (xmax(3)-xmin(3))/(n_th-1)

end subroutine init


subroutine evaluate(f, r, th_c, ph_c, mode_secders)

    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    integer, parameter :: reorder(3) = [1, 3, 2]  ! dr, dph, dth -> dr, dth, dph
    integer, parameter :: reorder2(6) = [1, 3, 2, 6, 5, 4]
    ! drdr, drdph, drdth, dphdph, dphdth, dthdth ->
    ! drdr, drdth, drdph, dthdth, dthdph, dphdph

    real(dp) :: x(3), a, da(3), d2a(6)

    x = [r, ph_c, th_c]

    call evaluate_splines_3d_der2(spl_A2, x, a, da, d2a)
    f%Aph = a
    f%dAph = da(reorder)
    f%d2Aph = d2a(reorder2)

    call evaluate_splines_3d_der2(spl_A3, x, a, da, d2a)
    f%Ath = a
    f%dAth = da(reorder)
    f%d2Ath = d2a(reorder2)

    call evaluate_splines_3d_der2(spl_h2, x, a, da, d2a)
    f%hph = a
    f%dhph = da(reorder)
    f%d2hph = d2a(reorder2)

    call evaluate_splines_3d_der2(spl_h3, x, a, da, d2a)
    f%hth = a
    f%dhth = da(reorder)
    f%d2hth = d2a(reorder2)

    call evaluate_splines_3d_der2(spl_Bmod, x, a, da, d2a)
    f%Bmod = a
    f%dBmod = da(reorder)
    f%d2Bmod = d2a(reorder2)

    n_field_evaluations = n_field_evaluations + 1

end subroutine evaluate


subroutine get_meiss_coordinates
    print *, 'field_can_meiss.init_transformation'
    call init_transformation

    print *, 'field_can_meiss.spline_transformation'
    call spline_transformation

    print *, 'field_can_meiss.init_canonical_field_components'
    call init_canonical_field_components
end subroutine get_meiss_coordinates


subroutine init_transformation
    integer :: i_r, i_phi, i_th, i_ctr

    allocate(lam_phi(n_r, n_phi, n_th), &
             chi_gauge(n_r, n_phi, n_th))

    i_ctr=0

    !$omp parallel private(i_r, i_phi, i_th)
    !$omp do
    do i_th=1,n_th
        !$omp critical
        i_ctr = i_ctr + 1
        call print_progress
        !$omp end critical

        do i_phi=1,n_phi
            lam_phi(1, i_phi, i_th) = 0.d0
            chi_gauge(1, i_phi, i_th) = 0.d0

            do i_r=2,n_r
                call integrate(i_r, i_phi, i_th)
            enddo
        enddo
    enddo
    !$omp end do
    !$omp end parallel

    contains

    subroutine print_progress
        write(*,'(A, I4, A, I4)',advance='no') 'integrate ODE: ', i_ctr, ' of ', n_th
        if (i_ctr < n_th) then
            write(*, '(A)', advance="no") char(13)
        else
            write(*, *)
        end if
    end subroutine print_progress
end subroutine init_transformation

subroutine integrate(i_r, i_phi, i_th)
    use odeint_sub, only: odeint_allroutines

    integer, intent(in) :: i_r, i_phi, i_th
    real(dp), dimension(2) :: y

    real(dp), parameter :: relerr=1d-11
    real(dp) :: r1, r2
    integer :: ndim=2

    y = 0d0

    r1 = xmin(1) + h_r*(i_r-2)
    r2 = xmin(1) + h_r*(i_r-1)

    call odeint_allroutines(y, ndim, r1, r2, relerr, rh_can_closure)

    lam_phi(i_r, i_phi, i_th) = y(1)
    chi_gauge(i_r, i_phi, i_th) = y(2)

    contains

    subroutine rh_can_closure(r_c, z, dz)
        real(dp), intent(in) :: r_c
        real(dp), dimension(2), intent(in) :: z
        real(dp), dimension(2), intent(inout) :: dz

        call rh_can(r_c, z, dz, i_phi, i_th)
    end subroutine rh_can_closure

end subroutine integrate

subroutine rh_can(r_c, z, dz, i_phi, i_th)
    real(dp), intent(in) :: r_c  ! plus threadprivate phi_c, th_c from module
    real(dp), dimension(2), intent(in) :: z  ! lam_phi, chi_gauge
    real(dp), dimension(2), intent(inout) :: dz
    integer, intent(in) :: i_phi, i_th
    real(dp) :: phi_c
    real(dp) :: Ar, Ap, hr, hp

    phi_c = xmin(2) + h_phi*(i_phi-1)
    call ah_cov_on_slice(r_c, modulo(phi_c + z(1), twopi), i_th, Ar, Ap, hr, hp)

    dz(1) = -hr/hp
    dz(2) = Ar + Ap*dz(1)
end subroutine rh_can


subroutine write_transformation(filename)
    character(*), intent(in) :: filename

    integer :: funit
    integer :: i_r, i_phi, i_th
    real(dp) :: r, phi, th

    open(newunit=funit, file=filename, status='unknown')
    write(funit, *) '#', ' r', ' phi', ' th', ' lam_phi', ' chi_gauge'

    do i_th=1,n_th
        th = xmin(3) + h_th*(i_th-1)
        do i_phi=1,n_phi
            phi = xmin(2) + h_phi*(i_phi-1)
            do i_r=1,n_r
                r = xmin(1) + h_r*(i_r-1)
                write(funit, *) r, phi, th, lam_phi(i_r, i_phi, i_th), &
                    chi_gauge(i_r, i_phi, i_th)
            enddo
        enddo
    enddo

    close(funit)
end subroutine write_transformation


subroutine spline_transformation
    call construct_splines_3d(xmin, xmax, lam_phi, order, periodic, spl_lam_phi)
    call construct_splines_3d(xmin, xmax, chi_gauge, order, periodic, spl_chi_gauge)
end subroutine spline_transformation


subroutine ah_cov_on_slice(r, phi, i_th, Ar, Ap, hr, hp)
    real(dp), intent(in) :: r, phi
    integer, intent(in) :: i_th
    real(dp), intent(out) :: Ar, Ap, hr, hp

    real(dp), dimension(3) :: Acov, hcov
    real(dp) :: Bmod
    real(dp) :: th

    ! TODO: Make this more efficient with slices
    ! TODO: Support not only VMEC field
    th = xmin(3) + h_th*(i_th-1)
    call magfie%evaluate([r, th, phi], Acov, hcov, Bmod)

    Ar = Acov(1)
    Ap = Acov(3)
    hr = hcov(1)
    hp = hcov(3)
end subroutine ah_cov_on_slice


subroutine init_canonical_field_components
    real(dp), dimension(:,:,:,:), allocatable :: xcan
    real(dp), dimension(:,:,:), allocatable :: A2, A3, h2, h3, Bmod
    real(dp) :: xbase(3), Acov(3), hcov(3), lam, dlam(3), chi, dchi(3), dummy(6)
    integer :: i_r, i_phi, i_th

    allocate(xcan(3,n_r,n_phi,n_th))
    call generate_regular_grid(xcan)

    allocate(A2(n_r,n_phi,n_th), A3(n_r,n_phi,n_th))
    allocate(h2(n_r,n_phi,n_th), h3(n_r,n_phi,n_th))
    allocate(Bmod(n_r,n_phi,n_th))

    do i_th=1,n_th
        do i_phi=1,n_phi
            do i_r=1,n_r
                call can_to_base(xcan(:,i_r,i_phi,i_th), xbase)
                call magfie%evaluate(xbase, Acov, hcov, Bmod(i_r, i_phi, i_th))

                call evaluate_splines_3d_der2(spl_lam_phi, &
                    xcan(:,i_r,i_phi,i_th), lam, dlam, dummy)
                call evaluate_splines_3d_der2(spl_chi_gauge, &
                    xcan(:,i_r,i_phi,i_th), chi, dchi, dummy)

                A2(i_r, i_phi, i_th) = Acov(3)*(1.0d0 + dlam(2)) - dchi(2)
                A3(i_r, i_phi, i_th) = Acov(2) + Acov(2)*dlam(3) - dchi(3)
                h2(i_r, i_phi, i_th) = hcov(3)*(1.0d0 + dlam(2))
                h3(i_r, i_phi, i_th) = hcov(2) + hcov(2)*dlam(3)
            end do
        end do
    end do

    call construct_splines_3d(xmin, xmax, A2(:,:,:), order, periodic, spl_A2)
    call construct_splines_3d(xmin, xmax, A3(:,:,:), order, periodic, spl_A3)
    call construct_splines_3d(xmin, xmax, h2(:,:,:), order, periodic, spl_h2)
    call construct_splines_3d(xmin, xmax, h3(:,:,:), order, periodic, spl_h3)
    call construct_splines_3d(xmin, xmax, Bmod, order, periodic, spl_Bmod)

    deallocate(Bmod)
    deallocate(h2, h3)
    deallocate(A2, A3)

    deallocate(xcan)
end subroutine init_canonical_field_components


subroutine can_to_base(xcan, xbase)
    real(dp), intent(in) :: xcan(3)
    real(dp), intent(out) :: xbase(3)
    real(dp) :: lam

    call evaluate_splines_3d(spl_lam_phi, xcan, lam)
    xbase(1) = xcan(1)
    xbase(2) = xcan(3)
    xbase(3) = modulo(xcan(2) + lam, twopi)
end subroutine can_to_base


pure subroutine generate_regular_grid(x)
real(dp), intent(inout) :: x(:,:,:,:)

integer :: i_r, i_phi, i_th

do i_th=1,n_th
    do i_phi=1,n_phi
        do i_r=1,n_r
            x(:, i_r, i_phi, i_th) = get_grid_point(i_r, i_phi, i_th)
        enddo
    enddo
enddo

end subroutine generate_regular_grid

pure function get_grid_point(i_r, i_phi, i_th)
    integer, intent(in) :: i_r, i_phi, i_th
    real(dp) :: get_grid_point(3)

    get_grid_point = [ &
        xmin(1) + h_r*dble(i_r-1), &
        xmin(2) + h_phi*dble(i_phi-1), &
        xmin(3) + h_th*dble(i_th-1) &
    ]
end function get_grid_point

end module field_can_meiss
