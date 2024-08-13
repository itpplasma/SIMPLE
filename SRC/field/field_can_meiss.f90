module field_can_meiss

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_base, only: FieldCan
use simple_magfie, only: MagneticField

implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

class(MagneticField), allocatable :: magfie
integer :: n_r=63, n_phi=64, n_th=65
real(dp) :: rmin=1d-12, rmax=1d0, thmin=0d0, thmax=twopi

real(dp) :: h_r, h_phi, h_th
real(dp), dimension(:,:,:), allocatable :: lam_phi, chi_gauge

contains

subroutine evaluate(f, r, th_c, ph_c, mode_secders)

    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    ! TODO
    f%Ath = 0d0
    f%Aph = 0d0

end subroutine evaluate

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

    r1 = rmin + h_r*(i_r-2)
    r2 = rmin + h_r*(i_r-1)

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

    phi_c = h_phi*(i_phi-1)
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
        th = thmin + h_th*(i_th-1)
        do i_phi=1,n_phi
            phi = h_phi*(i_phi-1)
            do i_r=1,n_r
                r = rmin + h_r*(i_r-1)
                write(funit, *) r, phi, th, lam_phi(i_r, i_phi, i_th), &
                    chi_gauge(i_r, i_phi, i_th)
            enddo
        enddo
    enddo

    close(funit)
end subroutine write_transformation


subroutine ah_cov_on_slice(r, phi, i_th, Ar, Ap, hr, hp)
    real(dp), intent(in) :: r, phi
    integer, intent(in) :: i_th
    real(dp), intent(out) :: Ar, Ap, hr, hp

    real(dp), dimension(3) :: Acov, hcov, sqgBctr
    real(dp) :: Bmod
    real(dp) :: th

    ! TODO: Make this more efficient with slices
    th = thmin + h_th*(i_th-1)
    call magfie%evaluate([r, th, phi], Acov, hcov, sqgBctr, Bmod)

    Ar = Acov(1)
    Ap = Acov(3)
    hr = hcov(1)
    hp = hcov(3)
end subroutine ah_cov_on_slice


subroutine init(magfie_, n_r_, n_phi_, n_th_, rmin_, rmax_, thmin_, thmax_)

  class(MagneticField), intent(in) :: magfie_
  integer, intent(in), optional :: n_r_, n_phi_, n_th_
  real(dp), intent(in), optional :: rmin_, rmax_, thmin_, thmax_

  if (allocated(magfie)) deallocate(magfie)
  allocate(magfie, source=magfie_)

  if (present(n_r_)) n_r = n_r_
  if (present(n_phi_)) n_phi = n_phi_
  if (present(n_th_)) n_th = n_th_

  if (present(rmin_)) rmin = rmin_
  if (present(rmax_)) rmax = rmax_
  if (present(thmin_)) thmin = thmin_
  if (present(thmax_)) thmax = thmax_

  h_r = (rmax-rmin)/(n_r-1)
  h_phi = twopi/(n_phi-1)
  h_th = (thmax-thmin)/(n_th-1)

end subroutine init

end module field_can_meiss
