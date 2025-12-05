module field_spline

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use interpolate, only : BatchSplineData3D, construct_batch_splines_3d
  use field_base, only : magnetic_field_t

  implicit none

contains

  subroutine build_field_splines(field, xmin, xmax, n_r, n_th, n_phi, spl)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: xmin(3)
    real(dp), intent(in) :: xmax(3)
    integer, intent(in) :: n_r, n_th, n_phi
    type(BatchSplineData3D), intent(out) :: spl

    real(dp) :: h_r, h_th, h_phi
    real(dp), dimension(:,:,:), allocatable :: Ar, Ath, Aphi, hr, hth, hphi, Bmod
    real(dp), dimension(:,:,:,:), allocatable :: y_batch
    real(dp) :: x(3)
    real(dp) :: Acov(3), hcov(3), Bmod_loc
    integer :: i_r, i_th, i_phi, dims(3)
    integer, parameter :: order(3) = [5, 5, 5]
    logical, parameter :: periodic(3) = [.false., .true., .true.]

    h_r = (xmax(1) - xmin(1)) / real(n_r - 1, dp)
    h_th = (xmax(2) - xmin(2)) / real(n_th - 1, dp)
    h_phi = (xmax(3) - xmin(3)) / real(n_phi - 1, dp)

    allocate(Ar(n_r, n_th, n_phi), Ath(n_r, n_th, n_phi), Aphi(n_r, n_th, n_phi))
    allocate(hr(n_r, n_th, n_phi), hth(n_r, n_th, n_phi), hphi(n_r, n_th, n_phi))
    allocate(Bmod(n_r, n_th, n_phi))

    do i_phi = 1, n_phi
      do i_th = 1, n_th
        do i_r = 1, n_r
          call grid_point(i_r, i_th, i_phi, xmin, h_r, h_th, h_phi, x)
          call field%evaluate(x, Acov, hcov, Bmod_loc)

          Ar(i_r, i_th, i_phi) = Acov(1)
          Ath(i_r, i_th, i_phi) = Acov(2)
          Aphi(i_r, i_th, i_phi) = Acov(3)

          hr(i_r, i_th, i_phi) = hcov(1)
          hth(i_r, i_th, i_phi) = hcov(2)
          hphi(i_r, i_th, i_phi) = hcov(3)

          Bmod(i_r, i_th, i_phi) = Bmod_loc
        end do
      end do
    end do

    Ar = Ar - Ar(1, 1, 1)
    Ath = Ath - Ath(1, 1, 1)
    Aphi = Aphi - Aphi(1, 1, 1)

    dims = shape(Ar)
    allocate(y_batch(dims(1), dims(2), dims(3), 7))

    y_batch(:,:,:,1) = Ar
    y_batch(:,:,:,2) = Ath
    y_batch(:,:,:,3) = Aphi
    y_batch(:,:,:,4) = hr
    y_batch(:,:,:,5) = hth
    y_batch(:,:,:,6) = hphi
    y_batch(:,:,:,7) = Bmod

    call construct_batch_splines_3d(xmin, xmax, y_batch, order, periodic, spl)
  end subroutine build_field_splines

  subroutine grid_point(i_r, i_th, i_phi, xmin, h_r, h_th, h_phi, x)
    integer, intent(in) :: i_r, i_th, i_phi
    real(dp), intent(in) :: xmin(3)
    real(dp), intent(in) :: h_r, h_th, h_phi
    real(dp), intent(out) :: x(3)

    x(1) = xmin(1) + h_r * real(i_r - 1, dp)
    x(2) = xmin(2) + h_th * real(i_th - 1, dp)
    x(3) = xmin(3) + h_phi * real(i_phi - 1, dp)
  end subroutine grid_point

end module field_spline

