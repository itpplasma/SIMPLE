module field_spline
  !> Generic field spline building utilities.
  !> Builds 7-component batch splines (Ar, Ath, Aphi, hr, hth, hphi, Bmod)
  !> over a reference coordinate grid for any magnetic_field_t implementation.

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use interpolate, only : BatchSplineData3D, construct_batch_splines_3d
  use field_base, only : magnetic_field_t

  implicit none

  abstract interface
    subroutine field_evaluator(field_, x_, Acov_, hcov_, Bmod_)
      import :: dp, magnetic_field_t
      class(magnetic_field_t), intent(in) :: field_
      real(dp), intent(in) :: x_(3)
      real(dp), intent(out) :: Acov_(3), hcov_(3), Bmod_
    end subroutine field_evaluator
  end interface

contains

  subroutine build_field_splines(field, xmin, xmax, n_r, n_th, n_phi, spl, &
                                  evaluator)
    !> Build 7-component batch splines for a magnetic field over the given grid.
    !> Uses OpenMP parallelization over the phi dimension.
    !>
    !> @param field The magnetic field to spline
    !> @param xmin Grid minimum bounds (r, theta, phi)
    !> @param xmax Grid maximum bounds (r, theta, phi)
    !> @param n_r Number of radial grid points
    !> @param n_th Number of theta grid points
    !> @param n_phi Number of phi grid points
    !> @param spl Output batch spline data
    !> @param evaluator Optional custom evaluation routine (for direct eval)
    class(magnetic_field_t), intent(in) :: field
    real(dp), intent(in) :: xmin(3)
    real(dp), intent(in) :: xmax(3)
    integer, intent(in) :: n_r, n_th, n_phi
    type(BatchSplineData3D), intent(out) :: spl
    procedure(field_evaluator), optional :: evaluator

    real(dp) :: h_r, h_th, h_phi
    real(dp), dimension(:,:,:), allocatable :: Ar, Ath, Aphi, hr, hth, hphi, Bmod
    real(dp), dimension(:,:,:,:), allocatable :: y_batch
    real(dp) :: x(3)
    real(dp) :: Acov(3), hcov(3), Bmod_loc
    integer :: i_r, i_th, i_phi, i_ctr, dims(3)
    integer, parameter :: order(3) = [5, 5, 5]
    logical, parameter :: periodic(3) = [.false., .true., .true.]

    h_r = (xmax(1) - xmin(1)) / (n_r - 1)
    h_th = (xmax(2) - xmin(2)) / (n_th - 1)
    h_phi = (xmax(3) - xmin(3)) / (n_phi - 1)

    allocate(Ar(n_r, n_th, n_phi), Ath(n_r, n_th, n_phi), Aphi(n_r, n_th, n_phi))
    allocate(hr(n_r, n_th, n_phi), hth(n_r, n_th, n_phi), hphi(n_r, n_th, n_phi))
    allocate(Bmod(n_r, n_th, n_phi))

    i_ctr = 0
    !$omp parallel private(i_r, i_th, i_phi, x, Acov, hcov, Bmod_loc)
    !$omp do
    do i_phi = 1, n_phi
      !$omp atomic
      i_ctr = i_ctr + 1
      !$omp critical
      if (mod(i_ctr, max(1, n_phi / 10)) == 0) then
        write(*, '(a,a,f6.2,a)') char(13), 'Splining field: ', &
            100.0_dp * i_ctr / n_phi, ' %'
        call flush(6)
      end if
      !$omp end critical
      do i_th = 1, n_th
        do i_r = 1, n_r
          x(1) = xmin(1) + h_r * dble(i_r - 1)
          x(2) = xmin(2) + h_th * dble(i_th - 1)
          x(3) = xmin(3) + h_phi * dble(i_phi - 1)

          if (present(evaluator)) then
            call evaluator(field, x, Acov, hcov, Bmod_loc)
          else
            call field%evaluate(x, Acov, hcov, Bmod_loc)
          end if

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
    !$omp end do
    !$omp end parallel
    write(*, *)  ! newline after progress

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

end module field_spline
