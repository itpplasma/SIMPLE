module field_spline
  !> Generic field spline building utilities.
  !> Builds 7-component batch splines (Ar, Ath, Aphi, hr, hth, hphi, Bmod)
  !> over a reference coordinate grid for any magnetic_field_t implementation.

  use, intrinsic :: iso_fortran_env, only: dp => real64
  use interpolate, only: BatchSplineData3D, construct_batch_splines_3d
  use field_base, only: magnetic_field_t

  implicit none

contains

  subroutine build_field_splines(field, xmin_in, xmax_in, n_r, n_th, n_phi, spl)
    !> Build 7-component batch splines for a magnetic field over the given grid.
    !> Uses OpenMP parallelization over the phi dimension.
    !> Calls field%evaluate_spline_input for each grid point.
    class(magnetic_field_t), intent(inout) :: field
    real(dp), intent(in) :: xmin_in(3), xmax_in(3)
    integer, intent(in) :: n_r, n_th, n_phi
    type(BatchSplineData3D), intent(out) :: spl

    ! Local copies to match original init_splines structure
    real(dp) :: xmin(3), xmax(3)
    real(dp) :: h_r, h_th, h_phi

    integer, parameter :: order(3) = [5, 5, 5]
    logical, parameter :: periodic(3) = [.False., .True., .True.]

    real(dp) :: x(3)
    real(dp), dimension(:,:,:), allocatable :: Ar, Ath, Aphi, hr, hth, hphi, Bmod
    real(dp), dimension(:,:,:,:), allocatable :: y_batch
    real(dp) :: Acov(3), hcov(3)
    integer :: i_r, i_th, i_phi, i_ctr, dims(3)

    xmin = xmin_in
    xmax = xmax_in

    h_r = (xmax(1)-xmin(1))/(n_r-1)
    h_th = (xmax(2)-xmin(2))/(n_th-1)
    h_phi = (xmax(3)-xmin(3))/(n_phi-1)

    allocate(Ar(n_r,n_th,n_phi), Ath(n_r,n_th,n_phi), Aphi(n_r,n_th,n_phi))
    allocate(hr(n_r,n_th,n_phi), hth(n_r,n_th,n_phi), hphi(n_r,n_th,n_phi))
    allocate(Bmod(n_r,n_th,n_phi))

    i_ctr = 0
    !$omp parallel private(i_r, i_th, i_phi, x, Acov, hcov)
    !$omp do
    do i_phi=1,n_phi
        !$omp atomic
        i_ctr = i_ctr + 1
        call print_progress(i_ctr)
        do i_th=1,n_th
            do i_r=1,n_r
                x = get_grid_point(i_r, i_th, i_phi)

                call field%evaluate_spline_input(x, Acov, hcov, Bmod(i_r, i_th, i_phi))

                Ar(i_r, i_th, i_phi) = Acov(1)
                Ath(i_r, i_th, i_phi) = Acov(2)
                Aphi(i_r, i_th, i_phi) = Acov(3)

                hr(i_r, i_th, i_phi) = hcov(1)
                hth(i_r, i_th, i_phi) = hcov(2)
                hphi(i_r, i_th, i_phi) = hcov(3)
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    Ar = Ar - Ar(1,1,1)
    Ath = Ath - Ath(1,1,1)
    Aphi = Aphi - Aphi(1,1,1)

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

    contains

    pure function get_grid_point(i_r_, i_th_, i_phi_)
    real(dp) :: get_grid_point(3)
    integer, intent(in) :: i_r_, i_th_, i_phi_

    get_grid_point = [ &
        xmin(1) + h_r*dble(i_r_-1), &
        xmin(2) + h_th*dble(i_th_-1), &
        xmin(3) + h_phi*dble(i_phi_-1) &
    ]
    end function get_grid_point

    subroutine print_progress(i)
        integer, intent(in) :: i
        !$omp critical
        if (mod(i, n_phi / 10) == 0) then
            write(*, '(a,f6.2,a)', advance='no') char(13), 100.0*i/n_phi, ' %'
            call flush(6)
        else
            write(*, *)
        end if
        !$omp end critical
    end subroutine print_progress
  end subroutine build_field_splines

end module field_spline
