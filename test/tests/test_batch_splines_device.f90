program test_batch_splines_device
    ! Verify that the batch spline single-point evaluators produce identical
    ! results whether called on the host or from inside an OpenACC device
    ! kernel (acc routine seq). This is the leaf of the per-particle GPU orbit
    ! tracing kernel: the field evaluation must run correctly on the device,
    ! reading the managed-memory coefficient array.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
        evaluate_batch_splines_3d_der2
    implicit none

    integer, parameter :: n1 = 17, n2 = 25, n3 = 25, nq = 2, npt = 4096
    type(BatchSplineData3D) :: spl
    real(dp) :: xmin(3), xmax(3), field(n1, n2, n3, nq)
    integer :: order(3)
    logical :: periodic(3)
    real(dp) :: xpts(3, npt)
    real(dp) :: y_h(nq, npt), dy_h(3, nq, npt), d2y_h(6, nq, npt)
    real(dp) :: y_d(nq, npt), dy_d(3, nq, npt), d2y_d(6, nq, npt)
    real(dp) :: yt(nq), dyt(3, nq), d2yt(6, nq)
    real(dp) :: maxdiff, dd
    integer :: i, j, k, q, ip

    xmin = [0.0_dp, 0.0_dp, 0.0_dp]
    xmax = [1.0_dp, 6.283185307179586_dp, 6.283185307179586_dp]
    order = [5, 5, 5]
    periodic = [.false., .true., .true.]

    ! Smooth synthetic field
    do q = 1, nq
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1
                    associate (r => xmin(1) + (xmax(1) - xmin(1))*(i - 1)/real(n1 - 1, dp), &
                               th => (j - 1)*xmax(2)/real(n2, dp), &
                               ph => (k - 1)*xmax(3)/real(n3, dp))
                        field(i, j, k, q) = real(q, dp) + cos(th)*sin(ph)*(1.0_dp + r*r) &
                                            + 0.3_dp*r*cos(2.0_dp*th)
                    end associate
                end do
            end do
        end do
    end do

    call construct_batch_splines_3d(xmin, xmax, field, order, periodic, spl)

    ! Evaluation points inside the domain (avoid exact edges in r)
    do ip = 1, npt
        xpts(1, ip) = 0.02_dp + 0.96_dp*real(modulo(ip*7, 1000), dp)/1000.0_dp
        xpts(2, ip) = xmax(2)*real(modulo(ip*13, 1000), dp)/1000.0_dp
        xpts(3, ip) = xmax(3)*real(modulo(ip*29, 1000), dp)/1000.0_dp
    end do

    ! Host reference
    do ip = 1, npt
        call evaluate_batch_splines_3d_der2(spl, xpts(:, ip), yt, dyt, d2yt)
        y_h(:, ip) = yt
        dy_h(:, :, ip) = dyt
        d2y_h(:, :, ip) = d2yt
    end do

    ! Device evaluation: one point per gang/vector lane
    !$acc parallel loop gang vector private(yt, dyt, d2yt) &
    !$acc&   copyin(xpts) copyout(y_d, dy_d, d2y_d)
    do ip = 1, npt
        call evaluate_batch_splines_3d_der2(spl, xpts(:, ip), yt, dyt, d2yt)
        y_d(:, ip) = yt
        dy_d(:, :, ip) = dyt
        d2y_d(:, :, ip) = d2yt
    end do

    maxdiff = 0.0_dp
    do ip = 1, npt
        do q = 1, nq
            dd = abs(y_h(q, ip) - y_d(q, ip));        maxdiff = max(maxdiff, dd)
            do k = 1, 3
                dd = abs(dy_h(k, q, ip) - dy_d(k, q, ip)); maxdiff = max(maxdiff, dd)
            end do
            do k = 1, 6
                dd = abs(d2y_h(k, q, ip) - d2y_d(k, q, ip)); maxdiff = max(maxdiff, dd)
            end do
        end do
    end do

    print '(a, i0, a)', 'Device batch eval over ', npt, ' points'
    print '(a, es12.4)', 'max |host - device| (value+der+der2) = ', maxdiff

    if (maxdiff > 1.0e-10_dp) then
        print *, 'FAILED: device evaluation differs from host'
        error stop 1
    end if
    print *, 'PASSED: device batch spline evaluation matches host'
end program test_batch_splines_device
