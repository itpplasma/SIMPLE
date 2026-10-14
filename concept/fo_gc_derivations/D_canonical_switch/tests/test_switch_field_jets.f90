program test_switch_field_jets
    ! Independent Cartesian differences of B check the logical-chart Hessian.
    use gc_fo_switch, only: dp, switch_jets_t, SW_OK
    use spectre_jets, only: sj_load, sj_cart, sj_field_at, sj_jets, sj_owner, sj_u
    implicit none
    type(switch_jets_t) :: jets
    real(dp) :: x(3), coarse(3, 3, 3), fine(3, 3, 3), scale, error, spacing
    character(len=1024) :: path
    integer :: owner, status, failures
    call get_command_argument(1, path)
    call sj_load(trim(path), status)
    if (status /= SW_OK) error stop 'switch jet fixture'
    failures = 0
    do owner = 1, 3
        sj_owner = owner
        sj_u = [real(owner, dp) - 0.5_dp, 0.4_dp, 0.3_dp]
        call sj_cart(sj_u, x)
        call sj_jets(x, jets, status)
        if (status /= SW_OK) error stop 'switch field jets'
        call cartesian_hessian(x, 3.0e-4_dp, coarse)
        call cartesian_hessian(x, 1.5e-4_dp, fine)
        scale = max(1.0_dp, maxval(abs(fine)))
        error = maxval(abs(jets%d2b - fine))/scale
        spacing = maxval(abs(coarse - fine))/scale
        print '(a,i0,a,2es14.5)', 'owner ', owner, &
            ' Hessian error / Cartesian spacing difference ', error, spacing
        if (error > 2.0e-4_dp .or. spacing > 2.0e-4_dp) failures = failures + 1
    end do
    if (failures /= 0) error stop 'Cartesian switch field Hessian oracle'
contains
    subroutine sample(x, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        integer :: ierr
        call sj_field_at(x, b, ierr)
        if (ierr /= SW_OK) error stop 'Cartesian field sample'
    end subroutine sample

    subroutine cartesian_hessian(x, h, hessian)
        real(dp), intent(in) :: x(3), h
        real(dp), intent(out) :: hessian(3, 3, 3)
        real(dp) :: e(3), d(3), b0(3), bp(3), bm(3), bpp(3), bmm(3)
        integer :: j, k
        call sample(x, b0)
        do j = 1, 3
            e = 0.0_dp
            e(j) = h
            call sample(x + e, bp)
            call sample(x - e, bm)
            call sample(x + 2*e, bpp)
            call sample(x - 2*e, bmm)
            hessian(:, j, j) = (-bpp + 16*bp - 30*b0 + 16*bm - bmm)/(12*h*h)
            do k = j + 1, 3
                d = 0.0_dp
                d(k) = h
                call sample(x + e + d, bpp)
                call sample(x + e - d, bp)
                call sample(x - e + d, bm)
                call sample(x - e - d, bmm)
                hessian(:, j, k) = (bpp - bp - bm + bmm)/(4*h*h)
                hessian(:, k, j) = hessian(:, j, k)
            end do
        end do
    end subroutine cartesian_hessian
end program test_switch_field_jets
