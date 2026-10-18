module spectre_jets
    ! Field provider for the item-D switch study: Cartesian B and grad B of a
    ! SPECTRE equilibrium (libneo field_spectre, SI units: m, T) in one owner
    ! volume. B_cart = e_i (sqrtg B^i)/sqrtg is evaluated in logical u; grad B
    ! uses an 8th-order central difference in u, mapped by the inverse
    ! covariant basis. The Hessian includes the geometry chain-rule term.
    ! The difference step is fixed, so the provider is a
    ! deterministic smooth function of X (round-trip checks rely on this).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_spectre, only: spectre_field_t, create_spectre_field
    use gc_fo_switch, only: switch_jets_t, SW_OK, SW_FIELD
    implicit none
    private
    public :: sj_load, sj_locate, sj_bcart, sj_jets, sj_field_at, sj_cart
    public :: sj_acov, sj_basis, sj_owner, sj_u
    type(spectre_field_t), save, target :: field
    integer, save :: sj_owner = 1
    real(dp), save :: sj_u(3) = 0.0_dp
    real(dp), parameter :: hstep(3) = [1.0e-3_dp, 1.0e-3_dp, 1.0e-3_dp]
    real(dp), parameter :: c8(4) = [4.0_dp/5.0_dp, -1.0_dp/5.0_dp, &
                                    4.0_dp/105.0_dp, -1.0_dp/280.0_dp]
contains
    subroutine sj_load(path, ierr)
        character(*), intent(in) :: path
        integer, intent(out) :: ierr
        call create_spectre_field(field, path, ierr)
    end subroutine sj_load

    subroutine sj_cart(u, x)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)
        call field%coords%evaluate_cart(u, x)
    end subroutine sj_cart

    subroutine sj_basis(u, ecov)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: ecov(3, 3)
        call field%coords%covariant_basis(u, ecov)
    end subroutine sj_basis

    subroutine sj_acov(u, acov)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: acov(3)
        real(dp) :: h(3), bm, sq(3)
        call field%evaluate(u, acov, h, bm, sq)
    end subroutine sj_acov

    !> Locate Cartesian x in the owner volume: Newton on evaluate_cart from the
    !> warm start sj_u (falls back to from_cyl), polished to roundoff.
    subroutine sj_locate(x, u, ierr)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr
        real(dp) :: xt(3), ecov(3, 3), du(3), xcyl(3), best, rn
        integer :: it, ipiv(3), info
        u = sj_u
        best = huge(1.0_dp)
        ierr = SW_FIELD
        do it = 1, 40
            call field%coords%evaluate_cart(u, xt)
            rn = norm2(xt - x)
            if (rn < 1.0e-14_dp*max(1.0_dp, norm2(x)) .or. &
                (rn >= 0.5_dp*best .and. rn < 1.0e-11_dp)) then
                ierr = SW_OK
                exit
            end if
            best = min(best, rn)
            call field%coords%covariant_basis(u, ecov)
            du = x - xt
            call dgesv(3, 1, ecov, 3, ipiv, du, 3, info)
            if (info /= 0) exit
            u = u + du
        end do
        if (ierr /= SW_OK) then
            xcyl = [hypot(x(1), x(2)), atan2(x(2), x(1)), x(3)]
            call field%coords%from_cyl(xcyl, u, info)
            if (info /= 0) return
            ierr = SW_OK
        end if
        if (u(1) <= real(sj_owner - 1, dp) .or. u(1) >= real(sj_owner, dp)) then
            ierr = SW_FIELD
            return
        end if
        sj_u = u
    end subroutine sj_locate

    subroutine sj_bcart(u, b)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: b(3)
        real(dp) :: acov(3), h(3), bm, sq(3), g(3, 3), ginv(3, 3), sqrtg, ecov(3, 3)
        call field%evaluate(u, acov, h, bm, sq)
        call field%coords%metric_tensor(u, g, ginv, sqrtg)
        call field%coords%covariant_basis(u, ecov)
        b = matmul(ecov, sq/sqrtg)
    end subroutine sj_bcart

    subroutine sj_field_at(x, b, ierr)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        integer, intent(out) :: ierr
        real(dp) :: u(3)
        b = 0.0_dp
        call sj_locate(x, u, ierr)
        if (ierr /= SW_OK) return
        call sj_bcart(u, b)
    end subroutine sj_field_at

    subroutine sj_jets(X, jets, ierr)
        real(dp), intent(in) :: X(3)
        type(switch_jets_t), intent(inout) :: jets
        integer, intent(out) :: ierr
        real(dp) :: u(3), up(3), um(3), bp(3), bm(3), dbdu(3, 3), ecov(3, 3)
        real(dp) :: einv(3, 3), d2bdu(3, 3, 3), d2xdu(3, 3, 3)
        real(dp) :: ep(3, 3), em(3, 3), corrected(3, 3), value(3)
        real(dp), parameter :: csecond(0:4) = [-205.0_dp/72, 8.0_dp/5, &
            -1.0_dp/5, 8.0_dp/315, -1.0_dp/560]
        integer :: i, j, k, l, m, ipiv(3), info
        call sj_locate(X, u, ierr)
        if (ierr /= SW_OK) return
        call sj_bcart(u, jets%b)
        do j = 1, 3
            dbdu(:, j) = 0.0_dp
            d2bdu(:, j, j) = csecond(0)*jets%b
            d2xdu(:, :, j) = 0.0_dp
            do k = 1, 4
                up = u
                um = u
                up(j) = u(j) + k*hstep(j)
                um(j) = u(j) - k*hstep(j)
                call sj_bcart(up, bp)
                call sj_bcart(um, bm)
                dbdu(:, j) = dbdu(:, j) + c8(k)*(bp - bm)
                d2bdu(:, j, j) = d2bdu(:, j, j) + csecond(k)*(bp + bm)
                call field%coords%covariant_basis(up, ep)
                call field%coords%covariant_basis(um, em)
                d2xdu(:, :, j) = d2xdu(:, :, j) + c8(k)*(ep - em)
            end do
            dbdu(:, j) = dbdu(:, j)/hstep(j)
            d2bdu(:, j, j) = d2bdu(:, j, j)/hstep(j)**2
            d2xdu(:, :, j) = d2xdu(:, :, j)/hstep(j)
        end do
        do j = 1, 3
            do k = j + 1, 3
                d2bdu(:, j, k) = 0.0_dp
                do l = -4, 4
                    if (l == 0) cycle
                    do m = -4, 4
                        if (m == 0) cycle
                        up = u
                        up(j) = up(j) + l*hstep(j)
                        up(k) = up(k) + m*hstep(k)
                        call sj_bcart(up, value)
                        d2bdu(:, j, k) = d2bdu(:, j, k) + &
                            sign(1.0_dp, real(l*m, dp))*c8(abs(l))*c8(abs(m))*value
                    end do
                end do
                d2bdu(:, j, k) = d2bdu(:, j, k)/(hstep(j)*hstep(k))
                d2bdu(:, k, j) = d2bdu(:, j, k)
                ! Both independent basis derivatives describe the same jet.
                value = 0.5_dp*(d2xdu(:, j, k) + d2xdu(:, k, j))
                d2xdu(:, j, k) = value
                d2xdu(:, k, j) = value
            end do
        end do
        call field%coords%covariant_basis(u, ecov)
        einv = 0.0_dp
        do j = 1, 3
            einv(j, j) = 1.0_dp
        end do
        call dgesv(3, 3, ecov, 3, ipiv, einv, 3, info)
        if (info /= 0) then
            ierr = SW_FIELD
            return
        end if
        jets%db = matmul(dbdu, einv)
        ! B_,u_j u_k = B_,x_a x_b E_a j E_b k + B_,x_a X_a,jk.
        ! Subtract the chart curvature before applying the inverse basis twice.
        do i = 1, 3
            corrected = d2bdu(i, :, :)
            do j = 1, 3
                corrected = corrected - jets%db(i, j)*d2xdu(j, :, :)
            end do
            jets%d2b(i, :, :) = matmul(transpose(einv), matmul(corrected, einv))
        end do
        sj_u = u
    end subroutine sj_jets
end module spectre_jets
