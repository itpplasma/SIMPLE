module minpack_interfaces
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    interface
        function enorm(n, x)
            import :: dp
            integer, intent(in) :: n
            real(dp), intent(in) :: x(n)
            real(dp) :: enorm
        end function enorm

        subroutine dogleg(n, r, lr, diag, qtb, delta, x)
            import :: dp
            integer, intent(in) :: n, lr
            real(dp), intent(in) :: r(lr), diag(n), qtb(n), delta
            real(dp), intent(out) :: x(n)
        end subroutine dogleg

        subroutine fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn)
            import :: dp
            external :: fcn
            integer, intent(in) :: n, ldfjac, ml, mu
            integer, intent(inout) :: iflag
            real(dp), intent(in) :: x(n), epsfcn
            real(dp), intent(inout) :: fvec(n), fjac(ldfjac, n)
        end subroutine fdjac1

        subroutine fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn, wa)
            import :: dp
            external :: fcn
            integer, intent(in) :: m, n, ldfjac
            integer, intent(inout) :: iflag
            real(dp), intent(in) :: x(n), epsfcn
            real(dp), intent(inout) :: fvec(m), fjac(ldfjac, n), wa(m)
        end subroutine fdjac2

        subroutine qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm)
            import :: dp
            integer, intent(in) :: m, n, lda, lipvt
            logical, intent(in) :: pivot
            real(dp), intent(inout) :: a(lda, n)
            integer, intent(out) :: ipvt(lipvt)
            real(dp), intent(out) :: rdiag(n), acnorm(n)
        end subroutine qrfac

        subroutine qform(m, n, q, ldq)
            import :: dp
            integer, intent(in) :: m, n, ldq
            real(dp), intent(inout) :: q(ldq, m)
        end subroutine qform

        subroutine qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)
            import :: dp
            integer, intent(in) :: n, ldr
            integer, intent(in) :: ipvt(n)
            real(dp), intent(inout) :: r(ldr, n), diag(n), qtb(n)
            real(dp), intent(out) :: x(n), sdiag(n)
            real(dp), intent(inout) :: wa(n)
        end subroutine qrsolv

        subroutine r1mpyq(m, n, a, lda, v, w)
            import :: dp
            integer, intent(in) :: m, n, lda
            real(dp), intent(inout) :: a(lda, n)
            real(dp), intent(in) :: v(n), w(n)
        end subroutine r1mpyq

        subroutine r1updt(m, n, s, ls, u, v, w, sing)
            import :: dp
            integer, intent(in) :: m, n, ls
            real(dp), intent(inout) :: s(ls), u(m), v(n)
            real(dp), intent(out) :: w(m)
            logical, intent(out) :: sing
        end subroutine r1updt

        subroutine rwupdt(n, r, ldr, w, b, alpha, cos, sin)
            import :: dp
            integer, intent(in) :: n, ldr
            real(dp), intent(inout) :: r(ldr, n), w(n), b(n)
            real(dp), intent(inout) :: alpha
            real(dp), intent(out) :: cos(n), sin(n)
        end subroutine rwupdt

        subroutine lmpar(n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag, wa1, wa2)
            import :: dp
            integer, intent(in) :: n, ldr
            integer, intent(in) :: ipvt(n)
            real(dp), intent(inout) :: r(ldr, n), diag(n), qtb(n), delta, par
            real(dp), intent(out) :: x(n), sdiag(n)
            real(dp), intent(inout) :: wa1(n), wa2(n)
        end subroutine lmpar

        subroutine hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
                         factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf, wa1, &
                         wa2, wa3, wa4)
            import :: dp
            external :: fcn
            integer, intent(in) :: n, maxfev, ml, mu, mode, nprint, ldfjac, lr
            integer, intent(out) :: info, nfev
            real(dp), intent(in) :: xtol, epsfcn, factor
            real(dp), intent(inout) :: x(n), diag(n)
            real(dp), intent(out) :: fvec(n), fjac(ldfjac, n), r(lr), qtf(n)
            real(dp), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(n)
        end subroutine hybrd

        subroutine hybrd1(fcn, n, x, fvec, tol, info)
            import :: dp
            external :: fcn
            integer, intent(in) :: n
            integer, intent(out) :: info
            real(dp), intent(in) :: tol
            real(dp), intent(inout) :: x(n)
            real(dp), intent(out) :: fvec(n)
        end subroutine hybrd1

        subroutine hybrj(fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
                         factor, nprint, info, nfev, njev, r, lr, qtf, wa1, wa2, &
                         wa3, wa4)
            import :: dp
            external :: fcn
            integer, intent(in) :: n, ldfjac, maxfev, mode, nprint, lr
            integer, intent(out) :: info, nfev, njev
            real(dp), intent(in) :: xtol, factor
            real(dp), intent(inout) :: x(n), diag(n)
            real(dp), intent(out) :: fvec(n), fjac(ldfjac, n), r(lr), qtf(n)
            real(dp), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(n)
        end subroutine hybrj

        subroutine lmder(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
                         diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, &
                         wa1, wa2, wa3, wa4)
            import :: dp
            external :: fcn
            integer, intent(in) :: m, n, ldfjac, maxfev, mode, nprint
            integer, intent(out) :: info, nfev, njev
            integer, intent(out) :: ipvt(n)
            real(dp), intent(in) :: ftol, xtol, gtol, factor
            real(dp), intent(inout) :: x(n), diag(n)
            real(dp), intent(out) :: fvec(m), fjac(ldfjac, n), qtf(n)
            real(dp), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(m)
        end subroutine lmder

        subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, &
                         mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf, &
                         wa1, wa2, wa3, wa4)
            import :: dp
            external :: fcn
            integer, intent(in) :: m, n, maxfev, mode, nprint, ldfjac
            integer, intent(out) :: info, nfev
            integer, intent(out) :: ipvt(n)
            real(dp), intent(in) :: ftol, xtol, gtol, epsfcn, factor
            real(dp), intent(inout) :: x(n), diag(n)
            real(dp), intent(out) :: fvec(m), fjac(ldfjac, n), qtf(n)
            real(dp), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(m)
        end subroutine lmdif

        subroutine lmstr(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
                         diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, &
                         wa1, wa2, wa3, wa4)
            import :: dp
            external :: fcn
            integer, intent(in) :: m, n, ldfjac, maxfev, mode, nprint
            integer, intent(out) :: info, nfev, njev
            integer, intent(out) :: ipvt(n)
            real(dp), intent(in) :: ftol, xtol, gtol, factor
            real(dp), intent(inout) :: x(n), diag(n)
            real(dp), intent(out) :: fvec(m), fjac(ldfjac, n), qtf(n)
            real(dp), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(m)
        end subroutine lmstr
    end interface

end module minpack_interfaces
