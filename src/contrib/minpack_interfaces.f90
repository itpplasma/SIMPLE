module minpack_interfaces
  implicit none

  interface
    function enorm(n, x)
      integer, intent(in) :: n
      real(8), intent(in) :: x(n)
      real(8) :: enorm
    end function enorm

    subroutine dogleg(n, r, lr, diag, qtb, delta, x)
      integer, intent(in) :: n, lr
      real(8), intent(in) :: r(lr), diag(n), qtb(n), delta
      real(8), intent(out) :: x(n)
    end subroutine dogleg

    subroutine fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn)
      external :: fcn
      integer, intent(in) :: n, ldfjac, ml, mu
      integer, intent(inout) :: iflag
      real(8), intent(in) :: x(n), epsfcn
      real(8), intent(inout) :: fvec(n), fjac(ldfjac,n)
    end subroutine fdjac1

    subroutine fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn, wa)
      external :: fcn
      integer, intent(in) :: m, n, ldfjac
      integer, intent(inout) :: iflag
      real(8), intent(in) :: x(n), epsfcn
      real(8), intent(inout) :: fvec(m), fjac(ldfjac,n), wa(m)
    end subroutine fdjac2

    subroutine qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm)
      integer, intent(in) :: m, n, lda, lipvt
      logical, intent(in) :: pivot
      real(8), intent(inout) :: a(lda,n)
      integer, intent(out) :: ipvt(lipvt)
      real(8), intent(out) :: rdiag(n), acnorm(n)
    end subroutine qrfac

    subroutine qform(m, n, q, ldq)
      integer, intent(in) :: m, n, ldq
      real(8), intent(inout) :: q(ldq,m)
    end subroutine qform

    subroutine qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)
      integer, intent(in) :: n, ldr
      integer, intent(in) :: ipvt(n)
      real(8), intent(inout) :: r(ldr,n), diag(n), qtb(n)
      real(8), intent(out) :: x(n), sdiag(n)
      real(8), intent(inout) :: wa(n)
    end subroutine qrsolv

    subroutine r1mpyq(m, n, a, lda, v, w)
      integer, intent(in) :: m, n, lda
      real(8), intent(inout) :: a(lda,n)
      real(8), intent(in) :: v(n), w(n)
    end subroutine r1mpyq

    subroutine r1updt(m, n, s, ls, u, v, w, sing)
      integer, intent(in) :: m, n, ls
      real(8), intent(inout) :: s(ls), u(m), v(n)
      real(8), intent(out) :: w(m)
      logical, intent(out) :: sing
    end subroutine r1updt

    subroutine rwupdt(n, r, ldr, w, b, alpha, cos, sin)
      integer, intent(in) :: n, ldr
      real(8), intent(inout) :: r(ldr,n), w(n), b(n)
      real(8), intent(inout) :: alpha
      real(8), intent(out) :: cos(n), sin(n)
    end subroutine rwupdt

    subroutine lmpar(n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag, wa1, wa2)
      integer, intent(in) :: n, ldr
      integer, intent(in) :: ipvt(n)
      real(8), intent(inout) :: r(ldr,n), diag(n), qtb(n), delta, par
      real(8), intent(out) :: x(n), sdiag(n)
      real(8), intent(inout) :: wa1(n), wa2(n)
    end subroutine lmpar

    subroutine hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
                     factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf, wa1, &
                     wa2, wa3, wa4)
      external :: fcn
      integer, intent(in) :: n, maxfev, ml, mu, mode, nprint, ldfjac, lr
      integer, intent(out) :: info, nfev
      real(8), intent(in) :: xtol, epsfcn, factor
      real(8), intent(inout) :: x(n), diag(n)
      real(8), intent(out) :: fvec(n), fjac(ldfjac,n), r(lr), qtf(n)
      real(8), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(n)
    end subroutine hybrd

    subroutine hybrd1(fcn, n, x, fvec, tol, info)
      external :: fcn
      integer, intent(in) :: n
      integer, intent(out) :: info
      real(8), intent(in) :: tol
      real(8), intent(inout) :: x(n)
      real(8), intent(out) :: fvec(n)
    end subroutine hybrd1

    subroutine hybrj(fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
                     factor, nprint, info, nfev, njev, r, lr, qtf, wa1, wa2, &
                     wa3, wa4)
      external :: fcn
      integer, intent(in) :: n, ldfjac, maxfev, mode, nprint, lr
      integer, intent(out) :: info, nfev, njev
      real(8), intent(in) :: xtol, factor
      real(8), intent(inout) :: x(n), diag(n)
      real(8), intent(out) :: fvec(n), fjac(ldfjac,n), r(lr), qtf(n)
      real(8), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(n)
    end subroutine hybrj

    subroutine lmder(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
                     diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, &
                     wa1, wa2, wa3, wa4)
      external :: fcn
      integer, intent(in) :: m, n, ldfjac, maxfev, mode, nprint
      integer, intent(out) :: info, nfev, njev
      integer, intent(out) :: ipvt(n)
      real(8), intent(in) :: ftol, xtol, gtol, factor
      real(8), intent(inout) :: x(n), diag(n)
      real(8), intent(out) :: fvec(m), fjac(ldfjac,n), qtf(n)
      real(8), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(m)
    end subroutine lmder

    subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, &
                     mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf, &
                     wa1, wa2, wa3, wa4)
      external :: fcn
      integer, intent(in) :: m, n, maxfev, mode, nprint, ldfjac
      integer, intent(out) :: info, nfev
      integer, intent(out) :: ipvt(n)
      real(8), intent(in) :: ftol, xtol, gtol, epsfcn, factor
      real(8), intent(inout) :: x(n), diag(n)
      real(8), intent(out) :: fvec(m), fjac(ldfjac,n), qtf(n)
      real(8), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(m)
    end subroutine lmdif

    subroutine lmstr(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
                     diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, &
                     wa1, wa2, wa3, wa4)
      external :: fcn
      integer, intent(in) :: m, n, ldfjac, maxfev, mode, nprint
      integer, intent(out) :: info, nfev, njev
      integer, intent(out) :: ipvt(n)
      real(8), intent(in) :: ftol, xtol, gtol, factor
      real(8), intent(inout) :: x(n), diag(n)
      real(8), intent(out) :: fvec(m), fjac(ldfjac,n), qtf(n)
      real(8), intent(out) :: wa1(n), wa2(n), wa3(n), wa4(m)
    end subroutine lmstr
  end interface

end module minpack_interfaces
