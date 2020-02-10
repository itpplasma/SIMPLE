subroutine odeint_allroutines(y, nvar, x1, x2, eps, derivs)
    implicit none

    external :: derivs
    integer, intent(in) :: nvar
    double precision, intent(in) :: x1, x2, eps
    double precision, dimension(nvar) :: y
    double precision, dimension(nvar) :: yp

    double precision :: epsabs = 1d-31
    integer :: flag = 1

    call r8_rkf45 ( derivs, nvar, y, yp, x1, x2, eps, epsabs, flag )
    ! if (flag > 2)  print *, flag
end subroutine odeint_allroutines
