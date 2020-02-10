subroutine odeint_allroutines(y, nvar, x1, x2, eps, derivs)
    implicit none

    external :: derivs
    integer, intent(in) :: nvar
    double precision, intent(in) :: x1, x2, eps
    double precision, dimension(nvar) :: y
    double precision, dimension(nvar) :: yp

    double precision :: epsrel, epsabs
    integer :: flag

    flag = 1
    epsrel = eps
    epsabs = 1d-31

    call r8_rkf45 ( derivs, nvar, y, yp, x1, x2, epsrel, epsabs, flag )

    if (flag == 6) then
        epsrel = 10*epsrel
        epsabs = 10*epsabs
        flag = 2
        call r8_rkf45 ( derivs, nvar, y, yp, x1, x2, epsrel, epsabs, flag )
    elseif (flag == 7) then
        flag = 2
        call r8_rkf45 ( derivs, nvar, y, yp, x1, x2, epsrel, epsabs, flag )
    endif

end subroutine odeint_allroutines
