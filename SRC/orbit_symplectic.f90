module orbit_symplectic

use parmot_mod, only : ro0_parmot => ro0
use field_can_mod, only: field_can, d_field_can, d2_field_can, eval_field, &
  get_val, get_derivatives, get_derivatives2, &
  H, pth, vpar, dvpar, dH, dpth, f, df, d2f, ro0, mu, d2pth, d2H

implicit none
save

double precision, parameter :: atol = 1e-15, rtol = 1e-7

! Current phase-space coordinates z and old pth
!$omp threadprivate(z, pthold)
double precision, dimension(4) :: z  ! z = (r, th, ph, pphi)
double precision :: pthold

! Buffer for Lagrange polynomial interpolation
integer, parameter :: nlag = 3 ! order
integer, parameter :: nbuf = 16 ! values to store back
!$omp threadprivate(kbuf, kt, k, bufind, zbuf, coef)
integer :: kbuf = 0
integer :: kt = 0
integer :: k = 0
integer :: bufind(0:nlag)
double precision, dimension(4, nbuf) :: zbuf
double precision, dimension(0:0, nlag+1) :: coef

! Timestep and variables from z0
!$omp threadprivate(ntau, dt, alambd, pabs)
integer :: ntau
double precision :: dt
double precision :: alambd, pabs

logical, parameter :: extrap_field = .true.

interface orbit_timestep_sympl
  module procedure orbit_timestep_sympl_euler1
end interface

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init(z0, ntau_init)
!
  double precision, intent(in) :: z0(5)
  integer, intent(in) :: ntau_init ! dtau/dtaumin, must be integer

  ntau = ntau_init

  call eval_field(z0(1), z0(2), z0(3), 0)

  pabs=z0(4)
  alambd=z0(5)

  mu = .5d0*pabs**2*(1.d0-alambd**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
  ro0 = ro0_parmot/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
  vpar = pabs*alambd*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

  z(1:3) = z0(1:3)  ! r, th, ph
  z(4) = vpar*f%hph + f%Aph/ro0 ! pphi
  pth = vpar*f%hth + f%Ath/ro0 ! pth

  call plag_coeff(nlag+1, 0, 1d0, 1d0*(/(k,k=-nlag,0)/), coef)
end subroutine orbit_sympl_init


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_sympl_euler1(n, x, fvec, iflag)
!
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(x(1), z(2), z(3), 2)
  call get_derivatives2(x(2))

  fvec(1) = dpth(1)*(pth - pthold) + dt*(dH(2)*dpth(1) - dH(1)*dpth(2))
  fvec(2) = dpth(1)*(x(2) - z(4))  + dt*(dH(3)*dpth(1) - dH(1)*dpth(3))

  !print *, fvec

end subroutine f_sympl_euler1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_sympl_euler1(x, jac)
!
  double precision, intent(in)  :: x(2)
  double precision, intent(out) :: jac(2, 2)

  jac(1,1) = d2pth(1)*(pth - pthold) + dpth(1)**2 &
    + dt*(d2H(2)*dpth(1) + dH(2)*d2pth(1) - d2H(1)*dpth(2) - dH(1)*d2pth(2))
  jac(1,2) = d2pth(7)*(pth - pthold) + dpth(1)*dpth(4) &
    + dt*(d2H(8)*dpth(1) + dH(2)*d2pth(7) - d2H(7)*dpth(2) - dH(1)*d2pth(8))
  jac(2,1) = d2pth(1)*(x(2) - z(4)) &
    + dt*(d2H(3)*dpth(1) + dH(3)*d2pth(1) - d2H(1)*dpth(3) - dH(1)*d2pth(3))
  jac(2,2) = d2pth(7)*(x(2) - z(4)) + dpth(1) &
    + dt*(d2H(9)*dpth(1) + dH(3)*d2pth(7) - d2H(7)*dpth(3) - dH(1)*d2pth(9))

end subroutine jac_sympl_euler1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine f_sympl_euler2(n, x, fvec, iflag)
!
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(x(1), x(2), x(3), 2)
  call get_derivatives2(z(4))

  fvec(1) = pth - pthold
  fvec(2) = dpth(1)*(x(2) - z(2)) - dt*dH(1)
  fvec(3) = dpth(1)*f%hph*(x(3) - z(3)) - dt*(dpth(1)*vpar - dH(1)*f%hth)

end subroutine f_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine jac_sympl_euler2(x, jac)
!
  double precision, intent(in)  :: x(3)
  double precision, intent(out) :: jac(3, 3)

  jac(1,:) = dpth(1:3)
  jac(2,:) = d2pth(1:3)*(x(2) - z(2)) - dt*d2H(1:3)
  jac(2,2) = jac(2,2) + dpth(1) 
  jac(3,:) = (d2pth(1:3)*f%hph + dpth(1)*df%dhph)*(x(3) - z(3)) &
    - dt*(d2pth(1:3)*vpar + dpth(1)*dvpar(1:3) - d2H(1)*f%hth - dH(1)*df%dhth)
  jac(3,3) = jac(3,3) + dpth(1)*f%hph 

end subroutine jac_sympl_euler2

subroutine newton1(x, atol, rtol, maxit, xlast)
  integer, parameter :: n = 2

  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(n)

  double precision :: fvec(n), fjac(n,n), ijac(n,n)
  integer :: kit

  do kit = 1, maxit
    if(x(1) < 0.0 .or. x(1) > 1.0) return
    call f_sympl_euler1(n, x, fvec, 1)
    call jac_sympl_euler1(x, fjac)
    ijac(1,1) = fjac(2,2)
    ijac(1,2) = -fjac(1,2)
    ijac(2,1) = -fjac(2,1)
    ijac(2,2) = fjac(1,1)
    ijac = ijac/(fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1))
    xlast = x
    x = x - matmul(ijac, fvec)
    if (all(abs(fvec) < atol) &
      .or. all(abs(x-xlast)/(abs(x)*(1d0+1d-30)) < rtol)) return
  enddo
  print *, 'Warning: maximum iterations reached in newton1: ', maxit
end subroutine

subroutine newton2(x, atol, rtol, maxit, xlast)
  integer, parameter :: n = 3
  integer :: kit

  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  double precision, intent(out) :: xlast(n)

  double precision :: fvec(n), fjac(n,n)
  integer :: pivot(n), info

  do kit = 1, maxit
    call f_sympl_euler2(n, x, fvec, 1)
    call jac_sympl_euler2(x, fjac)
    call dgesv(n, 1, fjac, n, pivot, fvec, 3, info) 
    xlast = x
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    if (all(abs(fvec) < atol) &
      .or. all(abs(x-xlast)/(abs(x)*(1d0+1d-30)) < rtol)) return
  enddo
  print *, 'Warning: maximum iterations reached in newton2: ', maxit
end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_euler1(z0, dtau, dtaumin, ierr)
!
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  integer, parameter :: n = 2
  integer, parameter :: maxit = 100

  double precision, dimension(n) :: x, xlast

  double precision :: tau2
  integer :: ktau
  
  ierr = 0

  if (z0(1) > 1.0) then
    ierr = 1
    return
  end if

  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0d0
  ktau = 0
  do while(ktau .lt. ntau)
    pthold = pth

    ! Initial guess with Lagrange extrapolation
    do k=0, nlag
      bufind(k) = kbuf-nlag+k
      if (bufind(k)<1) bufind(k) = bufind(k) + nbuf 
    end do
    
    if (nlag>0 .and. kt>nlag) then
      x(1)=sum(zbuf(1,bufind)*coef(0,:))
      x(2)=sum(zbuf(4,bufind)*coef(0,:))
    else
      x(1)=z(1)
      x(2)=z(4)
    end if

    ! correct if Lagrange messed up
    if (x(1) < 0.0 .or. x(1) > 1.0) then
      x(1) = z(1) 
      x(2) = z(4)
    end if
      
    call newton1(x, atol, rtol, maxit, xlast)
    if (x(1) < 0.0 .or. x(1) > 1.0) then
      ierr = 1
      return
    end if
    
    z(1) = x(1)
    z(4) = x(2)

    if (extrap_field) then
      dH(1) = dH(1) + d2H(1)*(x(1)-xlast(1)) + d2H(7)*(x(2)-xlast(2))
      dpth(1) = dpth(1) + d2pth(1)*(x(1)-xlast(1)) + d2pth(7)*(x(2)-xlast(2))
      vpar = vpar + dvpar(1)*(x(1)-xlast(1)) + dvpar(4)*(x(2)-xlast(2))
      f%hth = f%hth + df%dhth(1)*(x(1)-xlast(1))
      f%hph = f%hph + df%dhph(1)*(x(1)-xlast(1))
    else
      call eval_field(z(1), z(2), z(3), 0)
      call get_derivatives(z(4))
    endif

    z(2) = z(2) + dt*dH(1)/dpth(1)
    z(3) = z(3) + dt*(vpar - dH(1)/dpth(1)*f%hth)/f%hph

    kbuf = mod(kt, nbuf) + 1
    zbuf(:,kbuf) = z

    tau2 = tau2 + dtaumin
    kt = kt+1
    ktau = ktau+1
  enddo
  z0(1:3) = z(1:3)
  !z0(4) = H
  z0(5) = vpar/(pabs*dsqrt(2d0))  ! alambda

end subroutine orbit_timestep_sympl_euler1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_euler2(z0, dtau, dtaumin, ierr)
!
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  integer, parameter :: n = 3
  integer, parameter :: maxit = 100

  double precision, dimension(n) :: x, xlast

  double precision :: tau2

  integer :: ktau
  
  ierr = 0

  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0d0
  ktau = 0
  do while(ktau .lt. ntau)
    pthold = pth

    ! Initial guess with Lagrange extrapolation
    do k=0, nlag
      bufind(k) = kbuf-nlag+k
      if (bufind(k)<1) bufind(k) = bufind(k) + nbuf 
    end do
    
    if (nlag>0 .and. kt>nlag) then
      x(1)=sum(zbuf(1,bufind)*coef(0,:))
      x(2)=sum(zbuf(2,bufind)*coef(0,:))
      x(3)=sum(zbuf(3,bufind)*coef(0,:))
    else
      x = z(1:3)
    end if
    
    call newton2(x, atol, rtol, maxit, xlast)

    z(1:3) = x

    call eval_field(z(1), z(2), z(3), 0)
    call get_derivatives(z(4))

    pth = pthold - dt*(dH(2) - dH(1)*dpth(2)/dpth(1))
    z(4) = z(4) - dt*(dH(3) - dH(1)*dpth(3)/dpth(1))

    kbuf = mod(kt, nbuf) + 1
    zbuf(:,kbuf) = z

    tau2 = tau2 + dtaumin
    kt = kt+1
    ktau = ktau+1
  enddo
  z0(1:3) = z(1:3)
  !z0(4) = H
  z0(5) = vpar/(pabs*dsqrt(2d0))  ! alambda

end subroutine orbit_timestep_sympl_euler2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl_verlet(z0, dtau, dtaumin, ierr)
!
  integer, intent(out) :: ierr
  double precision, dimension(5), intent(inout) :: z0
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 2
  integer, parameter :: maxit = 100

  double precision, dimension(n1) :: x1, xlast1
  double precision, dimension(n2) :: x2, xlast2

  double precision :: tau2

  integer :: ktau

  ierr = 0

  dt = 0.5d0*dtaumin/dsqrt(2d0)
  tau2 = 0d0
  ktau = 0
  do while(ktau .lt. ntau)
    pthold = pth

    ! Verlet part one: impl/expl Euler

    ! Initial guess with Lagrange extrapolation
    
    if (nlag>0 .and. kt>nlag) then
      do k=0, nlag
        bufind(k) = kbuf-2*nlag+2*k-1
        if (bufind(k)<1) bufind(k) = bufind(k) + nbuf 
      end do
  !    print *, bufind
      x1(1)=sum(zbuf(1,bufind)*coef(0,:))
      x1(2)=sum(zbuf(2,bufind)*coef(0,:))
      x1(3)=sum(zbuf(3,bufind)*coef(0,:))
    else
      x1 = z(1:3)
    end if
    
    call newton2(x1, atol, rtol, maxit, xlast1)

    z(1:3) = x1

    call eval_field(z(1), z(2), z(3), 0)
    call get_derivatives(z(4))

    pthold = pthold - dt*(dH(2) - dH(1)*dpth(2)/dpth(1))
    z(4) = z(4) - dt*(dH(3) - dH(1)*dpth(3)/dpth(1))

    kbuf = mod(2*kt, nbuf) + 1
    zbuf(:,kbuf) = z

    ! Verlet part two: expl/impl Euler

    ! Initial guess with Lagrange extrapolation
    
    if (nlag>0 .and. kt>nlag) then
      do k=0, nlag
        bufind(k) = kbuf-2*nlag+2*k
        if (bufind(k)<1) bufind(k) = bufind(k) + nbuf 
      end do
      x2(1)=sum(zbuf(1,bufind)*coef(0,:))
      x2(2)=sum(zbuf(4,bufind)*coef(0,:))
    else
      x2(1)=z(1)
      x2(2)=z(4)
    end if
    
    call newton1(x2, atol, rtol, maxit, xlast2)

    z(1) = x2(1)
    z(4) = x2(2)

    call eval_field(z(1), z(2), z(3), 0)
    call get_derivatives(z(4))
    z(2) = z(2) + dt*dH(1)/dpth(1)
    z(3) = z(3) + dt*(vpar - dH(1)/dpth(1)*f%hth)/f%hph

    kbuf = mod(2*kt+1, nbuf) + 1
    zbuf(:,kbuf) = z
    
    tau2 = tau2 + dtaumin
    kt = kt+1
    ktau = ktau+1
  enddo
  z0(1:3) = z(1:3)
  !z0(4) = H
  z0(5) = vpar/(pabs*dsqrt(2d0))  ! alambda

end subroutine orbit_timestep_sympl_verlet

subroutine debug_root(x0)
  double precision :: x0(2), x(2)
  integer :: k, l, iflag
  integer, parameter :: n = 100
  double precision, parameter :: eps = 1e-15

  double precision :: fvec(2)

  do k = -n,n
    do l = -n,n
      x = x0 + l*eps/n*(/x0(1),0d0/) + k*eps/n*(/0d0,x0(2)/) 
      call f_sympl_euler1(2, x, fvec, iflag)
      write(5001,*) x(1), x(2), fvec
    end do
  end do

end subroutine debug_root

end module orbit_symplectic