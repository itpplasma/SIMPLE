module orbit_symplectic

use parmot_mod, only : ro0_parmot => ro0
use field_can_mod, only: field_can, d_field_can, d2_field_can, eval_field

implicit none

double precision, dimension(2) :: qold, wold
double precision, dimension(2) :: q, w  ! q = (th, ph), w = (r, pphi)
double precision :: pthold

double precision :: mu
double precision :: dt, ro0

double precision :: coala
double precision :: derphi(3)
double precision :: alambd, pabs

double precision :: H, pth, vpar
double precision, dimension(4) :: dvpar, dH, dpth

type(field_can) :: f
type(d_field_can) :: df
type(d2_field_can) :: d2f

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init(z)
  double precision, intent(in) :: z(5)

  call eval_field(z(1), z(2), z(3), 0, f, df, d2f)

  pabs=z(4)
  alambd=z(5)

  mu = .5d0*pabs**2*(1.d0-alambd**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
  ro0 = ro0_parmot/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
  vpar = pabs*alambd*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

  q = z(2:3)  ! theta, varphi
  w(1) = z(1) ! r
  w(2) = vpar*f%Bph/f%Bmod + f%Aph/ro0 ! p_ph

  pth = vpar*f%Bth/f%Bmod + f%Ath/ro0 ! p_th
end subroutine orbit_sympl_init

subroutine f_sympl_euler(n, x, fvec, iflag)
  integer, parameter :: mode = 2

  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(out) :: iflag

  w = x
  call get_derivatives()

  fvec(1) = dpth(1)*(pth - pthold)   + dt*(dH(2)*dpth(1) - dH(1)*dpth(2))
  fvec(2) = dpth(1)*(w(2) - wold(2)) + dt*(dH(3)*dpth(1) - dH(1)*dpth(3))

  !print *, pth, pthold
  !print *, dH
  !print *, dpth
  !print *, dvpar
  !stop

end subroutine f_sympl_euler


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine get_val()
!
! computes values of H, pth and vpar at z=(r, th, ph, pphi)
!
!
  call eval_field(w(1), q(1), q(2), 0, f, df, d2f)

  vpar = f%Bmod/f%Bph*(w(2) - f%Aph/ro0)
  H = vpar**2/2d0 + mu*f%Bmod
  pth = f%Bth/f%Bmod*vpar + f%Ath/ro0
  
end subroutine get_val

    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine get_derivatives()
!
! computes H, pth and vpar at z=(r, th, ph, pphi) and their derivatives  
!
!
  call get_val()

  dvpar(1:3) = df%dBmod*(w(2) - f%Aph/ro0) - f%Bmod*df%dAph/ro0
  dvpar(1:3) = (dvpar(1:3) - df%dBph*vpar)/f%Bph
  dvpar(4)   = f%Bmod/f%Bph

  dH(1:3) = vpar*dvpar(1:3) + mu*df%dBmod
  dH(4)   = vpar*f%Bmod/f%Bph
  
  dpth(1:3) = ((w(2)-f%Aph/ro0)*df%dBth - df%dAph/ro0*f%Bth            &
               -(pth-f%Ath/ro0)*df%dBph + df%dAth/ro0*f%Bph)/f%Bph
  dpth(4) = f%Bth/f%Bph

end subroutine get_derivatives

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl(z, dtau, dtaumin, ierr)

  implicit none
  
  integer, parameter :: ndim = 5
  integer, parameter :: n = 2
  double precision :: tol

  integer, intent(out) :: ierr
  double precision, dimension(ndim), intent(inout) :: z
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  double precision, dimension(2) :: fvec

  double precision tau2
  
  ierr = 0

  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0d0
  do while(tau2.lt.dtau)
    qold = q
    wold = w
    pthold = pth
    ! TODO: initial guess with Lagrange
    tol = 1d-10
    call hybrd1 (f_sympl_euler, n, w, fvec, tol, ierr)
    if(ierr > 1) stop 'error in root finding'
    call get_derivatives()
    q(1) = qold(1) + dt*dH(1)/dpth(1)
    q(2) = qold(2) + dt*(vpar*f%Bmod - dH(1)/dpth(1)*f%Bth)/f%Bph
    tau2 = tau2 + dtaumin
  enddo
  z(1) = w(1)
  z(2:3) = q
  z(4) = H
  z(5) = vpar/sqrt(H) ! alambda
end

end module orbit_symplectic
