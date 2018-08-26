! TODO: fix timestep in this file

module orbit_symplectic

use parmot_mod, only : ro0_parmot => ro0
use field_can_mod, only: field_can, d_field_can, d2_field_can, eval_field, &
  get_val, get_derivatives, H, pth, vpar, dvpar, dH, dpth, f, df, d2f, ro0, mu

implicit none

double precision :: pthold

double precision :: dt

double precision :: coala
double precision :: derphi(3)
double precision :: alambd, pabs

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_sympl_init(z)
!
  double precision, intent(in) :: z(5)

end subroutine orbit_sympl_init

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_sympl(z, dtau, dtaumin, ierr)

  implicit none
  
  integer, parameter :: ndim = 5
  integer, parameter :: n = 2
  double precision :: tol

  double precision, dimension(ndim), intent(inout) :: z
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin
  integer, intent(out) :: ierr

  double precision, dimension(2) :: fvec

  double precision :: tau2

  double precision :: Bctr(3), Bstarctr(3), curlh(3), vz(5), Bstarpar, pphi
  
  ierr = 0

  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0d0
  do while(tau2.lt.dtau)
    pabs=z(4)
    alambd=z(5)

    call eval_field(z(1), z(2), z(3), 0)

    mu = .5d0*pabs**2*(1.d0-alambd**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
    ro0 = ro0_parmot/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
    vpar = pabs*alambd*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules
    pphi = f%Bph/f%Bmod*vpar + f%Aph/ro0

    call get_derivatives(pphi)
    Bctr(1) = 0d0
    Bctr(2) = -df%dAph(1)
    Bctr(3) = df%dAth(1)

    curlh(1) = (df%dBph(2) - df%dBth(3))/f%Bmod &
             - (f%Bph*df%dBmod(2) - f%Bth*df%dBmod(3))/f%Bmod**2
    curlh(2) = -df%dBph(1)/f%Bmod + df%dBmod(1)*f%Bph/f%Bmod**2
    curlh(3) =  df%dBth(1)/f%Bmod - df%dBmod(1)*f%Bth/f%Bmod**2

    Bstarctr = Bctr + vpar*ro0*curlh
    Bstarpar = (Bstarctr(2) * f%Bth + Bstarctr(3) * f%Bph)/f%Bmod

    vz(1:3) = vpar*Bstarctr/Bstarpar
    vz(2) = vz(2) + ro0*f%Bph*mu*df%dBmod(2)/(f%Bmod*Bstarpar)
    vz(3) = vz(3) - ro0*f%Bth*mu*df%dBmod(3)/(f%Bmod*Bstarpar)
    vz(4) = 0d0 !... no electric field, momentum unchanged
    vz(5) = -(1d0-alambd**2)/alambd*sum(vz(1:3)*df%dBmod)/(2d0*f%Bmod)
    
    write(4004,*) z(1:3)
    write(4004,*) vpar/dsqrt(2d0)*Bstarctr/Bstarpar
    write(4004,*) vz/dsqrt(2d0)
    !write(4004,*) df%dAth
    !write(4004,*) df%dAph
    !write(4004,*) df%dBmod

    z = z + dt*vz
    tau2 = tau2 + dtaumin
  enddo
end

end module orbit_symplectic
