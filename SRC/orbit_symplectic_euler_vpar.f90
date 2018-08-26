! TODO: check if everything runs correctly here

module orbit_symplectic

use parmot_mod, only : ro0_parmot => ro0
use field_can_mod, only: field_can, d_field_can, d2_field_can, eval_field

implicit none

double precision, dimension(6) :: yold
double precision :: rmumag
double precision :: dt, ro0

double precision :: coala
double precision :: derphi(3)
double precision :: alambd, pabs

type(field_can) :: f
type(d_field_can) :: df
type(d2_field_can) :: d2f

contains

subroutine orbit_sympl_init(z)
  double precision, intent(in) :: z(5)

  call eval_field(z(1), z(2), z(3), 0, f, df, d2f)

  pabs=z(4)
  alambd=z(5)

  rmumag = .5d0*pabs**2*(1.d0-alambd**2)/f%Bmod*2d0 ! rmumag by factor 2 from other modules
  ro0 = ro0_parmot/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules

  yold(1:3) = z(1:3) ! r, theta, varphi
  yold(4) = pabs*alambd*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules
  yold(5) = yold(4)*f%Bth/f%Bmod + f%Ath/ro0
  yold(6) = yold(4)*f%Bph/f%Bmod + f%Aph/ro0

end subroutine orbit_sympl_init

subroutine f_sympl_euler(n, y, fvec, iflag)
  integer, parameter :: mode = 2

  integer, intent(in) :: n
  double precision, intent(in) :: y(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(out) :: iflag

  double precision, dimension(2) :: q, p, w
  double precision, dimension(2) :: dqdt, dpdt, dHdq, dHdw, pqw
  double precision, dimension(2,2) :: dwdq, dwdp

  if (mode==1) then
    q = y(2:3)
    p = yold(5:6)
    w = y((/1,4/))
  else if (mode==2) then ! TODO: debug mode 2
    q = yold(2:3)
    p = y(5:6)
    w = y((/1,4/))
  end if

  dqdt = (y(2:3)-yold(2:3))/dt
  dpdt = (y(5:6)-yold(5:6))/dt

  call get_derivatives(w(1), q(1), q(2), w(2), dHdq, dHdw, pqw, dwdq, dwdp)
  call step_sympl(q, p, w, dqdt, dpdt, dHdq, dHdw, pqw, dwdq, dwdp, fvec)
!  write(5001,*) y
!  write(5002,*) fvec
  !fvec = fvec + exp(-w(1)*1e2) ! to remove negative radius
end

subroutine step_sympl(q, p, w, dqdt, dpdt, dHdq, dHdw, pqw, dwdq, dwdp, ret)
  double precision, intent(out) :: ret(:)
  double precision, intent(in) :: q(:), p(:), w(:)
  double precision, intent(in) :: dqdt(:), dpdt(:)
  double precision, intent(in) :: dHdq(:), dHdw(:)
  double precision, intent(in) :: pqw(:)
  double precision, intent(in) :: dwdq(:,:), dwdp(:,:)

  integer :: k, nq

  nq = size(q)

  ret = 0d0
  do k = 1, nq
    ret(k) = sum(dHdw*dwdp(:,k)) - dqdt(k)
    ret(k + nq) = -sum(dHdw*dwdq(:,k)) - dHdq(k) - dpdt(k)
    ret(k + 2*nq) = p(k) - pqw(k)
  end do

end subroutine step_sympl

subroutine get_derivatives(r, th, ph, vpar, dHdq, dHdw, pqw, dwdq, dwdp)
  double precision, intent(in) :: r, th, ph, vpar
  double precision, dimension(2), intent(out) :: dHdq, dHdw, pqw
  double precision, dimension(2,2), intent(out) :: dwdq, dwdp
  double precision, dimension(3) :: x
  double precision, dimension(2,2) :: dpdw, dpdq

  call eval_field(r, th, ph, 0, f, df, d2f)

  x(1) = r
  x(2) = th
  x(3) = ph
  call elefie_can(x, derphi)

  dHdq(1) = rmumag*df%dBmod(2) + derphi(2)
  dHdq(2) = rmumag*df%dBmod(3) + derphi(3)
  dHdw(1) = rmumag*df%dBmod(1) + derphi(1)
  dHdw(2) = vpar

  pqw(1) = vpar*f%Bth/f%Bmod + f%Ath/ro0
  pqw(2) = vpar*f%Bph/f%Bmod + f%Aph/ro0

  dpdq(1,1) = vpar*(df%dBth(2) - f%Bth*df%dBmod(2)/f%Bmod)/f%Bmod + df%dAth(2)/ro0
  dpdq(1,2) = vpar*(df%dBth(3) - f%Bth*df%dBmod(3)/f%Bmod)/f%Bmod + df%dAth(3)/ro0
  dpdq(2,1) = vpar*(df%dBph(2) - f%Bph*df%dBmod(2)/f%Bmod)/f%Bmod + df%dAph(2)/ro0
  dpdq(2,2) = vpar*(df%dBph(3) - f%Bph*df%dBmod(3)/f%Bmod)/f%Bmod + df%dAph(3)/ro0

  dpdw(1,1) = vpar*(df%dBth(1) - f%Bth*df%dBmod(1)/f%Bmod)/f%Bmod + df%dAth(1)/ro0
  dpdw(1,2) = f%Bth/f%Bmod
  dpdw(2,1) = vpar*(df%dBph(1) - f%Bth*df%dBmod(1)/f%Bmod)/f%Bmod + df%dAph(1)/ro0
  dpdw(2,2) = f%Bph/f%Bmod

  ! inverse matrix
  dwdp(1,1) = dpdw(2,2)
  dwdp(1,2) = -dpdw(1,2)
  dwdp(2,1) = -dpdw(2,1)
  dwdp(2,2) = dpdw(1,1)
  dwdp = dwdp/(dpdw(1,1)*dpdw(2,2)-dpdw(1,2)*dpdw(2,1))

  ! matrix multiplication
  dwdq = -matmul(dwdp,dpdq)

end subroutine get_derivatives

subroutine step_forward(y)
  double precision, dimension(6), intent(inout) :: y
  double precision, dimension(2) :: dqdt, dpdt, dHdq, dHdw, pqw
  double precision, dimension(2,2) :: dwdq, dwdp

  double precision, dimension(2) :: dHdq_can, dHdp_can

  !call eval_field(y(1), y(2), y(3), 0, f, df, d2f)
  call get_derivatives(y(1), y(2), y(3), y(4), dHdq, dHdw, pqw, dwdq, dwdp)

  ! dH(q,p)/dq
  dHdq_can(1) = dHdq(1) + sum(dHdw*dwdq(:,1)) 
  dHdq_can(2) = dHdq(2) + sum(dHdw*dwdq(:,2))
  dHdp_can(1) = sum(dHdw*dwdp(:,1)) 
  dHdp_can(2) = sum(dHdw*dwdp(:,2))

  y(1) = y(1) + sum(dwdq(1,:)*dHdp_can - dwdp(1,:)*dHdq_can)*dt
  y(4) = y(4) + sum(dwdq(2,:)*dHdp_can - dwdp(2,:)*dHdq_can)*dt

  ! for using new r and vpar for angle timestep (pseudo-symplectic Euler)
  !call eval_field(y(1), y(2), y(3), 0, f, df, d2f)
  call get_derivatives(y(1), y(2), y(3), y(4), dHdq, dHdw, pqw, dwdq, dwdp)
  
  y(2) = y(2) + sum(dHdw*dwdp(:,1))*dt
  y(3) = y(3) + sum(dHdw*dwdp(:,2))*dt
  
end subroutine step_forward

subroutine orbit_timestep_sympl(z, dtau, dtaumin, ierr)

  implicit none
  
  integer, parameter :: ndim = 5
  integer, parameter :: n = 6
  double precision :: tol

  integer, intent(out) :: ierr
  double precision, dimension(ndim), intent(inout) :: z
  double precision, intent(in) :: dtau
  double precision, intent(in) :: dtaumin

  double precision, dimension(6) :: y, fvec

  double precision tau2
  
  ierr = 0

  ! for nleq1

  ! INTEGER NMAX, LIOPT, LIWK, LRWK, LUPRT
  ! PARAMETER ( NMAX=50 )
  ! PARAMETER ( LIOPT=50, LIWK=NMAX+50, LRWK=(NMAX+13)*NMAX+60 )
  ! PARAMETER ( LUPRT=6 )
  ! EXTERNAL FCN, NLEQ1
  ! DOUBLE PRECISION XSCAL(NMAX)
  ! INTEGER IOPT(LIOPT), IWK(LIWK)
  ! DOUBLE PRECISION RWK(LRWK)
  ! INTEGER I, NIW, NRW, DUMMY
  ! CHARACTER CHGDAT*20, PRODCT*8

  ! end for nleq1
  dt = dtaumin/dsqrt(2d0) ! factor 1/sqrt(2) due to velocity normalisation different from other modules
  tau2 = 0.0
  do while(tau2.lt.dtau)
    y = yold
    call step_forward(y)
    tol = 1d-10
    call hybrd1 (f_sympl_euler, n, y, fvec, tol, ierr)
    !IOPT = 0
    !iopt(11) = 1 ! print error messages
    !IWK = 0
    !iwk(31) = 1000 ! maximum iterations, default: 50
    !RWK = 0.0D0
    !XSCAL = 0.0D0
    !CALL nleq1(n,f_sympl_euler,dummy,y,XSCAL,tol,IOPT,ierr,LIWK,IWK,LRWK,RWK)
    !write(5001,*) '---------------------------------------------'
    !write(5002,*) '---------------------------------------------'
    if(ierr > 1) then
      print *, 'orbit_timestep_sympl: warning ', ierr, ', step at half accuracy'
      tol = 1d-5
      ierr = 0
      !y = yold
      !call step_forward(y)
      call hybrd1 (f_sympl_euler, n, y, fvec, tol, ierr)
      if(ierr > 1) stop
      !if(ierr /= 4 .and. ierr /= 5) stop
    end if
    yold = y
    tau2=tau2+dtaumin
  enddo
  z = 0.0
  z(1:4) = y(1:4)
  z(5) = y(4)**2/2d0 + rmumag*f%Bmod
end

end module orbit_symplectic
