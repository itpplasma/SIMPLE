module orbit_symplectic

use util, only: pi, twopi
use field_can_mod, only: field_can_t, get_val, get_derivatives, get_derivatives2, &
  eval_field => evaluate
use orbit_symplectic_base, only: symplectic_integrator_t, multistage_integrator_t, &
  RK45, EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT, GAUSS1, GAUSS2, GAUSS3, GAUSS4, &
  LOBATTO3, S_MAX, orbit_timestep_sympl_i, extrap_field, &
  coeff_rk_gauss, coeff_rk_lobatto, f_rk_lobatto, current_particle_index
use orbit_symplectic_quasi, only: orbit_timestep_quasi, timestep_expl_impl_euler_quasi, &
  timestep_impl_expl_euler_quasi, timestep_midpoint_quasi, orbit_timestep_rk45, &
  timestep_rk_gauss_quasi, timestep_rk_lobatto_quasi
use vector_potentail_mod, only: torflux
use lapack_interfaces, only: dgesv

implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)

procedure(orbit_timestep_sympl_i), pointer :: orbit_timestep_sympl => null()

contains

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init(si, f, z, dt, ntau, rtol_init, mode_init)
  !
  use plag_coeff_sub, only : plag_coeff

  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  real(dp), intent(in) :: z(:)
  real(dp), intent(in) :: dt
  integer, intent(in) :: ntau
  real(dp), intent(in) :: rtol_init
  integer, intent(in) :: mode_init

  integer :: k

  si%atol = 1d-15
  si%rtol = rtol_init

  si%ntau = ntau
  si%dt = dt

  si%z = z

  call eval_field(f, z(1), z(2), z(3), 0)
  call get_val(f, si%z(4)) ! for pth
  si%pthold = f%pth

  select case (mode_init)
    case (RK45)
      orbit_timestep_quasi => orbit_timestep_rk45
    case (EXPL_IMPL_EULER)
      orbit_timestep_sympl => orbit_timestep_sympl_expl_impl_euler
      orbit_timestep_quasi => timestep_expl_impl_euler_quasi
    case (IMPL_EXPL_EULER)
      orbit_timestep_sympl => orbit_timestep_sympl_impl_expl_euler
      orbit_timestep_quasi => timestep_impl_expl_euler_quasi
    case (MIDPOINT)
      orbit_timestep_sympl => orbit_timestep_sympl_midpoint
      orbit_timestep_quasi => timestep_midpoint_quasi
    case (GAUSS1)
      orbit_timestep_sympl => orbit_timestep_sympl_gauss1
      orbit_timestep_quasi => orbit_timestep_quasi_gauss1
    case (GAUSS2)
      orbit_timestep_sympl => orbit_timestep_sympl_gauss2
      orbit_timestep_quasi => orbit_timestep_quasi_gauss2
    case (GAUSS3)
      orbit_timestep_sympl => orbit_timestep_sympl_gauss3
      orbit_timestep_quasi => orbit_timestep_quasi_gauss3
    case (GAUSS4)
      orbit_timestep_sympl => orbit_timestep_sympl_gauss4
      orbit_timestep_quasi => orbit_timestep_quasi_gauss4
    case (LOBATTO3)
      orbit_timestep_sympl => orbit_timestep_sympl_lobatto3
      orbit_timestep_quasi => orbit_timestep_quasi_lobatto3
    case default
      print *, 'invalid mode for orbit_timestep_sympl: ', mode_init
      error stop
  end select
end subroutine orbit_sympl_init


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_gauss1(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_gauss(si, f, 1, ierr)
end subroutine orbit_timestep_sympl_gauss1

recursive subroutine orbit_timestep_sympl_gauss2(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_gauss(si, f, 2, ierr)
end subroutine orbit_timestep_sympl_gauss2

recursive subroutine orbit_timestep_sympl_gauss3(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_gauss(si, f, 3, ierr)
end subroutine orbit_timestep_sympl_gauss3

recursive subroutine orbit_timestep_sympl_gauss4(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_gauss(si, f, 4, ierr)
end subroutine orbit_timestep_sympl_gauss4

recursive subroutine orbit_timestep_sympl_lobatto3(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_lobatto(si, f, 3, ierr)
end subroutine orbit_timestep_sympl_lobatto3

recursive subroutine orbit_timestep_quasi_gauss1(ierr)
  integer, intent(out) :: ierr
  call timestep_rk_gauss_quasi(1, ierr)
end subroutine orbit_timestep_quasi_gauss1

recursive subroutine orbit_timestep_quasi_gauss2(ierr)
  integer, intent(out) :: ierr
  call timestep_rk_gauss_quasi(2, ierr)
end subroutine orbit_timestep_quasi_gauss2

recursive subroutine orbit_timestep_quasi_gauss3(ierr)
  integer, intent(out) :: ierr
  call timestep_rk_gauss_quasi(3, ierr)
end subroutine orbit_timestep_quasi_gauss3

recursive subroutine orbit_timestep_quasi_gauss4(ierr)
  integer, intent(out) :: ierr
  call timestep_rk_gauss_quasi(4, ierr)
end subroutine orbit_timestep_quasi_gauss4

recursive subroutine orbit_timestep_quasi_lobatto3(ierr)
  integer, intent(out) :: ierr
  call timestep_rk_lobatto_quasi(3, ierr)
end subroutine orbit_timestep_quasi_lobatto3


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine f_sympl_euler1(si, f, n, x, fvec, iflag)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(in) :: n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(f, x(1), si%z(2), si%z(3), 2)
  call get_derivatives2(f, x(2))

  fvec(1) = f%dpth(1)*(f%pth - si%pthold) + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
  fvec(2) = f%dpth(1)*(x(2) - si%z(4))  + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))

end subroutine f_sympl_euler1


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_sympl_euler1(si, f, x, jac)
  !
  type(symplectic_integrator_t), intent(in) :: si
  type(field_can_t), intent(inout) :: f

  real(dp), intent(in)  :: x(2)
  real(dp), intent(out) :: jac(2, 2)

  jac(1,1) = f%d2pth(1)*(f%pth - si%pthold) + f%dpth(1)**2 &
    + si%dt*(f%d2H(2)*f%dpth(1) + f%dH(2)*f%d2pth(1) - f%d2H(1)*f%dpth(2) - f%dH(1)*f%d2pth(2))
  jac(1,2) = f%d2pth(7)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(4) &
    + si%dt*(f%d2H(8)*f%dpth(1) + f%dH(2)*f%d2pth(7) - f%d2H(7)*f%dpth(2) - f%dH(1)*f%d2pth(8))
  jac(2,1) = f%d2pth(1)*(x(2) - si%z(4)) &
    + si%dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
  jac(2,2) = f%d2pth(7)*(x(2) - si%z(4)) + f%dpth(1) &
    + si%dt*(f%d2H(9)*f%dpth(1) + f%dH(3)*f%d2pth(7) - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9))

end subroutine jac_sympl_euler1

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine f_sympl_euler2(si, f, n, x, fvec, iflag)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(in) :: n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(f, x(1), x(2), x(3), 2)
  call get_derivatives2(f, si%z(4))

  fvec(1) = f%pth - si%pthold
  fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
  fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)

end subroutine f_sympl_euler2


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_sympl_euler2(si, f, x, jac)
  !
  type(symplectic_integrator_t), intent(in) :: si
  type(field_can_t), intent(inout) :: f
  real(dp), intent(in)  :: x(3)
  real(dp), intent(out) :: jac(3, 3)

  jac(1,:) = f%dpth(1:3)
  jac(2,:) = f%d2pth(1:3)*(x(2) - si%z(2)) - si%dt*f%d2H(1:3)
  jac(2,2) = jac(2,2) + f%dpth(1)
  jac(3,:) = (f%d2pth(1:3)*f%hph + f%dpth(1)*f%dhph)*(x(3) - si%z(3)) &
    - si%dt*(f%d2pth(1:3)*f%vpar + f%dpth(1)*f%dvpar(1:3) - f%d2H(1:3)*f%hth - f%dH(1)*f%dhth)
  jac(3,3) = jac(3,3) + f%dpth(1)*f%hph

end subroutine jac_sympl_euler2


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine f_midpoint_part1(si, f, n, x, fvec)
  !
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    real(dp), intent(out) :: fvec(n)

    ! evaluate at midpoint
    call eval_field(f, x(5), 0.5d0*(x(2) + si%z(2)), 0.5d0*(x(3) + si%z(3)), 2)
    call get_derivatives2(f, 0.5d0*(x(4) + si%z(4)))

    fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
    fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)
    fvec(4) = f%dpth(1)*(x(4) - si%z(4)) + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
    fvec(5) = f%dpth(1)*(f%pth - si%pthold) + 0.5d0*si%dt*(f%dpth(1)*f%dH(2)-f%dpth(2)*f%dH(1))

  end subroutine f_midpoint_part1

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine f_midpoint_part2(si, f, n, x, fvec)
  !
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    real(dp), intent(out) :: fvec(n)

    real(dp) :: dpthmid, pthdotbar

    ! save evaluation from midpoint
    dpthmid = f%dpth(1)
    pthdotbar = f%dpth(1)*f%dH(2) - f%dpth(2)*f%dH(1)

    ! evaluate at endpoint
    call eval_field(f, x(1), x(2), x(3), 2)
    call get_derivatives2(f, x(4))
    fvec(1) = dpthmid*(f%pth - si%pthold) + si%dt*pthdotbar

  end subroutine f_midpoint_part2

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_midpoint_part1(si, f, x, jac)
  !
    type(symplectic_integrator_t), intent(in) :: si
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in)  :: x(5)
    real(dp), intent(out) :: jac(5, 5)

    jac(2,1) = 0d0
    jac(2,5) = f%d2pth(1)*(x(2) - si%z(2)) - si%dt*f%d2H(1)
    jac(2,2:3) = 0.5d0*(f%d2pth(2:3)*(x(2) - si%z(2)) - si%dt*f%d2H(2:3))
    jac(2,2) = jac(2,2) + f%dpth(1)
    jac(2,4) = 0.5d0*(f%d2pth(7)*(x(2) - si%z(2)) - si%dt*f%d2H(7))

    jac(3,1) = 0d0
    jac(3,5) = (f%d2pth(1)*f%hph + f%dpth(1)*f%dhph(1))*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(1)*f%vpar + f%dpth(1)*f%dvpar(1) &
      - f%d2H(1)*f%hth - f%dH(1)*f%dhth(1))
    jac(3,2:3) = 0.5d0*((f%d2pth(2:3)*f%hph + f%dpth(1)*f%dhph(2:3))*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(2:3)*f%vpar + f%dpth(1)*f%dvpar(2:3) &
      - f%d2H(2:3)*f%hth - f%dH(1)*f%dhth(2:3)))
    jac(3,3) = jac(3,3) + f%dpth(1)*f%hph
    jac(3,4) = 0.5d0*(f%d2pth(7)*f%hph*(x(3) - si%z(3)) &
      - si%dt*(f%d2pth(7)*f%vpar + f%dpth(1)*f%dvpar(4) - f%d2H(7)*f%hth))

    jac(4,1) = 0d0
    jac(4,5) = f%d2pth(1)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) &
      - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
    jac(4,2) = 0.5d0*(f%d2pth(2)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(5)*f%dpth(1) + f%dH(3)*f%d2pth(2) &
      - f%d2H(2)*f%dpth(3) - f%dH(1)*f%d2pth(5)))
    jac(4,3) = 0.5d0*(f%d2pth(3)*(x(4) - si%z(4)) &
      + si%dt*(f%d2H(6)*f%dpth(1) + f%dH(3)*f%d2pth(3) &
      - f%d2H(3)*f%dpth(3) - f%dH(1)*f%d2pth(6)))
    jac(4,4) = 0.5d0*(f%d2pth(7)*(x(4) - si%z(4)) &
      + si%dt*(f%dH(3)*f%d2pth(7) + f%d2H(9)*f%dpth(1)&
      - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9)))
    jac(4,4) = jac(4,4) + f%dpth(1)

    jac(5,1) = 0d0
    jac(5,5) = f%d2pth(1)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(1)&
      + si%dt/2.0d0*(f%d2pth(1)*f%dH(2) + f%dpth(1)*f%d2H(2) &
      - f%d2pth(2)*f%dH(1) - f%dpth(2)*f%d2H(1))
    jac(5,2) = 0.5d0*(f%d2pth(2)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(2) &
      + si%dt/2.0d0*(f%d2pth(2)*f%dH(2) + f%dpth(1)*f%d2H(4) &
      - f%d2pth(4)*f%dH(1) - f%dpth(2)*f%d2H(2)))
    jac(5,3) = 0.5d0*(f%d2pth(3)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(3) &
      + si%dt/2.0d0*(f%d2pth(3)*f%dH(2) + f%dpth(1)*f%d2H(5) &
      - f%d2pth(5)*f%dH(1) - f%dpth(2)*f%d2H(3)))
    jac(5,4) = 0.5d0*(f%d2pth(7)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(4) &
      + si%dt/2.0d0*(f%d2pth(7)*f%dH(2) + f%dpth(1)*f%d2H(8) &
      - f%d2pth(8)*f%dH(1) - f%dpth(2)*f%d2H(7)))

end subroutine jac_midpoint_part1

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_midpoint_part2(si, f, fmid, x, jac)
  !
    type(symplectic_integrator_t), intent(in) :: si
    type(field_can_t), intent(inout) :: f
    type(field_can_t), intent(inout) :: fmid
    real(dp), intent(in)  :: x(5)
    real(dp), intent(out) :: jac(5, 5)

    ! fmid%dpth(1)*(f%pth - si%pthold) + si%dt*(fmid%dpth(1)*fmid%dH(2)-fmid%dpth(2)*fmid%dH(1))

    jac(1,1) = fmid%dpth(1)*f%dpth(1)
    jac(1,2) = 0.5d0*(fmid%d2pth(2)*(f%pth - si%pthold) &
      + si%dt*(fmid%d2pth(2)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(4) &
      - fmid%dpth(2)*fmid%d2H(2) - fmid%d2pth(4)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(2)
    jac(1,3) = 0.5d0*(fmid%d2pth(3)*(f%pth - si%pthold) &
      + si%dt*(fmid%d2pth(3)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(5) &
      - fmid%dpth(2)*fmid%d2H(3) - fmid%d2pth(5)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(3)
    jac(1,4) = 0.5d0*(fmid%d2pth(7)*(f%pth - si%pthold) &
        + si%dt*(fmid%d2pth(7)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(8) &
        - fmid%dpth(2)*fmid%d2H(7) - fmid%d2pth(8)*fmid%dH(1))) + fmid%dpth(1)*f%dpth(4)
    jac(1,5) = fmid%d2pth(1)*(f%pth - si%pthold) &
        + si%dt*(fmid%d2pth(1)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(2) &
        - fmid%dpth(2)*fmid%d2H(1) - fmid%d2pth(2)*fmid%dH(1))

end subroutine jac_midpoint_part2

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine newton1(si, f, x, maxit, xlast)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, parameter :: n = 2

  real(dp), intent(inout) :: x(n)
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(n)

  real(dp) :: fvec(n), fjac(n,n), ijac(n,n)
  real(dp) :: tolref(n)
  integer :: kit

  tolref(1) = 1d0
  tolref(2) = dabs(1d1*torflux/f%ro0)

  do kit = 1, maxit
    if(x(1) > 1d0) return
    if(x(1) < 0d0) x(1) = 0.01d0

    call f_sympl_euler1(si, f, n, x, fvec, 1)
    call jac_sympl_euler1(si, f, x, fjac)
    ijac(1,1) = 1d0/(fjac(1,1) - fjac(1,2)*fjac(2,1)/fjac(2,2))
    ijac(1,2) = -1d0/(fjac(1,1)*fjac(2,2)/fjac(1,2) - fjac(2,1))
    ijac(2,1) = -1d0/(fjac(1,1)*fjac(2,2)/fjac(2,1) - fjac(1,2))
    ijac(2,2) = 1d0/(fjac(2,2) - fjac(1,2)*fjac(2,1)/fjac(1,1))
    xlast = x
    x = x - matmul(ijac, fvec)

    ! Don't take too small values in pphi as tolerance reference
    tolref(2) = max(dabs(x(2)), tolref(2))

    if (all(dabs(fvec) < si%atol)) return
    if (all(dabs(x-xlast) < si%rtol*tolref)) return
  enddo
  print *, 'newton1: maximum iterations reached: ', maxit
  write(6601,*) x(1), x(2)
  write(6601,*) x-xlast
  write(6601,*) fvec
  write(6601,*) ''
  write(6601,*) fjac(1,1), fjac(1,2)
  write(6601,*) fjac(2,1), fjac(2,2)
  write(6601,*) ''
  write(6601,*) ijac(1,1), ijac(1,2)
  write(6601,*) ijac(2,1), ijac(2,2)
  write(6601,*) ''
  write(6601,*) si%z(2), si%z(3)
  write(6601,*) ''
  write(6601,*) ''
end subroutine

recursive subroutine newton2(si, f, x, atol, rtol, maxit, xlast)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f

  integer, parameter :: n = 3
  integer :: kit

  real(dp), intent(inout) :: x(n)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(n)

  real(dp) :: fvec(n), fjac(n,n), jinv(n,n)
  integer :: pivot(n), info

  real(dp) :: xabs(n), tolref(n), fabs(n)
  real(dp) :: det

  do kit = 1, maxit
    if(x(1) > 1.0) return
    if(x(1) < 0.0) x(1) = 0.01
    call f_sympl_euler2(si, f, n, x, fvec, 1)
    fabs = dabs(fvec)
    call jac_sympl_euler2(si, f, x, fjac)
    xlast = x

    det = fjac(1,1)*fjac(2,2)*fjac(3,3) + fjac(1,2)*fjac(2,3)*fjac(3,1) + fjac(1,3)*fjac(2,1)*fjac(3,2) &
        - fjac(1,3)*fjac(2,2)*fjac(3,1) - fjac(1,1)*fjac(2,3)*fjac(3,2) - fjac(1,2)*fjac(2,1)*fjac(3,3)

    jinv(1,1) = fjac(2,2)*fjac(3,3) - fjac(2,3)*fjac(3,2)
    jinv(1,2) = fjac(1,3)*fjac(3,2) - fjac(1,2)*fjac(3,3)
    jinv(1,3) = fjac(1,2)*fjac(2,3) - fjac(1,3)*fjac(2,2)

    jinv(2,1) = fjac(2,3)*fjac(3,1) - fjac(2,1)*fjac(3,3)
    jinv(2,2) = fjac(1,1)*fjac(3,3) - fjac(1,3)*fjac(3,1)
    jinv(2,3) = fjac(1,3)*fjac(2,1) - fjac(1,1)*fjac(2,3)

    jinv(3,1) = fjac(2,1)*fjac(3,2) - fjac(3,1)*fjac(2,2)
    jinv(3,2) = fjac(1,2)*fjac(3,1) - fjac(3,2)*fjac(1,1)
    jinv(3,3) = fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1)

    x = x - matmul(jinv, fvec)/det

    !call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    ! after solution: fvec = (xold-xnew)_Newton
    !x = x - fvec

    xabs = dabs(x-xlast)

    tolref(1) = 1d0
    tolref(2) = twopi
    tolref(3) = twopi

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton2: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6602,*) x(1), x(2), x(3), si%z(4), xabs
end subroutine

recursive subroutine newton_midpoint(si, f, x, atol, rtol, maxit, xlast)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f

  type(field_can_t) :: fmid

  integer, parameter :: n = 5
  integer :: kit

  real(dp), intent(inout) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(n)

  real(dp) :: fvec(n), fjac(n,n)
  integer :: pivot(n), info

  real(dp) :: xabs(n), tolref(n), fabs(n)

  tolref(1) = 1d0
  tolref(2) = twopi
  tolref(3) = twopi
  tolref(4) = max(dabs(f%Aph), dabs(1d1*torflux/f%ro0))
  tolref(5) = 1d0

  do kit = 1, maxit
    if(x(1) > 1.0) return
    if(x(1) < 0.0) x(1) = 0.01
    if(x(5) < 0.0) x(5) = 0.01
    call f_midpoint_part1(si, f, n, x, fvec)
    call jac_midpoint_part1(si, f, x, fjac)
    fmid = f
    call f_midpoint_part2(si, f, n, x, fvec)
    call jac_midpoint_part2(si, f, fmid, x, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)

    ! Don't take too small values in pphi as tolerance reference
    tolref(4) = max(dabs(x(4)), tolref(4))

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
end subroutine

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! Gauss-Legendre Runge-Kutta method with s internal stages (n=4*s variables)
  !
recursive subroutine f_rk_gauss(si, fs, s, x, fvec)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  integer, intent(in) :: s
  type(field_can_t), intent(inout) :: fs(:)
  real(dp), intent(in) :: x(4*s)  ! = (rend, thend, phend, pphend)
  real(dp), intent(out) :: fvec(4*s)

  real(dp) :: a(s,s), b(s), c(s), Hprime(s)
  integer :: k,l  ! counters

  call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

  ! evaluate stages
  do k = 1, s
    call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 2)
    call get_derivatives2(fs(k), x(4*k))
    Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
  end do

  do k = 1, s
    fvec(4*k-3) = fs(k)%pth - si%pthold
    fvec(4*k-2) = x(4*k-2)  - si%z(2)
    fvec(4*k-1) = x(4*k-1)  - si%z(3)
    fvec(4*k)   = x(4*k)    - si%z(4)
    do l = 1, s
      fvec(4*k-3) = fvec(4*k-3) + si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
      fvec(4*k-2) = fvec(4*k-2) - si%dt*a(k,l)*Hprime(l)                                      ! thdot
      fvec(4*k-1) = fvec(4*k-1) - si%dt*a(k,l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
      fvec(4*k)   = fvec(4*k)   + si%dt*a(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
    end do
  end do

  end subroutine f_rk_gauss



  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_rk_gauss(si, fs, s, jac)
  !
  type(symplectic_integrator_t), intent(in) :: si
  integer, intent(in) :: s
  type(field_can_t), intent(in) :: fs(:)
  real(dp), intent(out) :: jac(4*s, 4*s)

  real(dp) :: a(s,s), b(s), c(s), Hprime(s), dHprime(4*s)
  integer :: k,l,m  ! counters

  call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

  ! evaluate stages
  do k = 1, s
    m=4*k
    Hprime(k)      = fs(k)%dH(1)/fs(k)%dpth(1)
    dHprime(m-3) = (fs(k)%d2H(1) - Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)  ! d/dr
    dHprime(m-2) = (fs(k)%d2H(2) - Hprime(k)*fs(k)%d2pth(2))/fs(k)%dpth(1)  ! d/dth
    dHprime(m-1) = (fs(k)%d2H(3) - Hprime(k)*fs(k)%d2pth(3))/fs(k)%dpth(1)  ! d/dph
    dHprime(m)   = (fs(k)%d2H(7) - Hprime(k)*fs(k)%d2pth(7))/fs(k)%dpth(1)  ! d/dpph
  end do

  jac = 0d0

  do k = 1, s
    m=4*k
    jac(m-3, m-3) = fs(k)%dpth(1)
    jac(m-3, m-2) = fs(k)%dpth(2)
    jac(m-3, m-1) = fs(k)%dpth(3)
    jac(m-3, m)   = fs(k)%dpth(4)
    jac(m-2, m-2) = 1d0
    jac(m-1, m-1) = 1d0
    jac(m, m) = 1d0
  end do

  do l = 1, s
    do k = 1, s
        m=4*k
        jac(m-3, 4*l-3) = jac(m-3, 4*l-3) & ! d/dr
          + si%dt*a(k,l)*(fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-3))
        jac(m-3, 4*l-2) = jac(m-3, 4*l-2) & ! d/dth
          + si%dt*a(k,l)*(fs(l)%d2H(4) - fs(l)%d2pth(4)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-2))
        jac(m-3, 4*l-1) = jac(m-3, 4*l-1) & ! d/dph
          + si%dt*a(k,l)*(fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l-1))
        jac(m-3, 4*l) = jac(m-3, 4*l) & ! d/dpph
          + si%dt*a(k,l)*(fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(4*l))

      jac(m-2, 4*l-3)   = jac(m-2, 4*l-3)   - si%dt*a(k,l)*dHprime(4*l-3)   ! d/dr
      jac(m-2, 4*l-2) = jac(m-2, 4*l-2) - si%dt*a(k,l)*dHprime(4*l-2) ! d/dth
      jac(m-2, 4*l-1) = jac(m-2, 4*l-1) - si%dt*a(k,l)*dHprime(4*l-1) ! d/dph
      jac(m-2, 4*l) = jac(m-2, 4*l) - si%dt*a(k,l)*dHprime(4*l) ! d/dpph

      jac(m-1, 4*l-3) = jac(m-1, 4*l-3) & ! d/dr
        - si%dt*a(k,l)*(-fs(l)%dhph(1)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
          + (fs(l)%dvpar(1) - dHprime(4*l-3)*fs(l)%hth - Hprime(l)*fs(l)%dhth(1))/fs(l)%hph)
      jac(m-1, 4*l-2) = jac(m-1, 4*l-2) & ! d/dth
        - si%dt*a(k,l)*(-fs(l)%dhph(2)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
          + (fs(l)%dvpar(2) - dHprime(4*l-2)*fs(l)%hth - Hprime(l)*fs(l)%dhth(2))/fs(l)%hph)
      jac(m-1, 4*l-1) = jac(m-1, 4*l-1) & ! d/dph
        - si%dt*a(k,l)*(-fs(l)%dhph(3)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
          + (fs(l)%dvpar(3) - dHprime(4*l-1)*fs(l)%hth - Hprime(l)*fs(l)%dhth(3))/fs(l)%hph)
      jac(m-1, 4*l) = jac(m-1, 4*l) & ! d/dpph
        - si%dt*a(k,l)*((fs(l)%dvpar(4) - dHprime(4*l)*fs(l)%hth)/fs(l)%hph)

      jac(m, 4*l-3) = jac(m, 4*l-3) & ! d/dr
        + si%dt*a(k,l)*(fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-3))
      jac(m, 4*l-2) = jac(m, 4*l-2) & ! d/dth
        + si%dt*a(k,l)*(fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-2))
      jac(m, 4*l-1) = jac(m, 4*l-1) & ! d/dph
        + si%dt*a(k,l)*(fs(l)%d2H(6) - fs(l)%d2pth(6)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l-1))
      jac(m, 4*l) = jac(m, 4*l) & ! d/dpph
        + si%dt*a(k,l)*(fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(4*l))
    end do
  end do

end subroutine jac_rk_gauss

recursive subroutine newton_rk_gauss(si, fs, s, x, atol, rtol, maxit, xlast)
  type(symplectic_integrator_t), intent(inout) :: si

  integer, intent(in) :: s
  type(field_can_t), intent(inout) :: fs(:)
  integer :: kit, ks

  real(dp), intent(inout) :: x(4*s)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(4*s)

  real(dp) :: fvec(4*s), fjac(4*s, 4*s)
  integer :: pivot(4*s), info

  real(dp) :: xabs(4*s), tolref(4*s), fabs(4*s)

  do kit = 1, maxit

    ! Check if radius left the boundary
    do ks = 1, s
      if (x(4*ks-3) > 1d0) return
      if (x(4*ks-3) < 0.0) x(4*ks-3) = 0.01d0
    end do

    call f_rk_gauss(si, fs, s, x, fvec)
    call jac_rk_gauss(si, fs, s, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(4*s, 1, fjac, 4*s, pivot, fvec, 4*s, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)

    do ks = 1, s
      tolref(4*ks-3) = 1d0
      tolref(4*ks-2) = twopi
      tolref(4*ks-1) = twopi
      tolref(4*ks) = dabs(xlast(4*ks))
    end do

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton_rk_gauss: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6603,*) x, xabs, fvec
end subroutine newton_rk_gauss


recursive subroutine fixpoint_rk_gauss(si, fs, s, x, atol, rtol, maxit, xlast)
  ! TODO: this doesn't work well yet
  type(symplectic_integrator_t), intent(inout) :: si

  integer, intent(in) :: s
  type(field_can_t), intent(inout) :: fs(:)
  integer :: kit, ks

  real(dp), intent(inout) :: x(4*s)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(4*s)

  real(dp) :: fvec(4*s)

  real(dp) :: xabs(4*s), tolref(4*s), fabs(4*s)

  real(dp) :: a(s,s), b(s), c(s), Hprime(s), dHprimedr(s)
  integer :: k, l

  real(dp) :: pthnew, dpthnewdr, damp

  call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

  do kit = 1, 4096

    ! Check if radius left the boundary
    do ks = 1, s
      if (x(4*ks-3) > 1d0) return
      if (x(4*ks-3) < 0.0) x(4*ks-3) = 0.01d0
    end do

    call f_rk_gauss(si, fs, s, x, fvec)
    fabs = dabs(fvec)
    xlast = x

    do k = 1, s
      call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 2)
      call get_derivatives2(fs(k), x(4*k))
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
      dHprimedr(k) = (fs(k)%d2H(1) - Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)
    end do

    k = 1
    l = 1

    pthnew = si%pthold - si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
    dpthnewdr = -si%dt*a(k,l)*( &
    fs(l)%d2H(2) - Hprime(l)*fs(l)%d2pth(2) - dHprimedr(l)*fs(l)%dpth(2))
    !print *, pthnew, dpthnewdr, pthnew/dpthnewdr

    damp = 0.99d0

    x(1) = x(1) - (1d0-damp)*pthnew/dpthnewdr
    x(2) = damp*x(2) + (1d0-damp)*(si%z(2) + si%dt*a(k,l)*si%dt*a(k,l)*Hprime(l))
    x(3) = damp*x(3) + (1d0-damp)*(si%z(3) + si%dt*a(k,l)*(fs(l)%vpar &
                                    - Hprime(l)*fs(l)%hth)/fs(l)%hph)
    x(4) = damp*x(4) + (1d0-damp)*(si%z(4) - si%dt*a(k,l)*(fs(l)%dH(3) &
                                    - Hprime(l)*fs(l)%dpth(3)))

    print *, x

    xabs = dabs(x - xlast)

    do ks = 1, s
      xabs(4*ks-2) = modulo(xabs(4*ks-2), pi)
      xabs(4*ks-1) = modulo(xabs(4*ks-1), pi)
    end do

    do ks = 1, s
      tolref(4*ks-3) = 1d0
      tolref(4*ks-2) = twopi
      tolref(4*ks-1) = twopi
      tolref(4*ks) = dabs(x(4*ks))
    end do

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'fixpoint_rk_gauss: maximum iterations reached: ', maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6603,*) x, xabs, fvec
end subroutine fixpoint_rk_gauss


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_rk_lobatto(si, fs, s, jac)
  !
  type(symplectic_integrator_t), intent(in) :: si
  integer, intent(in) :: s
  type(field_can_t), intent(in) :: fs(:)
  real(dp), intent(out) :: jac(4*s-2, 4*s-2)

  real(dp) :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s), dHprime(4*s-2)
  integer :: k,l,m,n  ! counters

  call coeff_rk_lobatto(s, a, ahat, b, c)
  jac = 0d0

  Hprime = fs%dH(1)/fs%dpth(1)
  dHprime(1) = (fs(1)%d2H(1)-Hprime(1)*fs(1)%d2pth(1))/fs(1)%dpth(1)  ! d/dr
  dHprime(2) = (fs(1)%d2H(7)-Hprime(1)*fs(1)%d2pth(7))/fs(1)%dpth(1)  ! d/dpph
  do k = 2, s
    m = 4*k-2
    dHprime(m-3)=(fs(k)%d2H(1)-Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)  ! d/dr
    dHprime(m-2)=(fs(k)%d2H(2)-Hprime(k)*fs(k)%d2pth(2))/fs(k)%dpth(1)  ! d/dth
    dHprime(m-1)=(fs(k)%d2H(3)-Hprime(k)*fs(k)%d2pth(3))/fs(k)%dpth(1)  ! d/dph
    dHprime(m)  =(fs(k)%d2H(7)-Hprime(k)*fs(k)%d2pth(7))/fs(k)%dpth(1)  ! d/dpph
  end do

  ! evaluate first stage with only r and pph as unknowns
  jac(1, 1) = fs(1)%dpth(1)
  jac(1, 2) = fs(1)%dpth(4)
  jac(2, 2) = 1d0

  jac(1, 1) = jac(1, 1) + si%dt*ahat(1,1) * &
    (fs(1)%d2H(2) - fs(1)%d2pth(2)*Hprime(1) - fs(1)%dpth(2)*dHprime(1))
  jac(1, 2) = jac(1, 2) + si%dt*ahat(1,1)* &
    (fs(1)%d2H(8) - fs(1)%d2pth(8)*Hprime(1) - fs(1)%dpth(2)*dHprime(2))
  jac(2, 1) = jac(2, 1) + si%dt*ahat(1,1) * &
    (fs(1)%d2H(3) - fs(1)%d2pth(3)*Hprime(1) - fs(1)%dpth(3)*dHprime(1))
  jac(2, 2) = jac(2, 2) + si%dt*ahat(1,1) * &
    (fs(1)%d2H(9) - fs(1)%d2pth(9)*Hprime(1) - fs(1)%dpth(3)*dHprime(2))

  do l = 2, s
    n = 4*l-2
    jac(1, n-3) = jac(1, n-3) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-3))
    jac(1, n-2) = jac(1, n-2) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(4) - fs(l)%d2pth(4)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-2))
    jac(1, n-1) = jac(1, n-1) + si%dt*ahat(1,l)* &
      (fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-1))
    jac(1, n) = jac(1, n) + si%dt*ahat(1,l)* &
      (fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(n))

    jac(2, n-3) = jac(2, n-3) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-3))
    jac(2, n-2) = jac(2, n-2) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-2))
    jac(2, n-1) = jac(2, n-1) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(6) - fs(l)%d2pth(6)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-1))
    jac(2, n) = jac(2, n) + si%dt*ahat(1,l) * &
      (fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(n))
  end do

  ! evaluate remaining stages with r, th, ph, pph as unknowns

  do k = 2, s
    m = 4*k-2
    jac(m-3, m-3) = fs(k)%dpth(1)
    jac(m-3, m-2) = fs(k)%dpth(2)
    jac(m-3, m-1) = fs(k)%dpth(3)
    jac(m-3, m)   = fs(k)%dpth(4)
    jac(m-2, m-2) = 1d0
    jac(m-1, m-1) = 1d0
    jac(m, m) = 1d0
  end do

  l = 1
  do k = 2, s
    m = 4*k-2
    jac(m-3, 1) = jac(m-3, 1) + si%dt*ahat(k,l) * &              ! d/dr
      (fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(1))
    jac(m-3, 2) = jac(m-3, 2) + si%dt*ahat(k,l)* &                   ! d/dpph
      (fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(2))

    jac(m-2, 1) = jac(m-2, 1) - si%dt*a(k,l)*dHprime(1)  ! d/dr
    jac(m-2, 2) = jac(m-2, 2) - si%dt*a(k,l)*dHprime(2)  ! d/dpph

    jac(m-1, 1) = jac(m-1, 1) - si%dt*a(k,l) * & ! d/dr
      (-fs(l)%dhph(1)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
      +(fs(l)%dvpar(1)-dHprime(1)*fs(l)%hth-Hprime(l)*fs(l)%dhth(1))/fs(l)%hph)
    jac(m-1, 2) = jac(m-1, 2) - si%dt*a(k,l) * & ! d/dpph
      (fs(l)%dvpar(4)-dHprime(2)*fs(l)%hth)/fs(l)%hph

    jac(m, 1) = jac(m, 1) + si%dt*ahat(k,l) * & ! d/dr
      (fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(1))
    jac(m, 2) = jac(m, 2) + si%dt*ahat(k,l) * & ! d/dpph
      (fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(2))
  end do

  do l = 2, s
    do k = 2, s
      m = 4*k-2
      n = 4*l-2
      jac(m-3, n-3) = jac(m-3, n-3) + si%dt*ahat(k,l) * &              ! d/dr
        (fs(l)%d2H(2) - fs(l)%d2pth(2)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-3))
      jac(m-3, n-2) = jac(m-3, n-2) + si%dt*ahat(k,l) * &              ! d/dth
        (fs(l)%d2H(4) - fs(l)%d2pth(4)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-2))
      jac(m-3, n-1) = jac(m-3, n-1) + si%dt*ahat(k,l)* &               ! d/dph
        (fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(2)*dHprime(n-1))
      jac(m-3, n) = jac(m-3, n) + si%dt*ahat(k,l)* &                   ! d/dpph
        (fs(l)%d2H(8) - fs(l)%d2pth(8)*Hprime(l) - fs(l)%dpth(2)*dHprime(n))

      jac(m-2, n-3) = jac(m-2, n-3) - si%dt*a(k,l)*dHprime(n-3)
      jac(m-2, n-2) = jac(m-2, n-2) - si%dt*a(k,l)*dHprime(n-2)
      jac(m-2, n-1) = jac(m-2, n-1) - si%dt*a(k,l)*dHprime(n-1)
      jac(m-2, n)   = jac(m-2, n)   - si%dt*a(k,l)*dHprime(n)

      jac(m-1, n-3) = jac(m-1, n-3) - si%dt*a(k,l) * & ! d/dr
        (-fs(l)%dhph(1)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
        +(fs(l)%dvpar(1)-dHprime(n-3)*fs(l)%hth-Hprime(l)*fs(l)%dhth(1))/fs(l)%hph)
      jac(m-1, n-2) = jac(m-1, n-2) - si%dt*a(k,l) * & ! d/dth
        (-fs(l)%dhph(2)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
        +(fs(l)%dvpar(2)-dHprime(n-2)*fs(l)%hth-Hprime(l)*fs(l)%dhth(2))/fs(l)%hph)
      jac(m-1, n-1) = jac(m-1, n-1) - si%dt*a(k,l) * & ! d/dph
        (-fs(l)%dhph(3)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph**2 &
        +(fs(l)%dvpar(3)-dHprime(n-1)*fs(l)%hth-Hprime(l)*fs(l)%dhth(3))/fs(l)%hph)
      jac(m-1, n) = jac(m-1, n) - si%dt*a(k,l) * & ! d/dpph
        (fs(l)%dvpar(4)-dHprime(n)*fs(l)%hth)/fs(l)%hph

      jac(m, n-3) = jac(m, n-3) + si%dt*ahat(k,l) * & ! d/dr
        (fs(l)%d2H(3) - fs(l)%d2pth(3)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-3))
      jac(m, n-2) = jac(m, n-2) + si%dt*ahat(k,l) * & ! d/dth
        (fs(l)%d2H(5) - fs(l)%d2pth(5)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-2))
      jac(m, n-1) = jac(m, n-1) + si%dt*ahat(k,l) * & ! d/dph
        (fs(l)%d2H(6) - fs(l)%d2pth(6)*Hprime(l) - fs(l)%dpth(3)*dHprime(n-1))
      jac(m, n) = jac(m, n) + si%dt*ahat(k,l) * & ! d/dpph
        (fs(l)%d2H(9) - fs(l)%d2pth(9)*Hprime(l) - fs(l)%dpth(3)*dHprime(n))
    end do
  end do
end subroutine jac_rk_lobatto


recursive subroutine newton_rk_lobatto(si, fs, s, x, atol, rtol, maxit, xlast)
  type(symplectic_integrator_t), intent(inout) :: si

  integer, intent(in) :: s
  type(field_can_t), intent(inout) :: fs(:)
  integer :: kit, ks

  real(dp), intent(inout) :: x(4*s-2)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(4*s-2)

  real(dp) :: fvec(4*s-2), fjac(4*s-2, 4*s-2)
  integer :: pivot(4*s-2), info

  real(dp) :: xabs(4*s-2), tolref(4*s-2), fabs(4*s-2)

  do kit = 1, maxit

    ! Check if radius left the boundary
    if (x(1) > 1d0) return
    if (x(1) < 0.0) x(1) = 0.01d0
    do ks = 2, s
      if (x(4*ks-2-3) > 1d0) return
      if (x(4*ks-2-3) < 0.0) x(4*ks-3) = 0.01d0
    end do

    call f_rk_lobatto(si, fs, s, x, fvec, 2)
    call jac_rk_lobatto(si, fs, s, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(4*s-2, 1, fjac, 4*s-2, pivot, fvec, 4*s-2, info)
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)

    tolref(1) = 1d0
    tolref(2) = dabs(xlast(2))
    do ks = 2, s
      tolref(4*ks-2-3) = 1d0
      tolref(4*ks-2-2) = twopi
      tolref(4*ks-2-1) = twopi
      tolref(4*ks-2) = dabs(xlast(4*ks-2))
    end do

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  print *, 'newton_rk_lobatto: maximum iterations reached: ', &
           maxit, 'z = ', x(1), x(2), x(3), si%z(4)
  write(6603,*) x, xabs, fvec
end subroutine newton_rk_lobatto


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_multi(mi, f, ierr)
  !
  type(multistage_integrator_t), intent(inout) :: mi
  type(field_can_t), intent(inout) :: f

  integer, intent(out) :: ierr

  integer :: kstage

  call orbit_timestep_sympl(mi%stages(1), f, ierr)

  do kstage = 2, 2*mi%s
    mi%stages(kstage)%z = mi%stages(kstage-1)%z
    mi%stages(kstage)%pthold = f%pth
    call orbit_timestep_sympl(mi%stages(kstage), f, ierr)
  end do

  mi%stages(1)%z = mi%stages(2*mi%s)%z
  mi%stages(1)%pthold = f%pth

end subroutine orbit_timestep_sympl_multi


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  !
  type(multistage_integrator_t), intent(inout) :: mi
  type(field_can_t), intent(inout) :: f

  real(dp), intent(in) :: z(:)
  real(dp), intent(in) :: dtau
  integer, intent(in) :: ntau
  real(dp), intent(in) :: rtol_init

  real(dp), intent(in) :: alpha(:), beta(:)

  integer :: ks

  mi%s = size(alpha)

  do ks = 1, mi%s
    call orbit_sympl_init(mi%stages(2*ks-1), f, z, &
      alpha(ks)*dtau, ntau, rtol_init, 1)
    call orbit_sympl_init(mi%stages(2*ks), f, z, &
      beta(ks)*dtau, ntau, rtol_init, 2)
  end do
end subroutine orbit_sympl_init_multi


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_verlet(mi, f, z, dtau, ntau, rtol_init)
  !
  type(multistage_integrator_t), intent(inout) :: mi

  type(field_can_t), intent(inout) :: f

  real(dp), intent(in) :: z(:)
  real(dp), intent(in) :: dtau
  integer, intent(in) :: ntau
  real(dp), intent(in) :: rtol_init

  real(dp) :: alpha(1), beta(1)

  alpha(1) = 0.5d0
  beta(1)  = 0.5d0

  call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
end subroutine orbit_sympl_init_verlet


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_order4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=3
  !
  !
    type(multistage_integrator_t), intent(inout) :: mi
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dtau
    integer, intent(in) :: ntau
    real(dp), intent(in) :: rtol_init

    real(dp) :: alpha(5), beta(5)

    alpha(1) = 1d0/(2d0*(2d0 - 2d0**(1d0/3d0)))
    alpha(2) = 2d0**(1d0/3d0)/(2d0*(2d0 - 2d0**(1d0/3d0)))
    alpha(3) = alpha(1)

    beta = alpha

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_order4


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_mclachlan4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=5 by McLachlan (1995)
  !
  !
    type(multistage_integrator_t), intent(inout) :: mi
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dtau
    integer, intent(in) :: ntau
    real(dp), intent(in) :: rtol_init

    real(dp) :: alpha(5), beta(5)

    alpha(1) = (146d0 + 5d0*dsqrt(19d0))/540d0
    alpha(2) = (-2d0 + 10d0*dsqrt(19d0))/135d0
    alpha(3) = 1d0/5d0
    alpha(4) = (-23d0 + 20d0*dsqrt(19d0))/270d0
    alpha(5) = (14d0 - dsqrt(19d0))/108d0

    beta(5)  = alpha(1)
    beta(4) = alpha(2)
    beta(3) = alpha(3)
    beta(2) = alpha(4)
    beta(1) = alpha(5)

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_mclachlan4


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_blanes4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=6 by Blanes&Moan (2002)
  ! with coefficients in the form of Hairer (2002)
  !
  !
    type(multistage_integrator_t), intent(inout) :: mi
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dtau
    integer, intent(in) :: ntau
    real(dp), intent(in) :: rtol_init

    real(dp) :: alpha(6), beta(6)

    alpha(1) = 0.16231455076687d0
    alpha(2) = 0.37087741497958d0
    alpha(3) = 0.059762097006575d0
    alpha(4) = 0.40993371990193d0
    alpha(5) = 0.23399525073150d0
    alpha(6) = 0.082984406417405d0

    beta(6) = alpha(1)
    beta(5) = alpha(2)
    beta(4) = alpha(3)
    beta(3) = alpha(4)
    beta(2) = alpha(5)
    beta(1) = alpha(6)

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_blanes4


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_kahan6(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 6 with s=9 by Kahan&Li (1995)
  ! with coefficients in the form of Hairer (2002)
  !
  !
    type(multistage_integrator_t), intent(inout) :: mi
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dtau
    integer, intent(in) :: ntau
    real(dp), intent(in) :: rtol_init

    real(dp) :: gam(9)

    gam(1) = 0.39216144400731413927925056d0
    gam(2) = 0.33259913678935943859974864d0
    gam(3) = -0.70624617255763935980996482d0
    gam(4) = 0.08221359629355080023149045d0
    gam(5) = 0.79854399093482996339895035d0
    gam(6) = gam(4)
    gam(7) = gam(3)
    gam(8) = gam(2)
    gam(9) = gam(1)

    call orbit_sympl_init_multi(mi,f,z,dtau,ntau,rtol_init,gam/2d0,gam/2d0)
  end subroutine orbit_sympl_init_kahan6


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  recursive subroutine orbit_sympl_init_kahan8(mi, f, z, dtau, ntau, rtol_init)
    !
    ! Composition method of order 8 with s=17 by Kahan&Li (1995)
    ! with coefficients in the form of Hairer (2002)
    !
    !
      type(multistage_integrator_t), intent(inout) :: mi
      type(field_can_t), intent(inout) :: f

      real(dp), intent(in) :: z(:)
      real(dp), intent(in) :: dtau
      integer, intent(in) :: ntau
      real(dp), intent(in) :: rtol_init

      real(dp) :: gam(17)

      gam(1) =  0.13020248308889008087881763d0
      gam(2) =  0.56116298177510838456196441d0
      gam(3) = -0.38947496264484728640807860d0
      gam(4) =  0.15884190655515560089621075d0
      gam(5) = -0.39590389413323757733623154d0
      gam(6) =  0.18453964097831570709183254d0
      gam(7) =  0.25837438768632204729397911d0
      gam(8) =  0.29501172360931029887096624d0
      gam(9) = -0.60550853383003451169892108d0
      gam(10) = gam(8)
      gam(11) = gam(7)
      gam(12) = gam(6)
      gam(13) = gam(5)
      gam(14) = gam(4)
      gam(15) = gam(3)
      gam(16) = gam(2)
      gam(17) = gam(1)

      call orbit_sympl_init_multi(mi,f,z,dtau,ntau,rtol_init,gam/2d0,gam/2d0)
    end subroutine orbit_sympl_init_kahan8

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_expl_impl_euler(si, f, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  integer, parameter :: n = 2
  integer, parameter :: maxit = 32

  real(dp), dimension(n) :: x, xlast
  integer :: k, ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x(1)=si%z(1)
    x(2)=si%z(4)

    call newton1(si, f, x, maxit, xlast)

    if (x(1) > 1.0d0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0d0) then
      print *, 'r<0, ipart=', current_particle_index, ', z=', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01d0
    end if

    si%z(1) = x(1)
    si%z(4) = x(2)

    if (extrap_field) then
      f%pth = f%pth + f%dpth(1)*(x(1)-xlast(1))  + f%dpth(4)*(x(2)-xlast(2))
      f%dH(1) = f%dH(1) + f%d2H(1)*(x(1)-xlast(1)) + f%d2H(7)*(x(2)-xlast(2))
      f%dpth(1)=f%dpth(1)+f%d2pth(1)*(x(1)-xlast(1))+f%d2pth(7)*(x(2)-xlast(2))
      f%vpar = f%vpar + f%dvpar(1)*(x(1)-xlast(1)) + f%dvpar(4)*(x(2)-xlast(2))
      f%hth = f%hth + f%dhth(1)*(x(1)-xlast(1))
      f%hph = f%hph + f%dhph(1)*(x(1)-xlast(1))
    else
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    endif

    si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
    si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_expl_impl_euler


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_impl_expl_euler(si, f, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: n = 3
  integer, parameter :: maxit = 32

  real(dp), dimension(n) :: x, xlast, dz
  integer :: k, ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x = si%z(1:3)

    call newton2(si, f, x, si%atol, si%rtol, maxit, xlast)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, ipart=', current_particle_index, ', z=', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    si%z(1:3) = x

    if (extrap_field) then
      dz(1) = x(1)-xlast(1)
      dz(2) = x(2)-xlast(2)
      dz(3) = x(3)-xlast(3)

      f%dH(1) = f%dH(1) + f%d2H(1)*dz(1) + f%d2H(2)*dz(2) + f%d2H(3)*dz(3)
      f%dH(2) = f%dH(2) + f%d2H(2)*dz(1) + f%d2H(4)*dz(2) + f%d2H(5)*dz(3)
      f%dH(3) = f%dH(3) + f%d2H(3)*dz(1) + f%d2H(5)*dz(2) + f%d2H(6)*dz(3)

      f%dpth(1) = f%dpth(1)+f%d2pth(1)*dz(1)+f%d2pth(2)*dz(2)+f%d2pth(3)*dz(3)
      f%dpth(2) = f%dpth(2)+f%d2pth(2)*dz(1)+f%d2pth(4)*dz(2)+f%d2pth(5)*dz(3)
      f%dpth(3) = f%dpth(3)+f%d2pth(3)*dz(1)+f%d2pth(5)*dz(2)+f%d2pth(6)*dz(3)
    else
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    endif

    f%pth = si%pthold - si%dt*(f%dH(2) - f%dH(1)*f%dpth(2)/f%dpth(1))
    si%z(4) = si%z(4) - si%dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_impl_expl_euler


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_midpoint(si, f, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: n = 5
  integer, parameter :: maxit = 8

  real(dp), dimension(n) :: x, xlast
  integer :: k, ktau

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x(1:4) = si%z
    x(5) = si%z(1)

    call newton_midpoint(si, f, x, si%atol, si%rtol, maxit, xlast)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, ipart=', current_particle_index, ', z=', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    si%z = x(1:4)

    if (extrap_field) then
      f%pth = f%pth + f%dpth(1)*(x(1)-xlast(1) + x(5) - xlast(5)) &  ! d/dr
                    + f%dpth(2)*(x(2)-xlast(2)) &  ! d/dth
                    + f%dpth(3)*(x(3)-xlast(3)) &  ! d/dph
                    + f%dpth(4)*(x(4)-xlast(4))    ! d/dpph
    else
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_val(f, si%z(4))
    endif

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_midpoint


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_rk_gauss(si, f, s, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: maxit = 32

  integer, intent(in) :: s
  real(dp), dimension(4*s) :: x, xlast
  integer :: k, l, ktau

  type(field_can_t) :: fs(s)
  real(dp) :: a(s,s), b(s), c(s), Hprime(s), dz(4*s)

  do k = 1,s
    fs(k) = f
    fs(k) = f
  end do

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    do k = 1,s
      x((4*k-3):(4*k)) = si%z
    end do

    call newton_rk_gauss(si, fs, s, x, si%atol, si%rtol, maxit, xlast)
    !optionally try fixed point iterations, doesn't work yet
    !call fixpoint_rk_gauss(si, fs, s, x, si%atol, si%rtol, maxit, xlast)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, ipart=', current_particle_index, ', z=', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

    if (extrap_field) then
      do k = 1, s
        dz(1) = x(4*k-3)-xlast(4*k-3)
        dz(2) = x(4*k-2)-xlast(4*k-2)
        dz(3) = x(4*k-1)-xlast(4*k-1)
        dz(4) = x(4*k)-xlast(4*k)

        fs(k)%pth = fs(k)%pth + fs(k)%dpth(1)*dz(1) &  ! d/dr
                      + fs(k)%dpth(2)*dz(2) &          ! d/dth
                      + fs(k)%dpth(3)*dz(3) &          ! d/dph
                      + fs(k)%dpth(4)*dz(4)            ! d/dpph

        fs(k)%dpth(1) = fs(k)%dpth(1) + fs(k)%d2pth(1)*dz(1) &
                      + fs(k)%d2pth(2)*dz(2) &
                      + fs(k)%d2pth(3)*dz(3) &
                      + fs(k)%d2pth(7)*dz(4)

        fs(k)%dpth(2) = fs(k)%dpth(2) + fs(k)%d2pth(2)*dz(1)&
                      + fs(k)%d2pth(4)*dz(2) &
                      + fs(k)%d2pth(5)*dz(3) &
                      + fs(k)%d2pth(8)*dz(4)

        fs(k)%dpth(3) = fs(k)%dpth(3) + fs(k)%d2pth(2)*dz(1) &
                      + fs(k)%d2pth(5)*dz(2) &
                      + fs(k)%d2pth(6)*dz(3) &
                      + fs(k)%d2pth(9)*dz(4)

        fs(k)%dH(1) = fs(k)%dH(1) + fs(k)%d2H(1)*dz(1) &
                      + fs(k)%d2H(2)*dz(2) &
                      + fs(k)%d2H(3)*dz(3) &
                      + fs(k)%d2H(7)*dz(4)

        fs(k)%dH(2) = fs(k)%dH(2) + fs(k)%d2H(2)*dz(1) &
                      + fs(k)%d2H(4)*dz(2) &
                      + fs(k)%d2H(5)*dz(3) &
                      + fs(k)%d2H(8)*dz(4)

        fs(k)%dH(3) = fs(k)%dH(3) + fs(k)%d2H(3)*dz(1) &
                      + fs(k)%d2H(5)*dz(2) &
                      + fs(k)%d2H(6)*dz(3) &
                      + fs(k)%d2H(9)*dz(4)

        fs(k)%vpar = fs(k)%vpar + fs(k)%dvpar(1)*dz(1) &  ! d/dr
                   + fs(k)%dvpar(2)*dz(2) &
                   + fs(k)%dvpar(3)*dz(3)

        fs(k)%hth = fs(k)%hth + fs(k)%dhth(1)*dz(1) &  ! d/dr
                      + fs(k)%dhth(2)*dz(2) &
                      + fs(k)%dhth(3)*dz(3)

        fs(k)%hph = fs(k)%hph + fs(k)%dhph(1)*dz(1) &  ! d/dr
                      + fs(k)%dhph(2)*dz(2) &
                      + fs(k)%dhph(3)*dz(3)
      end do
    else
      do k = 1, s
        call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 0)
        call get_derivatives(fs(k), x(4*k))
      end do
    end if

    do k = 1, s
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
    end do

    f = fs(s)
    f%pth = si%pthold
    si%z(1) = x(4*s-3)

    do l = 1, s
      f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
      si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
      si%z(3) = si%z(3)+si%dt*b(l)*(fs(l)%vpar-Hprime(l)*fs(l)%hth)/fs(l)%hph
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
    end do

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_rk_gauss


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_rk_lobatto(si, f, s, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  integer, parameter :: maxit = 32

  integer, intent(in) :: s
  real(dp), dimension(4*s-2) :: x, xlast
  integer :: ktau, k, l

  type(field_can_t) :: fs(s)
  real(dp) :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s)

  do k = 1,s
    fs(k) = f
  end do

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x(1) = si%z(1)
    x(2) = si%z(4)
    do k = 2,s
      x((4*k-3-2):(4*k-2)) = si%z
    end do

    call newton_rk_lobatto(si, fs, s, x, si%atol, si%rtol, maxit, xlast)

    call coeff_rk_lobatto(s, a, ahat, b, c)

    call eval_field(fs(1), x(1), si%z(2), si%z(3), 0)
    call get_derivatives(fs(1), x(2))
    Hprime(1) = fs(1)%dH(1)/fs(1)%dpth(1)

    do k = 2, s
      call eval_field(fs(k), x(4*k-3-2), x(4*k-2-2), x(4*k-1-2), 0)
      call get_derivatives(fs(k), x(4*k-2))
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
    end do

    f = fs(s)
    f%pth   = si%pthold
    si%z(1) = x(4*s-3-2)

    do l = 1, s
      f%pth   = f%pth   - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
      si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
      si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar-Hprime(l)*fs(l)%hth)/fs(l)%hph
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
    end do

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_rk_lobatto


recursive subroutine debug_root(si, f, x0)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), pointer :: f

  real(dp) :: x0(2), x(2)
  integer :: k, l, iflag
  integer, parameter :: n = 100
  real(dp), parameter :: eps = 1d-15

  real(dp) :: fvec(2)

  do k = -n,n
    do l = -n,n
      x = x0 + l*eps/n*(/x0(1),0d0/) + k*eps/n*(/0d0,x0(2)/)
      call f_sympl_euler1(si, f, 2, x, fvec, iflag)
      write(5001,*) x(1), x(2), fvec
    end do
  end do

end subroutine debug_root

end module orbit_symplectic
