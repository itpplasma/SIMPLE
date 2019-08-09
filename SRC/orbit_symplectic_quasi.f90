module orbit_symplectic_quasi

use field_can_mod
use orbit_symplectic

implicit none
save

logical, parameter :: exact_steps = .False.
integer, parameter :: S_MAX_GAUSS = 3

type(FieldCan) :: f
type(FieldCan) :: fs(S_MAX_GAUSS)
type(SymplecticIntegrator) :: si
!$omp threadprivate(f, si)

contains

!
! Wrapper routines for ODEPACK
!

subroutine f_exact_quasi(n, x, fvec, iflag)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  type(FieldCan) :: f2

  f2 = f

  call eval_field(f2, x(1), si%z(2), si%z(3), 0)
  call get_derivatives(f2, si%z(4))

  fvec(1) = f%pth - f2%pth

  f%H = f2%H

end subroutine f_exact_quasi


subroutine f_euler1_quasi(n, x, fvec, iflag)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(f, x(1), si%z(2), si%z(3), 0)
  call get_derivatives(f, x(2))

  fvec(1) = f%dpth(1)*(f%pth - si%pthold) + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
  fvec(2) = f%dpth(1)*(x(2) - si%z(4))  + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))

end subroutine f_euler1_quasi


subroutine f_euler2_quasi(n, x, fvec, iflag)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(f, x(1), x(2), x(3), 0)
  call get_derivatives(f, si%z(4))

  fvec(1) = f%pth - si%pthold
  fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
  fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)

end subroutine f_euler2_quasi


subroutine f_midpoint_quasi(n, x, fvec, iflag)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: fvec(n)
    integer, intent(in) :: iflag
    
    double precision :: dpthmid, pthdotbar
    
    call eval_field(f, x(5), 0.5*(x(2) + si%z(2)), 0.5*(x(3) + si%z(3)), 0)
    call get_derivatives(f, 0.5*(x(4) + si%z(4)))

    fvec(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
    fvec(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)
    fvec(4) = f%dpth(1)*(x(4) - si%z(4)) + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
    fvec(5) = f%dpth(1)*(f%pth - si%pthold) + si%dt/2.0d0*(f%dpth(1)*f%dH(2)-f%dpth(2)*f%dH(1))
  
    ! save evaluation from midpoint
    dpthmid = f%dpth(1)
    pthdotbar = f%dpth(1)*f%dH(2) - f%dpth(2)*f%dH(1)

    ! evaluate at endpoint
    call eval_field(f, x(1), x(2), x(3), 0)
    call get_derivatives(f, x(4))
    fvec(1) = dpthmid*(f%pth - si%pthold) + si%dt*pthdotbar

end subroutine f_midpoint_quasi


subroutine f_rk_gauss_quasi(n, x, fvec, iflag)

  integer, intent(in) :: n
  double precision, intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
  double precision, intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  double precision :: a(n/4,n/4), b(n/4), c(n/4), Hprime(n/4)
  integer :: k,l  ! counters

  call coeff_rk_gauss(n/4, a, b, c)  ! TODO: move this to preprocessing

  ! evaluate stages
  do k = 1, n/4
    call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 0)
    call get_derivatives(fs(k), x(4*k))
    Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
  end do

  do k = 1, n/4
    fvec(4*k-3) = fs(k)%pth - si%pthold
    fvec(4*k-2) = x(4*k-2)  - si%z(2)
    fvec(4*k-1) = x(4*k-1)  - si%z(3)
    fvec(4*k)   = x(4*k)    - si%z(4)
    do l = 1, n/4
      fvec(4*k-3) = fvec(4*k-3) + si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
      fvec(4*k-2) = fvec(4*k-2) - si%dt*a(k,l)*Hprime(l)                                      ! thdot
      fvec(4*k-1) = fvec(4*k-1) - si%dt*a(k,l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
      fvec(4*k)   = fvec(4*k)   + si%dt*a(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
    end do
  end do

end subroutine f_rk_gauss_quasi

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Lobatto (IIIA)-(IIIB) Runge-Kutta method with s internal stages (n=4*s variables)
!
subroutine f_rk_lobatto_quasi(s, x, fvec)
  !
  integer, intent(in) :: s
  double precision, intent(in) :: x(4*s)  ! = (rend, thend, phend, pphend)
  double precision, intent(out) :: fvec(4*s)

  double precision :: a(s,s), ahat(s), b(s), c(s), Hprime(s)
  integer :: k,l  ! counters

  call coeff_rk_lobatto(s, a, ahat, b, c)  ! TODO: move this to preprocessing

  ! evaluate stages
  do k = 1, s
    call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 0)
    call get_derivatives(fs(k), x(4*k))
    Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
  end do

  do k = 1, s
    fvec(4*k-3) = fs(k)%pth - si%pthold
    fvec(4*k-2) = x(4*k-2)  - si%z(2)
    fvec(4*k-1) = x(4*k-1)  - si%z(3)
    fvec(4*k)   = x(4*k)    - si%z(4)
    do l = 1, s
      fvec(4*k-3) = fvec(4*k-3) + si%dt*ahat(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
      fvec(4*k-2) = fvec(4*k-2) - si%dt*a(k,l)*Hprime(l)                                         ! thdot
      fvec(4*k-1) = fvec(4*k-1) - si%dt*ahat(k,l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
      fvec(4*k)   = fvec(4*k)   + si%dt*a(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))           ! pphdot
    end do
  end do
  
  end subroutine f_rk_lobatto_quasi

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_quasi(ierr)
  
  integer, intent(out) :: ierr

  select case (si%mode)
    case (1)
      call timestep_euler1_quasi(ierr)
    case (2)
        call timestep_euler2_quasi(ierr)
    case (3)
      call timestep_midpoint_quasi(ierr)
    case (4)
      call timestep_rk_gauss_quasi(1, ierr)
    case (5)
      call timestep_rk_gauss_quasi(2, ierr)
    case (6)
      call timestep_rk_gauss_quasi(3, ierr)
    case (7)
      call timestep_rk_gauss_quasi(4, ierr)
    case (15)
      call timestep_rk_lobatto_quasi(si, f, 3, ierr)
    case default
      print *, 'invalid mode for orbit_timestep_quasi: ', si%mode
      stop
  end select

end subroutine orbit_timestep_quasi

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_multi_quasi(mi, ierr)
  !
    type(MultistageIntegrator), intent(inout) :: mi
  
    integer, intent(out) :: ierr
  
    integer :: kstage
  
    si = mi%stages(1)
    call orbit_timestep_quasi(ierr)
  
    do kstage = 2, 2*mi%s
      mi%stages(kstage)%z = si%z
      mi%stages(kstage)%pthold = f%pth
      si = mi%stages(kstage)
      call orbit_timestep_quasi(ierr)
    end do
    mi%stages(1)%z = si%z
    mi%stages(1)%pthold = f%pth
  
  end subroutine orbit_timestep_multi_quasi

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine timestep_midpoint_quasi(ierr)
    
  integer, intent(out) :: ierr

  integer, parameter :: n = 5
  integer, parameter :: maxit = 256

  double precision, dimension(n) :: x
  double precision :: fvec(n)
  integer :: ktau, info

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x(1:4) = si%z
    x(5) = si%z(1)

    call hybrd1(f_midpoint_quasi, n, x, fvec, si%rtol, info)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    si%z = x(1:4)

    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine timestep_midpoint_quasi


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine timestep_euler1_quasi(ierr)
  !  
  integer, intent(out) :: ierr

  integer, parameter :: n = 2
  integer, parameter :: maxit = 256

  double precision, dimension(n) :: x
  double precision :: fvec(n)
  integer :: ktau, info

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x(1)=si%z(1)
    x(2)=si%z(4)

    call hybrd1(f_euler1_quasi, n, x, fvec, si%rtol, info)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    si%z(1) = x(1)
    si%z(4) = x(2)

    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives(f, si%z(4))

    si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
    si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph

    if (exact_steps) then
      call hybrd1(f_exact_quasi, 1, x, fvec, si%rtol, info)
      si%z(1) = x(1)
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    end if

    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine timestep_euler1_quasi



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine timestep_euler2_quasi(ierr)
!
  integer, intent(out) :: ierr

  integer, parameter :: n = 3
  integer, parameter :: maxit = 256

  double precision, dimension(n) :: x
  double precision :: fvec(n)
  integer :: ktau, info

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth

    x = si%z(1:3)

    call hybrd1(f_euler2_quasi, n, x, fvec, si%rtol, info)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    si%z(1:3) = x

    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives(f, si%z(4))

    f%pth = si%pthold - si%dt*(f%dH(2) - f%dH(1)*f%dpth(2)/f%dpth(1))
    si%z(4) = si%z(4) - si%dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

    if (exact_steps) then
      call hybrd1(f_exact_quasi, 1, x, fvec, si%rtol, info)
      si%z(1) = x(1)
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    end if

    si%kt = si%kt+1
    ktau = ktau+1
  enddo
  
end subroutine timestep_euler2_quasi

subroutine timestep_rk_gauss_quasi(s, ierr)
!
  integer, intent(out) :: ierr

  integer, parameter :: maxit = 256

  integer, intent(in) :: s
  double precision, dimension(4*s) :: x
  double precision :: fvec(4*s)
  integer :: ktau, info, k, l

  double precision :: a(s,s), b(s), c(s), Hprime(s)

  do k = 1,s
    fs(k) = f
  end do

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth
    
    do k = 1,s
      x((4*k-3):(4*k)) = si%z
    end do

    call hybrd1(f_rk_gauss_quasi, 4*s, x, fvec, si%rtol, info)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

    f = fs(s)
    f%pth = si%pthold
    si%z(1) = x(4*s-3)

    do l = 1, s
      Hprime(l) = fs(l)%dH(1)/fs(l)%dpth(1)
      f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))            ! pthdot
      si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)                                      ! thdot
      si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
    end do
  
    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine timestep_rk_gauss_quasi


subroutine timestep_rk_lobatto_quasi(s, ierr)
!
  integer, intent(out) :: ierr

  integer, parameter :: maxit = 256

  integer, intent(in) :: s
  double precision, dimension(4*s) :: x
  double precision :: fvec(4*s)
  integer :: ktau, info, k, l

  double precision :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s)

  do k = 1,s
    fs(k) = f
  end do

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    si%pthold = f%pth
    
    do k = 1,s
      x((4*k-3):(4*k)) = si%z
    end do

    call hybrd1(f_rk_gauss_quasi, 4*s, x, fvec, si%rtol, info)

    if (x(1) > 1.0) then
      ierr = 1
      return
    end if

    if (x(1) < 0.0) then
      print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
      x(1) = 0.01
    end if

    call coeff_rk_lobatto(s, a, ahat, b, c)  ! TODO: move this to preprocessing

    f = fs(s)
    f%pth = si%pthold
    si%z(1) = x(4*s-3)

    do l = 1, s
      Hprime(l) = fs(l)%dH(1)/fs(l)%dpth(1)
      f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))            ! pthdot
      si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)                                      ! thdot
      si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph  ! phdot
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
    end do
  
    si%kt = si%kt+1
    ktau = ktau+1
  enddo

end subroutine timestep_rk_lobatto_quasi

end module orbit_symplectic_quasi
