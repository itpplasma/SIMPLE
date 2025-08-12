module orbit_symplectic_base
use field_can_mod, only: eval_field => evaluate, FieldCan, get_val, get_derivatives, &
  get_derivatives2

implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)

logical, parameter :: extrap_field = .True.  ! do extrapolation after final iteration

  ! Integration methods
integer, parameter :: RK45 = 0, EXPL_IMPL_EULER = 1, IMPL_EXPL_EULER = 2, &
  MIDPOINT = 3, GAUSS1 = 4, GAUSS2 = 5, GAUSS3 = 6, GAUSS4 = 7, LOBATTO3 = 15

type :: SymplecticIntegrator
  real(dp) :: atol
  real(dp) :: rtol

  ! Current phase-space coordinates z and old pth
  real(dp), dimension(4) :: z  ! z = (r, th, ph, pphi)
  real(dp) :: pthold

  ! Timestep and variables from z0
  integer :: ntau
  real(dp) :: dt
  real(dp) :: pabs
end type SymplecticIntegrator

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! Composition method with 2s internal stages according to Hairer, 2002 V.3.1
  !
integer, parameter :: S_MAX = 32
type :: MultistageIntegrator
  integer :: s
  real(dp) :: alpha(S_MAX), beta(S_MAX)
  type(SymplecticIntegrator) stages(2*S_MAX)
end type MultistageIntegrator

abstract interface
  subroutine orbit_timestep_sympl_i(si, f, ierr)
    import :: SymplecticIntegrator, FieldCan
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    integer, intent(out) :: ierr
  end subroutine orbit_timestep_sympl_i
end interface

abstract interface
  subroutine orbit_timestep_quasi_i(ierr)
    integer, intent(out) :: ierr
  end subroutine orbit_timestep_quasi_i
end interface

contains

subroutine coeff_rk_gauss(n, a, b, c)
  integer, intent(in) :: n
  real(dp), intent(inout) :: a(n,n), b(n), c(n)

  if (n == 1) then
    a(1,1) = 0.5d0
    b(1) = 1.0d0
    c(1) = 0.5d0
  elseif (n == 2) then
    a(1,1) =  0.25d0
    a(1,2) = -0.038675134594812d0
    a(2,1) =  0.538675134594812d0
    a(2,2) =  0.25d0

    b(1) = 0.5d0
    b(2) = 0.5d0

    c(1) = 0.211324865405187d0
    c(2) = 0.788675134594812d0
  elseif (n == 3) then
    a(1,1) =  0.1388888888888889d0
    a(1,2) = -0.03597666752493894d0
    a(1,3) =  0.009789444015308318d0
    a(2,1) =  0.3002631949808646d0
    a(2,2) =  0.2222222222222222d0
    a(2,3) = -0.022485417203086805d0
    a(3,1) = 0.26798833376246944d0
    a(3,2) = 0.48042111196938336d0
    a(3,3) = 0.1388888888888889d0

    b(1) = 0.2777777777777778d0
    b(2) = 0.4444444444444444d0
    b(3) = 0.2777777777777778d0

    c(1) = 0.1127016653792583d0
    c(2) = 0.5d0
    c(3) = 0.8872983346207417d0
  elseif (n == 4) then  ! with help of coefficients from GeometricIntegrators.jl of Michael Kraus
    a(1,1) = 0.086963711284363462428182d0
    a(1,2) = -0.026604180084998794303397d0
    a(1,3) = 0.012627462689404725035280d0
    a(1,4) = -0.003555149685795683332096d0

    a(2,1) = 0.188118117499868064967927d0
    a(2,2) = 0.163036288715636523694030d0
    a(2,3) = -0.027880428602470894855481d0
    a(2,4) = 0.006735500594538155853808d0

    a(3,1) = 0.167191921974188778543535d0
    a(3,2) = 0.353953006033743966529670d0
    a(3,3) = 0.163036288715636523694030d0
    a(3,4) = -0.014190694931141143581010d0

    a(4,1) = 0.177482572254522602550608d0
    a(4,2) = 0.313445114741868369190314d0
    a(4,3) = 0.352676757516271865977586d0
    a(4,4) = 0.086963711284363462428182d0

    b(1) = 0.173927422568726924856364d0
    b(2) = 0.326072577431273047388061d0
    b(3) = 0.326072577431273047388061d0
    b(4) = 0.173927422568726924856364d0

    c(1) = 0.069431844202973713731097d0
    c(2) = 0.330009478207571871344328d0
    c(3) = 0.669990521792428128655672d0
    c(4) = 0.930568155797026341780054d0
  else
    ! not implemented
    a = 0d0
    b = 0d0
    c = 0d0
  endif
end subroutine coeff_rk_gauss


subroutine coeff_rk_lobatto(n, a, ahat, b, c)
  integer, intent(in) :: n
  real(dp), intent(inout) :: a(n,n), ahat(n,n), b(n), c(n)

  if (n == 3) then
    a(1,1) =  0d0
    a(1,2) =  0d0
    a(1,3) =  0d0

    a(2,1) =  0.20833333333333334d0
    a(2,2) =  0.33333333333333333d0
    a(2,3) = -0.041666666666666664d0

    a(3,1) =  0.16666666666666667d0
    a(3,2) =  0.66666666666666667d0
    a(3,3) =  0.16666666666666667d0

    ahat(1,1) =  0.16666666666666667d0
    ahat(1,2) = -0.16666666666666667d0
    ahat(1,3) =  0d0

    ahat(2,1) =  0.16666666666666667d0
    ahat(2,2) =  0.33333333333333333d0
    ahat(2,3) =  0d0

    ahat(3,1) =  0.16666666666666667d0
    ahat(3,2) =  0.83333333333333333d0
    ahat(3,3) =  0d0

    b(1) =  0.16666666666666667d0
    b(2) =  0.66666666666666667d0
    b(3) =  0.16666666666666667d0

    c(1) = 0d0
    c(2) = 0.5d0
    c(3) = 1.0d0

  else
    ! not implemented
    a = 0d0
    b = 0d0
    c = 0d0
  endif
end subroutine coeff_rk_lobatto


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! Lobatto (IIIA)-(IIIB) Runge-Kutta method with s internal stages (n=4*s variables)
  !
subroutine f_rk_lobatto(si, fs, s, x, fvec, jactype)
  !
  type(SymplecticIntegrator), intent(inout) :: si
  integer, intent(in) :: s
  type(FieldCan), intent(inout) :: fs(:)
  real(dp), intent(in) :: x(4*s)  ! = (rend, thend, phend, pphend)
  real(dp), intent(out) :: fvec(4*s)
  integer, intent(in) :: jactype  ! 0 = no second derivatives, 2 = second derivatives

  real(dp) :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s)
  integer :: k,l  ! counters

  call coeff_rk_lobatto(s, a, ahat, b, c)

  call eval_field(fs(1), x(1), si%z(2), si%z(3), jactype)
  call get_derivatives(fs(1), x(2))

  do k = 2, s
    call eval_field(fs(k), x(4*k-3-2), x(4*k-2-2), x(4*k-1-2), jactype)
    call get_derivatives(fs(k), x(4*k-2))
  end do

  Hprime = fs%dH(1)/fs%dpth(1)

  fvec(1) = fs(1)%pth - si%pthold
  fvec(2) = x(2) - si%z(4)

  do l = 1, s
      fvec(1) = fvec(1) + si%dt*ahat(1,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
      fvec(2) = fvec(2) + si%dt*ahat(1,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
  end do

  do k = 2, s
    fvec(4*k-3-2) = fs(k)%pth - si%pthold
    fvec(4*k-2-2) = x(4*k-2-2)  - si%z(2)
    fvec(4*k-1-2) = x(4*k-1-2)  - si%z(3)
    fvec(4*k-2)   = x(4*k-2)    - si%z(4)
  end do

  do l = 1, s
    do k = 2, s
      fvec(4*k-3-2) = fvec(4*k-3-2) + si%dt*ahat(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))        ! pthdot
      fvec(4*k-2-2) = fvec(4*k-2-2) - si%dt*a(k,l)*Hprime(l)                                         ! thdot
      fvec(4*k-1-2) = fvec(4*k-1-2) - si%dt*a(k,l)*(fs(l)%vpar  - Hprime(l)*fs(l)%hth)/fs(l)%hph     ! phdot
      fvec(4*k-2)   = fvec(4*k-2)   + si%dt*ahat(k,l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))        ! pphdot
    end do
  end do

end subroutine f_rk_lobatto

end module orbit_symplectic_base
