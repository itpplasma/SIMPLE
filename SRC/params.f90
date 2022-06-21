module params
  use util
  use parmot_mod, only : ro0, rmu
  use neo_orb, only: NeoOrb
  implicit none

  integer          :: npoi,L1i,nper,i,ntestpart
  integer          :: notrace_passing,loopskip,iskip
  double precision :: dphi,phibeg,bmod00,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi,v0
  double precision :: RT0,R0i,cbfi,bz0i,bf0,rbig
  double precision :: sbeg,thetabeg
  double precision, dimension(:),   allocatable :: bstart,volstart
  double precision, dimension(:,:), allocatable :: xstart
  double precision, dimension(:,:), allocatable :: zstart, zend
  double precision, dimension(:), allocatable :: confpart_trap,confpart_pass
  double precision, dimension(:), allocatable :: times_lost
  double precision :: contr_pp
  integer          :: ibins
  integer          :: n_e,n_d
  integer          :: startmode

  integer :: ntau ! number of dtaumin in dtau
  integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet

  integer :: kpart = 0 ! progress counter for particles

  double precision :: relerr

  double precision, allocatable :: trap_par(:), perp_inv(:)
  integer,          allocatable :: iclass(:,:)

  integer, parameter :: n_tip_vars = 6  ! variables to evaluate at tip: z(1..5), par_inv
  integer :: nplagr,nder,npl_half
  integer :: norbper,nfp
  double precision :: fper, zerolam = 0d0

  double precision :: tcut = -1d0
  integer :: ntcut
  logical          :: class_plot     !<=AAA
  double precision :: cut_in_per     !<=AAA
  logical          :: local=.false.

  logical :: fast_class=.true.  !if .true. quit immeadiately after fast classification

! colliding with D-T reactor plasma. TODO: Make configurable
  logical :: swcoll
  double precision, parameter :: am1=2.0d0, am2=3.0d0, Z1=1.0d0, Z2=1.0d0, &
    densi1=0.5d14, densi2=0.5d14, tempi1=1.0d4, tempi2=1.0d4, tempe=1.0d4
  double precision :: dchichi,slowrate,dchichi_norm,slowrate_norm
  logical :: deterministic = .False.

contains

  subroutine init_params(E_alpha, bmod_ref, trace_time, ntimstep, npoiper, &
      npoiper2)
    double precision, intent(in) :: E_alpha     ! Particle energy in eV
    double precision, intent(in) :: bmod_ref    ! Reference magnetic field in G
    double precision, intent(in) :: trace_time  ! Tracing time in seconds
    integer, intent(in) :: ntimstep  ! Number of times to record loss fraction
    integer, intent(in) :: npoiper   ! Number of field-line integration steps
                                     ! per field period to init flux surface
    integer, intent(in) :: npoiper2  ! Minimum number of integration steps per
                                     ! field period (i.e. strongly passing)

  ! set alpha energy, velocity, and Larmor radius
    v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
    rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)

  ! Neglect relativistic effects by large inverse relativistic temperature
    rmu=1d8

  ! normalized slowing down time:
    tau=trace_time*v0
  ! normalized time step:
    dtau=tau/dble(ntimstep-1)
  ! parameters for the vacuum chamber:
    call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0) ! TODO: why again?
    rbig=rt0
  ! field line integration step step over phi (to check chamber wall crossing)
    dphi=2.d0*pi/(L1i*npoiper)
  ! orbit integration time step (to check chamber wall crossing)
    dtaumin=2.d0*pi*rbig/npoiper2  ! ntimstep =
    ntau=ceiling(dtau/dtaumin)
    dtaumin=dtau/ntau

    ntcut = ceiling(ntimstep*ntau*tcut/trace_time)

    norbper=ceiling(1d0*ntau*ntimstep/(L1i*npoiper2))
    nfp=L1i*norbper         !<= guess for footprint number

    zerolam=0.d0
    nplagr=4
    nder=0
    npl_half=nplagr/2

    fper = 2d0*pi/dble(L1i)   !<= field period
  end subroutine init_params
end module params
