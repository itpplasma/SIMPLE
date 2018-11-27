!
  use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
  use chamb_mod,  only : rnegflag
  use parmot_mod, only : rmu,ro0,eeff
  use velo_mod,   only : isw_field_type
use diag_mod, only : dodiag
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
  double precision,parameter  :: c=2.9979d10
  double precision,parameter  :: e_charge=4.8032d-10
  double precision,parameter  :: e_mass=9.1094d-28
  double precision,parameter  :: p_mass=1.6726d-24
  double precision,parameter  :: ev=1.6022d-12
!
  integer          :: npoi,ierr,L1i,nper,npoiper,i,ntimstep,ntestpart
  integer          :: ipart,notrace_passing,loopskip,iskip,ilost
  real             :: zzg
  double precision :: dphi,rbeg,phibeg,zbeg,bmod00,rcham,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: RT0,R0i,cbfi,bz0i,bf0,trap_par,rbig
  double precision :: sbeg,thetabeg
  double precision, dimension(5) :: z
  double precision, dimension(:),   allocatable :: bstart,volstart,confpart
  double precision, dimension(:,:), allocatable :: xstart
  double precision, dimension(:), allocatable :: confpart_trap,confpart_pass
! 24.03.2016
  integer          :: npoiper2
! 24.03.2016 end
! 13.02.2013
  double precision :: contr_pp
! 14.04.2013
  double precision :: facE_al
  integer          :: ibins
! 14.04.2013 end
! 22.09.2013
  integer          :: n_e,n_d,n_b
! 22.09.2013 end
! 13.02.2013 end
!
! inverse relativistic temperature
  rmu=1d5
! alpha particle energy, eV:
! 12.11.2011  E_alpha=3.5d6
!  E_alpha=3.5d6/16d0
! alpha particle velocity, cm/s
!  v0=sqrt(2.d0*E_alpha*ev/(4.d0*p_mass))
!
  open(1,file='alpha_lifetime_c.inp')
  read (1,*) notrace_passing   !skip tracing passing prts if notrace_passing=1
  read (1,*) nper              !number of periods for initial field line
  read (1,*) npoiper           !number of points per period on this field line
  read (1,*) ntimstep          !number of time steps per slowing down time
  read (1,*) ntestpart         !number of test particles
!<=2017  read (1,*) rcham             !vacuum chamber radius, cm
  read (1,*) bmod_ref          !reference field, G, for Boozer $B_{00}$
  read (1,*) trace_time        !slowing down time, s
  read (1,*) sbeg              !starting s for field line                       !<=2017
  read (1,*) phibeg            !starting phi for field line                     !<=2017
  read (1,*) thetabeg          !starting theta for field line                   !<=2017
  read (1,*) loopskip          !how many loops to skip to shift random numbers
! 13.02.2013
  read (1,*) contr_pp          !control of passing particle fraction
! 13.02.2013 end
! 14.04.2013
  read (1,*) facE_al           !facE_al test particle energy reduction factor
! 14.04.2013
! 24.03.2016
  read (1,*) npoiper2          !additional integration step split factor
! 24.03.2016 end
! 22.09.2013
  read (1,*) n_e               !test particle charge number (the same as Z)
  read (1,*) n_d               !test particle mass number (the same as A)
! 22.09.2013 end
! 26.09.2013
!<=2017  read (1,*) n_b               !determination bmax
  read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
  read (1,*) ns_s              !spline order for 3D quantities over s variable
  read (1,*) ns_tp             !spline order for 3D quantities over theta and phi
  read (1,*) multharm          !angular grid factor (n_grid=multharm*n_harm_max where n_harm_max - maximum Fourier index)
  close(1)
!
!
  call spline_vmec_data
!
  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0) 
!
  call get_canonical_coordinates
!
  isw_field_type=0
!
  z=0.d0                        !<=2017
  call chamb_can(z(1:2),z(3),ierr)  !<=2017
!
! 14.04.2013
! alpha particle energy, eV:
!  E_alpha=3.5d6/16d0
  E_alpha=3.5d6/facE_al
! alpha particle velocity, cm/s
! 22.09.2013  v0=sqrt(2.d0*E_alpha*ev/(4.d0*p_mass))
  v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
! 14.04.2013 end
!
! Larmor radius:
!  rlarm=v0*4.d0*p_mass*c/(e_charge*bmod_ref)
! 22.09.2013  rlarm=v0*2.d0*p_mass*c/(e_charge*bmod_ref)
  rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)
! normalized slowing down time:
  tau=trace_time*v0
! normalized time step:
  dtau=tau/dfloat(ntimstep-1)
!
! 14.11.2011  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)
  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)         !<=2017
!<=2017  call tj2vvo(RT0,R0i,L1i,cbfi,bz0i,bf0)
! 05.03.2014 for w7x
! 13.07.2016  rt0=145d0
! 13.07.2016  L1i=3
! 05.03.2014 end
!
! parameters for the vacuum chamber:
  rbig=rt0
! field line integration step step over phi (to check chamber wall crossing)
  dphi=2.d0*pi/(L1i*npoiper)
! orbit integration time step (to check chamber wall crossing)
! 24.03.2016  dtaumin=dphi*rbig
  dtaumin=dphi*rbig/npoiper2!
!  open(1,status='new',file='alpha_lifetime.log')
  open(1,file='alpha_lifetime.log')
  write (1,*) 'notrace_passing = ',notrace_passing
  write (1,*) 'nper = ',nper
  write (1,*) 'npoiper = ',npoiper
  write (1,*) 'ntimestep = ',ntimstep
  write (1,*) 'ntestpart = ',ntestpart
!<=2017  write (1,*) 'rcham = ',rcham
  write (1,*) 'bmod_ref = ',bmod_ref
  write (1,*) 'trace_time = ',trace_time
  write (1,*) 'sbeg = ',sbeg                !<=2017
  write (1,*) 'phibeg = ',phibeg            !<=2017
  write (1,*) 'thetabeg = ',thetabeg        !<=2017
  write (1,*) 'Rbig = ',rbig
  write (1,*) 'dphi = ',dphi
  write (1,*) 'v0 = ',v0
  write (1,*) 'rlarm = ',rlarm
  write (1,*) 'dtau = ',dtau
  write (1,*) 'dtaumin = ',dtaumin
! 16.11.2011
  write (1,*) 'E_alpha = ',E_alpha
! 13.02.2013
  write (1,*) 'contr_pp =', contr_pp
! 13.02.2013 end
! 22.09.2013 
  write (1,*) 'n_e = ',n_e
  write (1,*) 'n_d = ',n_d
! 22.09.2013 end
! 26.09.2013 
!<=2017  write (1,*) 'n_b = ',n_b
  write (1,*) 'ns_s = ',ns_s
  write (1,*) 'ns_tp = ',ns_tp
  write (1,*) 'multharm = ',multharm
  close(1)
!
! total number of starting points:
  npoi=nper*npoiper
  allocate(xstart(3,npoi),bstart(npoi),volstart(npoi))
  xstart=0.d0
  bstart=0.d0
  volstart=0.d0
!
! pre-compute starting points:
!
! 23.09.2013
! 05.03.2014  go to 110

  call integrate_mfl_can(npoi,dphi,sbeg,phibeg,thetabeg,          &
                         xstart,bstart,volstart,bmod00,ierr)
!
  if(ierr.ne.0) then
    print *,'starting field line has points outside the chamber'
    stop
  endif
!
! Larmor raidus corresponds to the field stregth egual to $B_{00}$ harmonic
! in Boozer coordinates:
! 14.11.2011  bmod00=bmod_ref  !<=deactivated, use value from the 'alpha_lifetime.inp'
  ro0=rlarm*bmod00  ! 23.09.2013
! maximum value of B module:
  bmax=maxval(bstart)
  bmin=minval(bstart)
!
! 23.09.2013
! 05.03.2014  110  continue
!  open(1,file='bminmax.d')
!  read(1,*)bmin,bmax,bmod00
!  close(1)
!  open(1,file='xstart.d')
!  open(2,file='line.d')
!  do i=1,npoi
!  read(1,*)xstart(1,i),xstart(2,i),xstart(3,i),bstart(i)
!  read(2,*)ibins,volstart(i)
!  end do
!  close(1)
!  close(2)
! 05.03.2014! 25.09.2013
! 05.03.2014! 26.09.2013  bmax=maxval(bstart)
! 05.03.2014! 26.09.2013  bmin=minval(bstart)
! 05.03.2014  if(n_b.eq.0) then
! 05.03.2014  bmax=maxval(bstart)
! 05.03.2014  bmin=minval(bstart)
! 05.03.2014  end if
  open(1,file='bminmax.dat')
  write(1,*)bmin,bmax,bmod00
  close(1)
  open(1,file='volstart.dat')
  do i=1,npoi
  write(1,*)i,volstart(i)
  end do
  close(1)
! 05.03.2014! 25.09.2013 end
! 05.03.2014  ro0=rlarm*bmod00
! 23.09.2013 end
  open(1,file='starting_surface.dat')
  do i=1,npoi,npoiper
    write (1,*) xstart(1,i),xstart(3,i)
  enddo
  close(1)
!<=2017  open(1,file='vacuum_chamber.dat')
!<=2017  do i=0,npoi
!<=2017    write (1,*) rbig+rcham*cos(2.d0*pi*dfloat(i)/dfloat(npoi)), &
!<=2017                rcham*sin(2.d0*pi*dfloat(i)/dfloat(npoi))
!<=2017  enddo
!<=2017  close(1)
!
  allocate(confpart_trap(ntimstep),confpart_pass(ntimstep))
  confpart_trap=0.d0
  confpart_pass=0.d0
!
  do iskip=1,loopskip
    do ipart=1,ntestpart
      xi=zzg()
      xi=zzg()
    enddo
  enddo
!
  open(2,file='start_p.dat')
  open(3,file='start_t.dat')
  do ipart=1,ntestpart
print *,ipart,' / ',ntestpart
! determine the starting point:
    xi=zzg()
!
    call binsrc(volstart,1,npoi,xi,i)
! 14.04.2013 
  ibins=i
! 14.04.2013 end
! coordinates: z(1) = R, z(2) = phi, z(3) = Z
    z(1:3)=xstart(:,i)
! normalized velocity module z(4) = v / v_0:
    z(4)=1.d0
! starting pitch z(5)=v_\parallel / v:
    xi=zzg()
    z(5)=2.d0*(xi-0.5d0)
!
    trap_par=((1.d0-z(5)**2)*bmax/bstart(i)-1.d0)*bmin/(bmax-bmin)
    if(z(5)**2.gt.1.d0-bstart(i)/bmax) then
! passing particle
      if(notrace_passing.eq.1) then
! no tracing of passing particles, assume that all are confined
        confpart_pass=confpart_pass+1.d0
        cycle
      endif
! trace passing particle
      confpart_pass(1)=confpart_pass(1)+1.d0
      ilost=ntimstep-1
      do i=2,ntimstep
!
! 13.02.2013
      if(trap_par.le.contr_pp) go to 111
! 13.02.2013 end
!
        call orbit_timestep_axis(z,dtau,dtaumin,ierr)
!
        if(ierr.ne.0) exit
! 13.02.2013
  111  continue
! 13.02.2013 end
        ilost=ntimstep-i
! 26.03.2016
print *,'passing particle ',ipart,' step ',i,' of ',ntimstep
! 26.03.2016 end
        confpart_pass(i)=confpart_pass(i)+1.d0
      enddo
! 14.04.2013      write(2,*) trap_par,ilost,xstart(:,i)
      write(2,*) trap_par,ilost,xstart(:,ibins),ierr
! 14.04.2013 end
      close(2)
      open(2,file='start_p.dat',position='append')
    else
! trapped particle (traced always)
      confpart_trap(1)=confpart_trap(1)+1.d0
      ilost=ntimstep-1
      do i=2,ntimstep
!
        call orbit_timestep_axis(z,dtau,dtaumin,ierr)
!
        if(ierr.ne.0) exit
        ilost=ntimstep-i
! 26.03.2016
print *,'trapped particle ',ipart,' step ',i,' of ',ntimstep
if(ipart.eq.15.and.i.eq.7633) dodiag=.true.
! 26.03.2016 end
        confpart_trap(i)=confpart_trap(i)+1.d0
      enddo
! 14.04.2013      write(3,*) trap_par,ilost,xstart(:,i)
      write(3,*) trap_par,ilost,xstart(:,ibins),ierr
! 14.04.2013 end
      close(3)
      open(3,file='start_t.dat',position='append')
    endif
    open(1,file='confined_fraction.dat')
    do i=1,ntimstep
!      write(1,*) dfloat(i-1)*dtau/v0,confpart_pass(i)/ntestpart, &
!                                     confpart_trap(i)/ntestpart
      write(1,*) dfloat(i-1)*dtau/v0,confpart_trap(i)/ntestpart, &
                                     confpart_pass(i)/ntestpart,ipart
    enddo
    close(1)
  enddo
  close(2)
  close(3)
  confpart_pass=confpart_pass/ntestpart
  confpart_trap=confpart_trap/ntestpart
! 
  open(1,file='confined_fraction.dat')
  do i=1,ntimstep
    write(1,*) dfloat(i-1)*dtau/v0,confpart_pass(i),confpart_trap(i),ntestpart
  enddo
  close(1)
!
  call deallocate_can_coord   !<=2017 NEW
!
  end
