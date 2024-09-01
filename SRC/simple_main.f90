module simple_main
  use omp_lib
  use util, only: pi, twopi, sqrt2
  use simple, only : init_sympl
  use diag_mod, only : icounter
  use collis_alp, only : loacol_alpha, stost
  use binsrc_sub, only : binsrc
  use params, only: swcoll, ntestpart, startmode, num_surf, dtau, dtaumin, ntau, v0, &
    kpart, confpart_pass, confpart_trap, times_lost, integmode, relerr, trace_time, &
    class_plot, ntcut, iclass, bmod00, loopskip, xi, idx, bmin, bmax, dphi, xstart, &
    zstart, zend, trap_par, perp_inv, volstart, sbeg, thetabeg, phibeg, npoiper, nper, &
    ntimstep, bstart, ibins, ierr, Tracer, should_skip, reset_seed_if_deterministic

  implicit none

  contains

  subroutine run(norb)
    type(Tracer), intent(inout) :: norb
    integer :: i

    call print_parameters
    if (swcoll) call init_collisions

    ! pre-compute starting flux surface
    call init_starting_surf

    ! local?
    if(1 == num_surf) then
      call init_starting_points
    else
      call init_starting_points_global
    endif
    if (startmode == 0) stop 'startmode == 0, stopping after generating start.dat'

    call init_counters

    !$omp parallel firstprivate(norb)
    !$omp do
    do i = 1, ntestpart
      !$omp critical
      kpart = kpart+1
      print *, kpart, ' / ', ntestpart, 'particle: ', i, 'thread: ', omp_get_thread_num()
      !$omp end critical
      call trace_orbit(norb, i)
    end do
    !$omp end do
    !$omp end parallel

    confpart_pass=confpart_pass/ntestpart
    confpart_trap=confpart_trap/ntestpart

  end subroutine run

  subroutine print_parameters
    print *, 'tau: ', dtau, dtaumin, min(dabs(mod(dtau, dtaumin)), &
                      dabs(mod(dtau, dtaumin)-dtaumin))/dtaumin, ntau
    print *, 'v0 = ', v0
  end subroutine print_parameters

  subroutine init_collisions
    use params, only: am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
    facE_al, dchichi, slowrate, dchichi_norm, slowrate_norm, v0

    real(8) :: v0_coll

    call loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe, &
      3.5d6/facE_al,v0_coll,dchichi,slowrate,dchichi_norm,slowrate_norm)

    if (abs(v0_coll - v0) > 1d-6) then
      error stop 'simple_main.init_collisions: v0_coll != v0'
    end if
  end subroutine init_collisions

  subroutine init_counters
    icounter=0 ! evaluation counter
    kpart=0

    ! initialize array of confined particle percentage
    confpart_trap=0.d0
    confpart_pass=0.d0

    ! initialize lost times when particles get lost
    times_lost = -1.d0
  end subroutine init_counters

  subroutine init_starting_surf
    use alpha_lifetime_sub, only : integrate_mfl_can

    xstart=0.d0
    bstart=0.d0
    volstart=0.d0 !ToDo add loop for all sbeg, check dimension

    call integrate_mfl_can( &
      npoiper*nper,dphi,sbeg(1),phibeg,thetabeg, &
      xstart,bstart,volstart,bmod00,ierr)

    if(ierr.ne.0) then
      print *,'starting field line has points outside the chamber'
      stop
    endif

  ! maximum value of B module:
    bmax=maxval(bstart)
    bmin=minval(bstart)

    print *, 'bmod00 = ', bmod00, 'bmin = ', bmin, 'bmax = ', bmax
  end subroutine init_starting_surf

  subroutine init_starting_points_ants(unit)
    use parse_ants, only : process_line
    use get_can_sub, only : vmec_to_can

    integer, intent(in) :: unit

    integer, parameter :: maxlen = 4096
    character(len=maxlen) :: line
    real(8) :: v_par, v_perp, u, v, s
    real(8) :: th, ph
    integer :: ipart

    do ipart=1,ntestpart
      read(unit, '(A)') line
      call process_line(line, v_par, v_perp, u, v, s)
      ! In the test case, u runs from 0 to 1 and v from 0 to 4
      th = 2d0*pi*u
      ph = 2d0*pi*v/4d0
      zstart(1, ipart) = s
      zstart(2, ipart) = th
      zstart(3, ipart) = ph
      zstart(4, ipart) = 1.d0
      zstart(5, ipart) = v_par / sqrt(v_par**2 + v_perp**2)
    enddo
  end subroutine

  subroutine init_starting_points
    use get_can_sub, only: can_to_vmec
    use samplers, only: START_FILE

    integer :: i, ipart, iskip

    ! skip random numbers according to configuration
      do iskip=1,loopskip
        do ipart=1,ntestpart
          call random_number(xi)
          call random_number(xi)
        enddo
      enddo

    ! files for storing starting coords
    open(1,file=START_FILE,recl=1024)
    ! determine the starting point:
    if (startmode == 0 .or. startmode == 1) then
      do ipart=1,ntestpart
        call random_number(xi)
        call binsrc(volstart,1,npoiper*nper,xi,i)
        ibins=i
        ! we store starting points in lab coordinates:
        call to_lab_coordinates(xstart(1:3,i), zstart(1:3,ipart))
        ! normalized velocity module z(4) = v / v_0:
        zstart(4,ipart)=1.d0
        ! starting pitch z(5)=v_\parallel / v:
        call random_number(xi)
        zstart(5,ipart)=2.d0*(xi-0.5d0)
        write(1,*) zstart(:,ipart)
      enddo
    else if (startmode == 2) then
      do ipart=1,ntestpart
        read(1,*) zstart(:,ipart)
      enddo
    else if (startmode == 3) then  ! ANTS input mode
      call init_starting_points_ants(1)
    else if (startmode == 4) then
      !select only the indices from batch and overwrite zstart.
      do ipart=idx(0),idx(ntestpart)
        read(1,*) zstart(:,ipart)
      enddo
      ! indices no longer needed
      deallocate(idx)
    endif

    close(1)
  end subroutine init_starting_points

  subroutine init_starting_points_global

    use find_bminmax_sub, only : get_bminmax
    use get_can_sub, only: can_to_vmec
    use samplers, only: START_FILE


    integer, parameter :: ns=1000
    integer :: ipart,iskip,is,s_idx,parts_per_s
    double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec
    double precision :: s,bmin,bmax
  !
    open(1,file='bminmax.dat',recl=1024)
    do is=0,ns
      s=dble(is)/dble(ns)
  !
      call get_bminmax(s,bmin,bmax)
  !
      write(1,*) s,bmin,bmax
    enddo
    close(1)

    ! skip random numbers according to configuration
      do iskip=1,loopskip
        do ipart=1,ntestpart
          call random_number(xi)
          call random_number(xi)
          call random_number(xi)
          call random_number(xi)
        enddo
      enddo

    ! files for storing starting coords
    open(1,file=START_FILE,recl=1024)
    ! determine the starting point:
    if (startmode == 0 .or. startmode == 1) then
      print *, "Initialising for", num_surf, "surfaces."
      do ipart=1,ntestpart
        if (0 == num_surf) then
          call random_number(xi)
          r = xi
        else if (1 < num_surf) then
          parts_per_s = int(ntestpart/num_surf)
          s_idx = (ipart/parts_per_s)+1
          r = sbeg(s_idx)
        else ! Should not happen (as we are not in "local mode"), however let's catch it anyway.
          r = sbeg(1)
        endif

        call random_number(xi)
        vartheta=twopi*xi
        call random_number(xi)
        varphi=twopi*xi
        ! we store starting points in lab coordinates:
        call to_lab_coordinates([r, vartheta, varphi], zstart(1:3,ipart))
        ! normalized velocity module z(4) = v / v_0:
        zstart(4,ipart)=1.d0
        ! starting pitch z(5)=v_\parallel / v:
        call random_number(xi)
        zstart(5,ipart)=2.d0*(xi-0.5d0)
        write(1,*) zstart(:,ipart)
      enddo
    else if (startmode == 2) then
      do ipart=1,ntestpart
        read(1,*) zstart(:,ipart)
      enddo
    else if (startmode == 3) then  ! ANTS input mode
      call init_starting_points_ants(1)
    else if (startmode == 4) then
      !select only the indices from batch and overwrite zstart.
      do ipart=idx(0),idx(ntestpart)
        read(1,*) zstart(:,ipart)
      enddo
      ! indices no longer needed
      deallocate(idx)
    endif

    close(1)
  end subroutine init_starting_points_global


  subroutine trace_orbit(anorb, ipart)
    use classification, only : trace_orbit_with_classifiers
    use callback, only : callbacks_macrostep

    type(Tracer), intent(inout) :: anorb
    integer, intent(in) :: ipart

    double precision, dimension(5) :: z
    integer :: it, ierr_orbit
    integer(8) :: kt
    logical :: passing

    ierr_orbit = 0

    call reset_seed_if_deterministic

    if (ntcut>0 .or. class_plot) then
      call trace_orbit_with_classifiers(anorb, ipart)
      return
    endif

    call from_lab_coordinates(zstart(1:3, ipart), z(1:3))
    z(4:5) = zstart(4:5, ipart)
    zend(:,ipart) = 0d0

    if (integmode>0) call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)

    call compute_pitch_angle_params(z, passing, trap_par(ipart), perp_inv(ipart))

    if(passing .and. should_skip(ipart)) then
      !$omp critical
      confpart_pass=confpart_pass+1.d0
      !$omp end critical
      return
    endif

    kt = 0
    do it = 1, ntimstep
      if (it >= 2) call macrostep(anorb, z, kt, ierr_orbit)
      call callbacks_macrostep(anorb, ipart, it, kt*dtaumin/v0, z, ierr_orbit)
      if(ierr_orbit .ne. 0) exit
      call increase_confined_count(it, passing)
    enddo

    !$omp critical
    call to_lab_coordinates(z(1:3), zend(1:3,ipart))
    zend(4:5, ipart) = z(4:5)
    times_lost(ipart) = kt*dtaumin/v0
    !$omp end critical
  end subroutine trace_orbit

  subroutine macrostep(anorb, z, kt, ierr_orbit)
    use alpha_lifetime_sub, only : orbit_timestep_axis
    use orbit_symplectic, only : orbit_timestep_sympl

    type(Tracer), intent(inout) :: anorb
    double precision, intent(inout) :: z(5)
    integer(8), intent(inout) :: kt
    integer, intent(out) :: ierr_orbit

    integer :: ktau

    do ktau=1,ntau
      if (integmode <= 0) then
        call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr_orbit)
      else
        if (swcoll) call update_momentum(anorb, z)
        call orbit_timestep_sympl(anorb%si, anorb%f, ierr_orbit)
        call to_standard_z_coordinates(anorb, z)
      endif
      if (swcoll) call collide(z, dtaumin) ! Collisions
      if (ierr_orbit .ne. 0) exit
      kt = kt+1
    enddo
  end subroutine macrostep

  subroutine to_standard_z_coordinates(anorb, z)
    type(Tracer), intent(in) :: anorb
    double precision, intent(inout) :: z(5)

    z(1:3) = anorb%si%z(1:3)
    z(4) = dsqrt(anorb%f%mu*anorb%f%Bmod+0.5d0*anorb%f%vpar**2)
    z(5) = anorb%f%vpar/(z(4)*sqrt2)
  end subroutine to_standard_z_coordinates

  subroutine from_lab_coordinates(xlab, x)
    use field_can_mod, only : ref_to_can
    double precision, intent(in) :: xlab(:)
    double precision, intent(inout) :: x(:)

    call ref_to_can(xlab, x)  ! TODO don't assume lab=ref
  end subroutine from_lab_coordinates

  subroutine to_lab_coordinates(x, xlab)
    use field_can_mod, only : can_to_ref
    double precision, intent(in) :: x(:)
    double precision, intent(inout) :: xlab(:)

    call can_to_ref(x, xlab)  ! TODO don't assume lab=ref
  end subroutine to_lab_coordinates

  subroutine increase_confined_count(it, passing)
    integer, intent(in) :: it
    logical, intent(in) :: passing

    !$omp critical
    if (passing) then
      confpart_pass(it) = confpart_pass(it) + 1.d0
    else
      confpart_trap(it) = confpart_trap(it) + 1.d0
    end if
    !$omp end critical
  end subroutine increase_confined_count

  subroutine compute_pitch_angle_params(z, passing, trap_par_, perp_inv_)
    use find_bminmax_sub, only : get_bminmax

    double precision, intent(in) :: z(5)
    logical, intent(out) :: passing
    double precision, intent(out) :: trap_par_, perp_inv_

    double precision :: bmod

    !$omp critical
    bmod = compute_bmod(z(1:3))
    if(num_surf > 1) then
      call get_bminmax(z(1),bmin,bmax)
    endif
    passing = z(5)**2.gt.1.d0-bmod/bmax
    trap_par_ = ((1.d0-z(5)**2)*bmax/bmod-1.d0)*bmin/(bmax-bmin)
    perp_inv_ = z(4)**2*(1.d0-z(5)**2)/bmod
    !$omp end critical
  end subroutine compute_pitch_angle_params

  function compute_bmod(z) result(bmod)
    use magfie_sub, only : magfie

    double precision :: bmod
    double precision, intent(in) :: z(3)

    double precision :: sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

    call magfie(z(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  end function compute_bmod

  subroutine update_momentum(anorb, z)
    use orbit_symplectic, only : get_val

    type(Tracer), intent(inout) :: anorb
    double precision, intent(in) :: z(5)

    anorb%si%pabs = z(4)
    anorb%f%vpar = z(4)*z(5)*sqrt2
    anorb%f%mu = z(4)**2*(1.d0-z(5)**2)/anorb%f%Bmod
    anorb%si%z(4) = anorb%f%vpar*anorb%f%hph + anorb%f%Aph/anorb%f%ro0
    call get_val(anorb%f, anorb%si%z(4)) ! for pth
    anorb%si%pthold = anorb%f%pth
  end subroutine update_momentum

  subroutine collide(z, dt)
    double precision, intent(in) :: z(5), dt
    integer :: ierr_coll

    call stost(z, dt, 1, ierr_coll)
    if (ierr_coll /= 0) then
      print *, 'Error in stost: ', ierr_coll, 'z = ', z, 'dtaumin = ', dtaumin
    endif
  end subroutine collide


  subroutine write_output

    integer :: i, num_lost
    double precision :: inverse_times_lost_sum

    open(1,file='times_lost.dat',recl=1024)
    num_lost = 0
    inverse_times_lost_sum = 0.0d0
    do i=1,ntestpart
      write(1,*) i, times_lost(i), trap_par(i), zstart(1,i), perp_inv(i), zend(:,i)
      if (times_lost(i) > 0.0d0 .and. times_lost(i) < trace_time) then
        num_lost = num_lost + 1
        inverse_times_lost_sum = inverse_times_lost_sum + 1.0/times_lost(i)
      end if
    enddo
    close(1)

    if (num_lost > 0) then
      open(1,file='avg_inverse_t_lost.dat',recl=1024) ! Write average loss time
      write(1,*) inverse_times_lost_sum/num_lost
      close(1)
    endif

    open(1,file='confined_fraction.dat',recl=1024)
    do i=1,ntimstep
      write(1,*) dble(i-1)*dtau/v0,confpart_pass(i),confpart_trap(i),ntestpart
    enddo
    close(1)

    if (ntcut>0 .or. class_plot) then
        open(1,file='class_parts.dat',recl=1024)
        do i=1,ntestpart
        write(1,*) i, zstart(1,i), perp_inv(i), iclass(:,i)
        enddo
        close(1)
    endif

  end subroutine write_output

end module simple_main
