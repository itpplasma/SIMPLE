module simple_main
    use omp_lib
    use util, only: pi, twopi, c, e_charge, e_mass, p_mass, ev, sqrt2
    use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp, &
      vmec_B_scale, vmec_RZ_scale

    use velo_mod,   only : isw_field_type
    use orbit_symplectic, only : orbit_timestep_sympl, get_val
    use simple, only : init_field, init_sympl, eval_field
    use cut_detector, only : fract_dimension
    use diag_mod, only : icounter
    use collis_alp, only : loacol_alpha, stost
    use params
    use binsrc_sub, only : binsrc
    use boozer_sub, only : vmec_to_boozer, boozer_to_vmec
    use check_orbit_type_sub, only : check_orbit_type

    implicit none

    public

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
    open(1,file='avg_inverse_t_lost.dat',recl=1024) ! Write average loss time
    write(1,*) inverse_times_lost_sum/num_lost
    close(1)

    open(1,file='confined_fraction.dat',recl=1024)
    do i=1,ntimstep
      write(1,*) dble(i-1)*dtau/v0,confpart_pass(i),confpart_trap(i),ntestpart
    enddo
    close(1)

    open(1,file='class_parts.dat',recl=1024)
    do i=1,ntestpart
      write(1,*) i, zstart(1,i), perp_inv(i), iclass(:,i)
    enddo
    close(1)

  end subroutine write_output

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
    real(8) :: th, ph, th_c, ph_c  ! Canonical flux coordinate angles
    integer :: ipart

    do ipart=1,ntestpart
      read(unit, '(A)') line
      call process_line(line, v_par, v_perp, u, v, s)
      ! In the test case, u runs from 0 to 1 and v from 0 to 4
      th = 2d0*pi*u
      ph = 2d0*pi*v/4d0
      call vmec_to_can(s, th, ph, th_c, ph_c)
      zstart(1, ipart) = s
      zstart(2, ipart) = ph_c
      zstart(3, ipart) = th_c
      zstart(4, ipart) = 1.d0
      zstart(5, ipart) = v_par / sqrt(v_par**2 + v_perp**2)
    enddo
  end subroutine

  subroutine init_starting_points
    use get_can_sub, only: can_to_vmec

    integer :: i, ipart, iskip
    double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec

    ! skip random numbers according to configuration
      do iskip=1,loopskip
        do ipart=1,ntestpart
          call random_number(xi)
          call random_number(xi)
        enddo
      enddo

    ! files for storing starting coords
    open(1,file='start.dat',recl=1024)
    ! determine the starting point:
    if (startmode == 0 .or. startmode == 1) then
      do ipart=1,ntestpart
        call random_number(xi)
        call binsrc(volstart,1,npoiper*nper,xi,i)
        ibins=i
        ! coordinates: z(1) = r, z(2) = vartheta, z(3) = varphi
        r=xstart(1,i)
        vartheta=xstart(2,i)
        varphi=xstart(3,i)
  !
  ! we store starting points in VMEC coordinates:
        if(isw_field_type.eq.0) then
          call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
        elseif(isw_field_type.eq.1) then
          theta_vmec=vartheta
          varphi_vmec=varphi
        elseif(isw_field_type.eq.2) then
          call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
        else
          print *,'init_starting_points: unknown field type'
        endif
  !
        zstart(1,ipart)=r
        zstart(2,ipart)=theta_vmec
        zstart(3,ipart)=varphi_vmec
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
    open(1,file='start.dat',recl=1024)
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
  !
  ! we store starting points in VMEC coordinates:
        if(isw_field_type.eq.0) then
          call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
        elseif(isw_field_type.eq.1) then
          theta_vmec=vartheta
          varphi_vmec=varphi
        elseif(isw_field_type.eq.2) then
          call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
        else
          print *,'init_starting_points: unknown field type'
        endif
  !
        zstart(1,ipart)=r
        zstart(2,ipart)=theta_vmec
        zstart(3,ipart)=varphi_vmec
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
    use find_bminmax_sub, only : get_bminmax
    use get_can_sub, only : vmec_to_can
    use magfie_sub, only : magfie
    use plag_coeff_sub, only : plag_coeff
    use alpha_lifetime_sub, only : orbit_timestep_axis
    use classification, only : trace_orbit_with_classifiers

    type(Tracer), intent(inout) :: anorb
    integer, intent(in) :: ipart
    integer :: ierr, ierr_coll
    double precision, dimension(5) :: z
    double precision :: bmod,sqrtg
    double precision, dimension(3) :: bder, hcovar, hctrvr, hcurl
    integer :: it, ktau
    integer(8) :: kt
    logical :: passing

    integer                                       :: ifp_tip,ifp_per
    integer,          dimension(:),   allocatable :: ipoi
    double precision, dimension(:),   allocatable :: xp
    double precision, dimension(:,:), allocatable :: coef,orb_sten
    double precision, dimension(:,:), allocatable :: zpoipl_tip,zpoipl_per,dummy2d
    double precision, dimension(n_tip_vars)       :: var_tip
    double precision :: phiper, alam_prev, par_inv
    integer :: iper, itip, kper, nfp_tip, nfp_per

    double precision :: fraction
    double precision :: r,theta_vmec,varphi_vmec
    logical :: regular

  ! Variables and settings for classification by J_parallel and ideal orbit condition:
    integer, parameter :: nfp_dim=3, nturns=8
    integer :: nfp_cot,ideal,ijpar,ierr_cot,iangvar
    double precision, dimension(nfp_dim) :: fpr_in

  ! for run with fixed random seed
    integer :: seedsize
    integer, allocatable :: seed(:)

    if (ntcut>0 .or. class_plot) then
      call trace_orbit_with_classifiers(anorb, ipart)
      return
    endif

    zend(:,ipart) = 0d0

    if (deterministic) then
      call random_seed(size = seedsize)
      if (.not. allocated(seed)) allocate(seed(seedsize))
      seed = 0
      call random_seed(put=seed)
    endif

    z = zstart(:, ipart)
    r=z(1)
    theta_vmec=z(2)
    varphi_vmec=z(3)

    if(isw_field_type.eq.0) then
        call vmec_to_can(r,theta_vmec,varphi_vmec,z(2),z(3))
    elseif(isw_field_type.eq.2) then
        call vmec_to_boozer(r,theta_vmec,varphi_vmec,z(2),z(3))
    endif

    if (integmode>0) call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)

    call magfie(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)

  !$omp critical
    if(num_surf > 1) then
      call get_bminmax(z(1),bmin,bmax)
    endif
    passing = z(5)**2.gt.1.d0-bmod/bmax
    trap_par(ipart) = ((1.d0-z(5)**2)*bmax/bmod-1.d0)*bmin/(bmax-bmin)
    perp_inv(ipart) = z(4)**2*(1.d0-z(5)**2)/bmod
  !$omp end critical

  ! Forced classification of passing as regular:
    if(passing .and. should_skip(ipart)) then
      !$omp critical
      confpart_pass=confpart_pass+1.d0
      !$omp end critical
      return
    endif

    kt = 0
    if (passing) then
      !$omp atomic
      confpart_pass(1)=confpart_pass(1)+1.d0
    else
      !$omp atomic
      confpart_trap(1)=confpart_trap(1)+1.d0
    end if
  !
    par_inv = 0d0
    regular = .False.
    do it=2,ntimstep
      if (regular) then  ! regular orbit, will not be lost
        if(passing) then
          !$omp atomic
          confpart_pass(it)=confpart_pass(it)+1.d0
        else
          !$omp atomic
          confpart_trap(it)=confpart_trap(it)+1.d0
        endif
        kt = kt+ntau
        cycle
      endif
      do ktau=1,ntau
        if (integmode <= 0) then
          call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
        else
          if (swcoll) then  ! TODO: move this inside modules
            anorb%si%pabs = z(4)
            anorb%f%vpar = z(4)*z(5)*sqrt2
            anorb%f%mu = z(4)**2*(1.d0-z(5)**2)/anorb%f%Bmod
            anorb%si%z(4) = anorb%f%vpar*anorb%f%hph + anorb%f%Aph/anorb%f%ro0
            call get_val(anorb%f, anorb%si%z(4)) ! for pth
            anorb%si%pthold = anorb%f%pth
          endif
          call orbit_timestep_sympl(anorb%si, anorb%f, ierr)
          z(1:3) = anorb%si%z(1:3)
          z(4) = dsqrt(anorb%f%mu*anorb%f%Bmod+0.5d0*anorb%f%vpar**2)
          z(5) = anorb%f%vpar/(z(4)*sqrt2)
        endif

        ! Collisions
        if (swcoll) then
          call stost(z, dtaumin, 1, ierr_coll)
          if (ierr_coll /= 0) then
            print *, 'Error in stost: ', ierr_coll, 'z = ', z, 'dtaumin = ', dtaumin
          endif
        endif

        if(ierr.ne.0) exit
        kt = kt+1
      enddo
      if(ierr.ne.0) exit
      if(passing) then
        !$omp atomic
        confpart_pass(it)=confpart_pass(it)+1.d0
      else
        !$omp atomic
        confpart_trap(it)=confpart_trap(it)+1.d0
      endif
    enddo

    !$omp critical
    zend(:,ipart) = z
    if(isw_field_type.eq.0) then
      ! TODO need to add can_to_vmec
      ! call can_to_vmec(z(1),z(2),z(3),zend(2,ipart),zend(3,ipart))
    elseif(isw_field_type.eq.2) then
      call boozer_to_vmec(z(1),z(2),z(3),zend(2,ipart),zend(3,ipart))
    endif
    times_lost(ipart) = kt*dtaumin/v0
    !$omp end critical
  end subroutine trace_orbit

  end module simple_main
