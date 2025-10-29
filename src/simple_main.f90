module simple_main
  use omp_lib
  use util, only: sqrt2
  use simple, only : init_vmec, init_sympl, Tracer
  use diag_mod, only : icounter
  use collis_alp, only : loacol_alpha, stost
  use samplers, only: sample
  use field_can_mod, only : can_to_ref, ref_to_can, init_field_can
  use callback, only : output_orbits_macrostep
  use params, only: swcoll, ntestpart, startmode, special_ants_file, num_surf, &
    grid_density, dtau, dtaumin, ntau, v0, &
    kpart, confpart_pass, confpart_trap, times_lost, integmode, relerr, trace_time, &
    class_plot, ntcut, iclass, bmin, bmax, &
    zstart, zend, trap_par, perp_inv, sbeg, &
    ntimstep, should_skip, reset_seed_if_deterministic, &
    field_input, isw_field_type, reuse_batch

  implicit none

  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)

  contains

  subroutine main
    use params, only : read_config, netcdffile, ns_s, ns_tp, multharm, &
      integmode, params_init, swcoll, generate_start_only, isw_field_type, &
      ntestpart, ntimstep
    use timing, only : init_timer, print_phase_time
    use magfie_sub, only : VMEC, init_magfie
    use samplers, only : init_starting_surf

    implicit none

    character(256) :: config_file
    type(Tracer) :: norb

    ! Initialize timing
    call init_timer()

    ! read configuration file name from command line arguments
    if (command_argument_count() == 0) then
      config_file = 'simple.in'
    else
      call get_command_argument(1, config_file)
    end if
    call print_phase_time('Command line parsing completed')

    ! Must be called in this order. TODO: Fix
    call read_config(config_file)
    call print_phase_time('Configuration reading completed')
    
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call print_phase_time('Field initialization completed')
    
    call params_init
    call print_phase_time('Parameter initialization completed')

    call print_parameters
    call print_phase_time('Parameter printing completed')
    
    if (swcoll) then
      call init_collisions
      call print_phase_time('Collision initialization completed')
    endif

    call init_magfie(VMEC)
    call print_phase_time('VMEC magnetic field initialization completed')
    
    call init_starting_surf
    call print_phase_time('Starting surface initialization completed')

    call sample_particles
    call print_phase_time('Particle sampling completed')

    if (generate_start_only) stop 'stopping after generating start.dat'

    call init_magfie(isw_field_type)
    call print_phase_time('Field type initialization completed')

    call init_counters
    call print_phase_time('Counter initialization completed')

    call trace_parallel(norb)
    call print_phase_time('Parallel particle tracing completed')

    confpart_pass=confpart_pass/ntestpart
    confpart_trap=confpart_trap/ntestpart
    call print_phase_time('Statistics normalization completed')
    
    call write_output
    call print_phase_time('Output writing completed')
  end subroutine main

  subroutine init_field(self, vmec_file, ans_s, ans_tp, amultharm, aintegmode)
    use field_base, only : MagneticField
    use field, only : field_from_file
    use timing, only : print_phase_time

    character(*), intent(in) :: vmec_file
    type(Tracer), intent(inout) :: self
    integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
    class(MagneticField), allocatable :: field_temp

    call init_vmec(vmec_file, ans_s, ans_tp, amultharm, self%fper)
    call print_phase_time('VMEC initialization completed')

    self%integmode = aintegmode
    if (self%integmode >= 0) then
      if(trim(field_input) == '') then
        call field_from_file(vmec_file, field_temp)
      else
        call field_from_file(field_input, field_temp)
      end if
      call print_phase_time('Field from file loading completed')
    end if

    if (isw_field_type == 0 .or. isw_field_type >= 2) then
      call init_field_can(isw_field_type, field_temp)
      call print_phase_time('Canonical field initialization completed')
    end if
  end subroutine init_field


  subroutine trace_parallel(norb)
    use netcdf_orbit_output, only : init_orbit_netcdf, close_orbit_netcdf, &
                                     write_orbit_to_netcdf

    type(Tracer), intent(inout) :: norb
    integer :: i
    real(dp), allocatable :: traj(:,:), times(:)

    if (output_orbits_macrostep) then
      call init_orbit_netcdf(ntestpart, ntimstep)
    endif

    !$omp parallel firstprivate(norb) private(traj, times, i)
    allocate(traj(5, ntimstep), times(ntimstep))

    !$omp do
    do i = 1, ntestpart
      !$omp critical
      kpart = kpart+1
      print *, kpart, ' / ', ntestpart, 'particle: ', i, 'thread: ', omp_get_thread_num()
      !$omp end critical

      call trace_orbit(norb, i, traj, times)

      if (output_orbits_macrostep) then
        !$omp critical
        call write_orbit_to_netcdf(i, traj, times)
        !$omp end critical
      endif
    end do
    !$omp end do
    !$omp end parallel

    if (output_orbits_macrostep) then
      call close_orbit_netcdf()
    endif
  end subroutine trace_parallel

  subroutine print_parameters
    print *, 'tau: ', dtau, dtaumin, min(dabs(mod(dtau, dtaumin)), &
                      dabs(mod(dtau, dtaumin)-dtaumin))/dtaumin, ntau
    print *, 'v0 = ', v0
  end subroutine print_parameters

  subroutine init_collisions
    use params, only: am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
    facE_al, dchichi, slowrate, dchichi_norm, slowrate_norm, v0

    real(dp) :: v0_coll

    call loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe, &
      3.5d6/facE_al,v0_coll,dchichi,slowrate,dchichi_norm,slowrate_norm)

    if (abs(v0_coll - v0) > 1d-6) then
      error stop 'simple_main.init_collisions: v0_coll != v0'
    end if
  end subroutine init_collisions

  subroutine sample_particles
    use samplers, only: sample, START_FILE

    if (1 == startmode) then
      if ((0d0 < grid_density) .and. (1d0 > grid_density)) then
        call sample(zstart, grid_density)
      else
        call sample(zstart)
      endif

    elseif (2 == startmode) then
      call sample(zstart, START_FILE)

    elseif (3 == startmode) then
      call sample(special_ants_file)

    elseif (4 == startmode) then
      call sample(zstart, reuse_batch)

    elseif (5 == startmode) then
      if (0 == num_surf) then
        call sample(zstart, 0.0d0, 1.0d0)
      elseif (1 == num_surf) then
        call sample(zstart, 0.0d0, sbeg(1))
      elseif (2 == num_surf) then
        call sample(zstart, sbeg(1), sbeg(num_surf))
      else
        print *, 'Invalid surface range for volume sample defined (2 < num_surf), stopping.'
        stop
      endif

    else
      print *, 'Unknown startmode: ', startmode
    endif
  end subroutine sample_particles

  subroutine init_counters
    icounter=0 ! evaluation counter
    kpart=0

    ! initialize array of confined particle percentage
    confpart_trap=0.d0
    confpart_pass=0.d0

    ! initialize lost times when particles get lost
    times_lost = -1.d0
  end subroutine init_counters

  subroutine trace_orbit(anorb, ipart, orbit_traj, orbit_times)
    use classification, only : trace_orbit_with_classifiers
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    type(Tracer), intent(inout) :: anorb
    integer, intent(in) :: ipart
    real(dp), intent(out) :: orbit_traj(:,:)  ! (5, ntimstep)
    real(dp), intent(out) :: orbit_times(:)   ! (ntimstep)

    real(dp), dimension(5) :: z
    integer :: it, ierr_orbit, it_final
    integer(8) :: kt
    logical :: passing

    ierr_orbit = 0

    call reset_seed_if_deterministic

    if (ntcut>0 .or. class_plot) then
      call trace_orbit_with_classifiers(anorb, ipart)
      return
    endif

    call ref_to_can(zstart(1:3, ipart), z(1:3))
    z(4:5) = zstart(4:5, ipart)
    zend(:,ipart) = 0d0

    if (integmode>0) call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)

    call compute_pitch_angle_params(z, passing, trap_par(ipart), perp_inv(ipart))

    if(passing .and. should_skip(ipart)) then
      ! Fill trajectory arrays with NaN since we're not tracing this particle
      orbit_traj = ieee_value(0.0d0, ieee_quiet_nan)
      orbit_times = ieee_value(0.0d0, ieee_quiet_nan)
      !$omp critical
      confpart_pass=confpart_pass+1.d0
      !$omp end critical
      return
    endif

    kt = 0
    it_final = 0
    do it = 1, ntimstep
      if (it >= 2) call macrostep(anorb, z, kt, ierr_orbit)
      if(ierr_orbit .ne. 0) then
        it_final = it
        exit
      endif

      ! Store trajectory data (after macrostep so time is correct)
      orbit_traj(:, it) = z
      orbit_times(it) = kt*dtaumin/v0

      call increase_confined_count(it, passing)
      it_final = it
    enddo

    ! Fill remaining timesteps with NaN if particle left domain early
    if (it_final < ntimstep) then
      do it = it_final + 1, ntimstep
        orbit_traj(:, it) = ieee_value(0.0d0, ieee_quiet_nan)
        orbit_times(it) = ieee_value(0.0d0, ieee_quiet_nan)
      enddo
    endif

    !$omp critical
    call can_to_ref(z(1:3), zend(1:3,ipart))
    zend(4:5, ipart) = z(4:5)
    times_lost(ipart) = kt*dtaumin/v0
    !$omp end critical
  end subroutine trace_orbit

  subroutine macrostep(anorb, z, kt, ierr_orbit)
    use alpha_lifetime_sub, only : orbit_timestep_axis
    use orbit_symplectic, only : orbit_timestep_sympl

    type(Tracer), intent(inout) :: anorb
    real(dp), intent(inout) :: z(5)
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
    real(dp), intent(inout) :: z(5)

    z(1:3) = anorb%si%z(1:3)
    z(4) = dsqrt(anorb%f%mu*anorb%f%Bmod+0.5d0*anorb%f%vpar**2)
    z(5) = anorb%f%vpar/(z(4)*sqrt2)
  end subroutine to_standard_z_coordinates

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

    real(dp), intent(in) :: z(5)
    logical, intent(out) :: passing
    real(dp), intent(out) :: trap_par_, perp_inv_

    real(dp) :: bmod

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

    real(dp) :: bmod
    real(dp), intent(in) :: z(3)

    real(dp) :: sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

    call magfie(z(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  end function compute_bmod

  subroutine update_momentum(anorb, z)
    use orbit_symplectic, only : get_val

    type(Tracer), intent(inout) :: anorb
    real(dp), intent(in) :: z(5)

    anorb%si%pabs = z(4)
    anorb%f%vpar = z(4)*z(5)*sqrt2
    anorb%f%mu = z(4)**2*(1.d0-z(5)**2)/anorb%f%Bmod
    anorb%si%z(4) = anorb%f%vpar*anorb%f%hph + anorb%f%Aph/anorb%f%ro0
    call get_val(anorb%f, anorb%si%z(4)) ! for pth
    anorb%si%pthold = anorb%f%pth
  end subroutine update_momentum

  subroutine collide(z, dt)
    real(dp), intent(in) :: z(5), dt
    integer :: ierr_coll

    call stost(z, dt, 1, ierr_coll)
    if (ierr_coll /= 0) then
      print *, 'Error in stost: ', ierr_coll, 'z = ', z, 'dtaumin = ', dtaumin
    endif
  end subroutine collide


  subroutine write_output
    use field_can_base, only: n_field_evaluations

    integer :: i, num_lost
    real(dp) :: inverse_times_lost_sum
    integer(8) :: total_field_evaluations

    ! Sum field evaluations across all threads
    total_field_evaluations = 0
    !$omp parallel reduction(+:total_field_evaluations)
    total_field_evaluations = total_field_evaluations + n_field_evaluations
    !$omp end parallel
    
    print *, "Total field evaluations: ", total_field_evaluations

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
