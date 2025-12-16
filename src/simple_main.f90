module simple_main
  use omp_lib
  use util, only : sqrt2
  use simple, only : init_vmec, init_sympl, tracer_t
  use diag_mod, only : icounter
  use collis_alp, only : loacol_alpha, stost
  use samplers, only: sample
  use field_can_mod, only : integ_to_ref, ref_to_integ, init_field_can
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
      ntestpart, ntimstep, coord_input
    use timing, only : init_timer, print_phase_time
    use magfie_sub, only : TEST, VMEC, init_magfie
    use samplers, only : init_starting_surf

    implicit none

    character(256) :: config_file
    type(tracer_t) :: norb

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

    if (isw_field_type == TEST) then
      ! TEST field uses analytic tokamak - no VMEC needed for sampling
      call init_magfie(TEST)
      call print_phase_time('TEST field initialization completed')
      call sample_particles_test_field
      call print_phase_time('TEST field particle sampling completed')
    else
      call init_magfie(VMEC)
      call print_phase_time('VMEC magnetic field initialization completed')

      call init_starting_surf
      call print_phase_time('Starting surface initialization completed')

      call sample_particles
      call print_phase_time('Particle sampling completed')

      if (generate_start_only) stop 'stopping after generating start.dat'

      call init_magfie(isw_field_type)
      call print_phase_time('Field type initialization completed')
    end if

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
    use field_base, only : magnetic_field_t
    use field, only : field_from_file
    use timing, only : print_phase_time
    use magfie_sub, only : TEST, CANFLUX, VMEC, BOOZER, MEISS, ALBERT, &
                           REFCOORDS, set_magfie_refcoords_field
    use field_splined, only : splined_field_t
    use util, only : twopi
    use reference_coordinates, only : init_reference_coordinates
    use params, only : coord_input, field_input

    character(*), intent(in) :: vmec_file
    type(tracer_t), intent(inout) :: self
    integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
    class(magnetic_field_t), allocatable :: field_temp

    self%integmode = aintegmode

    ! TEST field is analytic - no VMEC or field files needed
    if (isw_field_type == TEST) then
      self%fper = twopi  ! Full torus for analytic tokamak
      call print_phase_time('TEST field mode - no input files required')
    else
      call init_vmec(vmec_file, ans_s, ans_tp, amultharm, self%fper)
      call print_phase_time('VMEC initialization completed')

      call init_reference_coordinates(coord_input)
      call print_phase_time('Reference coordinate system initialization completed')

      if (self%integmode >= 0) then
        if (trim(field_input) == '') then
          print *, 'simple_main.init_field: field_input must be set (see params.apply_config_aliases)'
          error stop
        end if

        ! For VMEC fields we require that coord_input and field_input agree
        ! so that we can safely use the VMEC pass-through path.
        ! Pass-through optimization: When both coord_input and field_input point to
        ! the same VMEC file, vmec_field_t directly uses the libneo VMEC splines loaded
        ! by init_vmec rather than re-splining the field data. This preserves both
        ! performance and backward compatibility.
        block
          character(:), allocatable :: coord_file, field_file
          logical :: is_vmec_field

          coord_file = trim(coord_input)
          field_file = trim(field_input)

          is_vmec_field = .false.
          if (len_trim(field_file) >= 3) then
            if (field_file(len_trim(field_file)-2:len_trim(field_file)) == '.nc') then
              is_vmec_field = .true.
            end if
          end if

          if (is_vmec_field) then
            if (len_trim(coord_file) == 0) then
              print *, 'simple_main.init_field: coord_input must be set when using VMEC field_input'
              error stop
            end if
            if (coord_file /= field_file) then
              print *, 'simple_main.init_field: VMEC coord_input and field_input must match for VMEC pass-through'
              print *, '  coord_input = ', coord_file
              print *, '  field_input = ', field_file
              error stop
            end if
          end if
        end block

        call field_from_file(field_input, field_temp)
        call print_phase_time('Field from file loading completed')

        if (isw_field_type == REFCOORDS) then
          select type (field_temp)
          type is (splined_field_t)
            call set_magfie_refcoords_field(field_temp)
          class default
            print *, 'simple_main.init_field: REFCOORDS requires a splined field'
            print *, 'Use a coils or other cartesian field input that is splined'
            error stop
          end select
        end if
      end if
    end if

    if (self%integmode > 0) then
      select case (isw_field_type)
      case (VMEC)
        error stop 'Symplectic guiding-center integrators require canonical field ' &
          //'representation (set isw_field_type to TEST, CANFLUX, BOOZER, MEISS, or ALBERT)'
      case (TEST, CANFLUX, BOOZER, MEISS, ALBERT)
        continue
      case default
        error stop 'Unknown canonical field type for symplectic guiding-center integrator'
      end select
    end if

    if (isw_field_type == TEST) then
      ! TEST field is fully analytic - no field file needed
      call init_field_can(isw_field_type)
      call print_phase_time('Canonical field initialization completed')
    else if (isw_field_type == CANFLUX .or. isw_field_type == BOOZER .or. &
             isw_field_type == MEISS .or. isw_field_type == ALBERT) then
      call init_field_can(isw_field_type, field_temp)
      call print_phase_time('Canonical field initialization completed')
    end if
  end subroutine init_field


  subroutine trace_parallel(norb)
    use netcdf_orbit_output, only : init_orbit_netcdf, close_orbit_netcdf, &
                                     write_orbit_to_netcdf

    type(tracer_t), intent(inout) :: norb
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
      print *, kpart, ' / ', ntestpart, 'particle: ', i, 'thread: ', &
        omp_get_thread_num()
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

  subroutine classify_parallel(norb)
    use classification, only: trace_orbit_with_classifiers, classification_result_t
    use params, only: class_passing, class_lost

    type(tracer_t), intent(inout) :: norb
    integer :: i
    type(classification_result_t) :: class_result

    !$omp parallel firstprivate(norb) private(class_result, i)
    !$omp do
    do i = 1, ntestpart
      !$omp critical
      kpart = kpart+1
      print *, kpart, ' / ', ntestpart, 'particle: ', i, 'thread: ', &
        omp_get_thread_num()
      !$omp end critical

      call reset_seed_if_deterministic
      call trace_orbit_with_classifiers(norb, i, class_result)

      ! Store classification flags in params arrays
      class_passing(i) = class_result%passing
      class_lost(i) = class_result%lost
      ! iclass already populated by trace_orbit_with_classifiers
      ! Other results (zend, times_lost, trap_par, perp_inv) also already stored
    end do
    !$omp end do
    !$omp end parallel
  end subroutine classify_parallel

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
        print *, 'Invalid surface range for volume sample defined (2 < num_surf),', &
          ' stopping.'
        stop
      endif

    else
      print *, 'Unknown startmode: ', startmode
    endif
  end subroutine sample_particles

  subroutine sample_particles_test_field
    !> Sample particles for the analytic circular tokamak TEST field.
    !> TEST field uses (r, theta, phi) coordinates with B0=1, R0=1, a=0.5, iota=1.
    !> bmod = B0 * (1 - r/R0 * cos(theta))
    use util, only : twopi
    use params, only : ntestpart, sbeg, bmod00, bmin, bmax, zstart, &
      reset_seed_if_deterministic

    real(dp), parameter :: B0 = 1.0d0, R0 = 1.0d0, a = 0.5d0
    real(dp) :: r_start, tmp_rand
    integer :: ipart

    ! Use sbeg(1) as the starting minor radius (mapped to r for tokamak)
    ! sbeg(1) is input as a flux-like value; for TEST field, interpret as r/a
    r_start = sbeg(1) * a

    ! Set magnetic field bounds for this r value
    ! bmod = B0*(1 - r/R0*cos(theta))
    ! Maximum at theta=pi: bmax = B0*(1 + r/R0)
    ! Minimum at theta=0:  bmin = B0*(1 - r/R0)
    bmod00 = B0
    bmax = B0 * (1.0d0 + r_start / R0)
    bmin = B0 * (1.0d0 - r_start / R0)

    print *, 'TEST field: bmod00 = ', bmod00, 'bmin = ', bmin, 'bmax = ', bmax

    call reset_seed_if_deterministic

    ! Sample particles uniformly in (theta, phi) at fixed r
    do ipart = 1, ntestpart
      zstart(1, ipart) = r_start
      call random_number(tmp_rand)
      zstart(2, ipart) = twopi * tmp_rand
      call random_number(tmp_rand)
      zstart(3, ipart) = twopi * tmp_rand
      zstart(4, ipart) = 1.0d0
      call random_number(tmp_rand)
      zstart(5, ipart) = 2.0d0 * (tmp_rand - 0.5d0)
    end do

    call save_starting_points_test(zstart)

  contains
    subroutine save_starting_points_test(zs)
      real(dp), intent(in) :: zs(:,:)
      integer :: i, unit

      open(newunit=unit, file='start.dat', recl=1024)
      do i = 1, size(zs, 2)
        write(unit, *) zs(:, i)
      end do
      close(unit)
    end subroutine save_starting_points_test
  end subroutine sample_particles_test_field

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
    use classification, only : trace_orbit_with_classifiers, classification_result_t, &
      write_classification_results
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    type(tracer_t), intent(inout) :: anorb
    integer, intent(in) :: ipart
    real(dp), intent(out) :: orbit_traj(:,:)  ! (5, ntimstep)
    real(dp), intent(out) :: orbit_times(:)   ! (ntimstep)

    real(dp), dimension(5) :: z
    integer :: it, ierr_orbit, it_final
    integer(8) :: kt
    logical :: passing
    type(classification_result_t) :: class_result

    ierr_orbit = 0

    call reset_seed_if_deterministic

    if (ntcut>0 .or. class_plot) then
      call trace_orbit_with_classifiers(anorb, ipart, class_result)
      if(class_plot) then
        call write_classification_results(ipart, class_result)
      endif
      return
    endif

    call ref_to_integ(zstart(1:3, ipart), z(1:3))
    z(4:5) = zstart(4:5, ipart)
    zend(:,ipart) = 0d0

    if (integmode > 0) then
      call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)
    end if

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
    call integ_to_ref(z(1:3), zend(1:3,ipart))
    zend(4:5, ipart) = z(4:5)
    times_lost(ipart) = kt*dtaumin/v0
    !$omp end critical
  end subroutine trace_orbit

  subroutine macrostep(anorb, z, kt, ierr_orbit)
    use alpha_lifetime_sub, only : orbit_timestep_axis
    use orbit_symplectic, only : orbit_timestep_sympl

    type(tracer_t), intent(inout) :: anorb
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
    type(tracer_t), intent(in) :: anorb
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

    type(tracer_t), intent(inout) :: anorb
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
