module samplers
  use util

  implicit none

  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)

  character(len=*), parameter :: START_FILE = 'start.dat'
  character(len=*), parameter :: START_FILE_ANTS = 'start_ants.dat'
  character(len=*), parameter :: START_FILE_BATCH = 'batch.dat'

  ! Interface ################################
  INTERFACE sample
       MODULE PROCEDURE sample_read
       MODULE PROCEDURE sample_surface_fieldline
       MODULE PROCEDURE sample_grid
       MODULE PROCEDURE sample_volume_single
       MODULE PROCEDURE sample_random_batch
       MODULE PROCEDURE sample_points_ants
  END INTERFACE sample


  contains
  ! Functions #################################

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Generate random value in range
  function random_in_range(min_val, max_val) result(value)
    real(dp), intent(in) :: min_val, max_val
    real(dp) :: value, rand_val
    
    call random_number(rand_val)
    value = rand_val * (max_val - min_val) + min_val
  end function random_in_range

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Generate random pitch angle
  function random_pitch() result(pitch)
    real(dp) :: pitch, rand_val
    
    call random_number(rand_val)
    pitch = 2.0_dp * (rand_val - 0.5_dp)
  end function random_pitch

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Generate random angles (theta and phi)
  subroutine random_angles(theta, phi)
    real(dp), intent(out) :: theta, phi
    real(dp) :: rand_val
    
    call random_number(rand_val)
    theta = twopi * rand_val
    
    call random_number(rand_val)
    phi = twopi * rand_val
  end subroutine random_angles

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Initialize particle velocity
  subroutine initialize_particle_velocity(z, normalized_v, pitch)
    real(dp), dimension(:), intent(inout) :: z
    real(dp), intent(in), optional :: normalized_v, pitch
    
    ! Set normalized velocity (default to 1.0)
    if (present(normalized_v)) then
      z(4) = normalized_v
    else
      z(4) = 1.0_dp
    end if
    
    ! Set pitch angle (random if not provided)
    if (present(pitch)) then
      z(5) = pitch
    else
      z(5) = random_pitch()
    end if
  end subroutine initialize_particle_velocity

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Write particles to file
  subroutine write_particles_to_file(filename, particles)
    character(len=*), intent(in) :: filename
    real(dp), dimension(:,:), intent(in) :: particles
    integer :: ipart, unit_num, iostat
    
    ! Use newunit for automatic unit allocation
    open(newunit=unit_num, file=filename, recl=1024, iostat=iostat)
    if (iostat /= 0) then
      print *, 'Error opening file for writing: ', trim(filename)
      error stop 'File I/O error in write_particles_to_file'
    end if
    
    do ipart = 1, size(particles, 2)
      write(unit_num, *, iostat=iostat) particles(:, ipart)
      if (iostat /= 0) then
        print *, 'Error writing particle ', ipart
        close(unit_num)
        error stop 'Write error in write_particles_to_file'
      end if
    end do
    close(unit_num)
  end subroutine write_particles_to_file

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Read particles from file
  subroutine read_particles_from_file(filename, particles)
    character(len=*), intent(in) :: filename
    real(dp), dimension(:,:), intent(out) :: particles
    integer :: ipart, unit_num, iostat
    
    ! Use newunit for automatic unit allocation
    open(newunit=unit_num, file=filename, recl=1024, iostat=iostat)
    if (iostat /= 0) then
      print *, 'Error opening file for reading: ', trim(filename)
      error stop 'File I/O error in read_particles_from_file'
    end if
    
    do ipart = 1, size(particles, 2)
      read(unit_num, *, iostat=iostat) particles(:, ipart)
      if (iostat /= 0) then
        print *, 'Error reading particle ', ipart
        close(unit_num)
        error stop 'Read error in read_particles_from_file'
      end if
    end do
    close(unit_num)
  end subroutine read_particles_from_file

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute B field extrema
  subroutine compute_b_extrema(bstart, bmin, bmax)
    real(dp), dimension(:), intent(in) :: bstart
    real(dp), intent(out) :: bmin, bmax
    
    bmax = maxval(bstart)
    bmin = minval(bstart)
  end subroutine compute_b_extrema

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Calculate grid parameters
  subroutine calculate_grid_params(grid_density, xsize, ngrid, ntestpart)
    real(dp), intent(in) :: grid_density
    real(dp), intent(out) :: xsize
    integer, intent(out) :: ngrid, ntestpart
    
    xsize = 2.0_dp * pi * grid_density  ! angle density
    ngrid = int(1.0_dp / grid_density - 1.0_dp)
    ntestpart = ngrid * ngrid  ! total angle points
  end subroutine calculate_grid_params
  subroutine init_starting_surf
    use alpha_lifetime_sub, only : integrate_mfl_can
    use params, only: dphi, nper, npoiper, phibeg, thetabeg, volstart, &
        xstart, sbeg, bmin, bmax, bmod00

    integer :: ierr = 0
    real(dp), dimension(npoiper*nper) :: bstart

    ! Initialize arrays
    xstart = 0.0_dp
    bstart = 0.0_dp
    volstart = 0.0_dp

    ! Integrate along magnetic field line
    call integrate_mfl_can( &
      npoiper*nper, dphi, sbeg(1), phibeg, thetabeg, &
      xstart, bstart, volstart, bmod00, ierr)

    if (ierr /= 0) then
      error stop 'Starting field line has points outside the chamber'
    end if

    ! Compute B field extrema
    call compute_b_extrema(bstart, bmin, bmax)

    print *, 'bmod00 = ', bmod00, 'bmin = ', bmin, 'bmax = ', bmax
  end subroutine init_starting_surf

  subroutine load_starting_points(zstart, filename)
    real(dp), dimension(:,:), intent(inout) :: zstart
    character(len=*), intent(in) :: filename
    
    call read_particles_from_file(filename, zstart)
  end subroutine load_starting_points

  subroutine save_starting_points(zstart)
    real(dp), dimension(:,:), intent(in) :: zstart
    
    call write_particles_to_file(START_FILE, zstart)
  end subroutine save_starting_points

  subroutine sample_read(zstart, filename)
      real(dp), dimension(:,:), intent(inout) :: zstart
      character(len=*), intent(in) :: filename

      call load_starting_points(zstart, filename)
  end subroutine


  ! Samplers ################################
  subroutine sample_volume_single(zstart, s_inner, s_outer)
    use params, only: isw_field_type, num_surf
    use field_can_mod, only : can_to_ref

    real(dp), intent(in) :: s_inner, s_outer
    real(dp), dimension(:,:), intent(inout) :: zstart
    real(dp) :: r, vartheta, varphi
    integer :: ipart

    ! Ensure we have 2 surfaces for volume sampling
    if (num_surf /= 2) then
      num_surf = 2
    end if

    do ipart = 1, size(zstart, 2)
      ! Generate random radius within range
      r = random_in_range(s_inner, s_outer)
      
      ! Generate random angles
      call random_angles(vartheta, varphi)
      
      ! Convert to reference coordinates
      call can_to_ref([r, vartheta, varphi], zstart(1:3, ipart))
      
      ! Initialize particle velocity
      call initialize_particle_velocity(zstart(:, ipart))
    end do

    call save_starting_points(zstart)

  end subroutine sample_volume_single

  subroutine sample_surface_fieldline(zstart)
    use params, only: volstart, isw_field_type, ibins, xstart, npoiper, nper
    use binsrc_sub, only: binsrc

    real(dp), dimension(:,:), intent(inout) :: zstart
    real(dp) :: xi
    integer :: ipart, i

    do ipart = 1, size(zstart, 2)
      ! Select random point on field line surface
      call random_number(xi)
      call binsrc(volstart, 1, npoiper*nper, xi, i)
      ibins = i
      
      ! Copy coordinates from starting surface
      zstart(1:3, ipart) = xstart(:, i)
      
      ! Initialize particle velocity
      call initialize_particle_velocity(zstart(:, ipart))
    end do

    call save_starting_points(zstart)

  end subroutine sample_surface_fieldline

  subroutine sample_grid(zstart, grid_density)
    use params, only: ntestpart, zstart_dim1, zend, times_lost, &
        trap_par, perp_inv, iclass, xstart, sbeg
    use util, only: pi

    real(dp), dimension(:,:), allocatable, intent(inout) :: zstart
    real(dp), intent(in) :: grid_density
    real(dp) :: xsize
    integer :: ngrid, ipart, jpart, lidx

    ! Calculate grid parameters
    call calculate_grid_params(grid_density, xsize, ngrid, ntestpart)

    ! Reallocate arrays for new particle count
    call reallocate_grid_arrays(zstart, ntestpart)

    ! Generate grid of particles
    do ipart = 1, ngrid
      do jpart = 1, ngrid
        lidx = (jpart - 1) * ngrid + ipart
        zstart(1, lidx) = sbeg(1)
        zstart(2, lidx) = xsize * ipart
        zstart(3, lidx) = xsize * jpart
        call initialize_particle_velocity(zstart(:, lidx))
      end do
    end do

    call save_starting_points(zstart)

  end subroutine sample_grid

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Reallocate arrays for grid sampling
  subroutine reallocate_grid_arrays(zstart, ntestpart_new)
    use params, only: zstart_dim1, zend, times_lost, trap_par, perp_inv, iclass
    
    real(dp), dimension(:,:), allocatable, intent(inout) :: zstart
    integer, intent(in) :: ntestpart_new
    
    ! Deallocate existing arrays
    if (allocated(zstart)) deallocate(zstart)
    if (allocated(zend)) deallocate(zend)
    if (allocated(times_lost)) deallocate(times_lost)
    if (allocated(trap_par)) deallocate(trap_par)
    if (allocated(perp_inv)) deallocate(perp_inv)
    if (allocated(iclass)) deallocate(iclass)
    
    ! Allocate with new size
    allocate(zstart(zstart_dim1, ntestpart_new))
    allocate(zend(zstart_dim1, ntestpart_new))
    allocate(times_lost(ntestpart_new))
    allocate(trap_par(ntestpart_new))
    allocate(perp_inv(ntestpart_new))
    allocate(iclass(3, ntestpart_new))
  end subroutine reallocate_grid_arrays

  subroutine sample_random_batch(zstart, reuse_existing)
  ! Get random batch from preexisting zstart, allows reuse.
    use params, only: batch_size, ntestpart, zstart_dim1, idx

    integer :: ran_begin, ran_end, ipart
    real :: temp_ran
    real(dp), dimension(:,:), intent(inout) :: zstart
    real(dp), dimension(zstart_dim1,batch_size) :: zstart_batch
    logical, intent(in) :: reuse_existing

    if (reuse_existing .eqv. .True.) then
      call load_starting_points(zstart_batch, START_FILE_BATCH)
    else
      call load_starting_points(zstart_batch, START_FILE)
      call random_number(temp_ran)
      ran_begin = INT(temp_ran)
      ran_end = ran_begin+batch_size
      if ((ran_end).gt.(ntestpart)) then
        ran_begin = ran_begin - (ran_end-ntestpart)
      endif
      do ipart=0,batch_size
        zstart(:,ipart) = zstart_batch(:,(ipart+ran_begin))
      enddo
    endif

    do ipart=idx(0),idx(ntestpart)
      read(1,*) zstart(:,ipart)
    enddo

  end subroutine sample_random_batch

  subroutine sample_points_ants(use_special_ants_file)
    use parse_ants, only : process_line
    use get_can_sub, only : vmec_to_can
    use params, only: ntestpart, zstart ! ANTS sampler uses global zstart

    logical, intent(in) :: use_special_ants_file

    integer, parameter :: maxlen = 4096
    character(len=maxlen) :: line
    real(8) :: v_par, v_perp, u, v, s
    real(8) :: th, ph
    integer :: ipart

    do ipart=1,ntestpart
      if (use_special_ants_file) then
        open (1, file=START_FILE_ANTS, recl=1024)
        read(1, '(A)') line
        close(1)
      else
        open(1, file=START_FILE, recl=1024)
        read(1, '(A)') line
        close(1)
      endif

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
  end subroutine sample_points_ants


end module samplers
