module samplers
  use util

  implicit none

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
  subroutine init_starting_surf
    use alpha_lifetime_sub, only : integrate_mfl_can
    use params, only: dphi, nper, npoiper, phibeg, thetabeg, volstart, &
        xstart, sbeg, bmin, bmax, bmod00

    integer :: ierr=0
    double precision, dimension(npoiper*nper) :: bstart


    xstart=0.d0
    bstart=0.d0
    volstart=0.d0

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

  subroutine load_starting_points(zstart, filename)
    double precision, dimension(:,:), intent(inout) :: zstart
    character(len=*), intent(in) :: filename
    integer :: ipart

    open(1,file=filename,recl=1024)
    do ipart=1,size(zstart,2)
      read(1,*) zstart(:,ipart)
    enddo
    close(1)
  end subroutine load_starting_points

  subroutine save_starting_points(zstart)
    double precision, dimension(:,:), intent(in) :: zstart
    integer :: ipart

    open(1,file=START_FILE,recl=1024)
    do ipart=1,size(zstart,2)
      write(1,*) zstart(:,ipart)
    enddo
    close(1)
  end subroutine save_starting_points

  subroutine sample_read(zstart, filename)
      double precision, dimension(:,:), intent(inout) :: zstart
      character(len=*), intent(in) :: filename

      call load_starting_points(zstart, filename)
  end subroutine
  
  
  ! Samplers ################################
  subroutine sample_volume_single(zstart, s_inner, s_outer)
    use params, only: isw_field_type, num_surf
    use field_can_mod, only : can_to_ref

    double precision, intent(in) :: s_inner
    double precision, intent(in) :: s_outer
    double precision :: tmp_rand
    double precision :: r,vartheta,varphi
    double precision, dimension(:,:), intent(inout) :: zstart
    integer :: ipart

    ! If user wants to do volume with 0 or 1 surfaces,
    !   we "add" the constraints, therefore having 2 surfaces.
    if (2 /= num_surf) then
      num_surf = 2
    endif

    do ipart=1,size(zstart,2)
      call random_number(tmp_rand)
      r = tmp_rand * (s_outer - s_inner) + s_inner

      call random_number(tmp_rand)
      vartheta=twopi*tmp_rand
      call random_number(tmp_rand)
      varphi=twopi*tmp_rand
      ! we store starting points in reference coordinates:
      call can_to_ref([r, vartheta, varphi], zstart(1:3,ipart))
      ! normalized velocity module z(4) = v / v_0:
      zstart(4,ipart)=1.d0
      ! starting pitch z(5)=v_\parallel / v:
      call random_number(tmp_rand)
      zstart(5,ipart)=2.d0*(tmp_rand-0.5d0)
    enddo

    call save_starting_points(zstart)

  end subroutine sample_volume_single

  subroutine sample_surface_fieldline(zstart)
    use params, only: volstart, isw_field_type, ibins, xstart, npoiper, nper
    use binsrc_sub, only: binsrc

    double precision, dimension(:,:), intent(inout) :: zstart

    double precision :: r,vartheta,varphi
    double precision :: xi
    integer :: ipart, i

    call init_starting_surf
    do ipart=1,size(zstart,2)
      call random_number(xi)
      call binsrc(volstart,1,npoiper*nper,xi,i)
      ibins=i
      ! coordinates: z(1) = r, z(2) = vartheta, z(3) = varphi
      r=xstart(1,i)
      vartheta=xstart(2,i)
      varphi=xstart(3,i)

      zstart(1:3,ipart)=xstart(:,i)
      zstart(4,ipart)=1.d0  ! normalized velocity module z(4) = v / v_0
      call random_number(xi)
      zstart(5,ipart)=2.d0*(xi-0.5d0)  ! starting pitch z(5)=v_\parallel / v
    enddo

    call save_starting_points(zstart)

  end subroutine sample_surface_fieldline
  
  subroutine sample_grid(zstart, grid_density)
    use params, only: ntestpart, zstart_dim1, zend, times_lost, &
        trap_par, perp_inv, iclass, xstart, sbeg
    use util, only: pi

    double precision, dimension(:,:), allocatable, intent(inout) :: zstart
    double precision, intent(in) :: grid_density
    double precision :: ngrid, xi
    integer :: xsize, ipart, lidx
    
    xsize = (2*pi) * grid_density !angle density
    ngrid = (1 / grid_density) - 1
    ntestpart = ngrid ** 2 !number of total angle points

    ! Resize particle coord. arrays and result memory.
    if (allocated(zstart)) deallocate(zstart)
    if (allocated(zend)) deallocate(zend)
    allocate(zstart(zstart_dim1,ntestpart), zend(zstart_dim1,ntestpart))
    if (allocated(times_lost)) deallocate(times_lost)
    if (allocated(trap_par)) deallocate(trap_par)
    if (allocated(perp_inv)) deallocate(perp_inv)
    if (allocated(iclass)) deallocate(iclass)
    allocate(times_lost(ntestpart), trap_par(ntestpart), perp_inv(ntestpart), iclass(3,ntestpart))
    
    do ipart=1,ngrid
      zstart(1,ipart) = sbeg(1)
      zstart(2,ipart) = xsize * ipart
      zstart(3,ipart) = xsize * ipart
      zstart(4,ipart)=1.d0  ! normalized velocity module z(4) = v / v_0
      call random_number(xi)
      zstart(5,ipart)=2.d0*(xi-0.5d0)  ! starting pitch z(5)=v_\parallel / v
      do jpart=1,ngrid
        lidx = (jpart-1)*ntestpart+ipart
        zstart(1,lidx) = sbeg(1)
        zstart(2,lidx) = xsize * jpart
        zstart(3,lidx) = xsize * jpart
        zstart(4,lidx) = 1.d0  ! normalized velocity module z(4) = v / v_0
        call random_number(xi)
        zstart(5,lidx)=2.d0*(xi-0.5d0)  ! starting pitch z(5)=v_\parallel / v
      end do 
    enddo

    call save_starting_points(zstart)

  end subroutine sample_grid

  subroutine sample_random_batch(zstart, reuse_existing)
  ! Get random batch from preexisting zstart, allows reuse.
    use params, only: batch_size, ntestpart, zstart_dim1, idx

    integer :: ran_begin, ran_end, ipart
    real :: temp_ran
    double precision, dimension(:,:), intent(inout) :: zstart
    double precision, dimension(zstart_dim1,batch_size) :: zstart_batch
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
