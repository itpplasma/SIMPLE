module samplers
  use util

  implicit none

  character(len=*), parameter :: START_FILE = 'start.dat'
  character(len=*), parameter :: START_FILE_ANTS = 'start_ants.dat'
  character(len=*), parameter :: START_FILE_BATCH = 'batch.dat'
  contains
  ! Functions #################################
  subroutine init_starting_surf
    use alpha_lifetime_sub, only : integrate_mfl_can

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
    character(len=*), parameter, intent(in) :: filename = START_FILE
    integer :: ipart

    open(1,file=START_FILE,recl=1024)
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
      character(len=*), parameter, intent(in) :: filename = START_FILE
  
      call load_starting_points(zstart, filename)
  end subroutine

  subroutine sample_volume_single(zstart, s_inner, s_outer)
    use params, only: isw_field_type
    use boozer_sub, only: boozer_to_ref
    use get_can_sub, only: can_to_ref

    double precision, intent(in) :: s_inner
    double precision, intent(in) :: s_outer
    double precision :: tmp_rand
    double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec
    double precision, dimension(:,:), intent(inout) :: zstart
    integer :: ipart

    do ipart=1,size(zstart,2)
      call random_number(tmp_rand)
      if 0 == num_surf r = xi
      else r = tmp_rand * (s_outer - s_inner) + s_inner
      endif
     
      call random_number(tmp_rand)
      vartheta=twopi*tmp_rand
      call random_number(tmp_rand)
      varphi=twopi*tmp_rand
! we store starting points in VMEC coordinates:
      if(isw_field_type.eq.0) then
          call can_to_ref(r,vartheta,varphi,theta_vmec,varphi_vmec)
      elseif(isw_field_type.eq.1) then
          theta_vmec=vartheta
          varphi_vmec=varphi
      elseif(isw_field_type.eq.2) then
          call boozer_to_ref(r,vartheta,varphi,theta_vmec,varphi_vmec)
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
      call random_number(tmp_rand)
      zstart(5,ipart)=2.d0*(tmp_rand-0.5d0)
    enddo

    call save_starting_points(zstart)

  end subroutine sample_volume_single

  subroutine sample_surface_fieldline(zstart)
    use params, only: volstart, isw_field_type, ibins, xstart, npoiper, nper
    use boozer_sub, only: boozer_to_ref
    use get_can_sub, only: can_to_ref
    use binsrc_sub, only: binsrc

    double precision, dimension(:,:), intent(inout) :: zstart

    double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec
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

      ! we store starting points in VMEC coordinates:
      if(isw_field_type.eq.0) then
        call can_to_ref(r,vartheta,varphi,theta_vmec,varphi_vmec)
      elseif(isw_field_type.eq.1) then
        theta_vmec=vartheta
        varphi_vmec=varphi
      elseif(isw_field_type.eq.2) then
        call boozer_to_ref(r,vartheta,varphi,theta_vmec,varphi_vmec)
      else
        print *,'init_starting_points: unknown field type'
      endif

      zstart(1,ipart)=r
      zstart(2,ipart)=theta_vmec
      zstart(3,ipart)=varphi_vmec
      zstart(4,ipart)=1.d0  ! normalized velocity module z(4) = v / v_0
      call random_number(xi)
      zstart(5,ipart)=2.d0*(xi-0.5d0)  ! starting pitch z(5)=v_\parallel / v
    enddo

    call save_starting_points(zstart)

  end subroutine sample_surface_fieldline

  
  subroutine sample_random_batch(zstart, reuse_existing)
  ! Get random batch from preexisting zstart, allows reuse.
    use params, only: batch_size, ntestpart
    
    integer :: ran_begin, ran_end, ipart
    double precision, dimension(:,:) :: zstart_batch
    double precision, dimension(:,:), intent(inout) :: zstart
    logical, intent(in) :: reuse_existing

    if (reuse_batch.eq.1) then
      call load_starting_points(zstart_batch, START_FILE_BATCH)
    else
      call load_starting_points(zstart_batch, START_FILE)
      call random_number(ran_begin)
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
    
    deallocate(zstart_batch)
    
  end subroutine
  
  subroutine sample_points_ants(use_special_ants_file)
    use parse_ants, only : process_line
    use get_can_sub, only : vmec_to_can
    
    logical, intent(in) :: use_special_ants_file
    
    integer, parameter :: maxlen = 4096
    character(len=maxlen) :: line
    real(8) :: v_par, v_perp, u, v, s
    real(8) :: th, ph
    integer :: ipart

    do ipart=1,ntestpart
      if use_special_ants_file then
        read(START_FILE_ANTS, '(A)') line
      else
        read(START_FILE, '(A)') line
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
  end subroutine
  ! Interface #################################

  INTERFACE sample
    
    FUNCTION sample_read(zstart, filename)
      double precision, dimension(:,:), intent(inout) :: zstart
      character(len=*), parameter, intent(in) :: filename = START_FILE
    END FUNCTION sample_read

    FUNCTION sample_surface_fieldline(zstart)
      double precision, dimension(:,:), intent(inout) :: zstart
    END FUNCTION sample_surface_fieldline

    FUNCTION sample_volume_single(zstart, s_inner, s_outer)!global, mode 5
      double precision, dimension(:,:), intent(inout) :: zstart
      real, intent(in) :: s_inner
      real, intent(in) :: s_outer
    END FUNCTION sample_volume_single
    
    FUNCTION sample_random_batch(zstart, reuse_existing)
      double precision, dimension(:,:), intent(inout) :: zstart
      logical, intent(in) :: reuse_existing
    END FUNCTION
    
    FUNCTION sample_points_ants(use_special_ants_file)
      logical, intent(in) :: use_special_ants_file
    END FUNCTION
    
  END INTERFACE sample
end module samplers
