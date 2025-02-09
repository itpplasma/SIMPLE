module samplers
  use util

  implicit none

  character(len=*), parameter :: START_FILE = 'start.dat'

  INTERFACE sample
    
    MODULE PROCEDURE sample_volume_single
    MODULE PROCEDURE sample_surface_fieldline
    MODULE PROCEDURE sample_random_batch

  END INTERFACE sample

  contains
  ! Functions #################################
  subroutine load_starting_points(zstart)
    double precision, dimension(:,:), intent(inout) :: zstart
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

  subroutine sample_volume_single(zstart, s_inner, s_outer)
    use params, only: isw_field_type
    use boozer_sub, only: boozer_to_vmec
    use get_can_sub, only: can_to_vmec

    double precision, intent(in) :: s_inner
    double precision, intent(in) :: s_outer
    double precision :: tmp_rand
    double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec
    double precision, dimension(:,:), intent(inout) :: zstart
    integer :: ipart

    do ipart=1,size(zstart,2)
      call random_number(tmp_rand)
      r = tmp_rand * (s_outer - s_inner) + s_inner

      call random_number(tmp_rand)
      vartheta=twopi*tmp_rand
      call random_number(tmp_rand)
      varphi=twopi*tmp_rand

      !zstart is in VMEC-coordinates, a case construct on `isw_field_type` is here
      select case(isw_field_type)
        case(0,2,3)
          theta_vmec = vartheta
          varphi_vmec = varphi
        case(1,11) !boozer
          call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec) 
        case(21,31,32) !can
          call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec) 
      end select
      zstart(1,ipart)=r
      zstart(2,ipart)=theta_vmec
      zstart(3,ipart)=varphi_vmec
      zstart(4,ipart)=1.0d-3
      zstart(5,ipart)=1.0d0

    enddo

    !save to file
    call save_starting_points(zstart)

  end subroutine sample_volume_single

  subroutine sample_surface_fieldline(zstart)
    double precision, dimension(:,:), intent(inout) :: zstart
    integer :: ipart

    !to implement
    print *,"WARNING sample_surface_fieldline  not yet implemented!"
    !NEEDS: flux surface in canonical coordinates!!!!!!!!!!

    print *,"STUB: start points on a r=0.1 flux surface on field line in Boozer on phi-theta-grid"
    do ipart=1,size(zstart,2)
      ! this assumes particles on a canonical flux surface
      ! on Boozer this means constant Aphi and Ath
      ! where Aphi= psi_tor - aiota * psi_pol = Psi - aiota * Psi_pol;
      ! Psi is the toroidal flux devided by 2 pi normalized to 1 at the boundary
      zstart(1,ipart)=0.1d0
      zstart(2,ipart)=twopi/size(zstart,2)*dble(ipart)
      zstart(3,ipart)=twopi/size(zstart,2)*dble(ipart)
      zstart(4,ipart)=1.0d-3
      zstart(5,ipart)=1.0d0
    enddo

    call save_starting_points(zstart)

  end subroutine sample_surface_fieldline

  !FUNCTION sample_surface_regular_grid(phibeg, thetabeg)
  !  !TODO is the grid one above this one? then what is the grid one?

  !END FUNCTION sample_surface_regular_grid
  
  
  subroutine sample_random_batch(idx_begin, idx_end, zstart)
  ! Get random batch from preexisting zstart, allows reuse.
    use params, only: batch_size, reuse_batch, ntestpart
    
    integer :: ran_begin, ran_end, ipart
    integer, intent(in) :: idx_begin, idx_end
    double precision :: tmp_rand
    double precision, dimension(:,:), allocatable :: zstart_batch
    double precision, dimension(:,:), intent(inout) :: zstart
    
    allocate(zstart_batch(size(zstart,1), ntestpart))
    call load_starting_points(zstart_batch)
    if (reuse_batch) then
      do ipart=idx_begin, idx_end
        zstart(:,ipart-idx_begin+1) = zstart_batch(:,ipart)
      enddo
    else
      call random_number(tmp_rand)
      ran_begin = int(tmp_rand * (ntestpart - batch_size)) + 1
      ran_end = ran_begin+batch_size-1
      if (ran_end.gt.ntestpart) then
        ran_begin = ran_begin - (ran_end-ntestpart)
        ran_end = ntestpart
      endif
      do ipart=1,batch_size
        zstart(:,ipart) = zstart_batch(:,ran_begin+ipart-1)
      enddo
    endif 
    
    deallocate(zstart_batch)
    
  end subroutine sample_random_batch

end module samplers