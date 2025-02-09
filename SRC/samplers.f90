module samplers
  use util

  implicit none

  character(len=*), parameter :: START_FILE = 'start.dat'

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
      call random_number(tmp_rand)
      zstart(5,ipart)=2.d0*(tmp_rand-0.5d0)
    enddo

    call save_starting_points(zstart)

  end subroutine sample_volume_single

  subroutine sample_surface_fieldline(zstart)
    use params, only: volstart, isw_field_type, ibins, xstart, npoiper, nper
    use boozer_sub, only: boozer_to_vmec
    use get_can_sub, only: can_to_vmec
    use binsrc_sub, only: binsrc

    double precision, dimension(:,:), intent(inout) :: zstart

    double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec
    double precision :: xi
    integer :: ipart, i

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
        call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      elseif(isw_field_type.eq.1) then
        theta_vmec=vartheta
        varphi_vmec=varphi
      elseif(isw_field_type.eq.2) then
        call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
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

  !FUNCTION sample_surface_regular_grid(phibeg, thetabeg)
  !  !TODO is the grid one above this one? then what is the grid one?

  !END FUNCTION sample_surface_regular_grid
  ! Interface #################################

  INTERFACE sample
    
    !load_starting_points(zstart)
    FUNCTION sample_surface_read(zstart)
      double precision, dimension(:,:), intent(inout) :: zstart
    END FUNCTION sample_surface_read

    FUNCTION sample_surface_fieldline(n_start)
      integer, intent(in) :: n_start
    END FUNCTION sample_surface_fieldline

    FUNCTION sample_volume_single(n_start, s_inner, s_outer)
      integer, intent(in) :: n_start
      real, intent(in) :: s_inner
      real, intent(in) :: s_outer
    END FUNCTION sample_volume_single

  END INTERFACE sample
end module samplers
