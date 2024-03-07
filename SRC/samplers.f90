

! Functions #################################
FUNCTION sample_surface_read()
    ! Read from previous start.dat (startmode = 2)
    double precision, dimension(:,:) :: zstart ! Assume allocation in params_init
    open(1,file='start.dat',recl=1024)
    do ipart=1,ntestpart
        read(1,*) zstart(:,ipart)
    enddo

    sample_surface_read = zstart

END FUNCTION sample_surface_read

FUNCTION sample_volume_single(s_inner, s_outter)
    real, intent(in) :: s_inner
    real, intent(in) :: s_outter
    real :: s_rand, tmp_rand
    double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec
    double precision, dimension(:,:) :: zstart ! Assume allocation in params_init

    ! files for storing starting coords
    open(1,file='start.dat',recl=1024)

    do ipart=1,ntestpart
        s_rand = -0.1
        do while ((s_inner .gt. s_rand) .or. (s_outter .lt. s_rand))
            call random_number(s_rand)
        enddo

        call random_number(s_rand)
        vartheta=twopi*s_rand
        call random_number(s_rand)
        varphi=twopi*s_rand
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
        xi=zzg()
        zstart(5,ipart)=2.d0*(xi-0.5d0)
        write(1,*) zstart(:,ipart)
    enddo
    close(1)

    sample_volume_single = zstart

END FUNCTION sample_volume_single

FUNCTION sample_surface_regular_grid(n_start)
    integer, intent(in) :: n_start
    double precision, dimension(:,:) :: zstart ! Assume allocation in params_init

! files for storing starting coords
    open(1,file='start.dat',recl=1024)
    ! determine the starting point:
    if (startmode == 0 .or. startmode == 1) then
        do ipart=1,ntestpart
            xi=zzg()
            call binsrc(volstart,1,npoi,xi,i)
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
            xi=zzg()
            zstart(5,ipart)=2.d0*(xi-0.5d0)
            write(1,*) zstart(:,ipart)
        enddo
    endif
    close(1)

    sample_surface_regular_grid = zstart
END FUNCTION sample_surface_regular_grid


! Interface #################################

INTERFACE sample

    FUNCTION sample_surface_read()
    END FUNCTION sample_surface_read

    FUNCTION sample_surface_regular_grid(n_start)
        integer, intent(in) :: n_start
    END FUNCTION sample_surface_regular_grid

    FUNCTION sample_volume_single(s_inner, s_outter)
        real, intent(in) :: s_inner
        real, intent(in) :: s_outter
    END FUNCTION sample_volume_single

END INTERFACE sample