module classification
use omp_lib
use params, only: zstart, zend, times_lost, trap_par, perp_inv, iclass, &
    ntimstep, confpart_trap, confpart_pass, notrace_passing, contr_pp, &
    class_plot, ntcut, nturns, fast_class, n_tip_vars, nplagr, nder, npl_half, &
    nfp, fper, zerolam, num_surf, bmax, bmin, dtaumin, v0, cut_in_per, &
    integmode, relerr, ntau, should_skip
use util, only: twopi, sqrt2
use velo_mod,   only : isw_field_type
use orbit_symplectic, only : orbit_timestep_sympl, get_val
use simple, only : init_sympl, Tracer
use cut_detector, only : fract_dimension
use diag_mod, only : icounter
use get_can_sub, only : vmec_to_can, can_to_vmec
use boozer_sub, only : vmec_to_boozer, boozer_to_vmec
use magfie_sub, only : CANFLUX, BOOZER
use check_orbit_type_sub, only : check_orbit_type

  implicit none

  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)

  ! Classification result type - separates data from I/O
  ! Note: 0=unclassified means the classification was not computed
  ! This depends on orbit type (trapped/passing) and class_plot flag
  type :: classification_result_t
    logical :: passing           ! Trapped (false) or passing (true)
    logical :: lost              ! Orbit lost (true) or confined (false)
    integer :: fractal           ! Fractal: 0=unclassified, 1=regular, 2=chaotic
    integer :: jpar              ! J_parallel: 0=unclassified, 1=regular, 2=stochastic
    integer :: topology          ! Topology: 0=unclassified, 1=ideal, 2=non-ideal
  end type classification_result_t

  ! output files:
  ! iaaa_bou - trapped-passing boundary
  ! iaaa_pnt - forced regular passing
  ! iaaa_prp - lossed passing
  ! iaaa_prt - lossed trapped
  ! iaaa_rep - regular passing
  ! iaaa_ret - regular trapped
  ! iaaa_stp - stochastic passing
  ! iaaa_stt - stochastic trapped
integer, parameter :: iaaa_bou=20000, iaaa_pnt=10000, iaaa_prp=10001, iaaa_prt=10002, &
    iaaa_rep=10011, iaaa_ret=10012, iaaa_stp=10021, iaaa_stt=10022

  ! output files:
  ! iaaa_jre - regular trapped by J_parallel
  ! iaaa_jst - stochastic trapped by J_parallel
  ! iaaa_jer - non-classified trapped by J_parallel
  ! iaaa_ire - ideal trapped by recurrences and monotonicity
  ! iaaa_ist - non-ideal trapped by recurrences and monotonicity
  ! iaaa_ier - non-classified trapped by recurrences and monotonicity
integer, parameter :: iaaa_jre=40012, iaaa_jst=40022, iaaa_jer=40032, &
    iaaa_ire=50012, iaaa_ist=50022, iaaa_ier=50032

contains

subroutine trace_orbit_with_classifiers(anorb, ipart, class_result)
    use find_bminmax_sub, only : get_bminmax
    use magfie_sub, only : magfie
    use plag_coeff_sub, only : plag_coeff
  use alpha_lifetime_sub, only : orbit_timestep_axis, orbit_timestep_rk4, orbit_timestep_can

    type(Tracer), intent(inout) :: anorb
    integer, intent(in) :: ipart
    type(classification_result_t), intent(out) :: class_result
    integer :: ierr, ierr_coll
    real(dp), dimension(5) :: z
    real(dp) :: bmod,sqrtg
    real(dp), dimension(3) :: bder, hcovar, hctrvr, hcurl
    integer :: it, ktau
    integer(8) :: kt
    logical :: passing

    integer                                       :: ifp_tip,ifp_per
    integer,          dimension(:),   allocatable :: ipoi
    real(dp), dimension(:),   allocatable :: xp
    real(dp), dimension(:,:), allocatable :: coef,orb_sten
    real(dp), dimension(:,:), allocatable :: zpoipl_tip,zpoipl_per,dummy2d
    real(dp), dimension(n_tip_vars)       :: var_tip
    real(dp) :: phiper, alam_prev, par_inv
    integer :: iper, itip, kper, nfp_tip, nfp_per

    real(dp) :: fraction
    real(dp) :: r,theta_vmec,varphi_vmec
    logical :: regular

    ! Variables and settings for classification by J_parallel and ideal orbit condition:
    integer, parameter :: nfp_dim=3
    integer :: nfp_cot,ideal,ijpar,ierr_cot,iangvar
    real(dp), dimension(nfp_dim) :: fpr_in

    zend(:,ipart) = 0d0
    !
    iangvar=2
    ! End variables and settings for classification by J_parallel and ideal orbit condition
    !

    ! Initialize classification result - all unclassified
    class_result%passing = .false.
    class_result%lost = .false.
    class_result%fractal = 0
    class_result%jpar = 0
    class_result%topology = 0

    !  open(unit=10000+ipart, iostat=stat, status='old')
    !  if (stat == 0) close(10000+ipart, status='delete')
    !  open(unit=20000+ipart, iostat=stat, status='old')
    !  if (stat == 0) close(20000+ipart, status='delete')

    ! Write out trapped-passing boundary at the classification cut:
    if(class_plot) then
        if(ipart.eq.1) then
            z(1)=zstart(1,ipart)
            z(3)=cut_in_per*fper
            do kt=0,1000
                z(2)=1d-3*twopi*dble(kt)
                call magfie(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
                write(iaaa_bou,*) z(2),sqrt(1.d0-bmod/bmax)
            enddo
        endif
    endif
    ! End write out trapped-passing boundary at the classification cut
    !
    z = zstart(:, ipart)
    r=z(1)
    theta_vmec=z(2)
    varphi_vmec=z(3)
    !
    if(isw_field_type .eq. CANFLUX) then
        call vmec_to_can(r,theta_vmec,varphi_vmec,z(2),z(3))
    elseif(isw_field_type .eq. BOOZER) then
        call vmec_to_boozer(r,theta_vmec,varphi_vmec,z(2),z(3))
    endif

    ! In case of classification plot all starting points are moved to the classification cut:
    if(class_plot) then
        z(3)=cut_in_per*fper
        zstart(2,ipart)=modulo(zstart(2,ipart),twopi)
    endif
    ! End moving starting points to the classification cut

    if (integmode>0) call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)

    call magfie(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)

    !$omp critical
    if(num_surf > 1) then
        call get_bminmax(z(1),bmin,bmax)
    endif
    passing = z(5)**2.gt.1.d0-bmod/bmax
    trap_par(ipart) = ((1.d0-z(5)**2)*bmax/bmod-1.d0)*bmin/(bmax-bmin)
    perp_inv(ipart) = z(4)**2*(1.d0-z(5)**2)/bmod
    iclass(:,ipart) = 0
    !$omp end critical

    ! Store passing status in result
    class_result%passing = passing

    ! Forced classification of passing as regular:
    if(passing .and. should_skip(ipart)) then
        !$omp critical
        confpart_pass=confpart_pass+1.d0
        !$omp end critical
        if(class_plot) then
            !$omp critical
            write (iaaa_pnt,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            !$omp end critical
        endif
        iclass(:,ipart) = 1
        ! Mark as regular passing (fractal=1) and not lost
        class_result%fractal = 1
        class_result%lost = .false.
        return
    endif
    ! End forced classification of passing as regular

    !$omp critical
    if (.not. allocated(ipoi)) &
        allocate(ipoi(nplagr),coef(0:nder,nplagr),orb_sten(6,nplagr),xp(nplagr))
    !$omp end critical
    do it=1,nplagr
        ipoi(it)=it
    enddo

    nfp_tip=nfp             !<= initial array dimension for tips
    nfp_per=nfp             !<= initial array dimension for periods
    !$omp critical
    if (.not. allocated(zpoipl_tip)) &
        allocate(zpoipl_tip(2,nfp_tip),zpoipl_per(2,nfp_per))
    !$omp end critical

    !  open(unit=10000+ipart, recl=1024, position='append')
    !  open(unit=20000+ipart, recl=1024, position='append')

    ifp_tip=0               !<= initialize footprint counter on tips
    ifp_per=0               !<= initialize footprint counter on periods

    icounter=0
    phiper=0.0d0


    kt = 0
    if (passing) then
        !$omp atomic
        confpart_pass(1)=confpart_pass(1)+1.d0
    else
        !$omp atomic
        confpart_trap(1)=confpart_trap(1)+1.d0
    end if

    !--------------------------------
    ! Initialize tip detector

    itip=npl_half+1
    alam_prev=z(5)

    ! End initialize tip detector
    !--------------------------------
    ! Initialize period crossing detector

    iper=npl_half+1
    kper=int(z(3)/fper)

    ! End initialize period crossing detector
    !--------------------------------
    !
    ! Initialize classification by J_parallel and ideal orbit condition:
    nfp_cot=0
    ! End Initialize classification by J_parallel and ideal orbit condition
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
            if (integmode == -2) then
                call orbit_timestep_rk4(z, dtaumin, ierr)
            else if (integmode == -1) then
                call orbit_timestep_can(z, dtaumin, dtaumin, relerr, ierr)
            else if (integmode <= 0) then
                call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
            else
                call orbit_timestep_sympl(anorb%si, anorb%f, ierr)
                z(1:3) = anorb%si%z(1:3)
                z(4) = dsqrt(anorb%f%mu*anorb%f%Bmod+0.5d0*anorb%f%vpar**2)
                z(5) = anorb%f%vpar/(z(4)*sqrt2)
            endif

            ! Write starting data for orbits which were lost in case of classification plot
            ! Mark orbit as lost if integration failed
            if(ierr.ne.0) then
                class_result%lost = .true.
                if(class_plot) then
                    call output_lost_orbit_starting_data(ipart, passing)
                endif
                exit
            endif
            ! End handling of lost orbits
            kt = kt+1

            par_inv = par_inv+z(5)**2*dtaumin ! parallel adiabatic invariant
            if(kt.le.nplagr) then          !<=first nplagr points to initialize stencil
                orb_sten(1:5,kt)=z
                orb_sten(6,kt)=par_inv
            else                          !<=normal case, shift stencil
                orb_sten(1:5,ipoi(1))=z
                orb_sten(6,ipoi(1))=par_inv
                ipoi=cshift(ipoi,1)
            endif

            ! Tip detection and interpolation
            if(alam_prev.lt.0.d0.and.z(5).gt.0.d0) itip=0   !<=tip has been passed
            itip=itip+1
            alam_prev=z(5)
            if(kt.gt.nplagr) then          !<=use only initialized stencil
                if(itip.eq.npl_half) then   !<=stencil around tip is complete, interpolate
                    xp=orb_sten(5,ipoi)

                    call plag_coeff(nplagr,nder,zerolam,xp,coef)

                    var_tip=matmul(orb_sten(:,ipoi),coef(0,:))
                    var_tip(2)=modulo(var_tip(2),twopi)
                    var_tip(3)=modulo(var_tip(3),twopi)

                    !          write(10000+ipart,*) var_tip

                    ifp_tip=ifp_tip+1
                    if(ifp_tip.gt.nfp_tip) then   !<=increase the buffer for banana tips
                        !$omp critical
                        allocate(dummy2d(2,ifp_tip-1))
                        !$omp end critical
                        dummy2d=zpoipl_tip(:,1:ifp_tip-1)
                        !$omp critical
                        deallocate(zpoipl_tip)
                        !$omp end critical
                        nfp_tip=nfp_tip+nfp
                        !$omp critical
                        allocate(zpoipl_tip(2,nfp_tip))
                        !$omp end critical
                        zpoipl_tip(:,1:ifp_tip-1)=dummy2d
                        !$omp critical
                        deallocate(dummy2d)
                        !$omp end critical
                    endif
                    zpoipl_tip(:,ifp_tip)=var_tip(1:2)
                    par_inv = par_inv - var_tip(6)
                    !
                    ! Classification by J_parallel and ideal orbit conditions:
                    fpr_in(1)=var_tip(1)
                    fpr_in(2)=var_tip(iangvar)
                    fpr_in(3)=var_tip(6)
                    !
                    call check_orbit_type(nturns,nfp_cot,fpr_in,ideal,ijpar,ierr_cot)
                    !
                    iclass(1,ipart) = ijpar
                    iclass(2,ipart) = ideal
                    ! Store in classification result
                    class_result%jpar = ijpar
                    class_result%topology = ideal
                    if(fast_class) ierr=ierr_cot
                    !
                    ! End classification by J_parallel and ideal orbit conditions
                endif
            endif
            ! End tip detection and interpolation

            ! Periodic boundary footprint detection and interpolation
            if(z(3).gt.dble(kper+1)*fper) then
                iper=0   !<=periodic boundary has been passed
                phiper=dble(kper+1)*fper
                kper=kper+1
            elseif(z(3).lt.dble(kper)*fper) then
                iper=0   !<=periodic boundary has been passed
                phiper=dble(kper)*fper
                kper=kper-1
            endif
            iper=iper+1
            if(kt.gt.nplagr) then          !<=use only initialized stencil
                if(iper.eq.npl_half) then   !<=stencil around periodic boundary is complete, interpolate
                    xp=orb_sten(3,ipoi)-phiper

                    call plag_coeff(nplagr,nder,zerolam,xp,coef)

                    var_tip=matmul(orb_sten(:,ipoi),coef(0,:))
                    var_tip(2)=modulo(var_tip(2),twopi)
                    var_tip(3)=modulo(var_tip(3),twopi)
                    ! write(20000+ipart,*) var_tip
                    ifp_per=ifp_per+1
                    if(ifp_per.gt.nfp_per) then   !<=increase the buffer for periodic boundary footprints
                        !$omp critical
                        allocate(dummy2d(2,ifp_per-1))
                        !$omp end critical
                        dummy2d=zpoipl_per(:,1:ifp_per-1)
                        !$omp critical
                        deallocate(zpoipl_per)
                        !$omp end critical
                        nfp_per=nfp_per+nfp
                        !$omp critical
                        allocate(zpoipl_per(2,nfp_per))
                        !$omp end critical
                        zpoipl_per(:,1:ifp_per-1)=dummy2d
                        !$omp critical
                        deallocate(dummy2d)
                        !$omp end critical
                    endif
                    zpoipl_per(:,ifp_per)=var_tip(1:2)
                endif
            endif
            ! End periodic boundary footprint detection and interpolation

            ! Cut classification into regular or chaotic
            if (kt == ntcut) then
                regular = .True.

                if(ifp_per > 0) then

                    call fract_dimension(ifp_per,zpoipl_per(:,1:ifp_per),fraction)

                    if(fraction.gt.0.2d0) then
                        print *, ipart, ' chaotic per ', ifp_per
                        regular = .False.
                    else
                        print *, ipart, ' regular per', ifp_per
                    endif
                endif

                if(ifp_tip > 0) then

                    call fract_dimension(ifp_tip,zpoipl_tip(:,1:ifp_tip),fraction)

                    if(fraction.gt.0.2d0) then
                        print *, ipart, ' chaotic tip ', ifp_tip
                        regular = .False.
                        iclass(3,ipart) = 2
                    else
                        print *, ipart, ' regular tip ', ifp_tip
                        iclass(3,ipart) = 1
                    endif
                endif

                ! Store fractal classification in result
                if(regular) then
                    class_result%fractal = 1
                else
                    class_result%fractal = 2
                endif

                if(class_plot) then
                    call output_minkowsky_class(ipart, regular, passing)
                    ierr=1
                endif
            endif
            !
            if(ierr.ne.0) then
                if(class_plot .and. .not. passing) then
                    call output_jpar_class(ipart, ijpar)
                    call output_topological_class(ipart, ideal)
                    exit
                endif
            endif
            !    write(999, *) kt*dtaumin/v0, z
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
    if(isw_field_type .eq. CANFLUX) then
        call can_to_vmec(z(1),z(2),z(3),zend(2,ipart),zend(3,ipart))
    elseif(isw_field_type .eq. BOOZER) then
        call boozer_to_vmec(z(1),z(2),z(3),zend(2,ipart),zend(3,ipart))
    endif
    times_lost(ipart) = kt*dtaumin/v0
    deallocate(zpoipl_tip, zpoipl_per)
    !$omp end critical
    !  close(unit=10000+ipart)
    !  close(unit=10000+ipart)
end subroutine trace_orbit_with_classifiers


subroutine output_lost_orbit_starting_data(ipart, passing)
    integer, intent(in) :: ipart
    logical, intent(in) :: passing

    if(passing) then
        call write_output_line(iaaa_prp, ipart)
    else
        call write_output_line(iaaa_prt, ipart)
    endif
end subroutine output_lost_orbit_starting_data


subroutine output_minkowsky_class(ipart, regular, passing)
    integer, intent(in) :: ipart
    logical, intent(in) :: regular, passing

    if(regular) then
        if(passing) then
            call write_output_line(iaaa_rep, ipart)
        else
            call write_output_line(iaaa_ret, ipart)
        endif
    else
        if(passing) then
            call write_output_line(iaaa_stp, ipart)
        else
            call write_output_line(iaaa_stt, ipart)
        endif
    endif
end subroutine output_minkowsky_class


subroutine output_jpar_class(ipart, ijpar)
    integer, intent(in) :: ipart, ijpar

    select case(ijpar)
        case(0)
        call write_output_line(iaaa_jer, ipart)
        case(1)
        call write_output_line(iaaa_jre, ipart)
        case(2)
        call write_output_line(iaaa_jst, ipart)
    end select
end subroutine output_jpar_class


subroutine output_topological_class(ipart, ideal)
    integer, intent(in) :: ipart, ideal

    select case(ideal)
        case(0)
        call write_output_line(iaaa_ier, ipart)
        case(1)
        call write_output_line(iaaa_ire, ipart)
        case(2)
        call write_output_line(iaaa_ist, ipart)
    end select
end subroutine output_topological_class


subroutine write_output_line(iunit, ipart)
    integer, intent(in) :: iunit, ipart

    !$omp critical
    write (iunit, *) zstart(2,ipart), zstart(5,ipart), trap_par(ipart)
    !$omp end critical
end subroutine write_output_line


! Write classification results to fort.* files
! This subroutine centralizes all classification file I/O
! Note: Only writes classifications that were computed (non-zero values)
subroutine write_classification_results(ipart, class_result)
    integer, intent(in) :: ipart
    type(classification_result_t), intent(in) :: class_result
    logical :: regular

    ! Write lost orbits
    if(class_result%lost) then
        call output_lost_orbit_starting_data(ipart, class_result%passing)
        return
    endif

    ! Write fractal classification if computed
    if(class_result%fractal /= 0) then
        regular = (class_result%fractal == 1)
        call output_minkowsky_class(ipart, regular, class_result%passing)
    endif

    ! Write J_parallel and topological classification if computed
    ! These are only done for trapped orbits
    if(.not. class_result%passing) then
        if(class_result%jpar /= 0) then
            call output_jpar_class(ipart, class_result%jpar)
        endif
        if(class_result%topology /= 0) then
            call output_topological_class(ipart, class_result%topology)
        endif
    endif

end subroutine write_classification_results

end module classification
