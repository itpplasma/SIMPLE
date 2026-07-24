module alpha_lifetime_sub

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    integer, parameter :: axis_rhs_count = 5

    type :: toroidal_regularization_context_t
        real(dp) :: mu = 0.0_dp
        real(dp) :: major_radius_cm = 1.0_dp
    end type toroidal_regularization_context_t

    ! The Burby--Ellison position map is retained only through first order in
    ! gyroradius. Near a mirror point, its neglected second-order terms can put
    ! the mapped state just below the original constant-momentum energy shell.
    ! Clamp only that small mismatch; a larger negative radicand is a failed
    ! regularized step and must fall through to the next recovery method.
    real(dp), parameter :: toroidal_energy_shell_tolerance = 1.0e-5_dp

contains

    subroutine bstar_parallel_factor(z, factor)
        use parmot_mod, only: ro0
        use magfie_sub, only: magfie
        real(dp), intent(in) :: z(5)
        real(dp), intent(out) :: factor
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie(z(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        factor = 1.0_dp + z(4)*z(5)*(ro0/bmod)*dot_product(hcurl, hcovar)
    end subroutine bstar_parallel_factor

    subroutine toroidal_regularized_position_map(x, ppar, direction, x_mapped, &
            ierr)
        !> Leading near-identity map from Burby and Ellison, Phys. Plasmas 24,
        !> 110703 (2017), Eq. (13). direction=1 maps the conventional guiding
        !> centre to the toroidally regularized coordinate; direction=-1
        !> numerically inverts that first-order map. The toroidal component is
        !> unchanged, which preserves field-period closure without an angular
        !> seam.
        real(dp), intent(in) :: x(3), ppar
        integer, intent(in) :: direction
        real(dp), intent(out) :: x_mapped(3)
        integer, intent(out) :: ierr
        real(dp), parameter :: derivative_step = 1.0e-6_dp
        real(dp) :: shift(3), shift_plus(3), shift_minus(3)
        real(dp) :: residual(2), residual_trial(2), jacobian(2, 2), delta(2)
        real(dp) :: x_plus(3), x_minus(3), x_trial(3), determinant, damping
        real(dp) :: residual_norm, trial_norm
        integer :: iteration, line_search, shift_error

        ierr = 2
        x_mapped = x
        if (abs(direction) /= 1) return
        if (direction == 1) then
            call toroidal_position_shift(x, ppar, shift, shift_error)
            if (shift_error /= 0) return
            x_mapped = x + shift
            ierr = 0
            return
        end if

        ! Invert x_regular=x_standard+shift(x_standard) with a damped two-
        ! dimensional Newton solve. Fixed-point inversion is not contractive
        ! near the conventional B_parallel_star singularity.
        do iteration = 1, 20
            call toroidal_position_shift(x_mapped, ppar, shift, shift_error)
            if (shift_error /= 0) then
                x_mapped = x
                return
            end if
            residual = x_mapped(1:2) + shift(1:2) - x(1:2)
            residual_norm = maxval(abs(residual))
            if (residual_norm <= &
                    1.0e-12_dp*max(1.0_dp, maxval(abs(x)))) then
                ierr = 0
                return
            end if
            do line_search = 1, 2
                x_plus = x_mapped
                x_minus = x_mapped
                x_plus(line_search) = x_plus(line_search) + derivative_step
                x_minus(line_search) = x_minus(line_search) - derivative_step
                call toroidal_position_shift(x_plus, ppar, shift_plus, &
                    shift_error)
                if (shift_error /= 0) then
                    x_mapped = x
                    return
                end if
                call toroidal_position_shift(x_minus, ppar, shift_minus, &
                    shift_error)
                if (shift_error /= 0) then
                    x_mapped = x
                    return
                end if
                jacobian(:, line_search) = &
                    (shift_plus(1:2) - shift_minus(1:2))/ &
                    (2.0_dp*derivative_step)
                jacobian(line_search, line_search) = &
                    jacobian(line_search, line_search) + 1.0_dp
            end do
            determinant = jacobian(1, 1)*jacobian(2, 2) - &
                jacobian(1, 2)*jacobian(2, 1)
            if (abs(determinant) <= 100.0_dp*tiny(1.0_dp)) then
                x_mapped = x
                return
            end if
            delta(1) = (-jacobian(2, 2)*residual(1) + &
                jacobian(1, 2)*residual(2))/determinant
            delta(2) = (jacobian(2, 1)*residual(1) - &
                jacobian(1, 1)*residual(2))/determinant

            damping = 1.0_dp
            do line_search = 1, 12
                x_trial = x_mapped
                x_trial(1:2) = x_trial(1:2) + damping*delta
                if (x_trial(1) < 0.0_dp) then
                    damping = 0.5_dp*damping
                    cycle
                end if
                call toroidal_position_shift(x_trial, ppar, shift_plus, &
                    shift_error)
                if (shift_error /= 0) then
                    damping = 0.5_dp*damping
                    cycle
                end if
                residual_trial = x_trial(1:2) + shift_plus(1:2) - x(1:2)
                trial_norm = maxval(abs(residual_trial))
                if (trial_norm < residual_norm) then
                    x_mapped = x_trial
                    exit
                end if
                damping = 0.5_dp*damping
            end do
            if (line_search > 12) then
                x_mapped = x
                return
            end if
        end do
        x_mapped = x
    end subroutine toroidal_regularized_position_map

    subroutine toroidal_position_shift(x, ppar, shift, ierr)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use parmot_mod, only: ro0
        use magfie_sub, only: magfie_refcoords
        real(dp), intent(in) :: x(3), ppar
        real(dp), intent(out) :: shift(3)
        integer, intent(out) :: ierr
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: bphi

        ierr = 2
        shift = 0.0_dp
        call magfie_refcoords(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        bphi = bmod*hctrvr(3)
        if (abs(sqrtg*bphi) <= 100.0_dp*tiny(1.0_dp)) return
        shift(1) = ro0*ppar*hcovar(2)/(sqrtg*bphi)
        shift(2) = -ro0*ppar*hcovar(1)/(sqrtg*bphi)
        if (.not. all(ieee_is_finite(shift))) then
            shift = 0.0_dp
            return
        end if
        ierr = 0
    end subroutine toroidal_position_shift

    subroutine toroidal_regularized_rhs(tau, y, dy, context)
        use parmot_mod, only: rmu, ro0
        use magfie_sub, only: magfie_refcoords
        real(dp), intent(in) :: tau, y(:)
        real(dp), intent(out) :: dy(:)
        class(*), intent(in) :: context
        real(dp), parameter :: derivative_step = 1.0e-5_dp
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: bp, gp, bdp(3), hcp(3), htp(3), hxp(3)
        real(dp) :: bm, gm, bdm(3), hcm(3), htm(3), hxm(3)
        real(dp) :: xwork(3), dq(3), dkinetic(3), q, u, ppar, vpa
        real(dp) :: drift_scale, gamma
        real(dp) :: h
        integer :: i

        select type (ctx => context)
        type is (toroidal_regularization_context_t)
            call magfie_refcoords(y(1:3), bmod, sqrtg, bder, hcovar, &
                hctrvr, hcurl)
            q = ctx%major_radius_cm*hctrvr(3)
            if (abs(q) <= 100.0_dp*tiny(1.0_dp)) then
                dy = huge(1.0_dp)
                return
            end if

            do i = 1, 3
                h = derivative_step
                xwork = y(1:3)
                xwork(i) = xwork(i) + h
                call magfie_refcoords(xwork, bp, gp, bdp, hcp, htp, hxp)
                xwork(i) = y(i) - h
                call magfie_refcoords(xwork, bm, gm, bdm, hcm, htm, hxm)
                dq(i) = ctx%major_radius_cm*(htp(3) - htm(3))/(2.0_dp*h)
            end do

            u = y(4)
            ppar = q*u
            gamma = sqrt(1.0_dp + &
                2.0_dp*(ppar**2 + 2.0_dp*ctx%mu*bmod)/rmu)
            vpa = ppar/gamma
            dkinetic = (u*u*q*dq + ctx%mu*bmod*bder)/gamma
            drift_scale = ro0/(sqrtg*bmod*hctrvr(3))
            dy(1) = vpa*hctrvr(1) - drift_scale*dkinetic(2)
            dy(2) = vpa*hctrvr(2) + drift_scale*dkinetic(1)
            dy(3) = vpa*hctrvr(3)
            dy(4) = -dot_product(hctrvr, dkinetic)/q
        class default
            dy = huge(1.0_dp)
        end select
    end subroutine toroidal_regularized_rhs

    subroutine orbit_timestep_toroidal_regularized(z, dtau, relerr, ierr, &
            step_limit)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use odeint_allroutines_sub, only: odeint_allroutines, &
            odeint_has_failed
        use chamb_sub, only: chamb_can
        use magfie_sub, only: magfie_refcoords, has_magfie_refcoords_field
        use field_can_mod, only: integ_to_ref, ref_to_integ
        use new_vmec_stuff_mod, only: rmajor
        real(dp), intent(inout) :: z(5)
        real(dp), intent(in) :: dtau, relerr
        integer, intent(out) :: ierr
        integer, intent(in) :: step_limit
        type(toroidal_regularization_context_t) :: context
        real(dp) :: z_initial(5), y(4), chamber_y(2)
        real(dp) :: x_standard(3), x_regular(3)
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: q, ppar, ppar_raw, ppar_squared
        integer :: map_error

        z_initial = z
        if (.not. has_magfie_refcoords_field()) then
            ierr = 2
            return
        end if
        call integ_to_ref(z(1:3), x_standard)
        call magfie_refcoords(x_standard, bmod, sqrtg, bder, hcovar, hctrvr, &
            hcurl)
        context%major_radius_cm = 100.0_dp*rmajor
        context%mu = 0.5_dp*z(4)**2*(1.0_dp - z(5)**2)/bmod
        ppar = z(4)*z(5)
        call toroidal_regularized_position_map(x_standard, ppar, 1, x_regular, &
            map_error)
        if (map_error /= 0) then
            ierr = 2
            return
        end if
        y(1:3) = x_regular
        call magfie_refcoords(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        q = context%major_radius_cm*hctrvr(3)
        if (abs(q) <= 100.0_dp*tiny(1.0_dp)) then
            ierr = 2
            return
        end if
        y(4) = ppar/q

        call odeint_allroutines(y, 4, context, 0.0_dp, dtau, relerr, &
            toroidal_regularized_rhs, step_limit=step_limit)
        if (odeint_has_failed() .or. .not. all(ieee_is_finite(y))) then
            z = z_initial
            ierr = 2
            return
        end if

        call magfie_refcoords(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        q = context%major_radius_cm*hctrvr(3)
        ppar_raw = q*y(4)
        call toroidal_regularized_position_map(y(1:3), ppar_raw, -1, &
            x_standard, map_error)
        if (map_error /= 0) then
            z = z_initial
            ierr = 2
            return
        end if
        call magfie_refcoords(x_standard, bmod, sqrtg, bder, hcovar, hctrvr, &
            hcurl)
        ppar_squared = z_initial(4)**2 - 2.0_dp*context%mu*bmod
        if (ppar_squared < &
                -toroidal_energy_shell_tolerance*z_initial(4)**2) then
            z = z_initial
            ierr = 2
            return
        end if
        ppar = sign(sqrt(max(ppar_squared, 0.0_dp)), ppar_raw)
        call ref_to_integ(x_standard, z(1:3))
        z(4) = z_initial(4)
        z(5) = ppar/z(4)
        chamber_y = z(1:2)
        call chamb_can(chamber_y, z(3), ierr)
    end subroutine orbit_timestep_toroidal_regularized

    !ToDo make module from all global things like this one
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine elefie_can(x,derphi)
        !
        ! Computes the derivatives of the electrostatic potential over
        ! coordinates (covariant electric field). Potential is normalized
        ! to reference temperature : phinor = e phi / T.
        !
        !   Input parameters:
        !             formal:    x       -   array of coordinates
        !
        !   Output parameters:
        !             formal:    derphi  -   array of derivatives
        !
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: derphi(3)
        !
        derphi=0.d0
        !
    end subroutine elefie_can
    !
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine velo_can(tau,z,vz)
        !
        !
        !  Computes the components of the 5-D velocity          -  vz
        !  for given set of phase space coordinates             -  z.
        !
        !  Warning !!!  The dimensionless time is used (see below)
        !
        !  Order of the coordinates is the following:
        !  z(i) = x(i)  for i=1,2,3     - spatial coordinates with real
        !                                 dimension of the general covariant
        !                                 space coordinate system
        !  z(4) = p                     - momentum  module normalized to
        !                                 thermal momentum and sqrt(2);
        !  z(5) = alambd                - cosine of the pitch-angle
        !
        !  Input parameters:
        !            formal:  tau    -   dimensionless time: tau=sqrt(2*T/m)*t
        !                     z      -   see above
        !            common:  rmu    -   inverse relativistic temperature
        !                     ro0    -   Larmor radius for the reference
        !                                magnetic field and temperature:
        !                                ro0=sqrt(2*T/m)/(e*B_ref/m*c)
        !  Output parameters:
        !            formal:  vz     -   see above
        !
        !  Called routines: magfie, elefie_can
        !
        use parmot_mod, only : rmu,ro0
        use magfie_sub, only : magfie
        !
        implicit none
        !
        integer :: i
        !
        real(dp), intent(in) :: tau
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: vz(:)
        real(dp) x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
        real(dp) derphi
        real(dp) p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
        real(dp) rmumag,rovsqg,rosqgb,rovbm
        real(dp) a_phi,a_b,a_c,hstar
        real(dp) s_hc,hpstar,phidot,blodot,bra
        !
        dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
        dimension derphi(3)
        dimension a_phi(3),a_b(3),a_c(3),hstar(3)
        !
        do 1 i=1,3
            x(i)=z(i)
            1    continue
            !
            ! in magfie: x(i)   - set of 3 curvilinear space coordinates (input)
            !            bmod   - dimensionless magnetic field module: bmod=B/B_ref
            !            sqrtg  - Jacobian of space coordinates (square root of
            !                     metric tensor
            !            bder   - derivatives of logarithm of bmod over space coords
            !                     (covariant vector)
            !            hcovar - covariant components of the unit vector along
            !                     the magnetic field
            !            hctrvr - contravariant components of the unit vector along
            !                     the magnetic field
            !            hcurl  - contravariant components of the curl of this vector
            !
            call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
            !
            ! in elefie: x(i)   - space coords (input, see above)
            !            derphi - derivatives of the dimensionless electric potential
            !                     phihat=e*phi/T over space coords (covar. vector)
            !
            call elefie_can(x,derphi)
            !
            p=z(4)
            alambd=z(5)
            !
            p2=p*p
            ovmu=2.d0/rmu
            gamma2=p2*ovmu+1.d0
            gamma=dsqrt(gamma2)
            ppar=p*alambd
            ! vpa - dimensionless parallel velocity: vpa=v_parallel/sqrt(2*T/m)
            vpa=ppar/gamma
            coala=(1.d0-alambd**2)
            ! rmumag - magnetic moment
            rmumag=.5d0*p2*coala/bmod
            !
            rovsqg=ro0/sqrtg
            rosqgb=.5d0*rovsqg/bmod
            rovbm=ro0/bmod
            !
            a_phi(1)=(hcovar(2)*derphi(3)-hcovar(3)*derphi(2))*rosqgb
            a_b(1)=(hcovar(2)*bder(3)-hcovar(3)*bder(2))*rovsqg
            a_phi(2)=(hcovar(3)*derphi(1)-hcovar(1)*derphi(3))*rosqgb
            a_b(2)=(hcovar(3)*bder(1)-hcovar(1)*bder(3))*rovsqg
            a_phi(3)=(hcovar(1)*derphi(2)-hcovar(2)*derphi(1))*rosqgb
            a_b(3)=(hcovar(1)*bder(2)-hcovar(2)*bder(1))*rovsqg
            !
            s_hc=0.d0
            do i=1,3
                a_c(i)=hcurl(i)*rovbm
                s_hc=s_hc+a_c(i)*hcovar(i)
                hstar(i)=hctrvr(i)+ppar*a_c(i)
            enddo
            hpstar=1.d0+ppar*s_hc
            !
            ! velocities in the coordinate space
            !
            ! phidot - derivative of the dmls el. potential over dmls time
            ! blodot - derivative of the logarith of the mag. field module over dmls time
            phidot=0.d0
            blodot=0.d0
            do i=1,3
                bra=vpa*hstar(i)+a_phi(i)+a_b(i)*rmumag/gamma
                vz(i)=bra/hpstar
                phidot=phidot+vz(i)*derphi(i)
                blodot=blodot+vz(i)*bder(i)
            enddo
            !
            ! velocities in the phase space
            !
            vz(4)=-0.5d0*gamma*phidot/p
            !      vz(5)=coala/(alambd+dsign(1.d-32,alambd))*(vz(4)/p-0.5d0*blodot)
            vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p                 &
                + p*sum(hstar*bder)/gamma+alambd*sum(a_phi*bder))
            !
        end subroutine velo_can
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !
        subroutine orbit_timestep_can(z,dtau,dtaumin,relerr,ierr)
            use diag_mod, only : dodiag
            use odeint_allroutines_sub, only : odeint_allroutines
            use chamb_sub, only : chamb_can
            !
            implicit none
            !
            integer, parameter          :: ndim=5, nstepmax=1000000
            !
            integer :: ierr
            real(dp) :: dtau,dtaumin,phi,tau1,tau2
            !
            real(dp), dimension(2)    :: y
            real(dp), dimension(ndim) :: z
            real(dp) :: relerr
            !
            if(abs(dtaumin*nstepmax).le.abs(dtau)) then
                ierr=2
                print *,'orbit_timestep: number of steps exceeds nstepmax'
                return
            endif
            !
            ierr=0
            y(1)=z(1)
            y(2)=z(2)
            phi=z(3)
            !
            call chamb_can(y,phi,ierr)
            !
            if(ierr.eq.1) return
            tau1=0.d0
            tau2=dtaumin
            !
            do while((dtau>0 .and. (tau2 .lt. dtau)) .or. (dtau<0 .and. (tau2 .gt. dtau)))
                !
                call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
                !
                if(dodiag) write (123,*) tau2,z
                y(1)=z(1)
                y(2)=z(2)
                phi=z(3)
                !
                call chamb_can(y,phi,ierr)
                !
                if(ierr.eq.1) return
                tau1=tau2
                tau2=tau2+dtaumin
            enddo
            !
            tau2=dtau
            !
            call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
            !
            y(1)=z(1)
            y(2)=z(2)
            phi=z(3)
            !
            call chamb_can(y,phi,ierr)
            !
            if(ierr.eq.1) return
            !
        end subroutine orbit_timestep_can
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !
        subroutine rhs_mflint_can(phi,y,dery)
            !
            ! Computes the right hand side of the magnetic field line equation for
            ! the integration over the toroidal angle, subintegrand for the flux tube
            ! volume $1/B^\varphi$ and subintegrants for Boozer $B_{00}$ computation
            ! $B^2/B^\varphi$ and $B^3/B^\varphi$                       -   dery
            !
            ! Oder of coordinates in magfie is the following: x(1)=R (big radius),
            ! x(2)=phi (toroidal angle), x(3)=z (altitude).
            !
            !  Input parameters:
            !            formal:
            !                    y(1:2) - coordinates in the poloidal plane (phi=const)
            !                    y(3:5) - integrals
            !  Output parameters:
            !            formal:
            !                 dery(1:5) - vector of the right-hand side
            !  Called routines:  magfie_vmec, magfie_can, magfie_boozer
            !
            use magfie_sub, only : magfie
            !
            real(dp), intent(in) :: phi
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: dery(:)
            real(dp) x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
            dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
            !
            x(1)=y(1)
            x(2)=y(2)
            x(3)=phi
            !
            call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
            !
            dery(1)=hctrvr(1)/hctrvr(3)
            dery(2)=hctrvr(2)/hctrvr(3)
            dery(3)=1.d0/(bmod*hctrvr(3))
            dery(4)=bmod/hctrvr(3)
            dery(5)=bmod**2/hctrvr(3)
            !
        end subroutine rhs_mflint_can
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !Integrates along a magnetic field line to generate equidistant homogeneous volumetric distribution on flux surface
        subroutine integrate_mfl_can(npoi,dphi,rbeg,phibeg,zbeg,         &
                xstart,bstart,volstart,bmod00,ierr)
            !
            use odeint_allroutines_sub, only : odeint_allroutines
            use chamb_sub, only : chamb_can
            use magfie_sub, only : magfie
            !
            implicit none
            !
            real(dp), intent(out) :: bmod00
            !
            integer, parameter          :: ndim=5
            real(dp), parameter :: relerr=1d-10
            integer :: i
            integer, intent(in) :: npoi, ierr
            real(dp) :: phi,phiold
            real(dp) , intent(in):: dphi,rbeg,phibeg,zbeg
            real(dp), dimension(3,npoi), intent(out) :: xstart
            real(dp), dimension(npoi), intent(out)   :: bstart,volstart
            real(dp), dimension(ndim)   :: y
            !
            real(dp) x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
            dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
            !
            phi=phibeg
            y(1)=rbeg
            y(2)=zbeg
            y(3:5)=0.d0
            !
            call chamb_can(y(1:2),phi,ierr)
            !
            if(ierr.eq.1) return
            x(1)=y(1)
            x(2)=y(2)
            x(3)=phi
            !
            call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
            !
            xstart(:,1)=x
            bstart(1)=bmod
            volstart(1)=y(3)
            !
            do i=2,npoi
                phiold=phi
                phi=phiold+dphi
                !
                call odeint_allroutines(y,ndim,phiold,phi,relerr,rhs_mflint_can)
                !
                call chamb_can(y(1:2),phi,ierr)
                !
                if(ierr.eq.1) return
                x(1)=y(1)
                x(2)=y(2)
                x(3)=phi
                !
                call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
                !
                xstart(:,i)=x
                bstart(i)=bmod
                volstart(i)=y(3)
            enddo
            !
            volstart=volstart/volstart(npoi)
            bmod00=y(5)/y(4)
            !
        end subroutine integrate_mfl_can
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !
        subroutine velo_axis(tau,z_axis,vz_axis)
            real(dp), intent(in) :: tau
            real(dp), intent(in) :: z_axis(:)
            real(dp), intent(out) :: vz_axis(:)
            real(dp), parameter :: legacy_s_floor = 1.d-8

            call velo_axis_regularized(tau, z_axis, vz_axis, legacy_s_floor)
        end subroutine velo_axis

        subroutine velo_axis_regularized(tau,z_axis,vz_axis,s_floor)
            !
            ! Here variables z(1)=s and z(2)=theta are replaced with
            ! x_axis(1)=x=sqrt(s)*cos(theta) and x_axis(2)=y=sqrt(s)*sin(theta)
            !
            implicit none
            !
            real(dp), intent(in) :: tau
            real(dp), intent(in) :: z_axis(:)
            real(dp), intent(out) :: vz_axis(:)
            real(dp), intent(in) :: s_floor
            call evaluate_axis_rhs(tau, z_axis, s_floor, vz_axis)
        end subroutine velo_axis_regularized

        subroutine evaluate_axis_rhs(tau, z_axis, sample_s, vz_axis)
            real(dp), intent(in) :: tau, z_axis(:), sample_s
            real(dp), intent(out) :: vz_axis(:)
            integer, parameter :: nang = 8
            real(dp) :: values(axis_rhs_count, nang, 2)
            real(dp) :: mean(axis_rhs_count, 2)
            real(dp) :: m1c(axis_rhs_count, 2), m1s(axis_rhs_count, 2)
            real(dp) :: m2c(axis_rhs_count, 2), m2s(axis_rhs_count, 2)
            real(dp) :: coeff(axis_rhs_count), reconstructed(axis_rhs_count)
            real(dp) :: raw(axis_rhs_count), zsample(axis_rhs_count)
            real(dp) :: angle, radius, radius0, r2, blend, blend_arg
            integer :: iangle, iring

            r2=z_axis(1)*z_axis(1)+z_axis(2)*z_axis(2)
            if(r2.ge.sample_s) then
                call evaluate_axis_rhs_raw(tau, z_axis, vz_axis)
                return
            endif

            ! Reconstruct the physical guiding-centre vector field itself in
            ! regular Cartesian modes. Unlike independent tensor fits, this
            ! retains all algebraic identities among the magfie quantities.
            radius0=sqrt(sample_s)
            zsample(3:5)=z_axis(3:5)
            do iring=1,2
                radius=real(iring,dp)*radius0
                do iangle=1,nang
                    angle=2.d0*acos(-1.d0)*real(iangle-1,dp)/real(nang,dp)
                    zsample(1)=radius*cos(angle)
                    zsample(2)=radius*sin(angle)
                    call evaluate_axis_rhs_raw(tau, zsample, &
                        values(:,iangle,iring))
                enddo
            enddo

            mean=0.d0
            m1c=0.d0
            m1s=0.d0
            m2c=0.d0
            m2s=0.d0
            do iring=1,2
                do iangle=1,nang
                    angle=2.d0*acos(-1.d0)*real(iangle-1,dp)/real(nang,dp)
                    mean(:,iring)=mean(:,iring)+values(:,iangle,iring)
                    m1c(:,iring)=m1c(:,iring)+values(:,iangle,iring)*cos(angle)
                    m1s(:,iring)=m1s(:,iring)+values(:,iangle,iring)*sin(angle)
                    m2c(:,iring)=m2c(:,iring)+values(:,iangle,iring)*cos(2.d0*angle)
                    m2s(:,iring)=m2s(:,iring)+values(:,iangle,iring)*sin(2.d0*angle)
                enddo
            enddo
            mean=mean/real(nang,dp)
            m1c=2.d0*m1c/real(nang,dp)
            m1s=2.d0*m1s/real(nang,dp)
            m2c=2.d0*m2c/real(nang,dp)
            m2s=2.d0*m2s/real(nang,dp)

            reconstructed=(4.d0*mean(:,1)-mean(:,2))/3.d0
            coeff=(4.d0*m1c(:,1)/radius0 &
                -m1c(:,2)/(2.d0*radius0))/3.d0
            reconstructed=reconstructed+coeff*z_axis(1)
            coeff=(4.d0*m1s(:,1)/radius0 &
                -m1s(:,2)/(2.d0*radius0))/3.d0
            reconstructed=reconstructed+coeff*z_axis(2)
            coeff=(mean(:,2)-mean(:,1))/(3.d0*sample_s)
            reconstructed=reconstructed+coeff*r2
            coeff=(4.d0*m2c(:,1)/sample_s &
                -m2c(:,2)/(4.d0*sample_s))/3.d0
            reconstructed=reconstructed+coeff* &
                (z_axis(1)*z_axis(1)-z_axis(2)*z_axis(2))
            coeff=(4.d0*m2s(:,1)/sample_s &
                -m2s(:,2)/(4.d0*sample_s))/3.d0
            reconstructed=reconstructed+coeff* &
                (2.d0*z_axis(1)*z_axis(2))

            ! Join the reconstruction to the exact pullback with a C1
            ! smoothstep in the outer half of the disk.
            if(r2.gt.0.25d0*sample_s) then
                call evaluate_axis_rhs_raw(tau, z_axis, raw)
                blend_arg=2.d0*sqrt(r2/sample_s)-1.d0
                blend=blend_arg*blend_arg*(3.d0-2.d0*blend_arg)
                reconstructed=(1.d0-blend)*reconstructed+blend*raw
            endif
            vz_axis=reconstructed
        end subroutine evaluate_axis_rhs

        subroutine evaluate_axis_rhs_raw(tau, z_axis, vz_axis)
            real(dp), intent(in) :: tau, z_axis(:)
            real(dp), intent(out) :: vz_axis(:)
            real(dp) :: s, zpolar(5), vpolar(5)

            s=z_axis(1)*z_axis(1)+z_axis(2)*z_axis(2)
            if(.not.(s.gt.0.d0)) error stop &
                'raw axis-RHS evaluation requires nonzero radius'
            zpolar=[s,atan2(z_axis(2),z_axis(1)),z_axis(3:5)]
            call velo_can(tau,zpolar,vpolar)
            vz_axis(1)=0.5d0*vpolar(1)*z_axis(1)/s &
                -vpolar(2)*z_axis(2)
            vz_axis(2)=0.5d0*vpolar(1)*z_axis(2)/s &
                +vpolar(2)*z_axis(1)
            vz_axis(3:5)=vpolar(3:5)
        end subroutine evaluate_axis_rhs_raw
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !
        subroutine integrate_axis_segment(z,tau1,tau2,relerr,near_axis, &
                axis_floor,step_limit)
            use odeint_allroutines_sub, only : odeint_allroutines
            implicit none
            logical, intent(in) :: near_axis
            integer, intent(in), optional :: step_limit
            real(dp), intent(in) :: tau1,tau2,relerr,axis_floor
            real(dp), intent(inout) :: z(5)
            real(dp) :: context

            context=axis_floor
            if(present(step_limit)) then
                if(near_axis) then
                    call odeint_allroutines(z,5,context,tau1,tau2,relerr, &
                        velo_axis_context,step_limit=step_limit)
                else
                    call odeint_allroutines(z,5,context,tau1,tau2,relerr, &
                        velo_can_context,step_limit=step_limit)
                endif
            else if(near_axis) then
                call odeint_allroutines(z,5,context,tau1,tau2,relerr, &
                    velo_axis_context)
            else
                call odeint_allroutines(z,5,tau1,tau2,relerr,velo_can)
            endif
        end subroutine integrate_axis_segment

        subroutine velo_can_context(tau,z,vz,context)
            implicit none
            class(*), intent(in) :: context
            real(dp), intent(in) :: tau,z(:)
            real(dp), intent(out) :: vz(:)

            select type(context)
                type is(integer)
            end select
            call velo_can(tau,z,vz)
        end subroutine velo_can_context

        subroutine velo_axis_context(tau,z,vz,context)
            implicit none
            class(*), intent(in) :: context
            real(dp), intent(in) :: tau,z(:)
            real(dp), intent(out) :: vz(:)

            select type(context)
                type is(real(dp))
                call velo_axis_regularized(tau,z,vz,context)
            class default
                error stop 'invalid axis ODE context'
            end select
        end subroutine velo_axis_context

        subroutine orbit_timestep_axis(z,dtau,dtaumin,relerr,ierr,step_limit, &
                axis_floor,near_axis_limit)
            use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
            use odeint_allroutines_sub, only : odeint_has_failed
            use chamb_sub, only : chamb_can
            !
            implicit none
            !
            integer, parameter          :: ndim=5, nstepmax=1000000
            logical :: near_axis
            integer :: ierr
            integer, intent(in), optional :: step_limit
            real(dp), intent(in), optional :: axis_floor,near_axis_limit
            real(dp) :: dtau,dtaumin,phi,tau1,tau2,z1,z2,axis_floor_local
            real(dp) :: relerr,snear_axis
            !
            real(dp), dimension(2)    :: y
            real(dp), dimension(ndim) :: z, z_initial
            !
            if(abs(dtaumin*nstepmax).le.abs(dtau)) then
                ierr=2
                print *,'orbit_timestep: number of steps exceeds nstepmax'
                return
            endif
            !
            ierr=0
            axis_floor_local=1.d-8
            snear_axis=0.01d0
            if(present(axis_floor)) axis_floor_local=axis_floor
            if(present(near_axis_limit)) snear_axis=near_axis_limit
            if(.not.ieee_is_finite(axis_floor_local) .or. &
                axis_floor_local.le.0.d0 .or. &
                .not.ieee_is_finite(snear_axis) .or. &
                snear_axis.le.axis_floor_local .or. snear_axis.gt.1.d0) then
                ierr=2
                return
            endif
            z_initial=z
            y(1)=z(1)
            y(2)=z(2)
            phi=z(3)
            !
            call chamb_can(y,phi,ierr)
            !
            if(ierr.eq.1) return
            tau1=0.d0
            tau2=dtaumin
            !
            near_axis=.false.
            !
            do while((dtau>0 .and. (tau2 .lt. dtau)) .or. (dtau<0 .and. (tau2 .gt. dtau)))
                if(near_axis) then
                    if(z(1)**2+z(2)**2.gt.snear_axis) then
                        near_axis=.false.
                        z1=z(1)**2+z(2)**2
                        z2=atan2(z(2),z(1))
                        z(1)=z1
                        z(2)=z2
                        !
                        call integrate_axis_segment(z,tau1,tau2,relerr,.false., &
                            axis_floor_local,step_limit)
                        if (odeint_has_failed() .or. .not. all(ieee_is_finite(z))) then
                            z=z_initial
                            ierr=2
                            return
                        endif
                        !
                        y(1)=z(1)
                        y(2)=z(2)
                        phi=z(3)
                        !
                        call chamb_can(y,phi,ierr)
                        !
                        if(ierr.eq.1) return
                    else
                        !
                        call integrate_axis_segment(z,tau1,tau2,relerr,.true., &
                            axis_floor_local,step_limit)
                        if (odeint_has_failed() .or. .not. all(ieee_is_finite(z))) then
                            z=z_initial
                            ierr=2
                            return
                        endif
                        !
                    endif
                else
                    if(z(1).lt.snear_axis) then
                        near_axis=.true.
                        z1=sqrt(max(z(1),0.d0))*cos(z(2))
                        z2=sqrt(max(z(1),0.d0))*sin(z(2))
                        z(1)=z1
                        z(2)=z2
                        !
                        call integrate_axis_segment(z,tau1,tau2,relerr,.true., &
                            axis_floor_local,step_limit)
                        if (odeint_has_failed() .or. .not. all(ieee_is_finite(z))) then
                            z=z_initial
                            ierr=2
                            return
                        endif
                        !
                    else
                        !
                        call integrate_axis_segment(z,tau1,tau2,relerr,.false., &
                            axis_floor_local,step_limit)
                        if (odeint_has_failed() .or. .not. all(ieee_is_finite(z))) then
                            z=z_initial
                            ierr=2
                            return
                        endif
                        !
                        y(1)=z(1)
                        y(2)=z(2)
                        phi=z(3)
                        !
                        call chamb_can(y,phi,ierr)
                        !
                        if(ierr.eq.1) return
                        !
                    endif
                endif
                tau1=tau2
                tau2=tau2+dtaumin
            enddo
            !
            tau2=dtau
            !
            if(near_axis) then
                if(z(1)**2+z(2)**2.gt.snear_axis) then
                    near_axis=.false.
                    z1=z(1)**2+z(2)**2
                    z2=atan2(z(2),z(1))
                    z(1)=z1
                    z(2)=z2
                    !
                    call integrate_axis_segment(z,tau1,tau2,relerr,.false., &
                        axis_floor_local,step_limit)
                    if (odeint_has_failed() .or. .not. all(ieee_is_finite(z))) then
                        z=z_initial
                        ierr=2
                        return
                    endif
                    !
                    y(1)=z(1)
                    y(2)=z(2)
                    phi=z(3)
                    !
                    call chamb_can(y,phi,ierr)
                    !
                    if(ierr.eq.1) return
                else
                    !
                    call integrate_axis_segment(z,tau1,tau2,relerr,.true., &
                        axis_floor_local,step_limit)
                    if (odeint_has_failed() .or. .not. all(ieee_is_finite(z))) then
                        z=z_initial
                        ierr=2
                        return
                    endif
                    !
                endif
            else
                if(z(1).lt.snear_axis) then
                    near_axis=.true.
                    z1=sqrt(max(z(1),0.d0))*cos(z(2))
                    z2=sqrt(max(z(1),0.d0))*sin(z(2))
                    z(1)=z1
                    z(2)=z2
                    !
                    call integrate_axis_segment(z,tau1,tau2,relerr,.true., &
                        axis_floor_local,step_limit)
                    if (odeint_has_failed() .or. .not. all(ieee_is_finite(z))) then
                        z=z_initial
                        ierr=2
                        return
                    endif
                    !
                else
                    !
                    call integrate_axis_segment(z,tau1,tau2,relerr,.false., &
                        axis_floor_local,step_limit)
                    if (odeint_has_failed() .or. .not. all(ieee_is_finite(z))) then
                        z=z_initial
                        ierr=2
                        return
                    endif
                    !
                    y(1)=z(1)
                    y(2)=z(2)
                    phi=z(3)
                    !
                    call chamb_can(y,phi,ierr)
                    !
                    if(ierr.eq.1) return
                    !
                endif
            endif
            !
            if(near_axis) then
                z1=z(1)**2+z(2)**2
                z2=atan2(z(2),z(1))
                z(1)=z1
                z(2)=z2
            endif
            !
        end subroutine orbit_timestep_axis
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    end module alpha_lifetime_sub
