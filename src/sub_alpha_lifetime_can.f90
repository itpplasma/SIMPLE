module alpha_lifetime_sub

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    integer, parameter :: axis_tensor_count = 14

    type :: magfie_data_t
        real(dp) :: bmod, sqrtg
        real(dp) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
    end type magfie_data_t

contains

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
            real(dp) :: derphi(3)
            type(magfie_data_t) :: field

            ! Pull back the complete magfie tensor set before evaluating the
            ! coordinate-covariant guiding-centre equation. Transforming only
            ! the final polar velocity makes a radius floor select an arbitrary
            ! approach ray at the axis.
            call evaluate_axis_field(z_axis(1), z_axis(2), z_axis(3), &
                s_floor, field)
            derphi=0.d0
            call guiding_center_velocity(z_axis, field, derphi, vz_axis)
        end subroutine velo_axis_regularized

        subroutine guiding_center_velocity(z, field, derphi, vz)
            use parmot_mod, only : rmu,ro0
            real(dp), intent(in) :: z(:), derphi(3)
            type(magfie_data_t), intent(in) :: field
            real(dp), intent(out) :: vz(:)
            integer :: i
            real(dp) :: p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
            real(dp) :: rmumag,rovsqg,rosqgb,rovbm
            real(dp) :: a_phi(3),a_b(3),a_c(3),hstar(3)
            real(dp) :: s_hc,hpstar,phidot,blodot,bra

            p=z(4)
            alambd=z(5)
            p2=p*p
            ovmu=2.d0/rmu
            gamma2=p2*ovmu+1.d0
            gamma=dsqrt(gamma2)
            ppar=p*alambd
            vpa=ppar/gamma
            coala=(1.d0-alambd**2)
            rmumag=.5d0*p2*coala/field%bmod

            rovsqg=ro0/field%sqrtg
            rosqgb=.5d0*rovsqg/field%bmod
            rovbm=ro0/field%bmod

            a_phi(1)=(field%hcovar(2)*derphi(3) &
                -field%hcovar(3)*derphi(2))*rosqgb
            a_b(1)=(field%hcovar(2)*field%bder(3) &
                -field%hcovar(3)*field%bder(2))*rovsqg
            a_phi(2)=(field%hcovar(3)*derphi(1) &
                -field%hcovar(1)*derphi(3))*rosqgb
            a_b(2)=(field%hcovar(3)*field%bder(1) &
                -field%hcovar(1)*field%bder(3))*rovsqg
            a_phi(3)=(field%hcovar(1)*derphi(2) &
                -field%hcovar(2)*derphi(1))*rosqgb
            a_b(3)=(field%hcovar(1)*field%bder(2) &
                -field%hcovar(2)*field%bder(1))*rovsqg

            s_hc=0.d0
            do i=1,3
                a_c(i)=field%hcurl(i)*rovbm
                s_hc=s_hc+a_c(i)*field%hcovar(i)
                hstar(i)=field%hctrvr(i)+ppar*a_c(i)
            enddo
            hpstar=1.d0+ppar*s_hc

            phidot=0.d0
            blodot=0.d0
            do i=1,3
                bra=vpa*hstar(i)+a_phi(i)+a_b(i)*rmumag/gamma
                vz(i)=bra/hpstar
                phidot=phidot+vz(i)*derphi(i)
                blodot=blodot+vz(i)*field%bder(i)
            enddo

            vz(4)=-0.5d0*gamma*phidot/p
            vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p &
                +p*sum(hstar*field%bder)/gamma &
                +alambd*sum(a_phi*field%bder))
        end subroutine guiding_center_velocity

        subroutine evaluate_axis_field(x, y, phi, sample_s, field)
            real(dp), intent(in) :: x, y, phi, sample_s
            type(magfie_data_t), intent(out) :: field
            integer, parameter :: nang = 8
            real(dp) :: values(axis_tensor_count, nang, 2)
            real(dp) :: mean(axis_tensor_count, 2)
            real(dp) :: m1c(axis_tensor_count, 2), m1s(axis_tensor_count, 2)
            real(dp) :: m2c(axis_tensor_count, 2), m2s(axis_tensor_count, 2)
            real(dp) :: coeff(axis_tensor_count), packed(axis_tensor_count)
            real(dp) :: raw_packed(axis_tensor_count)
            real(dp) :: angle, radius, radius0, r2, blend, blend_arg
            integer :: iangle, iring
            type(magfie_data_t) :: sample

            r2=x*x+y*y
            if(r2.ge.sample_s) then
                call evaluate_axis_field_raw(x, y, phi, field)
                return
            endif

            ! Resolve the regular Cartesian modes through quadratic order on
            ! two circles. The second radius removes the leading radial
            ! contamination from the axis value and linear harmonics.
            radius0=sqrt(sample_s)
            do iring=1,2
                radius=real(iring,dp)*radius0
                do iangle=1,nang
                    angle=2.d0*acos(-1.d0)*real(iangle-1,dp)/real(nang,dp)
                    call evaluate_axis_field_raw(radius*cos(angle), &
                        radius*sin(angle), phi, sample)
                    call pack_axis_field(sample, values(:,iangle,iring))
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

            packed=(4.d0*mean(:,1)-mean(:,2))/3.d0
            coeff=(4.d0*m1c(:,1)/radius0 &
                -m1c(:,2)/(2.d0*radius0))/3.d0
            packed=packed+coeff*x
            coeff=(4.d0*m1s(:,1)/radius0 &
                -m1s(:,2)/(2.d0*radius0))/3.d0
            packed=packed+coeff*y
            coeff=(mean(:,2)-mean(:,1))/(3.d0*sample_s)
            packed=packed+coeff*r2
            coeff=(4.d0*m2c(:,1)/sample_s &
                -m2c(:,2)/(4.d0*sample_s))/3.d0
            packed=packed+coeff*(x*x-y*y)
            coeff=(4.d0*m2s(:,1)/sample_s &
                -m2s(:,2)/(4.d0*sample_s))/3.d0
            packed=packed+coeff*(2.d0*x*y)

            ! Join the reconstruction to the exact pullback with a C1
            ! smoothstep in the outer half of the disk.
            if(r2.gt.0.25d0*sample_s) then
                call evaluate_axis_field_raw(x, y, phi, sample)
                call pack_axis_field(sample, raw_packed)
                blend_arg=2.d0*sqrt(r2/sample_s)-1.d0
                blend=blend_arg*blend_arg*(3.d0-2.d0*blend_arg)
                packed=(1.d0-blend)*packed+blend*raw_packed
            endif
            call unpack_axis_field(packed, field)
        end subroutine evaluate_axis_field

        subroutine evaluate_axis_field_raw(x, y, phi, field)
            use magfie_sub, only : magfie
            real(dp), intent(in) :: x, y, phi
            type(magfie_data_t), intent(out) :: field
            real(dp) :: s, theta, xpolar(3)
            real(dp) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

            s=x*x+y*y
            if(.not.(s.gt.0.d0)) error stop &
                'raw axis-field evaluation requires nonzero radius'
            theta=atan2(y,x)
            xpolar=[s,theta,phi]
            call magfie(xpolar, field%bmod, field%sqrtg, bder, hcovar, &
                hctrvr, hcurl)

            ! ds = 2X dX + 2Y dY and
            ! dtheta = (-Y dX + X dY)/s; det(dx/dq) = 2 > 0.
            field%sqrtg=2.d0*field%sqrtg
            call polar_covector_to_axis(x, y, s, bder, field%bder)
            call polar_covector_to_axis(x, y, s, hcovar, field%hcovar)
            call polar_vector_to_axis(x, y, s, hctrvr, field%hctrvr)
            call polar_vector_to_axis(x, y, s, hcurl, field%hcurl)
        end subroutine evaluate_axis_field_raw

        pure subroutine polar_covector_to_axis(x, y, s, polar, axis)
            real(dp), intent(in) :: x, y, s, polar(3)
            real(dp), intent(out) :: axis(3)

            axis(1)=2.d0*x*polar(1)-y*polar(2)/s
            axis(2)=2.d0*y*polar(1)+x*polar(2)/s
            axis(3)=polar(3)
        end subroutine polar_covector_to_axis

        pure subroutine polar_vector_to_axis(x, y, s, polar, axis)
            real(dp), intent(in) :: x, y, s, polar(3)
            real(dp), intent(out) :: axis(3)

            axis(1)=x*polar(1)/(2.d0*s)-y*polar(2)
            axis(2)=y*polar(1)/(2.d0*s)+x*polar(2)
            axis(3)=polar(3)
        end subroutine polar_vector_to_axis

        pure subroutine pack_axis_field(field, packed)
            type(magfie_data_t), intent(in) :: field
            real(dp), intent(out) :: packed(axis_tensor_count)

            packed(1)=field%bmod
            packed(2)=field%sqrtg
            packed(3:5)=field%bder
            packed(6:8)=field%hcovar
            packed(9:11)=field%hctrvr
            packed(12:14)=field%hcurl
        end subroutine pack_axis_field

        pure subroutine unpack_axis_field(packed, field)
            real(dp), intent(in) :: packed(axis_tensor_count)
            type(magfie_data_t), intent(out) :: field

            field%bmod=packed(1)
            field%sqrtg=packed(2)
            field%bder=packed(3:5)
            field%hcovar=packed(6:8)
            field%hctrvr=packed(9:11)
            field%hcurl=packed(12:14)
        end subroutine unpack_axis_field
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
                axis_floor)
            use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
            use odeint_allroutines_sub, only : odeint_has_failed
            use chamb_sub, only : chamb_can
            !
            implicit none
            !
            integer, parameter          :: ndim=5, nstepmax=1000000
            real(dp), parameter :: snear_axis=0.01d0
            !
            logical :: near_axis
            integer :: ierr
            integer, intent(in), optional :: step_limit
            real(dp), intent(in), optional :: axis_floor
            real(dp) :: dtau,dtaumin,phi,tau1,tau2,z1,z2,axis_floor_local
            real(dp) :: relerr
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
            if(present(axis_floor)) axis_floor_local=axis_floor
            if(.not.ieee_is_finite(axis_floor_local) .or. &
                axis_floor_local.le.0.d0) then
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
