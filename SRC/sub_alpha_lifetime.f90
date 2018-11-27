!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine elefie(x,derphi)
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
      integer :: ierr
      double precision, dimension(3) :: x,derphi
      double precision :: r,phi,z,psi,phi_el,phi_el_pr,phi_el_prpr
!
      derphi=0.d0
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine velo(tau,z,vz)
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
!                     p0     -   dimensionless momentum module in the
!                                initial point
!                     alamb0 -   cos of pitch-angle in the initial point
!  Output parameters:
!            formal:  vz     -   see above
!
!  Called routines: magfie, elefie
!
      use parmot_mod, only : rmu,ro0,eeff
! 05.03.2014
      use magfie_mod
! 05.03.2014 end
! 26.03.2016
      use gbpi_mod
! 26.03.2016 end
!
      implicit none
!
      integer :: i
!
      double precision tau,z,vz
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      double precision derphi
      double precision p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
      double precision rmumag,rovsqg,rosqgb,rovbm
      double precision a_phi,a_b,a_c,hstar
      double precision s_hc,hpstar,phidot,blodot,bra
      double precision pardeb
! 05.03.2014
! 08.03.2014      double precision vscon,vtcon,vfcon
      double precision vxcon,vzcon,vfcon
! 05.03.2014 end
!
      dimension z(5),vz(5)
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
! 26.03.2016
      if(ierrfield.ne.0) then
        vz=0.d0
        return
      endif
! 26.03.2016 end
!
! in elefie: x(i)   - space coords (input, see above)
!            derphi - derivatives of the dimensionless electric potential
!                     phihat=e*phi/T over space coords (covar. vector)
!
      call elefie(x,derphi)
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
! 05.03.2014
!       vscon=vz(1)*gsr+vz(2)*gsf+vz(3)*gsz
!       vtcon=vz(1)*gtr+vz(2)*gtf+vz(3)*gtz
!       vfcon=vz(1)*gfr+vz(2)*gff+vz(3)*gfz
!       vz(1)=vscon
!       vz(2)=vfcon
!       vz(3)=vtcon
! 05.03.2014 end
! 25.03.2016 (...fco)
       vxcon=vz(1)*gxcr+vz(2)*gxcfco+vz(3)*gxcz
       vzcon=vz(1)*gzcr+vz(2)*gzcfco+vz(3)*gzcz
       vfcon=vz(1)*gfr+vz(2)*gffco+vz(3)*gfz
       vz(1)=vxcon
       vz(2)=vfcon
       vz(3)=vzcon
! 23.03.2016 end
      vz(4)=-0.5d0*gamma*phidot/p
!      vz(5)=coala/(alambd+dsign(1.d-32,alambd))*(vz(4)/p-0.5d0*blodot)
      vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p                 &
            + p*sum(hstar*bder)/gamma+alambd*sum(a_phi*bder))
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine orbit_timestep(z,dtau,dtaumin,ierr)
!
! 26.03.2016
      use gbpi_mod
! 26.03.2016 end
      implicit none
!
      integer, parameter          :: ndim=5, nstepmax=1000000
      double precision, parameter :: relerr=1d-9
!
      integer :: ierr,j
      double precision :: dtau,dtaumin,phi,tau1,tau2
!
      double precision, dimension(2)    :: y
      double precision, dimension(ndim) :: z
! 10.02.2014 ! 24.08.2011
! 17,07.20	6      double precision, parameter :: vdr_dv=0.03d0
! 10.02.2014 ! 24.08.2011 end
!
      external velo
!
      if(dtaumin*nstepmax.le.dtau) then
        ierr=2
        print *,'orbit_timestep: number of steps exceeds nstepmax'
        return
      endif
!
      ierr=0
      y(1)=z(1)
      phi=z(2)
      y(2)=z(3)
!
      call chamb(y,phi,ierr)
!
      if(ierr.ne.0) return
      tau1=0.d0
      tau2=dtaumin
! 17.07.2016      tau2=dtaumin/(dabs(z(5))+vdr_dv)
!
      do while(tau2.lt.dtau)
!
        call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo)
!
        y(1)=z(1)
        phi=z(2)
        y(2)=z(3)
!
        call chamb(y,phi,ierr)
! 26.03.2016
      ierr=max(ierr,ierrfield)
! 26.03.2016 end
!
        if(ierr.ne.0) return
        tau1=tau2
        tau2=tau2+dtaumin
! 17.07.2016      tau2=tau2+dtaumin/(dabs(z(5))+vdr_dv)
      enddo
!
      tau2=dtau
!
      call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo)
!
      y(1)=z(1)
      phi=z(2)
      y(2)=z(3)
!
      call chamb(y,phi,ierr)
! 26.03.2016
      ierr=max(ierr,ierrfield)
! 26.03.2016 end
!
      if(ierr.ne.0) return
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhs_mflint(phi,y,dery)
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
!  Called routines:  magfie
!
! 05.03.2014
      use magfie_mod
! 05.03.2014 end
! 26.03.2016
      use gbpi_mod
! 26.03.2016 end
      double precision :: phi
      double precision, dimension(5) :: y,dery
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
!
      x(1)=y(1)
      x(2)=phi
      x(3)=y(2)
!
      call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
! 26.03.2016
      if(ierrfield.ne.0) then
        dery=0.d0
        return
      endif
! 26.03.2016 end
! 05.03.2014      dery(1)=hctrvr(1)/hctrvr(2)
! 05.03.2014      dery(2)=hctrvr(3)/hctrvr(2)
! 05.03.2014      dery(3)=1.d0/(bmod*hctrvr(2))
! 05.03.2014      dery(4)=bmod/hctrvr(2)
! 05.03.2014      dery(5)=bmod**2/hctrvr(2)
! 05.03.2014
! 08.03.2014      dery(1)=bscon/bfcon
! 08.03.2014      dery(2)=btcon/bfcon
      dery(1)=bxcon/bfcon
      dery(2)=bzcon/bfcon
      dery(3)=1.d0/bfcon
      dery(4)=bmod**2/bfcon
      dery(5)=bmod**3/bfcon     
! 09.03.2014
!      print*,bxcon,bzcon,bfcon,'=bxcon,bzcon,bfcon'
! 05.03.2014 end
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrate_mfl(npoi,dphi,rbeg,phibeg,zbeg, &
! 17.07.2016      subroutine integrate_mfl(npoi,dphii,rbeg,phibeg,zbeg,        &
                               xstart,bstart,volstart,bmod00,ierr)
!
! 26.03.2016
      use gbpi_mod
! 26.03.2016 end
      implicit none
!
      integer, parameter          :: ndim=5
      double precision, parameter :: relerr=1d-9
      integer :: npoi,i,ierr
      double precision :: dphi,rbeg,phibeg,zbeg,bmod00,phi,phiold
      double precision, dimension(3,npoi) :: xstart
      double precision, dimension(npoi)   :: bstart,volstart
      double precision, dimension(ndim)   :: y
!
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
! 10.03.2014 ! 01.09.2011
! 03.09.2011      integer, parameter          :: ndphi=10
! 08.09.2011      integer, parameter          :: ndphi=4
! 08.09.2011      integer, parameter          :: ndphi=3
! 10.09.2011      integer, parameter          :: ndphi=5
! 17.07.2016      integer, parameter          :: ndphi=13
! 17.07.2016      integer          :: i2
! 17.07.2016      double precision :: dphii
! 10.03.2014 ! 01.09.2011 end
!
      external :: rhs_mflint
! 10.03.2014 ! 01.09.2011
! 17.07.2016      dphi=dphii/ndphi
! 10.03.2014 ! 01.09.2011 end

!
      ierr=0
!
      phi=phibeg
      y(1)=rbeg
      y(2)=zbeg
      y(3:5)=0.d0
!
      call chamb(y(1:2),phi,ierr)
!
      if(ierr.ne.0) return
      x(1)=y(1)
      x(2)=phi
      x(3)=y(2)
!
! 08.03.2014
!        print*,'start integrate_mfl'
!        print*,x(1),x(2),x(3),'=x123'
! 08.03.2014 end
      call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
! 26.03.2016
      ierr=max(ierr,ierrfield)
      if(ierr.ne.0) return
! 26.03.2016 end
      xstart(:,1)=x
      bstart(1)=bmod
      volstart(1)=y(3)
!
      do i=2,npoi
! 10.03.2014 ! 01.09.2011
! 17.07.2016      do i2=1,ndphi
! 10.03.2014 ! 01.09.2011 end
        phiold=phi
        phi=phiold+dphi
!
        call odeint_allroutines(y,ndim,phiold,phi,relerr,rhs_mflint)
! 26.03.2016
        if(ierrfield.ne.0) exit
! 26.03.2016 end
! 10.03.2014 ! 01.09.2011
! 17.07.2016      enddo
! 10.03.2014 ! 01.09.2011 end
!
        call chamb(y(1:2),phi,ierr)
!
        if(ierr.ne.0) return
        x(1)=y(1)
        x(2)=phi
        x(3)=y(2)
!
! 08.03.2014
!        print*,i,'=i'
!        print*,x(1),x(2),x(3),'=x123'
! 08.03.2014 end
        call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
! 26.03.2016
        ierr=max(ierr,ierrfield)
        if(ierr.ne.0) return
! 26.03.2016 end
!
        xstart(:,i)=x
        bstart(i)=bmod
        volstart(i)=y(3)
      enddo
!
      volstart=volstart/volstart(npoi)
      bmod00=y(5)/y(4)
!
      end subroutine integrate_mfl
!
