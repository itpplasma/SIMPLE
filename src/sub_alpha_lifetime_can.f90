module alpha_lifetime_sub

implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)

! Physical and numerical constants
integer, parameter :: ndim = 5
integer, parameter :: nstepmax = 1000000
real(dp), parameter :: snear_axis = 0.01d0
real(dp), parameter :: relerr_default = 1.0e-10_dp
real(dp), parameter :: s_min = 1.0e-8_dp

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
      integer :: ierr
      real(dp), dimension(3) :: x,derphi
      real(dp) :: r,phi,z,psi,phi_el,phi_el_pr,phi_el_prpr
  !
      derphi=0.d0
  !
      end subroutine elefie_can
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
      subroutine velo_can(tau,z,vz)
  !
  !  Computes the components of the 5-D velocity for given phase space coordinates
  !
      use parmot_mod, only : rmu,ro0
      use magfie_sub, only : magfie
  !
      implicit none
  !
      real(dp), intent(in) :: tau
      real(dp), dimension(5), intent(in) :: z
      real(dp), dimension(5), intent(out) :: vz
      
      real(dp), dimension(3) :: x, bder, hcovar, hctrvr, hcurl, derphi
      real(dp), dimension(3) :: a_phi, a_b, a_c, hstar
      real(dp) :: bmod, sqrtg
      real(dp) :: p, alambd, gamma, ppar, vpa, coala, rmumag
      real(dp) :: s_hc, hpstar, phidot, blodot
      
      ! Extract spatial coordinates
      x = z(1:3)
      
      ! Get magnetic field properties
      call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
      
      ! Get electric field
      call elefie_can(x,derphi)
      
      ! Extract momentum coordinates
      p = z(4)
      alambd = z(5)
      
      ! Compute particle properties
      call compute_particle_properties(p, alambd, rmu, bmod, &
                                       gamma, ppar, vpa, coala, rmumag)
      
      ! Compute drift velocities
      call compute_drift_velocities(hcovar, hctrvr, hcurl, derphi, bder, &
                                    ro0, sqrtg, bmod, ppar, &
                                    a_phi, a_b, a_c, hstar, s_hc, hpstar)
      
      ! Compute spatial velocities
      call compute_spatial_velocities(vpa, hstar, a_phi, a_b, rmumag, gamma, &
                                      hpstar, derphi, bder, vz(1:3), phidot, blodot)
      
      ! Compute phase space velocities
      call compute_phase_velocities(gamma, p, alambd, coala, hpstar, phidot, &
                                    hstar, derphi, bder, a_phi, vz(4:5))
  !
      end subroutine velo_can
      
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute particle properties
      subroutine compute_particle_properties(p, alambd, rmu, bmod, &
                                             gamma, ppar, vpa, coala, rmumag)
      implicit none
      
      real(dp), intent(in) :: p, alambd, rmu, bmod
      real(dp), intent(out) :: gamma, ppar, vpa, coala, rmumag
      
      real(dp) :: p2, ovmu, gamma2
      
      p2 = p * p
      ovmu = 2.0d0 / rmu
      gamma2 = p2 * ovmu + 1.0d0
      gamma = sqrt(gamma2)
      ppar = p * alambd
      vpa = ppar / gamma  ! dimensionless parallel velocity
      coala = 1.0d0 - alambd**2
      rmumag = 0.5d0 * p2 * coala / bmod  ! magnetic moment
      
      end subroutine compute_particle_properties
      
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute drift velocities
      subroutine compute_drift_velocities(hcovar, hctrvr, hcurl, derphi, bder, &
                                          ro0, sqrtg, bmod, ppar, &
                                          a_phi, a_b, a_c, hstar, s_hc, hpstar)
      implicit none
      
      real(dp), dimension(3), intent(in) :: hcovar, hctrvr, hcurl, derphi, bder
      real(dp), intent(in) :: ro0, sqrtg, bmod, ppar
      real(dp), dimension(3), intent(out) :: a_phi, a_b, a_c, hstar
      real(dp), intent(out) :: s_hc, hpstar
      
      real(dp) :: rovsqg, rosqgb, rovbm
      integer :: i
      
      rovsqg = ro0 / sqrtg
      rosqgb = 0.5d0 * rovsqg / bmod
      rovbm = ro0 / bmod
      
      ! Compute cross products for drift terms
      call compute_cross_products(hcovar, derphi, rosqgb, a_phi)
      call compute_cross_products(hcovar, bder, rovsqg, a_b)
      
      ! Compute curl-related terms
      s_hc = 0.0d0
      do i = 1, 3
        a_c(i) = hcurl(i) * rovbm
        s_hc = s_hc + a_c(i) * hcovar(i)
        hstar(i) = hctrvr(i) + ppar * a_c(i)
      end do
      hpstar = 1.0d0 + ppar * s_hc
      
      end subroutine compute_drift_velocities
      
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute cross products
      subroutine compute_cross_products(vec1, vec2, scale, result)
      implicit none
      
      real(dp), dimension(3), intent(in) :: vec1, vec2
      real(dp), intent(in) :: scale
      real(dp), dimension(3), intent(out) :: result
      
      result(1) = (vec1(2)*vec2(3) - vec1(3)*vec2(2)) * scale
      result(2) = (vec1(3)*vec2(1) - vec1(1)*vec2(3)) * scale
      result(3) = (vec1(1)*vec2(2) - vec1(2)*vec2(1)) * scale
      
      end subroutine compute_cross_products
      
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute spatial velocities
      subroutine compute_spatial_velocities(vpa, hstar, a_phi, a_b, rmumag, gamma, &
                                            hpstar, derphi, bder, vz, phidot, blodot)
      implicit none
      
      real(dp), intent(in) :: vpa, rmumag, gamma, hpstar
      real(dp), dimension(3), intent(in) :: hstar, a_phi, a_b, derphi, bder
      real(dp), dimension(3), intent(out) :: vz
      real(dp), intent(out) :: phidot, blodot
      
      real(dp) :: bra
      integer :: i
      
      phidot = 0.0d0
      blodot = 0.0d0
      
      do i = 1, 3
        bra = vpa*hstar(i) + a_phi(i) + a_b(i)*rmumag/gamma
        vz(i) = bra / hpstar
        phidot = phidot + vz(i)*derphi(i)
        blodot = blodot + vz(i)*bder(i)
      end do
      
      end subroutine compute_spatial_velocities
      
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute phase space velocities
      subroutine compute_phase_velocities(gamma, p, alambd, coala, hpstar, phidot, &
                                          hstar, derphi, bder, a_phi, vz_phase)
      implicit none
      
      real(dp), intent(in) :: gamma, p, alambd, coala, hpstar, phidot
      real(dp), dimension(3), intent(in) :: hstar, derphi, bder, a_phi
      real(dp), dimension(2), intent(out) :: vz_phase
      
      vz_phase(1) = -0.5d0 * gamma * phidot / p
      vz_phase(2) = -(0.5d0*coala/hpstar) * (sum(hstar*derphi)/p &
                    + p*sum(hstar*bder)/gamma + alambd*sum(a_phi*bder))
      
      end subroutine compute_phase_velocities
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
      subroutine orbit_timestep_can(z,dtau,dtaumin,relerr,ierr)
use diag_mod, only : dodiag
use odeint_sub, only : odeint_allroutines
use chamb_sub, only : chamb_can
  !
      implicit none
  !
      integer, parameter          :: ndim=5, nstepmax=1000000
  !
      integer :: ierr,j
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
      if(ierr.ne.0) return
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
        if(ierr.ne.0) return
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
      if(ierr.ne.0) return
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
      real(dp) :: phi
      real(dp), dimension(5) :: y,dery
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
      use odeint_sub, only : odeint_allroutines
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
      if(ierr.ne.0) return
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
        if(ierr.ne.0) return
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
  !
  ! Velocity computation near magnetic axis using transformed coordinates
  !
  implicit none
  !
  real(dp), intent(in) :: tau
  real(dp), dimension(5), intent(in) :: z_axis
  real(dp), dimension(5), intent(out) :: vz_axis
  
  real(dp), dimension(5) :: z, vz
  
  ! Transform from axis coordinates to standard coordinates
  call axis_to_standard_coords(z_axis, z)
  
  ! Compute velocity in standard coordinates
  call velo_can(tau, z, vz)
  
  ! Transform velocity back to axis coordinates
  call velocity_to_axis_coords(z_axis, z, vz, vz_axis)
  !
  end subroutine velo_axis
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Transform axis coordinates to standard coordinates
  subroutine axis_to_standard_coords(z_axis, z)
  implicit none
  
  real(dp), dimension(5), intent(in) :: z_axis
  real(dp), dimension(5), intent(out) :: z
  
  z(1) = sqrt(z_axis(1)**2 + z_axis(2)**2)
  z(1) = max(z(1), s_min)
  z(2) = atan2(z_axis(2), z_axis(1))
  z(3:5) = z_axis(3:5)
  
  end subroutine axis_to_standard_coords
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Transform velocity to axis coordinates
  subroutine velocity_to_axis_coords(z_axis, z, vz, vz_axis)
  implicit none
  
  real(dp), dimension(5), intent(in) :: z_axis, z, vz
  real(dp), dimension(5), intent(out) :: vz_axis
  
  real(dp) :: derlogsqs
  
  derlogsqs = vz(1) / z(1)
  vz_axis(1) = derlogsqs*z_axis(1) - vz(2)*z_axis(2)
  vz_axis(2) = derlogsqs*z_axis(2) + vz(2)*z_axis(1)
  vz_axis(3:5) = vz(3:5)
  
  end subroutine velocity_to_axis_coords
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
      subroutine orbit_timestep_axis(z,dtau,dtaumin,relerr,ierr)
      use odeint_sub, only : odeint_allroutines
      use chamb_sub, only : chamb_can
  !
      implicit none
  !
      integer, parameter          :: ndim=5, nstepmax=1000000
      real(dp), parameter :: snear_axis=0.01d0
  !
      logical :: near_axis
      integer :: ierr,j
      real(dp) :: dtau,dtaumin,phi,tau1,tau2,z1,z2
      real(dp) :: relerr
  !
      real(dp), dimension(2)    :: y
      real(dp), dimension(ndim) :: z
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
      if(ierr.ne.0) return
      tau1=0.d0
      tau2=dtaumin
  !
      near_axis=.false.
  !
      do while((dtau>0 .and. (tau2 .lt. dtau)) .or. (dtau<0 .and. (tau2 .gt. dtau)))
        if(near_axis) then
          if(z(1)**2+z(2)**2.gt.snear_axis**2) then
            near_axis=.false.
            z1=sqrt(z(1)**2+z(2)**2)
            z2=atan2(z(2),z(1))
            z(1)=z1
            z(2)=z2
  !
            call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
  !
            y(1)=z(1)
            y(2)=z(2)
            phi=z(3)
  !
            call chamb_can(y,phi,ierr)
  !
            if(ierr.ne.0) return
          else
  !
            call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_axis)
  !
          endif
        else
          if(z(1).lt.snear_axis) then
            near_axis=.true.
            z1=z(1)*cos(z(2))
            z2=z(1)*sin(z(2))
            z(1)=z1
            z(2)=z2
  !
            call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_axis)
  !
          else
  !
            call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
  !
            y(1)=z(1)
            y(2)=z(2)
            phi=z(3)
  !
            call chamb_can(y,phi,ierr)
  !
            if(ierr.ne.0) return
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
        if(z(1)**2+z(2)**2.gt.snear_axis**2) then
          near_axis=.false.
          z1=sqrt(z(1)**2+z(2)**2)
          z2=atan2(z(2),z(1))
          z(1)=z1
          z(2)=z2
  !
          call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
  !
          y(1)=z(1)
          y(2)=z(2)
          phi=z(3)
  !
          call chamb_can(y,phi,ierr)
  !
          if(ierr.ne.0) return
        else
  !
          call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_axis)
  !
        endif
      else
        if(z(1).lt.snear_axis) then
          near_axis=.true.
          z1=z(1)*cos(z(2))
          z2=z(1)*sin(z(2))
          z(1)=z1
          z(2)=z2
  !
          call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_axis)
  !
        else
  !
          call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
  !
          y(1)=z(1)
          y(2)=z(2)
          phi=z(3)
  !
          call chamb_can(y,phi,ierr)
  !
          if(ierr.ne.0) return
  !
        endif
      endif
  !
      if(near_axis) then
        z1=sqrt(z(1)**2+z(2)**2)
        z2=atan2(z(2),z(1))
        z(1)=z1
        z(2)=z2
      endif
  !
      end subroutine orbit_timestep_axis
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module alpha_lifetime_sub
