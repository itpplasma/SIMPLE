!> Module to replace splint_vmec_data with libneo interpolate
module spline_vmec_libneo
    use interpolate, only: SplineData3D, SplineData1D, &
                          construct_splines_3d, construct_splines_1d, &
                          evaluate_splines_3d_der, evaluate_splines_1d_der
    
    implicit none
    private
    
    integer, parameter :: dp = kind(1.0d0)
    
    ! Spline data structures for VMEC quantities
    type(SplineData3D), save :: spl_R, spl_Z, spl_lam
    type(SplineData1D), save :: spl_A_phi
    logical, save :: splines_initialized = .false.
    
    public :: init_vmec_splines_libneo, splint_vmec_data_libneo, cleanup_vmec_splines_libneo
    
contains
    
    !> Initialize spline structures from VMEC data
    subroutine init_vmec_splines_libneo()
        use new_vmec_stuff_mod, only: sR, sZ, slam, ns_s, ns_tp, n_theta, n_phi, h_theta, h_phi, nper, ns_A
        use vector_potentail_mod, only: sA_phi, ns, hs
        
        real(dp), parameter :: twopi = 2.0_dp * 3.14159265358979_dp
        real(dp) :: x_min(3), x_max(3)
        integer :: order(3), npoints(3)
        logical :: periodic(3)
        integer :: is, ith, iph, k1, k2, k3
        real(dp), allocatable :: data_3d(:,:,:)
        real(dp), allocatable :: data_1d(:)
        real(dp) :: x_min_1d, x_max_1d
        integer :: ns_actual
        
        ! Setup for 3D splines (R, Z, lambda)
        ns_actual = size(sR, 4)  ! Get actual radial grid size
        
        ! Grid parameters
        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max = [hs * (ns_actual-1), twopi, twopi/dble(nper)]
        order = [ns_s, ns_tp, ns_tp]
        npoints = [ns_actual, n_theta, n_phi]
        periodic = [.false., .true., .true.]
        
        ! Allocate temporary array for data
        allocate(data_3d(ns_actual, n_theta, n_phi))
        
        ! R spline - need to restructure data from sR(k1,k2,k3,is,ith,iph) to data_3d(is,ith,iph)
        ! We need the function values at grid points, which are the k1=0 coefficients
        do iph = 1, n_phi
            do ith = 1, n_theta
                do is = 1, ns_actual
                    data_3d(is, ith, iph) = sR(1, 1, 1, is, ith, iph)  ! k1=k2=k3=0 coefficients
                end do
            end do
        end do
        call construct_splines_3d(x_min, x_max, data_3d, order, periodic, spl_R)
        
        ! Z spline
        do iph = 1, n_phi
            do ith = 1, n_theta
                do is = 1, ns_actual
                    data_3d(is, ith, iph) = sZ(1, 1, 1, is, ith, iph)
                end do
            end do
        end do
        call construct_splines_3d(x_min, x_max, data_3d, order, periodic, spl_Z)
        
        ! Lambda spline
        do iph = 1, n_phi
            do ith = 1, n_theta
                do is = 1, ns_actual
                    data_3d(is, ith, iph) = slam(1, 1, 1, is, ith, iph)
                end do
            end do
        end do
        call construct_splines_3d(x_min, x_max, data_3d, order, periodic, spl_lam)
        
        deallocate(data_3d)
        
        ! 1D spline for A_phi (only if allocated)
        if (allocated(sA_phi)) then
            x_min_1d = 0.0_dp
            x_max_1d = hs * (ns-1)
            allocate(data_1d(ns))
            
            ! Extract function values (k=0 coefficients)
            do is = 1, ns
                data_1d(is) = sA_phi(1, is)
            end do
            
            call construct_splines_1d(x_min_1d, x_max_1d, data_1d, ns_A-1, .false., spl_A_phi)
            
            deallocate(data_1d)
        end if
        
        splines_initialized = .true.
        
    end subroutine init_vmec_splines_libneo
    
    
    !> Replacement for splint_vmec_data using libneo interpolate
    subroutine splint_vmec_data_libneo(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                      R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
        use vector_potentail_mod, only: torflux, hs
        use new_vmec_stuff_mod, only: nper, sR
        use spline_vmec_sub, only: splint_vmec_data
        
        real(dp), parameter :: twopi = 2.0_dp * 3.14159265358979_dp
        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
        real(dp), intent(out) :: R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
        
        real(dp) :: x_3d(3), dy_3d(3)
        real(dp) :: rho_tor, s_actual
        real(dp) :: dy_1d
        
        ! Check if VMEC data is loaded first
        if (.not. allocated(sR)) then
            ! VMEC data not loaded, call the original splint_vmec_data which will initialize it
            call splint_vmec_data(s, theta, varphi, &
                                 A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                 R, Z, alam, &
                                 dR_ds, dR_dt, dR_dp, &
                                 dZ_ds, dZ_dt, dZ_dp, &
                                 dl_ds, dl_dt, dl_dp)
            return
        end if
        
        ! Initialize splines if not done yet
        if (.not. splines_initialized) then
            call init_vmec_splines_libneo()
        end if
        
        ! A_theta and its derivative
        A_theta = torflux * s
        dA_theta_ds = torflux
        
        ! Evaluate A_phi using 1D spline
        s_actual = s * hs  ! Convert to actual coordinate
        call evaluate_splines_1d_der(spl_A_phi, s_actual, A_phi, dy_1d)
        dA_phi_ds = dy_1d / hs  ! Convert derivative back
        
        ! For R, Z, lambda: use sqrt(s) as radial coordinate
        rho_tor = sqrt(s)
        
        ! Setup 3D coordinates (rho_tor, theta, varphi)
        x_3d(1) = rho_tor * hs
        x_3d(2) = modulo(theta, twopi)
        x_3d(3) = modulo(varphi, twopi/dble(nper))
        
        ! Evaluate R and derivatives
        call evaluate_splines_3d_der(spl_R, x_3d, R, dy_3d)
        ! Convert derivatives: d/ds = d/d(rho) * d(rho)/ds = d/d(rho) / (2*sqrt(s))
        dR_ds = dy_3d(1) / (2.0_dp * rho_tor * hs)
        dR_dt = dy_3d(2)
        dR_dp = dy_3d(3)
        
        ! Evaluate Z and derivatives
        call evaluate_splines_3d_der(spl_Z, x_3d, Z, dy_3d)
        dZ_ds = dy_3d(1) / (2.0_dp * rho_tor * hs)
        dZ_dt = dy_3d(2)
        dZ_dp = dy_3d(3)
        
        ! Evaluate lambda and derivatives
        call evaluate_splines_3d_der(spl_lam, x_3d, alam, dy_3d)
        dl_ds = dy_3d(1) / (2.0_dp * rho_tor * hs)
        dl_dt = dy_3d(2)
        dl_dp = dy_3d(3)
        
        ! Rotational transform
        aiota = -dA_phi_ds / dA_theta_ds
        
    end subroutine splint_vmec_data_libneo
    
    
    !> Clean up spline structures
    subroutine cleanup_vmec_splines_libneo()
        if (allocated(spl_R%coeff)) deallocate(spl_R%coeff)
        if (allocated(spl_Z%coeff)) deallocate(spl_Z%coeff)
        if (allocated(spl_lam%coeff)) deallocate(spl_lam%coeff)
        if (allocated(spl_A_phi%coeff)) deallocate(spl_A_phi%coeff)
        splines_initialized = .false.
    end subroutine cleanup_vmec_splines_libneo
    
end module spline_vmec_libneo