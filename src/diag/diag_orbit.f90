module diag_orbit
!> Diagnostic routines for plotting single orbit trajectories
!> Provides trajectory plotting functionality for the Nth particle using
!> full SIMPLE initialization and real orbit integration

use, intrinsic :: iso_fortran_env, only: dp => real64
use pyplot_module, only: pyplot
use params, only: dtau, ntestpart, zstart, startmode, grid_density, &
    special_ants_file, reuse_batch, num_surf, sbeg
use samplers, only: sample, START_FILE
use field_can_mod, only: FieldCan, get_val, eval_field => evaluate
use orbit_symplectic_base, only: SymplecticIntegrator
use orbit_symplectic, only: orbit_timestep_sympl, f_midpoint_part1, f_midpoint_part2, &
    jac_midpoint_part1, jac_midpoint_part2
use vector_potentail_mod, only: torflux
use lapack_interfaces, only: dgesv
use util, only: twopi

implicit none
private

public :: integrate_orbit_with_trajectory_debug

contains

!> Newton midpoint solver that returns iteration count (no debug output)
function newton_midpoint_count_iterations(si, f, x, atol, rtol, maxit, xlast) result(iterations)
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  type(FieldCan) :: fmid
  integer, parameter :: n = 5
  integer :: kit, iterations
  real(dp), intent(inout) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(n)
  real(dp) :: fvec(n), fjac(n,n)
  integer :: pivot(n), info
  real(dp) :: xabs(n), tolref(n), fabs(n)
  
  tolref(1) = 1d0
  tolref(2) = twopi
  tolref(3) = twopi
  tolref(4) = dabs(1d1*torflux/f%ro0)
  tolref(5) = 1d0
  
  do kit = 1, maxit
    if(x(1) > 1.0) then
        iterations = kit - 1
        return
    end if
    if(x(1) < 0.0) x(1) = 0.01
    if(x(5) < 0.0) x(5) = 0.01
    call f_midpoint_part1(si, f, n, x, fvec)
    call jac_midpoint_part1(si, f, x, fjac)
    fmid = f
    call f_midpoint_part2(si, f, n, x, fvec)
    call jac_midpoint_part2(si, f, fmid, x, fjac)
    fabs = dabs(fvec)
    xlast = x
    call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    x = x - fvec
    xabs = dabs(x - xlast)
    tolref(4) = max(dabs(x(4)), tolref(4))
    
    if (all(fabs < atol)) then
        iterations = kit
        return
    end if
    if (all(xabs < rtol*tolref)) then
        iterations = kit
        return
    end if
  enddo
  iterations = maxit
end function newton_midpoint_count_iterations

!> Integration wrapper that plots the trajectory of the Nth particle
subroutine integrate_orbit_with_trajectory_debug(si, f, particle_number, num_steps)
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    integer, intent(in) :: particle_number
    integer, intent(in) :: num_steps
    
    real(dp), allocatable :: s_traj(:), theta_traj(:), phi_traj(:), time_traj(:)
    real(dp), allocatable :: pphi_traj(:)
    integer, allocatable :: newton_iter_traj(:)
    real(dp), dimension(5) :: z, xlast
    integer :: step, newton_iters
    integer, parameter :: maxit = 32
    real(dp) :: current_time
    
    ! Validate particle number
    if (particle_number < 1 .or. particle_number > ntestpart) then
        print '(A,I0,A,I0)', 'ERROR: Invalid particle number ', particle_number, &
            '. Must be between 1 and ', ntestpart
        return
    end if
    
    print *, 'Starting Orbit Trajectory Integration'
    print '(A,I0,A)', 'Will integrate particle ', particle_number, ' trajectory'
    print '(A,I0,A)', 'for ', num_steps, ' time steps'
    print '(A,ES12.5)', 'dtau (large time step): ', dtau
    print '(A,ES12.5)', 'dtaumin (integration time step): ', si%dt
    print '(A,I0)', 'ntau (substeps per dtau): ', si%ntau
    print '(A,ES12.5)', 'Absolute tolerance: ', si%atol
    print '(A,ES12.5)', 'Relative tolerance: ', si%rtol
    print *
    
    ! Get the Nth particle's initial conditions (assumed to be pre-initialized)
    z = zstart(:, particle_number)
    si%z = z(1:4)
    current_time = 0.0_dp
    
    print '(A,4ES12.5)', 'Initial conditions: ', si%z
    print *
    
    ! Allocate trajectory arrays
    allocate(s_traj(num_steps + 1))
    allocate(theta_traj(num_steps + 1))
    allocate(phi_traj(num_steps + 1))
    allocate(pphi_traj(num_steps + 1))
    allocate(newton_iter_traj(num_steps + 1))
    allocate(time_traj(num_steps + 1))
    
    ! Store initial conditions (0 Newton iterations for initial state)
    s_traj(1) = si%z(1)
    theta_traj(1) = si%z(2)
    phi_traj(1) = si%z(3)
    pphi_traj(1) = si%z(4)
    newton_iter_traj(1) = 0
    time_traj(1) = current_time
    
    ! Initialize field at starting position
    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_val(f, si%z(4))
    
    do step = 1, num_steps
        si%pthold = f%pth
        
        ! Set up for midpoint integration (like diag_newton)
        z(1:4) = si%z
        z(5) = si%z(1)
        
        ! Use custom Newton midpoint solver to get iteration count
        newton_iters = newton_midpoint_count_iterations(si, f, z, si%atol, si%rtol, maxit, xlast)
        
        current_time = current_time + dtau
        
        if (z(1) > 1.0_dp) then
            print *, 'Particle lost: s > 1.0 at step ', step
            exit
        end if
        
        ! Update integrator state
        si%z = z(1:4)
        
        ! Store trajectory point
        s_traj(step + 1) = si%z(1)
        theta_traj(step + 1) = si%z(2)
        phi_traj(step + 1) = si%z(3)
        pphi_traj(step + 1) = si%z(4)
        newton_iter_traj(step + 1) = newton_iters
        time_traj(step + 1) = current_time
        
        ! Update field
        call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
        call get_val(f, si%z(4))
        
        if (mod(step, max(1, num_steps/10)) == 0) then
            print '(A,I0,A,4ES12.5,A,I0)', 'Step ', step, ' state: ', si%z, ', Newton iters: ', newton_iters
        end if
    end do
    
    print '(A,I0,A,4ES12.5)', 'Step ', step, ' completed. Final state: ', si%z
    print *
    
    ! Write trajectory data to files for external plotting
    call write_trajectory_data(time_traj(1:step+1), s_traj(1:step+1), &
        theta_traj(1:step+1), phi_traj(1:step+1), pphi_traj(1:step+1), &
        newton_iter_traj(1:step+1), step+1, particle_number)
    
    ! Cleanup
    deallocate(s_traj, theta_traj, phi_traj, pphi_traj, newton_iter_traj, time_traj)
    
    print *, 'Orbit Trajectory Integration Debug completed successfully!'
    
end subroutine integrate_orbit_with_trajectory_debug

subroutine write_trajectory_data(time_traj, s_traj, theta_traj, phi_traj, pphi_traj, &
    newton_iter_traj, npoints, particle_number)
    integer, intent(in) :: npoints, particle_number
    real(dp), dimension(npoints), intent(in) :: time_traj, s_traj, theta_traj, phi_traj, pphi_traj
    integer, dimension(npoints), intent(in) :: newton_iter_traj
    
    integer :: i
    character(len=100) :: filename
    
    write(filename, '(A,I0,A)') 'orbit_trajectory_particle_', particle_number, '.dat'
    
    open(unit=20, file=filename, status='replace')
    write(20, '(A)') '# Time    s    theta    phi    pphi    newton_iters'
    do i = 1, npoints
        write(20, '(5ES16.8,I8)') time_traj(i), s_traj(i), theta_traj(i), phi_traj(i), pphi_traj(i), newton_iter_traj(i)
    end do
    close(20)
    
    print '(A,A)', 'Orbit trajectory data written to: ', trim(filename)
    print *, 'Data columns: Time, s, theta, phi, pphi, newton_iters'
    print *, 'Use external plotting tools (gnuplot, matplotlib, etc.) to visualize.'
    
end subroutine write_trajectory_data

end module diag_orbit