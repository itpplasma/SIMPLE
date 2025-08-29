module diag_orbit
!> Diagnostic routines for plotting single orbit trajectories
!> Provides trajectory plotting functionality for the Nth particle using
!> full SIMPLE initialization and real orbit integration

use, intrinsic :: iso_fortran_env, only: dp => real64
use params, only: dtau, dtaumin, ntestpart, ntimstep, ntau, zstart, startmode, grid_density, &
    special_ants_file, reuse_batch, num_surf, sbeg, integmode, relerr, reset_seed_if_deterministic
use samplers, only: sample, START_FILE
use field_can_mod, only: FieldCan, get_val, eval_field => evaluate, ref_to_can
use orbit_symplectic_base, only: SymplecticIntegrator, extrap_field
use orbit_symplectic, only: orbit_timestep_sympl, f_midpoint_part1, f_midpoint_part2, &
    jac_midpoint_part1, jac_midpoint_part2
use simple, only: init_sympl
use vector_potentail_mod, only: torflux
use lapack_interfaces, only: dgesv
use util, only: twopi

implicit none
private

public :: integrate_orbit_with_trajectory_debug

contains

!> Newton midpoint solver that returns iteration count (no debug output)
function newton_midpoint_count_iterations(si, f, x, atol, rtol, maxit, xlast, field_evals) result(iterations)
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  type(FieldCan) :: fmid
  integer, parameter :: n = 5
  integer :: kit, iterations
  real(dp), intent(inout) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(n)
  integer, intent(inout) :: field_evals
  real(dp) :: fvec(n), fjac(n,n)
  integer :: pivot(n), info
  real(dp) :: xabs(n), tolref(n), fabs(n)
  
  ! Buffers to store all iteration data (only printed if max iterations reached)
  real(dp) :: x_buffer(n,maxit), fabs_buffer(n,maxit), xabs_buffer(n,maxit)
  real(dp) :: x_initial(n)
  integer :: k
  
  tolref(1) = 1d0
  tolref(2) = twopi
  tolref(3) = twopi
  tolref(4) = max(dabs(f%Aph), dabs(1d1*torflux/f%ro0))  ! Use actual Aph from field
  tolref(5) = 1d0
  
  ! Store initial conditions
  x_initial = x
  
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
    ! Each Newton iteration involves multiple field evaluations
    ! f_midpoint_part1 and f_midpoint_part2 each do field evaluations
    field_evals = field_evals + 2
    fabs = dabs(fvec)
    xlast = x
    call dgesv(n, 1, fjac, n, pivot, fvec, n, info)
    x = x - fvec
    xabs = dabs(x - xlast)
    tolref(4) = max(dabs(x(4)), tolref(4))
    
    ! Store iteration data in buffers
    x_buffer(:,kit) = x
    fabs_buffer(:,kit) = fabs
    xabs_buffer(:,kit) = xabs
    
    if (all(fabs < atol)) then
        iterations = kit
        return
    end if
    if (all(xabs < rtol*tolref)) then
        iterations = kit
        return
    end if
  enddo
  
  ! Maximum iterations reached - print complete iteration history
  write(*,'(A)') '=== NEWTON SOLVER FAILURE: MAXIMUM ITERATIONS REACHED ==='
  write(*,'(A,I0)') 'Maximum iterations: ', maxit
  write(*,'(A,5ES12.5)') 'Initial x = [', x_initial, ']'
  write(*,*)
  write(*,'(A)') 'Complete iteration history:'
  write(*,'(A)') 'Iter |    max(fabs)    |    max(xabs)    | Result'
  write(*,'(A)') '-----|----------------|----------------|-------'
  
  do k = 1, maxit
    write(*,'(I4,A,ES12.5,A,ES12.5,A)',advance='no') k, ' | ', maxval(fabs_buffer(:,k)), &
        ' | ', maxval(xabs_buffer(:,k)), ' | '
    
    if (all(fabs_buffer(:,k) < atol)) then
        write(*,'(A)') 'fabs < atol'
    elseif (all(xabs_buffer(:,k) < rtol*tolref)) then
        write(*,'(A)') 'xabs < rtol*tolref'
    else
        write(*,'(A)') 'continuing...'
    end if
  enddo
  
  write(*,*)
  write(*,'(A,5ES12.5)') 'Final fabs = [', fabs_buffer(:,maxit), ']'
  write(*,'(A,5ES12.5)') 'Final xabs = [', xabs_buffer(:,maxit), ']'
  write(*,'(A,5ES12.5)') 'rtol*tolref= [', rtol*tolref, ']'
  write(*,'(A,5ES12.5)') 'Final x    = [', x_buffer(:,maxit), ']'
  write(*,*)
  
  iterations = maxit
  error stop 'Newton solver failed to converge within maximum iterations'
end function newton_midpoint_count_iterations

!> Integration wrapper that plots the trajectory of the Nth particle
subroutine integrate_orbit_with_trajectory_debug(si, f, particle_number)
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    integer, intent(in) :: particle_number
    
    real(dp), allocatable :: s_traj(:), theta_traj(:), phi_traj(:), time_traj(:)
    real(dp), allocatable :: pphi_traj(:)
    integer, allocatable :: newton_iter_traj(:)
    real(dp), dimension(5) :: z, xlast
    integer :: it, ktau, point_idx, newton_iters, ierr_orbit
    integer, parameter :: maxit = 32
    real(dp) :: current_time
    integer :: total_points, field_eval_count
    
    ! Validate particle number
    if (particle_number < 1 .or. particle_number > ntestpart) then
        print '(A,I0,A,I0)', 'ERROR: Invalid particle number ', particle_number, &
            '. Must be between 1 and ', ntestpart
        return
    end if
    
    ! CRITICAL: Follow simple_main.f90 trace_orbit EXACTLY
    ! 1. Reset random seed if deterministic
    call reset_seed_if_deterministic
    
    ! 2. Get particle coordinates and transform ref -> can (CRITICAL STEP MISSING!)
    call ref_to_can(zstart(1:3, particle_number), z(1:3))
    z(4:5) = zstart(4:5, particle_number)
    
    ! 3. Initialize symplectic integrator with TRANSFORMED coordinates
    if (integmode > 0) then
        call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)
    end if
    
    current_time = 0.0_dp
    
    ! Calculate total number of timesteps (macrosteps * substeps + initial)
    total_points = ntimstep * ntau + 1
    field_eval_count = 0
    
    ! Allocate trajectory arrays
    allocate(s_traj(total_points))
    allocate(theta_traj(total_points))
    allocate(phi_traj(total_points))
    allocate(pphi_traj(total_points))
    allocate(newton_iter_traj(total_points))
    allocate(time_traj(total_points))
    
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
    field_eval_count = field_eval_count + 1
    
    point_idx = 1
    
    ! Initialize xlast for field extrapolation to current coordinates
    xlast(1:4) = si%z
    xlast(5) = si%z(1)
    
    ! Use our custom Newton solver to get iteration counts (but with proper initialization now)
    do it = 1, ntimstep
        do ktau = 1, ntau
            si%pthold = f%pth
            
            ! Set up for midpoint integration (like diag_newton)
            z(1:4) = si%z
            z(5) = si%z(1)
            
            ! Use custom Newton midpoint solver to get iteration count
            newton_iters = newton_midpoint_count_iterations(si, f, z, si%atol, si%rtol, maxit, xlast, field_eval_count)
            
            current_time = current_time + dtaumin
            
            if (z(1) > 1.0_dp) then
                exit
            end if
            
            ! Update integrator state
            si%z = z(1:4)
            
            ! Store trajectory point
            point_idx = point_idx + 1
            s_traj(point_idx) = si%z(1)
            theta_traj(point_idx) = si%z(2)
            phi_traj(point_idx) = si%z(3)
            pphi_traj(point_idx) = si%z(4)
            newton_iter_traj(point_idx) = newton_iters
            time_traj(point_idx) = current_time
            
            ! Update field with extrapolation like production integrator
            if (extrap_field) then
                f%pth = f%pth + f%dpth(1)*(z(1)-xlast(1) + z(5) - xlast(5)) &  ! d/dr
                              + f%dpth(2)*(z(2)-xlast(2)) &  ! d/dth
                              + f%dpth(3)*(z(3)-xlast(3)) &  ! d/dph
                              + f%dpth(4)*(z(4)-xlast(4))    ! d/dpph
            else
                call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
                call get_val(f, si%z(4))
                field_eval_count = field_eval_count + 1
            endif
        end do
        if (z(1) > 1.0_dp) exit
    end do
    
    ! Write trajectory data to files for external plotting
    call write_trajectory_data(time_traj(1:point_idx), s_traj(1:point_idx), &
        theta_traj(1:point_idx), phi_traj(1:point_idx), pphi_traj(1:point_idx), &
        newton_iter_traj(1:point_idx), point_idx, particle_number)
    
    ! Output field evaluation count
    print '(A,I0)', 'Total field evaluations: ', field_eval_count
    
    ! Cleanup
    deallocate(s_traj, theta_traj, phi_traj, pphi_traj, newton_iter_traj, time_traj)
    
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
    
end subroutine write_trajectory_data

end module diag_orbit