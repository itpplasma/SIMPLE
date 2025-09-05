module diag_newton
!> Diagnostic routines for newton_midpoint function
!> Provides detailed analysis of Newton iteration convergence behavior
!> during orbit integration using midpoint rule

use, intrinsic :: iso_fortran_env, only: dp => real64
use util, only: pi, twopi
use field_can_mod, only: FieldCan, get_val, eval_field => evaluate
use orbit_symplectic_base, only: SymplecticIntegrator
use vector_potentail_mod, only: torflux
use lapack_interfaces, only: dgesv
use params, only: dtau
use orbit_symplectic, only: f_midpoint_part1, f_midpoint_part2, &
    jac_midpoint_part1, jac_midpoint_part2

implicit none
private

public :: newton_midpoint_debug, integrate_orbit_with_newton_debug

contains

!> EXACT copy of newton_midpoint from orbit_symplectic.f90 with debug output
subroutine newton_midpoint_debug(si, f, x, atol, rtol, maxit, xlast, step_num)
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  type(FieldCan) :: fmid
  integer, parameter :: n = 5
  integer :: kit
  real(dp), intent(inout) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(n)
  integer, intent(in) :: step_num
  real(dp) :: fvec(n), fjac(n,n)
  integer :: pivot(n), info
  real(dp) :: xabs(n), tolref(n), fabs(n)
  
  ! DEBUG OUTPUT
  print '(A,I0)', 'Newton Midpoint Debug - Step ', step_num
  print '(A,ES12.5,A,ES12.5)', 'Tolerances: atol = ', atol, ', rtol = ', rtol
  
  tolref(1) = 1d0
  tolref(2) = twopi
  tolref(3) = twopi
  tolref(4) = dabs(1d1*torflux/f%ro0)
  tolref(5) = 1d0
  
  print '(A,5ES12.5)', 'Initial Tolref = [', tolref, ']'
  print '(A,5ES12.5)', 'Initial x = [', x, ']'
  print *
  
  do kit = 1, maxit
    if(x(1) > 1.0) return
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
    ! after solution: fvec = (xold-xnew)_Newton
    x = x - fvec
    xabs = dabs(x - xlast)
    ! Don't take too small values in pphi as tolerance reference
    tolref(4) = max(dabs(x(4)), tolref(4))
    
    ! DEBUG OUTPUT
    print '(A,I0,A,5ES12.5)', 'Iteration ', kit, ': fabs = [', fabs, ']'
    print '(A,I0,A,5ES12.5)', 'Iteration ', kit, ': xabs = [', xabs, ']'
    print '(A,5ES12.5)', 'Tolref = [', tolref, ']'
    print '(A,5ES12.5)', 'rtol*tolref = [', rtol*tolref, ']'
    
    if (all(fabs < atol)) then
        print '(A,I0,A)', 'Iteration ', kit, ': Convergence achieved (fabs < atol)'
        return
    end if
    if (all(xabs < rtol*tolref)) then
        print '(A,I0,A)', 'Iteration ', kit, ': Convergence achieved (xabs < rtol*tolref)'
        return
    end if
    
    print '(A,I0,A,5ES12.5)', 'Iteration ', kit, ': Updated x = [', x, ']'
    print *
  enddo
  print '(A,I0)', 'newton_midpoint: maximum iterations reached: ', maxit
  !write(6603,*) x(1), x(2), x(3), x(4), x(5), xabs, fvec
  ! TODO fix criterion for convergence
end subroutine newton_midpoint_debug

!> Integration wrapper that calls debug newton_midpoint for specified number of steps
subroutine integrate_orbit_with_newton_debug(si, f, num_steps)
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    integer, intent(in) :: num_steps
    
    integer, parameter :: n = 5, maxit = 32
    real(dp), dimension(n) :: x, xlast
    integer :: step
    
    print *, 'Starting Newton Midpoint Integration Debug'
    print '(A,I0,A)', 'Will integrate for ', num_steps, ' time steps'
    print '(A,ES12.5)', 'dtau (large time step): ', dtau
    print '(A,ES12.5)', 'dtaumin (integration time step): ', si%dt
    print '(A,I0)', 'ntau (substeps per dtau): ', si%ntau
    print '(A,ES12.5)', 'Absolute tolerance: ', si%atol
    print '(A,ES12.5)', 'Relative tolerance: ', si%rtol
    print '(A,4ES12.5)', 'Initial conditions: ', si%z
    print *
    
    do step = 1, num_steps
        si%pthold = f%pth
        
        x(1:4) = si%z
        x(5) = si%z(1)
        
        call newton_midpoint_debug(si, f, x, si%atol, si%rtol, maxit, xlast, step)
        
        if (x(1) > 1.0_dp) then
            print *, 'Particle lost: s > 1.0 at step ', step
            exit
        end if
        
        si%z = x(1:4)
        
        ! Update field
        call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
        call get_val(f, si%z(4))
        
        print '(A,I0,A,4ES12.5)', 'Step ', step, ' completed. Final state: ', si%z
        print *
    end do
    
    print *, 'Newton Midpoint Integration Debug completed successfully!'
    
end subroutine integrate_orbit_with_newton_debug

end module diag_newton