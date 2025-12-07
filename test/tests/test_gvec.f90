program test_gvec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_gvec, only: gvec_field_t, create_gvec_field

implicit none

class(gvec_field_t), allocatable :: gvec_field
character(len=256) :: test_file

! Test data file path - using relative path from build/test/tests/ directory
test_file = '../../../test/test_data/GVEC_elliptok_State_final.dat'

print *, 'Testing gvec_field_t creation...'

! Test: Create gvec_field_t using constructor
print *, 'Creating gvec_field_t with create_gvec_field...'
call create_gvec_field(test_file, gvec_field)

if (allocated(gvec_field)) then
    print *, 'SUCCESS: gvec_field_t created successfully'
    print *, 'Filename stored: ', trim(gvec_field%filename)
else
    print *, 'FAILED: gvec_field_t creation failed'
    error stop 1
end if

! Test the evaluate method
print *, ''
print *, 'Testing gvec_field_t evaluate method...'

block
    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    real(dp) :: s_test, theta_test, phi_test
    
    ! Test at a few different flux surfaces
    s_test = 0.5_dp        ! Half-radius
    theta_test = 0.0_dp    ! Poloidal angle
    phi_test = 0.0_dp      ! Toroidal angle
    
    x(1) = sqrt(s_test)    ! r = sqrt(s)
    x(2) = theta_test      ! theta
    x(3) = phi_test        ! phi
    
    call gvec_field%evaluate(x, Acov, hcov, Bmod)
    
    print *, 'Field evaluation at s=0.5:'
    print *, '  Acov = ', Acov
    print *, '  hcov = ', hcov  
    print *, '  Bmod = ', Bmod
    
    if (Bmod > 0.0_dp) then
        print *, 'SUCCESS: evaluate() method working'
    else
        print *, 'ERROR: evaluate() returned invalid Bmod'
        error stop 1
    end if
end block

print *, ''
print *, 'All tests passed!'
print *, 'gvec_field_t constructor and evaluate() working correctly.'

end program test_gvec
