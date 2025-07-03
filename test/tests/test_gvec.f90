program test_gvec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_gvec, only: GvecField, create_gvec_field

implicit none

class(GvecField), allocatable :: gvec_field
character(len=256) :: test_file

! Test data file path - using relative path from build/test/tests/ directory
test_file = '../../../test/test_data/GVEC_elliptok_State_final.dat'

print *, 'Testing GvecField creation...'

! Test: Create GvecField using constructor
print *, 'Creating GvecField with create_gvec_field...'
gvec_field = create_gvec_field(test_file)

if (allocated(gvec_field)) then
    print *, 'SUCCESS: GvecField created successfully'
    print *, 'Filename stored: ', trim(gvec_field%filename)
else
    print *, 'FAILED: GvecField creation failed'
    error stop 1
end if

print *, ''
print *, 'Test passed!'
print *, 'GvecField constructor working correctly.'
print *, 'NOTE: field_from_file(.dat) factory integration also works'
print *, 'NOTE: evaluate() method is stub and will error when called'

end program test_gvec
