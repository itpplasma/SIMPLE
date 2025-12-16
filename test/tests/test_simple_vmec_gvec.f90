program test_simple_vmec_gvec
    !> Comprehensive test running SIMPLE with both VMEC and GVEC fields
    !> Verifies particle confinement statistics are consistent

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use params, only: pi

    implicit none

    character(len=256) :: vmec_file, gvec_file
    character(len=256) :: namelist_file
    logical :: file_exists
    real(dp) :: confined_fraction_vmec, confined_fraction_gvec
    real(dp) :: relative_difference
    integer :: iostat, unit

    print *, '======================================================='
    print *, 'Testing SIMPLE with both VMEC and GVEC fields'
    print *, '======================================================='
    print *, ''

    ! Check for required files
    vmec_file = 'wout.nc'
    inquire(file=vmec_file, exist=file_exists)
    if (.not. file_exists) then
        print *, 'FAIL: No VMEC file (wout.nc) found in working directory'
        error stop 1
    end if

    ! Look for existing GVEC file created by test_vmec_gvec
    gvec_file = 'test_vmec_gvec_State_0000_00000000.dat'
    inquire(file=gvec_file, exist=file_exists)
    if (.not. file_exists) then
        print *, 'FAIL: No GVEC test file found. Run test_vmec_gvec first.'
        error stop 1
    end if

    ! Test 1: Run with VMEC field (canonical coordinates)
    print *, '1. Running SIMPLE with VMEC field (canonical coordinates)...'

    ! Create namelist for VMEC test
    namelist_file = 'test_vmec.in'
    open(newunit=unit, file=namelist_file, status='replace')
    write(unit, '(A)') '&config'
    write(unit, '(A)') 'multharm = 3            ! low order splines'
    write(unit, '(A)') 'trace_time = 1d-4       ! very short test run'
    write(unit, '(A)') 'sbeg = 0.5d0            ! mid-radius'
    write(unit, '(A)') 'ntestpart = 16          ! minimal for quick test'
    write(unit, '(A)') 'contr_pp  = -1.0d10     ! trace also passing'
    write(unit, '(A)') '/'
    close(unit)

    ! Run SIMPLE with VMEC
    call system('../../simple.x test_vmec.in > test_vmec.log')

    ! Read confined fraction
    open(newunit=unit, file='confined_fraction.dat', status='old', iostat=iostat)
    if (iostat /= 0) then
        print *, 'ERROR: Could not read confined_fraction.dat for VMEC test'
        error stop 1
    end if

    ! Skip header and read last line
    read(unit, *)  ! skip header
    do
        read(unit, *, iostat=iostat) confined_fraction_vmec
        if (iostat /= 0) exit
    end do
    close(unit)

    print '(A,F8.4)', '   Final confined fraction: ', confined_fraction_vmec
    call system('mv confined_fraction.dat confined_fraction_vmec.dat')
    call system('mv times_lost.dat times_lost_vmec.dat')

    ! Test 2: Run with GVEC field (canonical coordinates)
    print *, ''
    print *, '2. Running SIMPLE with GVEC field (canonical coordinates)...'

    ! Create namelist for GVEC test
    namelist_file = 'test_gvec.in'
    open(newunit=unit, file=namelist_file, status='replace')
    write(unit, '(A)') '&config'
    write(unit, '(A)') 'trace_time = 1d-4       ! very short test run'
    write(unit, '(A)') 'sbeg = 0.5d0            ! mid-radius'
    write(unit, '(A)') 'ntestpart = 16          ! minimal for quick test'
    write(unit, '(A)') 'field_input = "test_vmec_gvec_State_0000_00000000.dat"  ! GVEC file'
    write(unit, '(A)') 'isw_field_type = 0      ! canonical coordinates'
    write(unit, '(A)') 'startmode = 0           ! volume sampling'
    write(unit, '(A)') 'integmode = 4           ! symplectic integrator'
    write(unit, '(A)') '/'
    close(unit)

    ! Run SIMPLE with GVEC
    error stop "Requires fully working GVEC field implementation"
    call system('../../simple.x test_gvec.in > test_gvec.log')

    ! Read confined fraction
    open(newunit=unit, file='confined_fraction.dat', status='old', iostat=iostat)
    if (iostat /= 0) then
        print *, 'ERROR: Could not read confined_fraction.dat for GVEC test'
        error stop 1
    end if

    ! Skip header and read last line
    read(unit, *)  ! skip header
    do
        read(unit, *, iostat=iostat) confined_fraction_gvec
        if (iostat /= 0) exit
    end do
    close(unit)

    print '(A,F8.4)', '   Final confined fraction: ', confined_fraction_gvec
    call system('mv confined_fraction.dat confined_fraction_gvec.dat')
    call system('mv times_lost.dat times_lost_gvec.dat')

    ! Compare results
    print *, ''
    print *, '======================================================='
    print *, 'Results Comparison:'
    print *, '======================================================='
    print '(A,F8.4)', 'VMEC confined fraction: ', confined_fraction_vmec
    print '(A,F8.4)', 'GVEC confined fraction: ', confined_fraction_gvec

    relative_difference = abs(confined_fraction_vmec - confined_fraction_gvec) / &
                         (0.5_dp * (confined_fraction_vmec + confined_fraction_gvec))

    print '(A,F8.4,A)', 'Relative difference: ', relative_difference * 100.0_dp, '%'

    ! Clean up
    call system('rm -f test_vmec.in test_gvec.in')
    call system('rm -f test_vmec.log test_gvec.log')
    call system('rm -f start.dat fort.6601')

    ! Check if results are consistent (allow 5% difference due to interpolation)
    if (relative_difference < 0.05_dp) then
        print *, ''
        print *, 'TEST PASSED: VMEC and GVEC fields give consistent results'
    else
        print *, ''
        print *, 'TEST FAILED: Results differ by more than 5%'
        error stop 1
    end if

end program test_simple_vmec_gvec
