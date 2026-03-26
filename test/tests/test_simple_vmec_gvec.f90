program test_simple_vmec_gvec
   use, intrinsic :: iso_fortran_env, only: dp => real64

   implicit none

   real(dp), parameter :: logical_start_points(5, 4) = reshape([ &
                                                               0.35_dp, 0.2_dp, 0.1_dp, 1.0_dp, 0.2_dp, &
                                                               0.42_dp, 1.0_dp, 0.4_dp, 1.0_dp, -0.1_dp, &
                                                               0.58_dp, 2.2_dp, 0.7_dp, 1.0_dp, 0.4_dp, &
                                                               0.71_dp, 4.1_dp, 1.0_dp, 1.0_dp, -0.3_dp &
                                                               ], [5, 4])
   real(dp) :: vmec_confined, gvec_confined, boozer_confined
   real(dp) :: rel_diff_gvec
   real(dp) :: rel_diff_boozer
   real(dp) :: rel_diff_internal

   call write_logical_start_file()
   call write_vmec_config()
   call run_case('test_vmec_refcoords.in', 'test_vmec_refcoords.log', vmec_confined)

   call write_logical_start_file()
   call write_gvec_config()
   call run_case('test_gvec_refcoords.in', 'test_gvec_refcoords.log', gvec_confined)

   call copy_file('start_boozer.dat', 'start.dat')
   call write_gvec_boozer_config()
   call run_case('test_gvec_boozer_refcoords.in', 'test_gvec_boozer_refcoords.log', &
                 boozer_confined)

   print '(A,F8.4)', 'vmec confined fraction = ', vmec_confined
   print '(A,F8.4)', 'gvec confined fraction = ', gvec_confined
   print '(A,F8.4)', 'boozer confined fraction = ', boozer_confined

   rel_diff_gvec = abs(vmec_confined - gvec_confined)/ &
                   max(1.0e-12_dp, abs(vmec_confined))
   rel_diff_boozer = abs(vmec_confined - boozer_confined)/ &
                     max(1.0e-12_dp, abs(vmec_confined))
   rel_diff_internal = abs(gvec_confined - boozer_confined)/ &
                       max(1.0e-12_dp, abs(gvec_confined))

   print '(A,ES12.4)', 'relative diff logical-vmec = ', rel_diff_gvec
   print '(A,ES12.4)', 'relative diff boozer-vmec  = ', rel_diff_boozer
   print '(A,ES12.4)', 'relative diff boozer-gvec  = ', rel_diff_internal

   if (rel_diff_gvec > 5.0e-2_dp) error stop 'test_simple_vmec_gvec: logical mismatch'
   if (rel_diff_boozer > 5.0e-2_dp) error stop 'test_simple_vmec_gvec: Boozer mismatch'
   if (rel_diff_internal > 5.0e-2_dp) &
      error stop 'test_simple_vmec_gvec: GVEC coordinate-family mismatch'

contains

   subroutine write_logical_start_file()
      integer :: unit

      open (newunit=unit, file='start.dat', status='replace')
      write (unit, '(5ES24.16)') logical_start_points(:, 1)
      write (unit, '(5ES24.16)') logical_start_points(:, 2)
      write (unit, '(5ES24.16)') logical_start_points(:, 3)
      write (unit, '(5ES24.16)') logical_start_points(:, 4)
      close (unit)
   end subroutine write_logical_start_file

   subroutine write_vmec_config()
      integer :: unit

      open (newunit=unit, file='test_vmec_refcoords.in', status='replace')
      write (unit, '(A)') '&config'
      write (unit, '(A)') 'coord_input = "wout.nc"'
      write (unit, '(A)') 'field_input = "wout.nc"'
      write (unit, '(A)') 'isw_field_type = 5'
      write (unit, '(A)') 'integmode = 0'
      write (unit, '(A)') 'startmode = 2'
      write (unit, '(A)') 'trace_time = 1.0d-5'
      write (unit, '(A)') 'ntestpart = 4'
      write (unit, '(A)') 'output_results_netcdf = .false.'
      write (unit, '(A)') '/'
      close (unit)
   end subroutine write_vmec_config

   subroutine write_gvec_config()
      integer :: unit

      open (newunit=unit, file='test_gvec_refcoords.in', status='replace')
      write (unit, '(A)') '&config'
      write (unit, '(A)') 'coord_input = "wout.gvec_export.nc"'
      write (unit, '(A)') 'field_input = "wout.gvec_export.nc"'
      write (unit, '(A)') 'isw_field_type = 5'
      write (unit, '(A)') 'integmode = 0'
      write (unit, '(A)') 'startmode = 2'
      write (unit, '(A)') 'trace_time = 1.0d-5'
      write (unit, '(A)') 'ntestpart = 4'
      write (unit, '(A)') 'output_results_netcdf = .false.'
      write (unit, '(A)') '/'
      close (unit)
   end subroutine write_gvec_config

   subroutine write_gvec_boozer_config()
      integer :: unit

      open (newunit=unit, file='test_gvec_boozer_refcoords.in', status='replace')
      write (unit, '(A)') '&config'
      write (unit, '(A)') 'coord_input = "wout.gvec_boozer_export.nc"'
      write (unit, '(A)') 'field_input = "wout.gvec_boozer_export.nc"'
      write (unit, '(A)') 'isw_field_type = 5'
      write (unit, '(A)') 'integmode = 0'
      write (unit, '(A)') 'startmode = 2'
      write (unit, '(A)') 'trace_time = 1.0d-5'
      write (unit, '(A)') 'ntestpart = 4'
      write (unit, '(A)') 'output_results_netcdf = .false.'
      write (unit, '(A)') '/'
      close (unit)
   end subroutine write_gvec_boozer_config

   subroutine run_case(input_file, log_file, confined_fraction)
      character(*), intent(in) :: input_file
      character(*), intent(in) :: log_file
      real(dp), intent(out) :: confined_fraction

      integer :: unit
      integer :: iostat
      character(len=1024) :: command

      write (command, '(A)') '../../simple.x '//trim(input_file)//' > '//trim(log_file)
      call execute_command_line(trim(command), exitstat=iostat)
      if (iostat /= 0) error stop 'test_simple_vmec_gvec: SIMPLE execution failed'

      open (newunit=unit, file='confined_fraction.dat', status='old', iostat=iostat)
      if (iostat /= 0) error stop 'test_simple_vmec_gvec: missing confined_fraction.dat'
      read (unit, *)
      do
         read (unit, *, iostat=iostat) confined_fraction
         if (iostat /= 0) exit
      end do
      close (unit)
   end subroutine run_case

   subroutine copy_file(source, destination)
      character(*), intent(in) :: source
      character(*), intent(in) :: destination

      integer :: in_unit
      integer :: out_unit
      integer :: iostat
      character(len=256) :: line

      open (newunit=in_unit, file=source, status='old', iostat=iostat)
      if (iostat /= 0) error stop 'test_simple_vmec_gvec: missing source start file'
      open (newunit=out_unit, file=destination, status='replace', iostat=iostat)
      if (iostat /= 0) error stop 'test_simple_vmec_gvec: failed to open destination'
      do
         read (in_unit, '(A)', iostat=iostat) line
         if (iostat /= 0) exit
         write (out_unit, '(A)') trim(line)
      end do
      close (in_unit)
      close (out_unit)
   end subroutine copy_file

end program test_simple_vmec_gvec
