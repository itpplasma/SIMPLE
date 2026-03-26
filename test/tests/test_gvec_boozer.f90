program test_gvec_boozer
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use field_gvec, only: gvec_field_t, create_gvec_field
   use field_vmec, only: vmec_field_t, create_vmec_field
   use new_vmec_stuff_mod, only: netcdffile, multharm
   use spline_vmec_sub, only: spline_vmec_data

   implicit none

   class(gvec_field_t), allocatable :: boozer_field
   class(vmec_field_t), allocatable :: vmec_field
   real(dp), parameter :: pos_tol = 5.0e-4_dp
   real(dp), parameter :: h_tol = 5.0e-2_dp
   real(dp), parameter :: bmod_tol = 4.0e-2_dp

   real(dp) :: x_boozer(3)
   real(dp) :: x_vmec(3)
   real(dp) :: xcart_boozer(3)
   real(dp) :: xcart_vmec(3)
   real(dp) :: xcyl(3)
   real(dp) :: Acov_boozer(3)
   real(dp) :: Acov_vmec(3)
   real(dp) :: hcov_boozer(3)
   real(dp) :: hcov_vmec(3)
   real(dp) :: h_cart_boozer(3)
   real(dp) :: h_cart_vmec(3)
   real(dp) :: Bmod_boozer
   real(dp) :: Bmod_vmec
   real(dp) :: start_row(5)
   real(dp) :: max_pos_err
   real(dp) :: max_h_err
   real(dp) :: max_bmod_err
   integer :: ierr
   integer :: i
   integer :: unit

   netcdffile = 'wout.nc'
   multharm = 7
   call spline_vmec_data
   allocate (vmec_field_t :: vmec_field)
   call create_vmec_field(vmec_field)
   call create_gvec_field('wout.gvec_boozer_export.nc', boozer_field)

   max_pos_err = 0.0_dp
   max_h_err = 0.0_dp
   max_bmod_err = 0.0_dp

   open (newunit=unit, file='start_boozer.dat', status='old')
   do i = 1, 4
      read (unit, *) start_row
      x_boozer = start_row(1:3)

      call boozer_field%coords%evaluate_cyl(x_boozer, xcyl)
      call vmec_field%coords%from_cyl(xcyl, x_vmec, ierr)
      if (ierr /= 0) error stop 'test_gvec_boozer: VMEC inverse mapping failed'

      call boozer_field%coords%evaluate_cart(x_boozer, xcart_boozer)
      call vmec_field%coords%evaluate_cart(x_vmec, xcart_vmec)
      max_pos_err = max(max_pos_err, maxval(abs(xcart_boozer - xcart_vmec)))

      call boozer_field%evaluate(x_boozer, Acov_boozer, hcov_boozer, Bmod_boozer)
      call vmec_field%evaluate(x_vmec, Acov_vmec, hcov_vmec, Bmod_vmec)
      call boozer_field%coords%cov_to_cart(x_boozer, hcov_boozer, h_cart_boozer)
      call vmec_field%coords%cov_to_cart(x_vmec, hcov_vmec, h_cart_vmec)

      max_h_err = max(max_h_err, maxval(abs(h_cart_boozer - h_cart_vmec))/ &
                      max(1.0_dp, maxval(abs(h_cart_vmec))))
      max_bmod_err = max(max_bmod_err, abs(Bmod_boozer - Bmod_vmec)/ &
                         max(1.0_dp, abs(Bmod_vmec)))
   end do
   close (unit)

   print '(A,ES12.4)', 'max Boozer position error = ', max_pos_err
   print '(A,ES12.4)', 'max Boozer h-cart error   = ', max_h_err
   print '(A,ES12.4)', 'max Boozer |B| error      = ', max_bmod_err

   if (max_pos_err > pos_tol) error stop 'test_gvec_boozer: position mismatch'
   if (max_h_err > h_tol) error stop 'test_gvec_boozer: field direction mismatch'
   if (max_bmod_err > bmod_tol) error stop 'test_gvec_boozer: field strength mismatch'
end program test_gvec_boozer
