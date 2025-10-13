program test_vmec_meiss_transform

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use params, only : read_config, isw_field_type
  use field, only : field_from_file, MagneticField, VmecField
  use field_can_meiss, only : init_meiss, n_r, n_th, n_phi, xmin, h_r, h_th, &
    h_phi, ah_cov_on_slice
  use pyplot_module, only : pyplot, pyplot_wp

  implicit none

  character(len=256), parameter :: config_path = '../../../simple.in'
  character(len=256), parameter :: vmec_path = '../../wout.nc'
  integer, parameter :: num_slices = 4
  real(dp), parameter :: hp_min_threshold = 1.0d-6
  real(dp), parameter :: integrand_limit = 5.0d0

  class(MagneticField), allocatable :: field_obj
  integer :: i_phi_idx, i_phi, i_th, i_r
  integer, dimension(num_slices) :: slice_indices
  real(dp) :: phi_val, th_val, r_val
  real(dp) :: Ar, Ap, hr, hp, integrand
  real(dp) :: max_abs_integrand, min_abs_hp
  character(len=64) :: filename_png, filename_py
  type(pyplot) :: plt
  real(pyplot_wp), allocatable :: theta_vals(:), radial_vals(:)
  real(pyplot_wp), allocatable :: slice_values(:,:)
  character(len=80) :: plot_title

  call read_config(config_path)
  isw_field_type = 3

  call field_from_file(trim(vmec_path), field_obj)

  select type(field_obj)
  type is (VmecField)
    continue
  class default
    error stop 'field_from_file did not return VmecField for VMEC setup'
  end select

  call init_meiss(field_obj)

  max_abs_integrand = 0.0_dp
  min_abs_hp = huge(1.0_dp)

  allocate(theta_vals(n_th), radial_vals(n_r), slice_values(n_th, n_r))
  slice_indices = [1, max(1, n_phi/4 + 1), max(1, n_phi/2 + 1), n_phi]
  do i_th = 1, n_th
    theta_vals(i_th) = real(xmin(2) + h_th*real(i_th - 1, dp), pyplot_wp)
  end do
  do i_r = 1, n_r
    radial_vals(i_r) = real(xmin(1) + h_r*real(i_r - 1, dp), pyplot_wp)
  end do

  do i_phi_idx = 1, num_slices
    i_phi = slice_indices(i_phi_idx)
    phi_val = xmin(3) + h_phi*real(i_phi - 1, dp)

    do i_th = 1, n_th
      th_val = xmin(2) + h_th*real(i_th - 1, dp)
      do i_r = 1, n_r
        r_val = xmin(1) + h_r*real(i_r - 1, dp)
        call ah_cov_on_slice(r_val, phi_val, i_th, Ar, Ap, hr, hp)

        if (abs(hp) <= hp_min_threshold) then
          write(*, '(A,1X,ES12.4,1X,ES12.4,1X,I0)') &
            'hp below threshold at r, theta, phi index:', r_val, th_val, i_phi
          error stop 'Meiss integrand encountered |hp| below threshold'
        end if

        integrand = -hr/hp
        if (ieee_is_nan(integrand) .or. ieee_is_nan(hr) .or. ieee_is_nan(hp)) then
          integrand = 0.0_dp
        else
          max_abs_integrand = max(max_abs_integrand, abs(integrand))
          min_abs_hp = min(min_abs_hp, abs(hp))
        end if

        slice_values(i_th, i_r) = real(integrand, pyplot_wp)
      end do
    end do

    write(filename_png, '(A,I2.2,A)') 'vmec_meiss_integrand_phi', i_phi, '.png'
    write(filename_py, '(A,I2.2,A)') 'vmec_meiss_integrand_phi', i_phi, '.py'
    write(plot_title, '(A,I0,A)') 'VMEC Meiss integrand slice at phi index ', &
      i_phi, ' (stellarator)'

    call plt%initialize(grid=.true., xlabel='theta', ylabel='s', &
      title=trim(plot_title), figsize=[10,8], tight_layout=.true.)
    call plt%add_contour(theta_vals, radial_vals, slice_values, linestyle='-', &
      filled=.true., cmap='viridis', colorbar=.true.)
    call plt%savefig(trim(filename_png), pyfile=trim(filename_py))
    call plt%destroy()
  end do

  if (max_abs_integrand > integrand_limit) then
    write(*, '(A,1X,ES12.4)') 'Max |integrand| observed:', max_abs_integrand
    error stop 'Meiss integrand exceeds expected bound'
  end if

  if (min_abs_hp <= hp_min_threshold) then
    write(*, '(A,1X,ES12.4)') 'Minimum |hp| observed:', min_abs_hp
    error stop 'hp magnitude fell below threshold'
  end if

  deallocate(theta_vals, radial_vals, slice_values)

end program test_vmec_meiss_transform
