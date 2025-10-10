module tokamak_config_mod

use, intrinsic :: iso_fortran_env, only: dp => real64

implicit none

real(dp) :: tok_R0 = 6.2d0
real(dp) :: tok_epsilon = 0.32d0
real(dp) :: tok_kappa = 1.0d0
real(dp) :: tok_delta = 0.0d0
real(dp) :: tok_A_param = -0.142d0
real(dp) :: tok_B0 = 5.3d0
integer :: tok_Nripple = 0
real(dp) :: tok_a0 = 1.984d0
real(dp) :: tok_alpha0 = 2.0d0
real(dp) :: tok_delta0 = 0.0d0
real(dp) :: tok_z0 = 0.0d0

end module tokamak_config_mod
