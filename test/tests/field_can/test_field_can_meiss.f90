program test_field_can_meiss

use, intrinsic :: iso_fortran_env, only: dp => real64
use params, only: read_config
use simple, only: Tracer, init_field
use velo_mod, only: isw_field_type
use simple_magfie, only: VmecField
use field_can_mod, only: eval_field => evaluate, FieldCan, FieldCan_init
use magfie_sub, only: CANFLUX, BOOZER, VMEC, MEISS
use field_can_meiss, only: init_meiss => init, init_transformation, write_transformation
use new_vmec_stuff_mod, only : old_axis_healing, old_axis_healing_boundary
implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

type(Tracer) :: norb

isw_field_type = MEISS

print *, 'init_field'
call init_field(norb, 'wout.nc', 5, 5, 3, 0)
call init_meiss(VmecField(), 128, 4, 4, 0.01d0, 1.0d0, 0.0d0, twopi)

call test_covar_components

stop

print *, 'FieldCanMeiss%init_transformation'
call init_transformation

print *, 'FieldCanMeiss%write_transformation'
call write_transformation('lam_chi.out')

contains

subroutine test_covar_components
    use field_can_meiss, only: thmin, h_th, rmin, h_r, h_phi, ah_cov_on_slice, &
        n_r, n_phi, n_th

    real(dp) :: r, phi, th
    real(dp) :: Ar, Ap, hr, hp
    integer :: i_r, i_phi, i_th
    integer :: funit

    open(newunit=funit, file='covar_components.out')
    write(funit, *) '#', ' r', ' phi', ' th', ' Arcov', ' Apcov', ' hrcov', ' hpcov'
    do i_th = 1, n_th
        th = thmin + h_th*(i_th-1)
        do i_phi = 1, n_phi
            phi = h_phi*(i_phi-1)
            do i_r = 1, n_r
                r = rmin + h_r*(i_r-1)
                call ah_cov_on_slice(r, phi, i_th, Ar, Ap, hr, hp)
                write(funit, *) r, phi, th, Ar, Ap, hr, hp
            end do
        end do
    end do
    close(funit)
end subroutine test_covar_components

end program test_field_can_meiss
