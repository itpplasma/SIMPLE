program test_field_can_meiss

use, intrinsic :: iso_fortran_env, only: dp => real64
use params, only: read_config
use simple, only: tracer_t
use simple_main, only: init_field
use velo_mod, only: isw_field_type
use field, only: vmec_field_t, create_vmec_field
use field_can_mod, only: eval_field => evaluate, field_can_t, field_can_init
use magfie_sub, only: MEISS
use field_can_meiss, only: init_meiss, init_transformation, &
    spline_transformation, init_canonical_field_components, &
    xmin, h_r, h_phi, h_th, ah_cov_on_slice, n_r, n_phi, n_th, lam_phi, chi_gauge
use new_vmec_stuff_mod, only : old_axis_healing, old_axis_healing_boundary
implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

type(tracer_t) :: norb
type(vmec_field_t) :: magfie

isw_field_type = MEISS
call create_vmec_field(magfie)

print *, 'init_field'
call init_field(norb, 'wout.nc', 5, 5, 3, 0)
call init_meiss(magfie, 128, 4, 4, 0.01d0, 1.0d0, 0.0d0, twopi)

print *, 'test_covar_components'
call test_covar_components

print *, 'field_can_meiss.write_transformation'
call write_transformation('lam_chi.out')

print *, 'test_evaluate_vmec'
call test_evaluate_vmec

print *, 'test_evaluate_meiss'
call test_evaluate_meiss

contains

subroutine test_covar_components
    real(dp) :: r, phi, th
    real(dp) :: Ar, Ap, hr, hp
    integer :: i_r, i_phi, i_th
    integer :: funit

    open(newunit=funit, file='covar_components.out')
    write(funit, *) '#', ' r', ' phi', ' th', ' Arcov', ' Apcov', ' Atcov', &
                    ' hrcov', ' hpcov', ' htcov', ' Bmod'
    do i_phi = 1, n_phi
        phi = xmin(3) + h_phi*(i_phi-1)
        do i_th = 1, n_th
            th = xmin(2) + h_th*(i_th-1)
                do i_r = 1, n_r
                    r = xmin(1) + h_r*(i_r-1)
                    call ah_cov_on_slice(r, phi, i_th, Ar, Ap, hr, hp)
                    write(funit, *) r, phi, th, Ar, Ap, 0d0, hr, hp, 0d0, 0d0
                end do
        end do
    end do
    close(funit)
end subroutine test_covar_components


subroutine write_transformation(filename)
    character(*), intent(in) :: filename

    integer :: funit
    integer :: i_r, i_th, i_phi
    real(dp) :: r, th, phi

    open(newunit=funit, file=filename, status='unknown')
    write(funit, *) '#', ' r', ' phi', ' th', ' lam_phi', ' chi_gauge'

    do i_th=1,n_th
        th = xmin(2) + h_th*(i_th-1)
        do i_phi=1,n_phi
            phi = xmin(3) + h_phi*(i_phi-1)
            do i_r=1,n_r
                r = xmin(1) + h_r*(i_r-1)
                write(funit, *) r, phi, th, lam_phi(i_r, i_th, i_phi), &
                    chi_gauge(i_r, i_th, i_phi)
            enddo
        enddo
    enddo

    close(funit)
end subroutine write_transformation


subroutine test_evaluate_vmec
    real(dp) :: r, phi, th
    real(dp) :: Acov(3), hcov(3), Bmod
    integer :: i_r, i_phi, i_th
    integer :: funit

    open(newunit=funit, file='field_vmec.out')
    write(funit, *) '#', ' r', ' phi', ' th', ' Arcov', ' Apcov', ' Atcov', &
                    ' hrcov', ' hpcov', ' htcov', ' Bmod'
    do i_th = 1, n_th
        th = xmin(2) + h_th*(i_th-1)
        do i_phi = 1, n_phi
            phi = xmin(3) + h_phi*(i_phi-1)
            do i_r = 1, n_r
                r = xmin(1) + h_r*(i_r-1)
                call magfie%evaluate([r, th, phi], Acov, hcov, Bmod)
                write(funit, *) r, phi, th, Acov(1), Acov(3), Acov(2), &
                    hcov(1), hcov(3), hcov(2), Bmod
            end do
        end do
    end do
    close(funit)
end subroutine test_evaluate_vmec


subroutine test_evaluate_meiss
    real(dp) :: r, phi, th
    type(field_can_t) :: f
    integer :: i_r, i_phi, i_th
    integer :: funit

    open(newunit=funit, file='field_can_meiss.out')
    write(funit, *) '#', ' r', ' phi', ' th', ' Arcov', ' Apcov', ' Atcov', &
                    ' hrcov', ' hpcov', ' htcov', ' Bmod'
    do i_th = 1, n_th
        th = xmin(3) + h_th*(i_th-1)
        do i_phi = 1, n_phi
            phi = xmin(2) + h_phi*(i_phi-1)
            do i_r = 1, n_r
                r = xmin(1) + h_r*(i_r-1)
                call eval_field(f, r, th, phi, 0)
                write(funit, *) r, phi, th, 0d0, f%Aph, f%Ath, 0d0, f%hph, f%hth, f%Bmod
            end do
        end do
    end do
    close(funit)
end subroutine test_evaluate_meiss

end program test_field_can_meiss
