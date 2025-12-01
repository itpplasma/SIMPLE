program test_magfie_coils

use, intrinsic :: iso_fortran_env, only: dp => real64
use simple, only : init_vmec
use magfie_sub, only : VMEC
use velo_mod, only: isw_field_type
use field, only: VmecField, CoilsField, create_coils_field
use magfie_sub, only: magfie_vmec
use magfie_coils_sub, only: init_magfie_coils_from_file, magfie_coils
use util, only: twopi

implicit none

class(VmecField), allocatable :: vmec_field
class(CoilsField), allocatable :: coils_field
real(dp) :: dummy, x(3), Acov(3), hcov(3), Bmod

isw_field_type = VMEC

call init_vmec('wout.nc', 5, 5, 5, dummy)
allocate(vmec_field)
call create_coils_field('coils.5C', coils_field)
call init_magfie_coils_from_file('coils.5C')

x = [0.3d0, 0.2d0, 0.1d0]

print *, 'x = ', x

call vmec_field%evaluate(x, Acov, hcov, Bmod)
print *, 'vmec_field%evaluate'
print *, 'A = ', Acov
print *, 'h = ', hcov
print *, 'B = ', Bmod

call coils_field%evaluate_direct(x, Acov, hcov, Bmod)
print *, 'coils_field%evaluate_direct'
print *, 'A = ', Acov
print *, 'h = ', hcov
print *, 'B = ', Bmod

call coils_field%evaluate(x, Acov, hcov, Bmod)
print *, 'coils_field%evaluate'
print *, 'A = ', Acov
print *, 'h = ', hcov
print *, 'B = ', Bmod

call test_curve
call test_magfie
call test_magfie_curve
call test_can
call test_can_curve

contains

subroutine test_curve
    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    integer :: i, N=1000

    x = [0.3d0, 0.2d0, 0.0d0]

    do i = 0, N
        x(3) = x(3) + twopi/N
        call vmec_field%evaluate(x, Acov, hcov, Bmod)
        write(1, *) x, Acov, hcov, Bmod
        call coils_field%evaluate(x, Acov, hcov, Bmod)
        write(2, *) x, Acov, hcov, Bmod
    end do
end subroutine test_curve


subroutine test_magfie

    real(dp) :: bmod, sqrtg
    real(dp), dimension(3) :: bder, hcovar, hctrvr, hcurl

    call magfie_vmec(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    print *, 'magfie_vmec'
    print *, 'B = ', bmod
    print *, 'sqrtg = ', sqrtg
    print *, 'Bder = ', bder
    print *, 'hcovar = ', hcovar
    print *, 'hctrvr = ', hctrvr
    print *, 'hcurl = ', hcurl

    call magfie_coils(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    print *, 'magfie_coils'
    print *, 'B = ', bmod
    print *, 'sqrtg = ', sqrtg
    print *, 'Bder = ', bder
    print *, 'hcovar = ', hcovar
    print *, 'hctrvr = ', hctrvr
    print *, 'hcurl = ', hcurl

end subroutine test_magfie


subroutine test_magfie_curve
    real(dp) :: bmod, sqrtg
    real(dp), dimension(3) :: x, bder, hcovar, hctrvr, hcurl
    integer :: i, N=1000

    x = [0.3d0, 0.2d0, 0.0d0]

    do i = 0, N
        x(3) = x(3) + twopi/N
        call magfie_vmec(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        write(11, *) x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl
        call magfie_coils(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        write(12, *) x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl
    end do
end subroutine test_magfie_curve



subroutine test_can
    use field_can_mod, only: FieldCan
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, evaluate_meiss

    type(FieldCan) :: f
    real(dp) :: r, th, ph

    r = 0.3d0
    th = 0.2d0
    ph = 0.1d0

    call init_meiss(coils_field)
    call get_meiss_coordinates
    call evaluate_meiss(f, r, th, ph, 0)
    print *, 'field_can_meiss(coils_field)'
    print *, 'A = ', f%Ath, f%Aph
    print *, 'h = ', f%hth, f%hph
    print *, 'B = ', f%Bmod
    print *, 'sqrtgBctr = ', f%dAph(2) - f%dAth(3), -f%dAph(1), f%dAth(1)

    call init_meiss(vmec_field)
    call get_meiss_coordinates
    call evaluate_meiss(f, r, th, ph, 0)
    print *, 'field_can_meiss(vmec_field)'
    print *, 'A = ', f%Ath, f%Aph
    print *, 'h = ', f%hth, f%hph
    print *, 'B = ', f%Bmod
    print *, 'sqrtgBctr = ', f%dAph(2) - f%dAth(3), -f%dAph(1), f%dAth(1)
end subroutine test_can


subroutine test_can_curve
    use field_can_mod, only: FieldCan
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, evaluate_meiss

    type(FieldCan) :: f
    real(dp) :: r, th, ph
    integer :: i, N=1000

    r = 0.1d0
    th = 0.2d0

    call init_meiss(coils_field)
    call get_meiss_coordinates

    do i = 0, N
        ph = i*twopi/N
        call evaluate_meiss(f, r, th, ph, 0)
        write(21, *) r, th, ph, f%Ath, f%Aph, f%hth, f%hph, f%Bmod
        write(31, *) r, th, ph, f%dAph(2) - f%dAth(3), -f%dAph(1), f%dAth(1)
    end do

    call init_meiss(vmec_field)
    call get_meiss_coordinates

    do i = 0, N
        ph = i*twopi/N
        call evaluate_meiss(f, r, th, ph, 0)
        write(22, *) r, th, ph, f%Ath, f%Aph, f%hth, f%hph, f%Bmod
        write(32, *) r, th, ph, f%dAph(2) - f%dAth(3), -f%dAph(1), f%dAth(1)
    end do
end subroutine test_can_curve


end program test_magfie_coils
