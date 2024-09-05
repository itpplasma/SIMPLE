program test_magfie_coils

use, intrinsic :: iso_fortran_env, only: dp => real64
use simple, only : init_vmec
use magfie_sub, only : VMEC
use velo_mod, only: isw_field_type
use simple_magfie, only: VmecField, CoilsField, create_coils_field
use magfie_sub, only: magfie_vmec
use util, only: twopi

implicit none

class(VmecField), allocatable :: vmec_field
class(CoilsField), allocatable :: coils_field
real(dp) :: dummy, x(3), Acov(3), hcov(3), Bmod

isw_field_type = VMEC

call init_vmec('wout.nc', 5, 5, 5, dummy)
allocate(vmec_field)
coils_field = create_coils_field('coils.5C')

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


subroutine magfie_coils(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    use interpolate, only: evaluate_splines_3d_der

    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

    real(dp) :: Ar, Ath, Aphi, hr, hth, hphi
    real(dp), dimension(3) :: dAr, dAth, dAphi, dhr, dhth, dhphi

    real(dp) :: sqrtg_bmod

    call evaluate_splines_3d_der(coils_field%spl_Ar, x, Ar, dAr)
    call evaluate_splines_3d_der(coils_field%spl_Ath, x, Ath, dAth)
    call evaluate_splines_3d_der(coils_field%spl_Aphi, x, Aphi, dAphi)

    call evaluate_splines_3d_der(coils_field%spl_hr, x, hr, dhr)
    call evaluate_splines_3d_der(coils_field%spl_hth, x, hth, dhth)
    call evaluate_splines_3d_der(coils_field%spl_hphi, x, hphi, dhphi)

    call evaluate_splines_3d_der(coils_field%spl_Bmod, x, bmod, bder)
    bder = bder/bmod

    sqrtg_bmod = dAth(1)*hphi - dAphi(1)*hth + &
                 dAphi(2)*hr - dAr(2)*hphi + &
                 dAr(3)*hth - dAth(3)*hr

    sqrtg = sqrtg_bmod/bmod

    hcovar(1) = hr
    hcovar(2) = hth
    hcovar(3) = hphi

    hctrvr(1) = (dAphi(2) - dAth(3))/sqrtg_bmod
    hctrvr(2) = (dAr(3) - dAphi(1))/sqrtg_bmod
    hctrvr(3) = (dAth(1) - dAr(2))/sqrtg_bmod

    hcurl(1) = (dhphi(2) - dhth(3))/sqrtg
    hcurl(2) = (dhr(3) - dhphi(1))/sqrtg
    hcurl(3) = (dhth(1) - dhr(2))/sqrtg
end subroutine magfie_coils


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

    call init_meiss(vmec_field)
    call get_meiss_coordinates
    call evaluate_meiss(f, r, th, ph, 0)
    print *, 'field_can_meiss(vmec_field)'
    print *, 'A = ', f%Ath, f%Aph
    print *, 'h = ', f%hth, f%hph
    print *, 'B = ', f%Bmod
end subroutine test_can


subroutine test_can_curve
    use field_can_mod, only: FieldCan
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, evaluate_meiss

    type(FieldCan) :: f
    real(dp) :: r, th, ph
    integer :: i, N=1000

    r = 0.3d0
    th = 0.2d0

    call init_meiss(coils_field)
    call get_meiss_coordinates

    do i = 0, N
        ph = i*twopi/N
        call evaluate_meiss(f, r, th, ph, 0)
        write(21, *) r, th, ph, f%Ath, f%Aph, f%hth, f%hph, f%Bmod
    end do

    call init_meiss(vmec_field)
    call get_meiss_coordinates

    do i = 0, N
        ph = i*twopi/N
        call evaluate_meiss(f, r, th, ph, 0)
        write(22, *) r, th, ph, f%Ath, f%Aph, f%hth, f%hph, f%Bmod
    end do
end subroutine test_can_curve


end program test_magfie_coils
