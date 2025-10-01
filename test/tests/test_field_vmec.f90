program test_field_vmec

    use, intrinsic :: iso_fortran_env, only : dp => real64
    use field, only : field_from_file, MagneticField
    use field_vmec, only : VmecField

    implicit none

    class(MagneticField), allocatable :: field_obj
    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    real(dp) :: norm_sq

    call field_from_file('wout.nc', field_obj)

    select type(field_obj)
    type is (VmecField)
        x = [sqrt(0.25_dp), 0.1_dp, 0.2_dp]
        call field_obj%evaluate(x, Acov, hcov, Bmod)

        if (Bmod <= 0.0_dp) then
            error stop 'VmecField: Bmod must be positive'
        end if

        norm_sq = sum(hcov**2)
        if (abs(norm_sq - 1.0_dp) > 5.0d-10) then
            error stop 'VmecField: hcov not normalized'
        end if

    class default
        error stop 'field_from_file did not return VmecField for VMEC input'
    end select

end program test_field_vmec
