program test_params_energy_config
    use params, only: alpha_energy_ev, effective_alpha_energy_ev, facE_al, &
                      read_config
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer :: errors

    errors = 0

    call test_explicit_energy(errors)
    call test_scaled_energy(errors)

    if (errors /= 0) then
        print *, "ERROR:", errors, "params energy config test(s) failed"
        stop 1
    end if

    print *, "All params energy config tests passed!"

contains

    subroutine test_explicit_energy(errors)
        integer, intent(inout) :: errors
        character(len=256) :: path

        path = "test_alpha_energy_explicit.in"
        call reset_energy_defaults
        call write_config(path, 10.0d0, 1.25d6)
        call read_config(path)

        call expect_close("facE_al explicit case", facE_al, 10.0d0, errors)
        call expect_close("alpha_energy_ev", alpha_energy_ev, 1.25d6, errors)
        call expect_close("explicit effective energy", effective_alpha_energy_ev(), &
                          1.25d6, errors)
    end subroutine test_explicit_energy

    subroutine test_scaled_energy(errors)
        integer, intent(inout) :: errors
        character(len=256) :: path

        path = "test_alpha_energy_scaled.in"
        call reset_energy_defaults
        call write_config(path, 100.0d0, -1.0d0)
        call read_config(path)

        call expect_close("facE_al scaled case", facE_al, 100.0d0, errors)
        call expect_close("unset alpha_energy_ev", alpha_energy_ev, -1.0d0, errors)
        call expect_close("scaled effective energy", effective_alpha_energy_ev(), &
                          3.5d4, errors)
    end subroutine test_scaled_energy

    subroutine write_config(path, fac, energy)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: fac, energy
        integer :: unit

        open (newunit=unit, file=path, status="replace", action="write")
        write (unit, '(a)') "&config"
        write (unit, '(a, es24.16)') "facE_al = ", fac
        write (unit, '(a, es24.16)') "alpha_energy_ev = ", energy
        write (unit, '(a)') "/"
        close (unit)
    end subroutine write_config

    subroutine reset_energy_defaults
        facE_al = 1d0
        alpha_energy_ev = -1d0
    end subroutine reset_energy_defaults

    subroutine expect_close(label, got, expected, errors)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: got, expected
        integer, intent(inout) :: errors

        if (abs(got - expected) > 1d-12*max(1d0, abs(expected))) then
            print *, "ERROR:", trim(label), "got", got, "expected", expected
            errors = errors + 1
        end if
    end subroutine expect_close
end program test_params_energy_config
