module diag_meiss
!> Diagnostic routines for field_can_meiss.f90
!> Outputs data to CSV for plotting with matplotlib

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_meiss, only: rh_can, ah_cov_on_slice, n_r, xmin, xmax, h_r

implicit none
private

public :: plot_rh_can_vs_rc

contains

subroutine plot_rh_can_vs_rc(i_th_in, i_phi_in, filename)
    !> Output rh_can data to CSV for plotting with matplotlib
    !>
    !> Outputs: r_c, dz1 (-hr/hp), dz2 (Ar + Ap*dz1), hr, hp, Ar, Ap
    !>
    !> @param i_th_in  Theta grid index (default: 1)
    !> @param i_phi_in Phi grid index (default: 1)
    !> @param filename Output CSV filename (default: "diag_meiss.csv")

    integer, intent(in), optional :: i_th_in, i_phi_in
    character(len=*), intent(in), optional :: filename

    integer :: i_th, i_phi
    integer :: i, n_points, unit_num
    real(dp) :: r_c, dz1, dz2, hr, hp, Ar, Ap, phi_c
    real(dp), dimension(2) :: z, dz
    character(len=256) :: output_file

    i_th = 1
    i_phi = 1
    if (present(i_th_in)) i_th = i_th_in
    if (present(i_phi_in)) i_phi = i_phi_in

    output_file = "diag_meiss.csv"
    if (present(filename)) then
        output_file = trim(filename)
        if (index(output_file, ".pdf") > 0) then
            output_file = output_file(1:len_trim(output_file)-4) // ".csv"
        end if
    end if

    n_points = n_r
    z = [0.0_dp, 0.0_dp]
    phi_c = xmin(3) + (xmax(3) - xmin(3)) * real(i_phi-1, dp) / real(n_r-1, dp)

    open(newunit=unit_num, file=trim(output_file), status='replace', action='write')
    write(unit_num, '(A)') "# i_th=" // trim(itoa(i_th)) // " i_phi=" // trim(itoa(i_phi))
    write(unit_num, '(A)') "r_c,dz1,dz2,hr,hp,Ar,Ap"

    do i = 1, n_points
        r_c = xmin(1) + (xmax(1) - xmin(1)) * real(i-1, dp) / real(n_points-1, dp)

        call rh_can(r_c, z, dz, i_th, i_phi)
        dz1 = dz(1)
        dz2 = dz(2)

        call ah_cov_on_slice(r_c, phi_c, i_th, Ar, Ap, hr, hp)

        write(unit_num, '(7(ES18.10,A))') r_c, ",", dz1, ",", dz2, ",", &
            hr, ",", hp, ",", Ar, ",", Ap, ""
    end do

    close(unit_num)
    print *, "Diagnostic data saved to: ", trim(output_file)

contains
    function itoa(i) result(str)
        integer, intent(in) :: i
        character(len=20) :: str
        write(str, '(I0)') i
    end function itoa

end subroutine plot_rh_can_vs_rc

end module diag_meiss
