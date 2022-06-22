program test_parse_ants
use iso_fortran_env
use parse_ants, only: process_line
use util, only: pi
implicit none

integer, parameter :: maxlen = 4096
character(len=maxlen) :: line
integer :: unit, iostat
real(8) :: v_par, v_perp, u, v, s

open(newunit=unit, file='launch_0.3')

read(unit, '(A)', iostat=iostat) line  ! Skip header
do  ! Read and parse every line
    read(unit, '(A)', iostat=iostat) line
    call process_line(line, v_par, v_perp, u, v, s)
    if(is_iostat_end(iostat)) exit
    print *, s, 2d0*pi*u, 2d0*pi*v, v_par, v_perp
end do

close(unit)
end program test_parse_ants
