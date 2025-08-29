module diag_albert
!> Diagnostic routines for Albert canonical coordinate system
!> Provides line plots of vector potential and magnetic field strength

use, intrinsic :: iso_fortran_env, only: dp => real64
use fortplot, only: figure, plot, savefig, xlabel, ylabel, title
use field_can_albert, only: Aph_of_xc, Bmod_of_xc, n_r, n_th, n_phi, &
    xmin, xmax

implicit none
private

public :: plot_albert_contours

contains

subroutine plot_albert_contours()
    !> Generate line plots of Aph_of_xc and Bmod_of_xc vs theta and phi
    !> for three radial slices (inner, middle, outer)

    integer :: i_r_inner, i_r_middle, i_r_outer
    integer :: i_th, i_ph
    real(dp), dimension(:), allocatable :: th_array, ph_array
    real(dp), dimension(:), allocatable :: profile_values
    character(len=100) :: filename
    character(len=80) :: plot_title_str
    real(dp) :: s_inner, s_middle, s_outer
    integer :: i_ph_mid, i_th_mid

    ! Define radial slice indices - inner (25%), middle (50%), outer (75%)
    i_r_inner = max(1, n_r / 4)
    i_r_middle = n_r / 2
    i_r_outer = max(1, 3 * n_r / 4)
    
    ! Middle indices for profiles
    i_th_mid = n_th / 2
    i_ph_mid = n_phi / 2
    
    ! Calculate corresponding s values
    s_inner = xmin(1) + (xmax(1) - xmin(1)) * real(i_r_inner - 1, dp) / &
        real(n_r - 1, dp)
    s_middle = xmin(1) + (xmax(1) - xmin(1)) * real(i_r_middle - 1, dp) / &
        real(n_r - 1, dp)
    s_outer = xmin(1) + (xmax(1) - xmin(1)) * real(i_r_outer - 1, dp) / &
        real(n_r - 1, dp)

    ! Allocate arrays
    allocate(th_array(n_th))
    allocate(ph_array(n_phi))
    allocate(profile_values(max(n_th, n_phi)))

    ! Create coordinate arrays
    do i_th = 1, n_th
        th_array(i_th) = xmin(2) + (xmax(2) - xmin(2)) * &
            real(i_th - 1, dp) / real(n_th - 1, dp)
    end do

    do i_ph = 1, n_phi
        ph_array(i_ph) = xmin(3) + (xmax(3) - xmin(3)) * &
            real(i_ph - 1, dp) / real(n_phi - 1, dp)
    end do

    ! Generate plots for inner slice
    
    ! Aph_of_xc vs theta (inner slice)
    profile_values(1:n_th) = Aph_of_xc(i_r_inner, :, i_ph_mid)
    filename = "albert_Aph_inner_theta.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc vs theta at s=", s_inner, " (Albert)"
    
    call figure()
    call plot(th_array, profile_values(1:n_th))
    call xlabel("theta")
    call ylabel("Aph_of_xc") 
    call title(trim(plot_title_str))
    call savefig(trim(filename))

    ! Aph_of_xc vs phi (inner slice)
    profile_values(1:n_phi) = Aph_of_xc(i_r_inner, i_th_mid, :)
    filename = "albert_Aph_inner_phi.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc vs phi at s=", s_inner, " (Albert)"
    
    call figure()
    call plot(ph_array, profile_values(1:n_phi))
    call xlabel("phi")
    call ylabel("Aph_of_xc") 
    call title(trim(plot_title_str))
    call savefig(trim(filename))

    ! Bmod_of_xc vs theta (inner slice)
    profile_values(1:n_th) = Bmod_of_xc(i_r_inner, :, i_ph_mid)
    filename = "albert_Bmod_inner_theta.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc vs theta at s=", s_inner, " (Albert)"
    
    call figure()
    call plot(th_array, profile_values(1:n_th))
    call xlabel("theta")
    call ylabel("Bmod_of_xc") 
    call title(trim(plot_title_str))
    call savefig(trim(filename))

    ! Bmod_of_xc vs phi (inner slice)
    profile_values(1:n_phi) = Bmod_of_xc(i_r_inner, i_th_mid, :)
    filename = "albert_Bmod_inner_phi.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc vs phi at s=", s_inner, " (Albert)"
    
    call figure()
    call plot(ph_array, profile_values(1:n_phi))
    call xlabel("phi")
    call ylabel("Bmod_of_xc") 
    call title(trim(plot_title_str))
    call savefig(trim(filename))

    ! Generate plots for middle slice
    
    ! Aph_of_xc vs theta (middle slice)
    profile_values(1:n_th) = Aph_of_xc(i_r_middle, :, i_ph_mid)
    filename = "albert_Aph_middle_theta.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc vs theta at s=", s_middle, " (Albert)"
    
    call figure()
    call plot(th_array, profile_values(1:n_th))
    call xlabel("theta")
    call ylabel("Aph_of_xc") 
    call title(trim(plot_title_str))
    call savefig(trim(filename))

    ! Aph_of_xc vs phi (middle slice)
    profile_values(1:n_phi) = Aph_of_xc(i_r_middle, i_th_mid, :)
    filename = "albert_Aph_middle_phi.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc vs phi at s=", s_middle, " (Albert)"
    
    call figure()
    call plot(ph_array, profile_values(1:n_phi))
    call xlabel("phi")
    call ylabel("Aph_of_xc") 
    call title(trim(plot_title_str))
    call savefig(trim(filename))

    ! Generate plots for outer slice
    
    ! Aph_of_xc vs theta (outer slice)
    profile_values(1:n_th) = Aph_of_xc(i_r_outer, :, i_ph_mid)
    filename = "albert_Aph_outer_theta.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Aph_of_xc vs theta at s=", s_outer, " (Albert)"
    
    call figure()
    call plot(th_array, profile_values(1:n_th))
    call xlabel("theta")
    call ylabel("Aph_of_xc") 
    call title(trim(plot_title_str))
    call savefig(trim(filename))

    ! Bmod_of_xc vs theta (outer slice)
    profile_values(1:n_th) = Bmod_of_xc(i_r_outer, :, i_ph_mid)
    filename = "albert_Bmod_outer_theta.png"
    write(plot_title_str, '(A,F5.3,A)') &
        "Bmod_of_xc vs theta at s=", s_outer, " (Albert)"
    
    call figure()
    call plot(th_array, profile_values(1:n_th))
    call xlabel("theta")
    call ylabel("Bmod_of_xc") 
    call title(trim(plot_title_str))
    call savefig(trim(filename))

    ! Cleanup
    deallocate(th_array, ph_array, profile_values)

    print *, "Albert coordinate diagnostic plots generated successfully!"
    print *, "Files created:"
    print *, "  Inner slice: albert_Aph_inner_theta.png, albert_Aph_inner_phi.png"
    print *, "              albert_Bmod_inner_theta.png, albert_Bmod_inner_phi.png"
    print *, "  Middle slice: albert_Aph_middle_theta.png, albert_Aph_middle_phi.png"
    print *, "  Outer slice: albert_Aph_outer_theta.png, albert_Bmod_outer_theta.png"

end subroutine plot_albert_contours

end module diag_albert