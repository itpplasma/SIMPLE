program simple_diag_chartmap
    !> Diagnostic for chartmap vs VMEC coordinate overlays in R-Z plane.
    !>
    !> Generates visual plots of poloidal cross-sections showing coordinate
    !> lines from both VMEC and chartmap coordinate systems for manual verification.
    !>
    !> Usage:
    !>   diag_chartmap.x wout.nc chartmap.nc [zeta_index]
    !>
    !> Output:
    !>   chartmap_overlay_zeta<N>.pdf - R-Z cross-section at toroidal slice N

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use libneo_coordinates, only: coordinate_system_t, &
        make_vmec_coordinate_system, make_chartmap_coordinate_system
    use timing, only: init_timer
    use params, only: coord_input, field_input
    use util, only: twopi
    use fortplot, only: figure, plot, savefig, xlabel, ylabel, title, legend

    implicit none

    character(256) :: vmec_file, chartmap_file, arg
    integer :: nargs, izeta_plot, nfp
    real(dp) :: dummy

    nargs = command_argument_count()
    if (nargs < 2) then
        print *, 'Usage: diag_chartmap.x <wout.nc> <chartmap.nc> [zeta_index]'
        print *, ''
        print *, 'Generates R-Z overlay plots of VMEC and chartmap coordinate lines.'
        stop 1
    end if

    call get_command_argument(1, vmec_file)
    call get_command_argument(2, chartmap_file)

    izeta_plot = 1
    if (nargs >= 3) then
        call get_command_argument(3, arg)
        read(arg, *) izeta_plot
    end if

    call init_timer()

    coord_input = vmec_file
    field_input = vmec_file
    call init_vmec(trim(vmec_file), 5, 5, 5, dummy)

    call get_nfp_from_vmec(vmec_file, nfp)
    call plot_coordinate_overlay(vmec_file, chartmap_file, izeta_plot, nfp)

    print *, 'Diagnostic complete.'

contains

    subroutine get_nfp_from_vmec(filename, nfp)
        use netcdf
        character(*), intent(in) :: filename
        integer, intent(out) :: nfp
        integer :: ncid, varid, ierr

        ierr = nf90_open(filename, NF90_NOWRITE, ncid)
        ierr = nf90_inq_varid(ncid, 'nfp', varid)
        ierr = nf90_get_var(ncid, varid, nfp)
        ierr = nf90_close(ncid)
    end subroutine get_nfp_from_vmec


    subroutine plot_coordinate_overlay(vmec_file, chartmap_file, izeta, nfp)
        character(*), intent(in) :: vmec_file, chartmap_file
        integer, intent(in) :: izeta, nfp

        class(coordinate_system_t), allocatable :: vmec_coords, chartmap_coords
        real(dp), allocatable :: R_vmec(:), Z_vmec(:)
        real(dp), allocatable :: R_chart(:), Z_chart(:)
        real(dp) :: u(3), x_cyl(3), x_cart(3)
        real(dp) :: zeta_val
        integer :: nth, nsurf, i, j
        real(dp) :: s_vals(5)
        character(100) :: outfile, plot_title

        call make_vmec_coordinate_system(vmec_coords)
        call make_chartmap_coordinate_system(chartmap_coords, chartmap_file)

        nth = 128
        nsurf = 5
        s_vals = [0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 1.0_dp]

        zeta_val = twopi / nfp * real(izeta - 1, dp) / 4.0_dp

        allocate(R_vmec(nth + 1), Z_vmec(nth + 1))
        allocate(R_chart(nth + 1), Z_chart(nth + 1))

        write(outfile, '(A,I0,A)') 'chartmap_overlay_zeta', izeta, '.pdf'
        write(plot_title, '(A,F6.3,A)') 'VMEC vs Chartmap at zeta = ', zeta_val, ' rad'

        call figure()

        do j = 1, nsurf
            do i = 1, nth + 1
                u(1) = s_vals(j)
                u(2) = twopi * real(i - 1, dp) / real(nth, dp)
                u(3) = zeta_val

                call vmec_coords%evaluate_point(u, x_cyl)
                R_vmec(i) = x_cyl(1)
                Z_vmec(i) = x_cyl(3)

                call chartmap_coords%evaluate_point(u, x_cart)
                R_chart(i) = sqrt(x_cart(1)**2 + x_cart(2)**2)
                Z_chart(i) = x_cart(3)
            end do

            if (j == 1) then
                call plot(R_vmec, Z_vmec, linestyle='b-', label='VMEC')
                call plot(R_chart, Z_chart, linestyle='r--', label='chartmap')
            else
                call plot(R_vmec, Z_vmec, linestyle='b-')
                call plot(R_chart, Z_chart, linestyle='r--')
            end if
        end do

        call xlabel('R [cm]')
        call ylabel('Z [cm]')
        call title(trim(plot_title))
        call legend()
        call savefig(trim(outfile))

        print *, 'Saved: ', trim(outfile)

        deallocate(R_vmec, Z_vmec, R_chart, Z_chart)
        deallocate(vmec_coords, chartmap_coords)
    end subroutine plot_coordinate_overlay

end program simple_diag_chartmap
