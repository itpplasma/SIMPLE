module meiss_line_cut_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    implicit none

    real(dp), parameter :: r0 = 300.0_dp, b0 = 2.0e3_dp, bz = 1.5e4_dp

    type, extends(magnetic_field_t) :: tilt_field_t
    contains
        procedure :: evaluate => tilt_evaluate
    end type tilt_field_t

contains

    ! A = Bz R/2 e_phi - B0 R0 ln(R/R0) e_z, so B = B0 R0/R e_phi + Bz e_z.
    subroutine tilt_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        class(tilt_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)
        real(dp) :: r, ephi(3), b(3)

        r = norm2(x(1:2))
        ephi = [-x(2), x(1), 0.0_dp]/r
        Acov = 0.5_dp*bz*r*ephi - [0.0_dp, 0.0_dp, b0*r0*log(r/r0)]
        b = b0*r0/r*ephi + [0.0_dp, 0.0_dp, bz]
        Bmod = norm2(b)
        hcov = b/Bmod
        if (present(sqgBctr)) sqgBctr = b
    end subroutine tilt_evaluate

end module meiss_line_cut_field

program test_meiss_line_cut
    ! A chart whose outer surfaces tilt toroidally: Z = a rho sin(theta) +
    ! d rho^6 cos(phi). In B = B0 R0/R e_phi + Bz e_z, h_phi = b . e_phi is
    ! proportional to B0 R0 - Bz d rho^6 sin(phi) and changes sign near the
    ! chart edge, where the radial Meiss ODE has no solution. |h_phi| falls
    ! below 0.1 of its axis value only for rho^6 sin(phi) > 0.6, rho > 0.9, so
    ! some lines, not all, must end, none before rho = 0.85, the canonical
    ! field must stay finite on them, and points past an end leave the domain.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use field_splined, only: splined_field_t, create_splined_field
    use cartesian_coordinates, only: cartesian_coordinate_system_t
    use libneo_coordinates, only: coordinate_system_t, make_chartmap_coordinate_system
    use libneo_chartmap_vmec_generator, only: write_chartmap_file
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, set_meiss_line_cut, &
                               meiss_cut_lines, evaluate_meiss, beyond_meiss_cut
    use field_can_base, only: field_can_t
    use meiss_line_cut_field, only: tilt_field_t, r0
    implicit none

    real(dp), parameter :: a = 100.0_dp, d = 60.0_dp

    type(tilt_field_t) :: source
    type(splined_field_t) :: splined
    class(coordinate_system_t), allocatable :: chart
    type(field_can_t) :: f
    integer :: ierr, n_lines, k, i
    integer, allocatable :: n_nodes(:)
    real(dp), allocatable :: nodes(:, :, :)
    real(dp) :: r_last

    call write_chart('meiss_line_cut.chartmap.nc')
    call make_chartmap_coordinate_system(chart, 'meiss_line_cut.chartmap.nc')
    nper = 1
    allocate (cartesian_coordinate_system_t :: source%coords)
    call create_splined_field(source, chart, splined)

    call set_meiss_line_cut(.true.)
    call init_meiss(splined, n_r_=24, n_th_=17, n_phi_=17, rmin=0.1_dp)
    call get_meiss_coordinates(ierr)
    if (ierr /= 0) error stop 'construction with line cuts failed'
    call meiss_cut_lines(n_lines, n_nodes, nodes)
    print '(a,i0,a,i0)', ' cut lines ', n_lines, ' of ', 17*17
    if (n_lines == 0) error stop 'no line was cut'
    if (n_lines == 17*17) error stop 'every line was cut'
    r_last = huge(1.0_dp)
    do k = 1, n_lines
        r_last = min(r_last, nodes(1, n_nodes(k), k))
        do i = 1, n_nodes(k)
            call evaluate_meiss(f, nodes(1, i, k), nodes(2, i, k), nodes(3, i, k), 0)
            if (.not. (abs(f%hph) > 0 .and. abs(f%Bmod) < huge(1.0_dp))) &
                error stop 'canonical field not finite on a cut line'
        end do
    end do
    print '(a,f7.3)', ' smallest last valid radius ', r_last
    k = 1
    if (.not. beyond_meiss_cut(nodes(:, n_nodes(k), k) + [0.05_dp, 0.0_dp, 0.0_dp])) &
        error stop 'a point beyond a cut must lie outside the Meiss domain'
    if (beyond_meiss_cut([0.5_dp, nodes(2, 1, k), nodes(3, 1, k)])) &
        error stop 'a point inside every line must lie in the Meiss domain'
    if (r_last < 0.85_dp) error stop 'a line was cut inside the regular region'
    call set_meiss_line_cut(.false.)
    print *, 'PASSED: Meiss lines end where h_phi vanishes'

contains

    subroutine write_chart(path)
        character(*), intent(in) :: path
        integer, parameter :: nrho = 25, nth = 33, nzeta = 33
        real(dp) :: rho(nrho), th(nth), ze(nzeta), radius
        real(dp), dimension(nrho, nth, nzeta) :: x, y, z
        integer :: i, j, l, ierr
        character(len=256) :: message

        rho = [(real(i - 1, dp)/(nrho - 1), i=1, nrho)]
        th = [(twopi*(j - 1)/nth, j=1, nth)]
        ze = [(twopi*(l - 1)/nzeta, l=1, nzeta)]
        do l = 1, nzeta
            do j = 1, nth
                do i = 1, nrho
                    radius = r0 + a*rho(i)*cos(th(j))
                    x(i, j, l) = radius*cos(ze(l))
                    y(i, j, l) = radius*sin(ze(l))
                    z(i, j, l) = a*rho(i)*sin(th(j)) + d*rho(i)**6*cos(ze(l))
                end do
            end do
        end do
        call write_chartmap_file(path, rho, th, ze, x, y, z, 1, ierr, message)
        if (ierr /= 0) error stop trim(message)
    end subroutine write_chart

end program test_meiss_line_cut
