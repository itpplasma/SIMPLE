program test_frozen_current_meiss
    ! A frozen-current field file (coils plus plasma potential, libneo format)
    ! as field_input: the splined field must reproduce A and B of the file, agree
    ! with the guiding-centre tables of the same file, and carry Meiss
    ! coordinates on which the implicit midpoint rule keeps H and p_phi, and
    ! RK45 in the reference coordinates of the same field must follow the
    ! same orbit.
    ! The field is axisymmetric (a current on the z axis and an n = 0 plasma
    ! potential with arbitrary coefficients), so it has no nested flux surfaces
    ! by construction and p_phi is an exact invariant of the orbit.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper, rmajor
    use field_base, only: magnetic_field_t
    use field, only: field_from_file
    use field_splined, only: splined_field_t
    use field_gc_tables, only: gc_tables_t, evaluate_gc_tables
    use field_gc_tables_input, only: build_gc_tables_from_file
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use libneo_chartmap_vmec_generator, only: write_chartmap_file
    use neo_biotsavart, only: coils_init
    use neo_frozen_current_field, only: frozen_current_field_t, frozen_current_patch_t
    use fortfem_meridian_bspline_n, only: meridian_bspline_build
    use neo_frozen_current_io, only: save_frozen_current_field
    use field_can_meiss, only: init_meiss, get_meiss_coordinates
    use field_can_mod, only: field_can_t, field_can_from_id, ref_to_integ, &
                             integ_to_ref, eval_field => evaluate, get_val
    use magfie_sub, only: init_magfie, MEISS, REFCOORDS, set_magfie_refcoords_field
    use alpha_lifetime_sub, only: orbit_timestep_axis
    use parmot_mod, only: rk_ro0 => ro0, rmu
    use orbit_symplectic_base, only: symplectic_integrator_t, MIDPOINT
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl
    implicit none

    real(dp), parameter :: r0 = 300.0_dp, a = 100.0_dp !! chart torus, cm
    real(dp), parameter :: b0 = 2.0e4_dp, bpol = 2.0e3_dp !! G at r0
    real(dp), parameter :: c_light = 2.99792458e10_dp
    character(len=*), parameter :: chart_file = 'frozen_current_torus.chartmap.nc', &
                                   field_file = 'frozen_current_torus.h5'
    integer, parameter :: npts = 200, nsteps = 400
    type(frozen_current_field_t) :: exact
    class(magnetic_field_t), allocatable :: loaded
    type(gc_tables_t) :: tables
    real(dp) :: u(3, npts), err_b, err_h, err_a, err_tab
    integer :: k

    call write_torus_chart
    call write_field_file
    call init_reference_coordinates(chart_file)
    call field_from_file(field_file, loaded)
    select type (loaded)
    type is (splined_field_t)
    class default
        error stop 'a frozen-current file must give a splined field'
    end select

    call random_seed(put=[(20261109 + k, k=1, 64)])
    call random_number(u)
    u(1, :) = 0.05_dp + 0.9_dp*u(1, :)
    u(2:3, :) = twopi*u(2:3, :)
    call field_errors(err_b, err_h, err_a)
    call build_gc_tables_from_file(field_file, ref_coords, 3, tables)
    err_tab = table_difference()
    print '(a,4es11.3)', ' |B|, h, A, tables relative errors', err_b, err_h, &
        err_a, err_tab
    if (err_b > 1.0e-3_dp) error stop '|B| differs from the field file'
    if (err_h > 1.0e-3_dp) error stop 'field direction differs from the field file'
    if (err_a > 1.0e-3_dp) error stop 'A differs from the field file'
    if (err_tab > 2.0e-3_dp) error stop 'splined field and tables differ'

    call check_midpoint_invariants

contains

    subroutine write_torus_chart
        integer, parameter :: nrho = 33, nth = 65, nzeta = 9
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
                    z(i, j, l) = a*rho(i)*sin(th(j))
                end do
            end do
        end do
        call write_chartmap_file(chart_file, rho, th, ze, x, y, z, 1, ierr, message)
        if (ierr /= 0) error stop trim(message)
        nper = 1
        rmajor = r0/100.0_dp
    end subroutine write_torus_chart

    ! Axis current for b0 at r0 (G, statA) and an n = 0 plasma potential with
    ! random coefficients, scaled to a largest poloidal field of bpol on the
    ! chart.
    subroutine write_field_file
        type(frozen_current_patch_t) :: patch
        real(dp) :: current, re(2), bc(3), ac(3), bmax, rz(2)
        integer :: status, i

        current = b0*c_light*r0/2.0_dp
        call coils_init([0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], [-1.0e6_dp, 1.0e6_dp], &
                        [current, 0.0_dp], exact%coils)
        patch%r0 = r0 - 1.5_dp*a
        patch%r1 = r0 + 1.5_dp*a
        patch%z0 = -1.5_dp*a
        patch%z1 = 1.5_dp*a
        patch%nr = 6
        patch%nz = 6
        patch%degree = 3
        patch%kmax = 0
        call meridian_bspline_build(patch%r0, patch%r1, patch%z0, patch%z1, patch%nr, &
                                    patch%nz, patch%degree, exact%plasma%sp, status)
        if (status /= 0) error stop 'meridian spline build'
        exact%plasma%nfp = 1
        exact%plasma%kmax = 0
        allocate (exact%plasma%u(exact%plasma%sp%n1, 0:0))
        call random_seed(put=[(4711 + i, i=1, 64)])
        do i = 1, exact%plasma%sp%n1
            call random_number(re)
            exact%plasma%u(i, 0) = cmplx(re(1) - 0.5_dp, 0.0_dp, dp)
        end do
        block
            type(frozen_current_field_t) :: plasma_only
            plasma_only%plasma = exact%plasma
            bmax = 0
            do i = 1, 1000
                call random_number(rz)
                rz = (2*rz - 1)*a
                call plasma_only%compute_abfield([r0 + rz(1), 0.0_dp, rz(2)], ac, bc)
                bmax = max(bmax, norm2(bc))
            end do
        end block
        exact%plasma%u = exact%plasma%u*(bpol/bmax)
        call save_frozen_current_field(field_file, patch, exact)
    end subroutine write_field_file

    subroutine exact_covariant(x, acov, bcov, bmod)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: acov(3), bcov(3), bmod
        real(dp) :: xc(3), e(3, 3), ac(3), bc(3)

        call ref_coords%evaluate_cart(x, xc)
        call ref_coords%covariant_basis(x, e)
        call exact%compute_abfield(xc, ac, bc)
        acov = matmul(ac, e)
        bcov = matmul(bc, e)
        bmod = norm2(bc)
    end subroutine exact_covariant

    ! The splines hold A up to a constant shift per component, so A is
    ! compared through differences to a reference point.
    subroutine field_errors(err_b, err_h, err_a)
        real(dp), intent(out) :: err_b, err_h, err_a
        real(dp) :: ae(3), be(3), bme, as(3), hs(3), bms, ae0(3), as0(3)
        real(dp) :: unused(3), bmod0, ascale, hscale
        integer :: i

        call exact_covariant(u(:, 1), ae0, unused, bmod0)
        call loaded%evaluate(u(:, 1), as0, unused, bmod0)
        err_b = 0
        err_h = 0
        err_a = 0
        ascale = 0
        hscale = 0
        do i = 2, npts
            call exact_covariant(u(:, i), ae, be, bme)
            call loaded%evaluate(u(:, i), as, hs, bms)
            ascale = max(ascale, maxval(abs(ae - ae0)))
            hscale = max(hscale, maxval(abs(be/bme)))
            err_b = max(err_b, abs(bms - bme)/bme)
            err_h = max(err_h, maxval(abs(hs - be/bme)))
            err_a = max(err_a, maxval(abs((as - as0) - (ae - ae0))))
        end do
        err_h = err_h/hscale
        err_a = err_a/ascale
    end subroutine field_errors

    real(dp) function table_difference() result(err)
        real(dp) :: bmod, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3), as(3), hs(3), &
                    bms
        integer :: i

        err = 0
        do i = 1, npts
            call evaluate_gc_tables(tables, u(:, i), bmod, sqrtg, bder, hcov, hctr, &
                                    hcurl)
            call loaded%evaluate(u(:, i), as, hs, bms)
            err = max(err, abs(bms - bmod)/bmod, maxval(abs(hs - hcov))/norm2(hcov))
        end do
    end function table_difference

    ! A 55 keV proton started across the chart, traced with the implicit
    ! midpoint rule in Meiss coordinates of the loaded field.
    subroutine check_midpoint_invariants
        real(dp), parameter :: ev = 1.602176634e-12_dp, m_p = 1.67262192369e-24_dp, &
                               e_charge = 4.80320471e-10_dp
        type(field_can_t) :: f
        type(symplectic_integrator_t) :: integ
        real(dp) :: v0, ro0, x_ref(3), z(4), lambda, dtau, h0, p0, dh, dp_phi
        real(dp) :: x_end(3), z_rk(5), x_sym(3), x_rk(3), gap
        integer :: j, ierr, taken

        call init_meiss(loaded, n_r_=32, n_th_=33, n_phi_=8, rmin=0.05_dp, &
                        rmax=0.95_dp)
        call get_meiss_coordinates
        call field_can_from_id(MEISS)
        call init_magfie(MEISS)

        v0 = sqrt(2.0_dp*55.0e3_dp*ev/m_p)
        ro0 = v0*m_p*c_light/e_charge
        x_ref = [0.5_dp, 1.0_dp, 0.3_dp]
        lambda = 0.3_dp
        call ref_to_integ(x_ref, z(1:3))
        call eval_field(f, z(1), z(2), z(3), 0)
        f%mu = (1.0_dp - lambda**2)/f%Bmod
        f%ro0 = ro0/sqrt(2.0_dp)
        f%vpar = lambda*sqrt(2.0_dp)
        z(4) = f%vpar*f%hph + f%Aph/f%ro0
        call get_val(f, z(4))
        h0 = f%H
        p0 = z(4)
        dtau = twopi*r0/256.0_dp/sqrt(2.0_dp)
        call orbit_sympl_init(integ, f, z, dtau, 1, 1.0e-12_dp, MIDPOINT)
        dh = 0
        dp_phi = 0
        taken = 0
        do j = 1, nsteps
            call orbit_timestep_sympl(integ, f, ierr)
            if (ierr /= 0) exit
            if (integ%z(1) > 0.95_dp .or. integ%z(1) < 0.05_dp) exit
            call eval_field(f, integ%z(1), integ%z(2), integ%z(3), 0)
            call get_val(f, integ%z(4))
            dh = max(dh, abs(f%H - h0)/h0)
            dp_phi = max(dp_phi, abs(integ%z(4) - p0)/abs(p0))
            taken = j
        end do
        print '(a,i0,a,2es11.3)', ' midpoint steps ', taken, &
            ', max dH/H, dp_phi/p_phi', dh, dp_phi
        if (taken < nsteps) error stop 'orbit left the chart'
        if (dh > 1.0e-3_dp) error stop 'energy not conserved'
        if (dp_phi > 1.0e-9_dp) error stop 'p_phi not conserved'

        ! RK45 time is v0 t; the symplectic step dtau is in units of v0 t/sqrt(2).
        call integ_to_ref(integ%z(1:3), x_end)
        rk_ro0 = ro0
        rmu = 1.0e8_dp
        call set_magfie_refcoords_field(loaded)
        call init_magfie(REFCOORDS)
        z_rk = [x_ref, 1.0_dp, lambda]
        do j = 1, nsteps
            call orbit_timestep_axis(z_rk, sqrt(2.0_dp)*dtau, sqrt(2.0_dp)*dtau, &
                                     1.0e-10_dp, ierr)
            if (ierr /= 0) error stop 'RK45 orbit failed'
        end do
        call ref_coords%evaluate_cart(x_end, x_sym)
        call ref_coords%evaluate_cart(z_rk(1:3), x_rk)
        gap = norm2(x_sym - x_rk)
        print '(a,es11.3,a)', ' midpoint and RK45 end points differ by ', gap, ' cm'
        if (gap > 1.0e-3_dp*a) error stop 'RK45 and midpoint orbits differ'
        print *, 'PASSED: frozen-current field file in Meiss coordinates'
    end subroutine check_midpoint_invariants

end program test_frozen_current_meiss
