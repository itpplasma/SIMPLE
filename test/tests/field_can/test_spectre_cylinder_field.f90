program test_spectre_cylinder_field
    ! Large-aspect-ratio SPECTRE torus against the straight-cylinder Lundquist
    ! field built from the same fluxes and mu: in volume k, curl B = mu_k B,
    ! B_s = cJ J0(|mu| r) + cY Y0(|mu| r), B_theta = sign(mu) (cJ J1 + cY Y1),
    ! with cY = 0 in the core. Toroidicity makes them differ by O(a/R0); the
    ! interface carries a sheet current with a field-direction jump.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_sub, only: spectre_field, set_magfie_spectre_field, &
                          M_TO_CM, TESLA_TO_GAUSS
    use field_spectre, only: spectre_field_t, create_spectre_field
    use spectre_fo_field, only: spectre_fo_field_local, spectre_fo_local_t
    use spectre_fo_hybrid, only: spectre_fo_configure
    implicit none

    real(dp), parameter :: PI = 3.14159265358979323846_dp
    integer, parameter :: NTH = 16, NRHO = 48
    type(spectre_field_t) :: sf
    type(spectre_fo_local_t) :: loc
    character(len=1024) :: h5file
    real(dp) :: r0, a, r1, length, k1, cj(2), cy(2), tor(2), pol(2), m(2, 2)
    real(dp) :: u(3), x(3), bvec(3), braw(3), bmod_raw, r, eps
    real(dp) :: bt_an, bs_an, b_an, pitch_an, pitch_fo, pitch_raw
    real(dp) :: dev_b, dev_p, dev_b_raw, dev_p_raw, pitch_edge(2), jump
    integer :: status, ith, ir, owner, unit, im1, im0

    call get_command_argument(1, h5file)
    call create_spectre_field(sf, trim(h5file), status)
    if (status /= 0) error stop 'fixture initialization failed'
    call set_magfie_spectre_field(sf)
    call spectre_fo_configure(1.0e-9_dp)

    associate (d => spectre_field%data)
        if (d%Mvol /= 2 .or. any(d%in /= 0)) error stop 'need 2 axisymmetric volumes'
        im0 = findloc(d%im, 0, dim=1)
        im1 = findloc(d%im, 1, dim=1)
        r0 = d%Rbc(im0, 2)
        a = d%Rbc(im1, 2)
        r1 = d%Rbc(im1, 1)
        length = 2*PI*r0
        tor = d%phiedge*[d%tflux(1), d%tflux(2) - d%tflux(1)]
        pol = d%phiedge*[d%pflux(1), d%pflux(2) - d%pflux(1)]
        cj(1) = tor(1)*abs(d%mu(1))/(2*PI*r1*bessel_j1(abs(d%mu(1))*r1))
        cy(1) = 0.0_dp
        k1 = abs(d%mu(2))
        m(1, :) = 2*PI/k1*[a*bessel_j1(k1*a) - r1*bessel_j1(k1*r1), &
                           a*bessel_y1(k1*a) - r1*bessel_y1(k1*r1)]
        m(2, :) = sign(1.0_dp, d%mu(2))*length/k1* &
                [bessel_j0(k1*r1) - bessel_j0(k1*a), bessel_y0(k1*r1) - bessel_y0(k1*a)]
        call solve2(m, tor(2), pol(2), cj(2), cy(2))
    end associate
    print '(A,3F10.5)', 'R0, a, r1 [m]: ', r0, a, r1
    print '(A,4ES14.6)', 'Lundquist cJ1, cJ2, cY2: ', cj(1), cj(2), cy(2)

    dev_b = 0; dev_p = 0; dev_b_raw = 0; dev_p_raw = 0; jump = huge(1.0_dp)
    open (newunit=unit, file='cylinder_field.dat', status='replace', action='write')
    write (unit, '(A)') '# theta r_m owner B_fo_T B_raw_T B_an_T '// &
        'pitch_fo pitch_raw pitch_an'
    do ith = 1, NTH
        do ir = 1, 2*NRHO
            owner = merge(1, 2, ir <= NRHO)
            eps = real(mod(ir - 1, NRHO), dp)/real(NRHO - 1, dp)
            ! Clustered at the interface, 1e-9 away from it and the wall.
            eps = 0.02_dp + 0.98_dp*sin(0.5_dp*PI*eps)**2
            if (owner == 1) then
                u = [min(eps, 1 - 1.0e-9_dp), 2*PI*(ith - 1)/NTH, 0.0_dp]
            else
                u = [max(2 - eps, 1 + 1.0e-9_dp), 2*PI*(ith - 1)/NTH, 0.0_dp]
            end if
            call spectre_fo_field_local(u, owner, loc)
            x = loc%x/M_TO_CM
            bvec = loc%b*loc%bmod/TESLA_TO_GAUSS
            call raw_field(u, braw, bmod_raw)
            r = hypot(hypot(x(1), x(2)) - r0, x(3))
            call lundquist(owner, r, bt_an, bs_an)
            b_an = hypot(bt_an, bs_an)
            pitch_an = atan2(bt_an, bs_an)
            pitch_fo = cylinder_pitch(x, bvec)
            pitch_raw = cylinder_pitch(x, braw)
            dev_b = max(dev_b, abs(norm2(bvec)/b_an - 1))
            dev_p = max(dev_p, abs(pitch_fo - pitch_an))
            dev_b_raw = max(dev_b_raw, abs(bmod_raw/b_an - 1))
            dev_p_raw = max(dev_p_raw, abs(pitch_raw - pitch_an))
            if (mod(ir, NRHO) == 0) pitch_edge(owner) = pitch_fo
            write (unit, '(F9.5,ES16.8,I3,6ES16.8)') u(2), r, owner, &
                norm2(bvec), bmod_raw, b_an, pitch_fo, pitch_raw, pitch_an
        end do
        jump = min(jump, abs(pitch_edge(1) - pitch_edge(2)))
    end do
    close (unit)
    print '(A,ES11.3)', 'a/R0:                              ', a/r0
    print '(A,2ES11.3)', 'max |B| rel. deviation (FO, raw):  ', dev_b, dev_b_raw
    print '(A,2ES11.3)', 'max pitch deviation (FO, raw) rad: ', dev_p, dev_p_raw
    print '(A,F8.4)', 'min pitch jump at the sheet (rad): ', jump

    ! Toroidal corrections are below 2 a/R0; the pitch jump of the fixture
    ! (0.284 rad in the cylinder) must survive.
    if (max(dev_b, dev_b_raw) > 2*a/r0) error stop '|B| differs from Lundquist'
    if (max(dev_p, dev_p_raw) > 2*a/r0) error stop 'pitch differs from Lundquist'
    if (jump < 0.25_dp) error stop 'no field-direction jump at the sheet'

contains

    subroutine solve2(mat, b1, b2, x1, x2)
        real(dp), intent(in) :: mat(2, 2), b1, b2
        real(dp), intent(out) :: x1, x2
        real(dp) :: det
        det = mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1)
        x1 = (b1*mat(2, 2) - mat(1, 2)*b2)/det
        x2 = (mat(1, 1)*b2 - mat(2, 1)*b1)/det
    end subroutine solve2

    subroutine lundquist(k, rr, bt, bs)
        integer, intent(in) :: k
        real(dp), intent(in) :: rr
        real(dp), intent(out) :: bt, bs
        real(dp) :: q
        q = abs(spectre_field%data%mu(k))*rr
        bs = cj(k)*bessel_j0(q)
        bt = cj(k)*bessel_j1(q)
        if (k > 1) then
            bs = bs + cy(k)*bessel_y0(q)
            bt = bt + cy(k)*bessel_y1(q)
        end if
        bt = sign(1.0_dp, spectre_field%data%mu(k))*bt
    end subroutine lundquist

    real(dp) function cylinder_pitch(xx, bb) result(p)
        ! Cylinder axis e_phi, theta about it right-handed: e_theta = e_phi x e_r.
        real(dp), intent(in) :: xx(3), bb(3)
        real(dp) :: big_r, e_phi(3), e_r(3), e_th(3)
        big_r = hypot(xx(1), xx(2))
        e_phi = [-xx(2), xx(1), 0.0_dp]/big_r
        e_r = [xx(1)/big_r*(big_r - r0), xx(2)/big_r*(big_r - r0), xx(3)]
        e_r = e_r/norm2(e_r)
        e_th = [e_phi(2)*e_r(3) - e_phi(3)*e_r(2), e_phi(3)*e_r(1) - &
                e_phi(1)*e_r(3), e_phi(1)*e_r(2) - e_phi(2)*e_r(1)]
        p = atan2(dot_product(bb, e_th), dot_product(bb, e_phi))
    end function cylinder_pitch

    subroutine raw_field(uu, bb, bmod)
        real(dp), intent(in) :: uu(3)
        real(dp), intent(out) :: bb(3), bmod
        real(dp) :: acov(3), hcov(3), sqgb(3), g(3, 3), gi(3, 3), sg, e(3, 3)
        call spectre_field%coords%covariant_basis(uu, e)
        call spectre_field%coords%metric_tensor(uu, g, gi, sg)
        call spectre_field%evaluate(uu, acov, hcov, bmod, sqgb)
        bb = matmul(e, sqgb/sg)
    end subroutine raw_field
end program test_spectre_cylinder_field
