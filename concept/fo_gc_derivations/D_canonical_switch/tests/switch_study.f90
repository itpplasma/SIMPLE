program switch_study
    ! Item D prototype check on SPECTRE fixtures (read-only):
    !  (1) round trips T^-1(T(z)) and T(T^-1(w)) for the first-order switch,
    !      in Cartesian GC variables, SIMPLE GC variables (rho,theta,zeta,
    !      y4=|v|/sqrt2, y5=v_par/|v|) and the canonical FO chart
    !      (q=u, P_i = e_i.v + A_i/eps);
    !  (2) mu defect: max over an exact (RK4-resolved) particle orbit of
    !      |mu(T^-1(x(t),v(t))) - mu(0)|/mu(0) for the first-order inverse,
    !      the zeroth-order inverse and the production-style energy match.
    ! Input lines (stdin): label path owner rho theta zeta lambda phi eps
    !                      nperiods steps_per_period samples_per_period
    ! Output: one JSON object per line.
    use gc_fo_switch
    use spectre_jets
    implicit none
    character(512) :: line, label, path, loaded
    integer :: owner, nper, spp, sper, ios, ierr
    real(dp) :: useed(3), lam, phi, eps
    loaded = ''
    do
        read (*, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        read (line, *) label, path, owner, useed, lam, phi, eps, nper, spp, sper
        if (trim(path) /= trim(loaded)) then
            call sj_load(trim(path), ierr)
            if (ierr /= 0) error stop 'field load failed'
            loaded = path
        end if
        call run_case()
    end do
contains
    subroutine run_case()
        type(switch_jets_t) :: jets
        real(dp) :: X0(3), upar, mu, bn, speed, z0(6), z(6), zs(6), xp(3), vp(3)
        real(dp) :: x2(3), v2(3), err_z(0:1), err_w(0:1), err_simple(0:1)
        real(dp) :: err_chart, err_chart_rel, rho_l, dt, t_g, mu_ref(0:2), dmu(0:2), state(6)
        real(dp) :: zprev(0:1, 6), zz(6), resid, q(3), pc(3), u(3)
        real(dp) :: em(4), y_in(5), y_out(5)
        integer :: order, it, n, iters(0:1), maxit(0:1), nsteps
        sj_owner = owner
        sj_u = useed
        call sj_cart(useed, X0)
        jets%epsilon = eps
        call sj_jets(X0, jets, ierr)
        if (ierr /= SW_OK) error stop 'jets at seed failed'
        bn = norm2(jets%b)
        speed = sqrt(2.0_dp)
        upar = lam*speed
        mu = (speed**2 - upar**2)/(2.0_dp*bn)
        z0 = [X0, upar, mu, phi]
        rho_l = eps*sqrt(2.0_dp*mu*bn)/bn
        y_in = [useed, speed/sqrt(2.0_dp), lam]
        ! (1) round trips
        do order = 0, 1
            call sj_jets(z0(1:3), jets, ierr)
            call switch_forward(order, jets, z0(1:3), z0(4), z0(5), z0(6), xp, vp, ierr)
            if (ierr /= SW_OK) error stop 'forward failed'
            z = 0.0_dp
            call switch_inverse(order, sj_jets, jets, xp, vp, z, ierr, it, resid)
            if (ierr /= SW_OK) error stop 'inverse failed'
            iters(order) = it
            err_z(order) = zerr(z, z0, rho_l, speed, mu)
            ! SIMPLE variables of the recovered GC
            call sj_locate(z(1:3), u, ierr)
            call sj_jets(z(1:3), jets, ierr)
            y_out(1:3) = u
            y_out(4) = sqrt((z(4)**2 + 2.0_dp*z(5)*norm2(jets%b))/2.0_dp)
            y_out(5) = z(4)/(sqrt(2.0_dp)*y_out(4))
            err_simple(order) = maxval(abs([y_out(1) - y_in(1), &
                angdiff(y_out(2), y_in(2)), angdiff(y_out(3), y_in(3)), &
                y_out(4) - y_in(4), y_out(5) - y_in(5)]))
            ! w -> T^-1 -> T, starting from a particle state not produced by T
            x2 = xp + 0.3_dp*rho_l*[0.2_dp, -0.1_dp, 0.3_dp]
            v2 = vp + 0.05_dp*speed*[0.1_dp, 0.2_dp, -0.15_dp]
            zs = 0.0_dp
            call switch_inverse(order, sj_jets, jets, x2, v2, zs, ierr)
            if (ierr /= SW_OK) error stop 'inverse (w) failed'
            call sj_jets(zs(1:3), jets, ierr)
            call switch_forward(order, jets, zs(1:3), zs(4), zs(5), zs(6), xp, vp, ierr)
            err_w(order) = max(maxval(abs(xp - x2))/rho_l, maxval(abs(vp - v2))/speed)
        end do
        ! canonical FO chart round trip at the first-order particle state
        call sj_jets(z0(1:3), jets, ierr)
        call switch_forward(1, jets, z0(1:3), z0(4), z0(5), z0(6), xp, vp, ierr)
        call fo_encode(xp, vp, eps, q, pc)
        call fo_decode(q, pc, eps, x2, v2)
        err_chart = max(maxval(abs(x2 - xp))/rho_l, maxval(abs(v2 - vp))/speed)
        err_chart_rel = max(maxval(abs(x2 - xp))/norm2(xp), maxval(abs(v2 - vp))/speed)
        ! (2) mu defect along the exact particle orbit
        t_g = 8.0_dp*atan(1.0_dp)*eps/bn
        dt = t_g/spp
        nsteps = nper*spp
        state = [xp, vp]
        dmu = 0.0_dp
        maxit = 0
        do n = 0, nsteps
            if (mod(n, spp/sper) == 0) then
                do order = 0, 1
                    if (n == 0) then
                        zz = z0
                    else
                        zz = zprev(order, :)
                    end if
                    call switch_inverse(order, sj_jets, jets, state(1:3), state(4:6), &
                                        zz, ierr, it, warm=.true.)
                    if (ierr /= SW_OK) error stop 'orbit inverse failed'
                    maxit(order) = max(maxit(order), it)
                    zprev(order, :) = zz
                    if (n == 0) mu_ref(order) = zz(5)
                    dmu(order) = max(dmu(order), abs(zz(5) - mu_ref(order))/mu_ref(order))
                end do
                call energy_matching_inverse(sj_jets, jets, state(1:3), state(4:6), &
                                             em(1:3), em(4), zz(5), zz(6), ierr)
                if (ierr /= SW_OK) error stop 'energy matching failed'
                if (n == 0) mu_ref(2) = zz(5)
                dmu(2) = max(dmu(2), abs(zz(5) - mu_ref(2))/mu_ref(2))
            end if
            if (n < nsteps) call rk4(state, dt, eps)
        end do
        write (*, '(a)') '{"label":"'//trim(label)//'","owner":'//itoa(owner)// &
            ',"seed":['//r(useed(1))//','//r(useed(2))//','//r(useed(3))//'],"lambda":'// &
            r(lam)//',"phi":'//r(phi)//',"eps":'//r(eps)//',"B":'//r(bn)// &
            ',"rho_larmor":'//r(rho_l)//',"nperiods":'//itoa(nper)//',"steps_per_period":'// &
            itoa(spp)//',"samples_per_period":'//itoa(sper)// &
            ',"roundtrip_z":['//r(err_z(0))//','//r(err_z(1))//']'// &
            ',"roundtrip_w":['//r(err_w(0))//','//r(err_w(1))//']'// &
            ',"roundtrip_simple":['//r(err_simple(0))//','//r(err_simple(1))//']'// &
            ',"roundtrip_fo_chart":'//r(err_chart)//',"roundtrip_fo_chart_rel":'//r(err_chart_rel)// &
            ',"newton_iterations":['//itoa(iters(0))//','//itoa(iters(1))//']'// &
            ',"orbit_newton_max_iterations":['//itoa(maxit(0))//','//itoa(maxit(1))//']'// &
            ',"mu_defect":{"order0":'//r(dmu(0))//',"order1":'//r(dmu(1))// &
            ',"energy_matching":'//r(dmu(2))//'}}'
        flush (6)
    end subroutine run_case

    real(dp) function zerr(z, z0, lscale, vscale, mu0)
        real(dp), intent(in) :: z(6), z0(6), lscale, vscale, mu0
        zerr = max(maxval(abs(z(1:3) - z0(1:3)))/lscale, abs(z(4) - z0(4))/vscale, &
                   abs(z(5) - z0(5))/mu0, abs(angdiff(z(6), z0(6))))
    end function zerr

    real(dp) function angdiff(a, b)
        real(dp), intent(in) :: a, b
        real(dp) :: twopi
        twopi = 8.0_dp*atan(1.0_dp)
        angdiff = modulo(a - b + 0.5_dp*twopi, twopi) - 0.5_dp*twopi
    end function angdiff

    !> Canonical FO chart: q = logical u of x, P_i = e_i(q).v + A_i(q)/eps.
    subroutine fo_encode(x, v, eps, q, pc)
        real(dp), intent(in) :: x(3), v(3), eps
        real(dp), intent(out) :: q(3), pc(3)
        real(dp) :: ecov(3, 3), acov(3)
        integer :: ierr
        call sj_locate(x, q, ierr)
        call sj_basis(q, ecov)
        call sj_acov(q, acov)
        pc = matmul(transpose(ecov), v) + acov/eps
    end subroutine fo_encode

    subroutine fo_decode(q, pc, eps, x, v)
        real(dp), intent(in) :: q(3), pc(3), eps
        real(dp), intent(out) :: x(3), v(3)
        real(dp) :: ecov(3, 3), acov(3), et(3, 3)
        integer :: ipiv(3), info
        call sj_cart(q, x)
        call sj_basis(q, ecov)
        call sj_acov(q, acov)
        et = transpose(ecov)
        v = pc - acov/eps
        call dgesv(3, 1, et, 3, ipiv, v, 3, info)
    end subroutine fo_decode

    subroutine lorentz(s, eps, f)
        real(dp), intent(in) :: s(6), eps
        real(dp), intent(out) :: f(6)
        real(dp) :: b(3)
        integer :: ierr
        call sj_field_at(s(1:3), b, ierr)
        if (ierr /= SW_OK) error stop 'particle left the owner volume'
        f(1:3) = s(4:6)
        f(4:6) = [s(5)*b(3) - s(6)*b(2), s(6)*b(1) - s(4)*b(3), &
                  s(4)*b(2) - s(5)*b(1)]/eps
    end subroutine lorentz

    subroutine rk4(s, h, eps)
        real(dp), intent(inout) :: s(6)
        real(dp), intent(in) :: h, eps
        real(dp) :: k1(6), k2(6), k3(6), k4(6)
        call lorentz(s, eps, k1)
        call lorentz(s + 0.5_dp*h*k1, eps, k2)
        call lorentz(s + 0.5_dp*h*k2, eps, k3)
        call lorentz(s + h*k3, eps, k4)
        s = s + h*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
    end subroutine rk4

    function r(x) result(s)
        real(dp), intent(in) :: x
        character(:), allocatable :: s
        character(32) :: buf
        write (buf, '(es24.16e3)') x
        s = trim(adjustl(buf))
    end function r

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(16) :: buf
        write (buf, '(i0)') i
        s = trim(buf)
    end function itoa
end program switch_study
