program test_switch_pullback
    ! Independent dGamma oracle, in Z=(X,u,mu,phi) and particle (x,v).
    ! Gamma_GC=(A/eps+u*b-eps*mu*Rstar).dX-eps*mu*dphi,
    ! Rstar=R+curl(b)/2 for selector 1; R for selector 0.
    use gc_fo_switch, only: dp, switch_jets_t, switch_jacobian, switch_forward, SW_OK
    use spectre_jets, only: sj_load, sj_cart, sj_jets, sj_owner, sj_u
    implicit none
    real(dp), parameter :: epsilons(4) = [0.02_dp, 0.01_dp, 0.005_dp, 0.0025_dp]
    real(dp) :: z(6)
    integer :: fixture, owner, status, failures
    logical :: spectral
    character(len=1024) :: path
    failures = 0
    spectral = command_argument_count() > 0
    if (spectral) then
        call get_command_argument(1, path)
        call sj_load(trim(path), status)
        if (status /= SW_OK) error stop 'pullback fixture'
        do owner = 1, 3
            sj_owner = owner
            sj_u = [real(owner, dp) - 0.5_dp, 0.4_dp, 0.3_dp]
            call sj_cart(sj_u, z(1:3))
            z(4:6) = [0.7_dp, 0.6_dp, 0.9_dp]
            call check_case(owner, z)
        end do
    else
        do fixture = 0, 3
            z = [0.2_dp, 0.1_dp, -0.15_dp, 0.7_dp, 0.6_dp, 0.9_dp]
            call check_case(fixture, z)
        end do
    end if
    if (failures /= 0) error stop 'switch two-form pullback oracle'
contains
    subroutine sample_jets(x, jets)
        real(dp), intent(in) :: x(3)
        type(switch_jets_t), intent(inout) :: jets
        integer :: ierr
        if (spectral) then
            jets%axis = [0.0_dp, 0.0_dp, 1.0_dp]
            call sj_jets(x, jets, ierr)
            if (ierr /= SW_OK) error stop 'pullback field jets'
        else
            jets%axis = [1.0_dp, 0.0_dp, 0.0_dp]
            jets%b = [0.0_dp, 0.0_dp, 1.0_dp]
            jets%db = 0.0_dp
            jets%d2b = 0.0_dp
            if (fixture >= 1) then
                jets%b(3) = 1 + 0.3_dp*x(1)
                jets%db(3, 1) = 0.3_dp
            end if
            if (fixture >= 2) then
                jets%b(1) = -0.2_dp*x(3)
                jets%b(2) = 0.4_dp*x(1)
                jets%db(1, 3) = -0.2_dp
                jets%db(2, 1) = 0.4_dp
            end if
            if (fixture >= 3) then
                ! Nonlinear divergence-free field with nonzero B Hessians.
                jets%b = jets%b + [0.2_dp*x(2)**2, 0.15_dp*x(3)**2, 0.1_dp*x(2)**2]
                jets%db(1, 2) = 0.4_dp*x(2)
                jets%db(2, 3) = 0.3_dp*x(3)
                jets%db(3, 2) = 0.2_dp*x(2)
                jets%d2b(1, 2, 2) = 0.4_dp
                jets%d2b(2, 3, 3) = 0.3_dp
                jets%d2b(3, 2, 2) = 0.2_dp
            end if
        end if
    end subroutine sample_jets

    subroutine check_case(label, state)
        integer, intent(in) :: label
        real(dp), intent(in) :: state(6)
        type(switch_jets_t) :: jets, particle
        real(dp) :: xp(3), vp(3), jac(6, 6), op(6, 6), og(6, 6)
        real(dp) :: defect(4, 0:1), rate(0:1), jac_error
        integer :: order, k, i, ierr
        call sample_jets(state(1:3), jets)
        do order = 0, 1
            do k = 1, 4
                jets%epsilon = epsilons(k)
                call switch_jacobian(order, jets, state(1:3), state(4), state(5), &
                    state(6), xp, vp, jac, ierr)
                if (ierr /= SW_OK) error stop 'pullback map'
                call sample_jets(xp, particle)
                op = 0.0_dp
                op(1:3, 1:3) = -skew(particle%b/jets%epsilon)
                do i = 1, 3
                    op(i, i + 3) = -1
                    op(i + 3, i) = 1
                end do
                call gc_form(jets, state, order, og)
                defect(k, order) = maxval(abs(matmul(transpose(jac), &
                    matmul(op, jac)) - og))
                print '(a,i0,a,i0,a,2es18.9)', 'case ', label, ' order ', order, &
                    ' epsilon / absolute two-form defect ', epsilons(k), defect(k, order)
            end do
            rate(order) = log(max(defect(2, order), tiny(1.0_dp))/ &
                max(defect(4, order), tiny(1.0_dp)))/log(4.0_dp)
        end do
        if (.not. spectral .and. label == 0) then
            if (maxval(defect) > 1.0e-11_dp) failures = failures + 1
        else
            ! Selector 0 is an order control, not a first-order canonical map.
            if (rate(1) < 1.7_dp .or. rate(0) > 1.5_dp) failures = failures + 1
        end if
        call check_jacobian(state, epsilons(4), jac_error)
        print '(a,i0,a,3es16.7)', 'case ', label, ' rates(0,1) / Jacobian error ', &
            rate, jac_error
        if (jac_error > 2.0e-6_dp) failures = failures + 1
    end subroutine check_case

    subroutine gc_form(jets, state, order, omega)
        type(switch_jets_t), intent(in) :: jets
        real(dp), intent(in) :: state(6)
        integer, intent(in) :: order
        real(dp), intent(out) :: omega(6, 6)
        real(dp) :: b(3), db(3, 3), d2b(3, 3, 3), bn, p(3), axis(3), c
        real(dp) :: rstar(3), curlr(3), curlb(3), cc(3), curvature(3)
        integer :: i, j, k
        bn = norm2(jets%b)
        b = jets%b/bn
        do j = 1, 3
            p(j) = dot_product(b, jets%db(:, j))
            db(:, j) = (jets%db(:, j) - b*p(j))/bn
        end do
        do j = 1, 3
            do k = 1, 3
                d2b(:, j, k) = (jets%d2b(:, j, k) - db(:, k)*p(j) - db(:, j)*p(k) &
                    - b*(dot_product(db(:, k), jets%db(:, j)) + &
                    dot_product(b, jets%d2b(:, j, k))))/bn
            end do
        end do
        axis = jets%axis/norm2(jets%axis)
        c = dot_product(axis, b)
        do j = 1, 3
            rstar(j) = -c*dot_product(cross(b, axis), db(:, j))/(1 - c*c)
        end do
        curlb = [db(3, 2) - db(2, 3), db(1, 3) - db(3, 1), db(2, 1) - db(1, 2)]
        curlr = -[dot_product(b, cross(db(:, 2), db(:, 3))), &
            dot_product(b, cross(db(:, 3), db(:, 1))), &
            dot_product(b, cross(db(:, 1), db(:, 2)))]
        if (order == 1) then
            cc = 0.0_dp
            do i = 1, 3
                do j = 1, 3
                    cc(i) = cc(i) + d2b(j, i, j) - d2b(i, j, j)
                end do
            end do
            rstar = rstar + 0.5_dp*curlb
            curlr = curlr + 0.5_dp*cc
        end if
        curvature = jets%b/jets%epsilon + state(4)*curlb - &
            jets%epsilon*state(5)*curlr
        omega = 0.0_dp
        omega(1:3, 1:3) = -skew(curvature)
        omega(1:3, 4) = -b
        omega(4, 1:3) = b
        omega(1:3, 5) = jets%epsilon*rstar
        omega(5, 1:3) = -jets%epsilon*rstar
        omega(5, 6) = -jets%epsilon
        omega(6, 5) = jets%epsilon
    end subroutine gc_form

    subroutine check_jacobian(state, eps, error)
        real(dp), intent(in) :: state(6), eps
        real(dp), intent(out) :: error
        type(switch_jets_t) :: jets
        real(dp) :: jac(6, 6), fd(6, 6), points(6, 4), trial(6), xp(3), vp(3)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp), parameter :: offsets(4) = [-2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp]
        real(dp), parameter :: weights(4) = [1.0_dp, -8.0_dp, 8.0_dp, -1.0_dp]
        integer :: i, k, ierr
        call sample_jets(state(1:3), jets)
        jets%epsilon = eps
        call switch_jacobian(1, jets, state(1:3), state(4), state(5), state(6), &
            xp, vp, jac, ierr)
        if (ierr /= SW_OK) error stop 'analytic switch Jacobian'
        do i = 1, 6
            do k = 1, 4
                trial = state
                trial(i) = trial(i) + offsets(k)*h
                call sample_jets(trial(1:3), jets)
                jets%epsilon = eps
                call switch_forward(1, jets, trial(1:3), trial(4), trial(5), trial(6), &
                    xp, vp, ierr)
                if (ierr /= SW_OK) error stop 'finite-difference switch map'
                points(:, k) = [xp, vp]
            end do
            fd(:, i) = matmul(points, weights)/(12*h)
        end do
        error = maxval(abs(jac - fd))/max(1.0_dp, maxval(abs(fd)))
    end subroutine check_jacobian

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

    pure function skew(b) result(s)
        real(dp), intent(in) :: b(3)
        real(dp) :: s(3, 3)
        s(1, :) = [0.0_dp, -b(3), b(2)]
        s(2, :) = [b(3), 0.0_dp, -b(1)]
        s(3, :) = [-b(2), b(1), 0.0_dp]
    end function skew
end program test_switch_pullback
