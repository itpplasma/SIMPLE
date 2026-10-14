module c_bridge
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use gc_api, only: dp, field_derivs_t, gc_transform, gc_inverse, gc_rhs
    implicit none
contains
    subroutine gc_eval_c(mode, n, epsilon, axis, origin, b, db, h, z, out, ierr) &
        bind(C, name='gc_eval_c')
        integer(c_int), value, intent(in) :: mode, n, axis
        real(c_double), value, intent(in) :: epsilon
        real(c_double), intent(in) :: origin(3), b(3), db(3, 3), h(3, 3, 3), z(6)
        real(c_double), intent(out) :: out(6)
        integer(c_int), intent(out) :: ierr
        type(field_derivs_t) :: f
        real(dp) :: bs(3), bp, drift(3)
        integer :: status
        f%origin = origin
        f%b = b
        f%db = db
        f%d2b = h
        f%epsilon = epsilon
        f%gauge_axis = axis
        select case (mode)
        case (0)
            call gc_transform(n, f, z(1:3), z(4), z(5), z(6), &
                out(1:3), out(4:6), status)
        case (1)
            call gc_inverse(n, f, z(1:3), z(4:6), out(1:3), &
                out(4), out(5), out(6), status)
        case (2)
            call gc_rhs(n, f, z(1:3), z(4), z(5), z(6), out(1:3), out(4), out(5), &
                        out(6), bs, bp, drift, status)
        case default
            out = 0.0_dp
            status = 1
        end select
        ierr = status
    end subroutine gc_eval_c
end module c_bridge
