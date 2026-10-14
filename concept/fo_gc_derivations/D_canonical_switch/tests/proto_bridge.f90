module proto_bridge
    ! C bridge to the SIMPLE-gcho prototype (598f041, compiled read-only from its
    ! sources): GC->particle map and its analytic 6x6 Jacobian for selectors
    ! 0, 1, 11 (1+), 2. Columns (X1,X2,X3,u,mu,phi); rows (x1,x2,x3,v1,v2,v3).
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use gc_api, only: dp, field_derivs_t, gc_transform, gc_transform_jacobian
    implicit none
contains
    subroutine gc_jac_c(n, epsilon, axis, origin, b, db, h, third, z, xv, jac, ierr) &
        bind(C, name='gc_jac_c')
        integer(c_int), value, intent(in) :: n, axis
        real(c_double), value, intent(in) :: epsilon
        real(c_double), intent(in) :: origin(3), b(3), db(3, 3), h(3, 3, 3), third(3, 3, 3, 3), z(6)
        real(c_double), intent(out) :: xv(6), jac(6, 6)
        integer(c_int), intent(out) :: ierr
        type(field_derivs_t) :: f
        integer :: status
        f%origin = origin
        f%b = b
        f%db = db
        f%d2b = h
        f%d3b = third
        f%epsilon = epsilon
        f%gauge_axis = axis
        call gc_transform(n, f, z(1:3), z(4), z(5), z(6), xv(1:3), xv(4:6), status)
        if (status == 0) call gc_transform_jacobian(n, f, z(1:3), z(4), z(5), z(6), jac, status)
        ierr = status
    end subroutine gc_jac_c
end module proto_bridge
