program test_spectre_rescale
    !> A SPECTRE file rescaled by tools/spectre_rescale.py in length L and field
    !> B must give, at the same logical point, positions scaled by L, |B| by B,
    !> the same unit field direction and A scaled by B L^2.
    !>
    !> usage: test_spectre_rescale.x <original.h5> <rescaled.h5> <L> <B>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_spectre, only: spectre_field_t, create_spectre_field

    implicit none

    real(dp), parameter :: TOL = 1.0e-11_dp
    type(spectre_field_t) :: base, scaled
    character(len=1024) :: arg
    real(dp) :: length, field, u(3), x0(3), x1(3), a0(3), a1(3), h0(3), h1(3)
    real(dp) :: b0, b1, err
    integer :: ierr, i, j, nfail

    if (command_argument_count() /= 4) error stop 'usage: <orig> <rescaled> <L> <B>'
    call get_command_argument(3, arg)
    read (arg, *) length
    call get_command_argument(4, arg)
    read (arg, *) field
    call get_command_argument(1, arg)
    call create_spectre_field(base, trim(arg), ierr)
    if (ierr /= 0) error stop 'cannot load original file'
    call get_command_argument(2, arg)
    call create_spectre_field(scaled, trim(arg), ierr)
    if (ierr /= 0) error stop 'cannot load rescaled file'

    nfail = 0
    do i = 1, 7
        do j = 1, 5
            u = [real(base%data%Mvol, dp)*(real(i, dp) - 0.5_dp)/7.0_dp, &
                 0.9_dp*j, 0.37_dp*j]
            call base%coords%evaluate_cart(u, x0)
            call scaled%coords%evaluate_cart(u, x1)
            call base%evaluate(u, a0, h0, b0)
            call scaled%evaluate(u, a1, h1, b1)
            ! Covariant h_i = e_i . b scales with the basis, i.e. with L.
            err = max(maxval(abs(x1 - length*x0))/maxval(abs(length*x0)), &
                      abs(b1 - field*b0)/(field*b0), &
                      maxval(abs(h1 - length*h0))/maxval(abs(length*h0)), &
                      maxval(abs(a1 - field*length**2*a0)) &
                      /max(maxval(abs(field*length**2*a0)), tiny(1.0_dp)))
            if (err > TOL) then
                nfail = nfail + 1
                print '(A,3F8.4,A,ES10.3)', 'FAIL at u = ', u, ': relative error ', err
            end if
        end do
    end do
    if (nfail /= 0) error stop 'rescaled SPECTRE field does not scale'
    print '(A)', 'PASS'
end program test_spectre_rescale
