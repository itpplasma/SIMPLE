program spectre_interface_jumps
    !> Direction and magnitude jump of B across each SPECTRE interface.
    !>
    !> Build against a SIMPLE build directory B (link line of any SPECTRE test
    !> in B/build.ninja), e.g.
    !>   gfortran -O2 -fopenmp -I$B/src -I$B/include \
    !>       -I$B/_deps/libneo-build/include spectre_interface_jumps.f90 <libs>
    !> usage: spectre_interface_jumps.x FILE.h5
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_spectre, only: spectre_field_t, create_spectre_field
    implicit none
    integer, parameter :: NTH = 128, NZE = 64
    real(dp), parameter :: EPS = 1.0e-10_dp, TWOPI = 6.283185307179586_dp
    character(1024) :: h5
    type(spectre_field_t) :: field
    real(dp) :: xm(3), xp(3), am(3), ap(3), hm(3), hp(3), bm, bp, sgm(3), sgp(3)
    real(dp) :: sqrtg, cosa, ang, angmax, angrms, rmax, rmin, rrms, dA, normal
    integer :: ierr, k, i, j, mvol, n
    call get_command_argument(1, h5)
    call create_spectre_field(field, trim(h5), ierr)
    if (ierr /= 0) error stop 'load'
    mvol = field%data%Mvol
    print '(a,a,a,i0)', 'file ', trim(h5), '  Mvol ', mvol
    do k = 1, mvol - 1
        angmax = 0; angrms = 0; rmax = -huge(1.0_dp); rmin = huge(1.0_dp); rrms = 0
        dA = 0; normal = 0; n = 0
        do j = 0, NZE - 1
            do i = 0, NTH - 1
                xm = [real(k, dp) - EPS, TWOPI*i/NTH, TWOPI*j/(NZE*field%data%Nfp)]
                xp = xm; xp(1) = real(k, dp) + EPS
                call field%evaluate(xm, am, hm, bm, sgm)
                call field%evaluate(xp, ap, hp, bp, sgp)
                sqrtg = dot_product(hm, sgm)/bm
                cosa = dot_product(sgm, hp)/(sqrtg*bm)
                ang = acos(max(-1.0_dp, min(1.0_dp, cosa)))
                angmax = max(angmax, ang); angrms = angrms + ang**2
                rmax = max(rmax, bp/bm - 1); rmin = min(rmin, bp/bm - 1)
                rrms = rrms + (bp/bm - 1)**2
                dA = max(dA, maxval(abs(ap(2:3) - am(2:3)))/max(maxval(abs(am(2:3))), 1e-30_dp))
                normal = max(normal, abs(sgm(1))/maxval(abs(sgm)))
                n = n + 1
            end do
        end do
        print '(a,i0,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,es9.2,a,es9.2)', 'interface ', k, &
            ':  dalpha rad max ', angmax, ' rms ', sqrt(angrms/n), &
            '   |B+|/|B-|-1 min ', rmin, ' max ', rmax, ' rms ', sqrt(rrms/n), &
            '   max rel [[A_t]] ', dA, '  max B^rho/|B| ', normal
    end do
end program spectre_interface_jumps
