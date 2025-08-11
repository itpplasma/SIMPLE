program test_sympl

use orbit_symplectic
use orbit_symplectic_quasi, only: f_quasi => f, si_quasi => si, &
    orbit_timestep_quasi, orbit_timestep_multi_quasi
use field_can_mod
use diag_mod, only : icounter

implicit none
save

double precision, parameter :: qe = 1d0, m = 1d0, c = 1d0, mu = 1d-5
integer :: steps_per_bounce, nbounce

integer :: ierr, kt

double precision :: z0(4), vpar0, dt, taub

type(FieldCan) :: f

type(SymplecticIntegrator) :: integ_euler1, integ_euler2, integ_midpoint, integ_gauss2, &
    integ_gauss4, integ_lobatto4, integ_gauss6, integ_gauss8
type(MultistageIntegrator) :: verlet, mclachlan4, blanes4, kahan6, kahan8

! Initial conditions
z0(1) = 0.1d0  ! r
z0(2) = 1.5d0  ! theta
z0(3) = 0.0d0  ! phi
vpar0 = 0.0d0  ! parallel velocity

call field_can_from_name('test')
call FieldCan_init(f, mu, c*m/qe, vpar0)

! Compute toroidal momentum from initial conditions
call eval_field(f, z0(1), z0(2), z0(3), 0)
z0(4) = m*vpar0*f%hph + qe/c*f%Aph  ! p_phi

taub = 7800d0  ! estimated bounce time
nbounce = 1000


steps_per_bounce = 8
dt = taub/steps_per_bounce

print *, 'timesteps: ', nbounce*steps_per_bounce

call orbit_sympl_init(integ_euler1, f, z0, dt, 1, 1d-12, 1)
call test_single(integ_euler1, 'euler1.out')

stop
call orbit_sympl_init(integ_euler2, f, z0, dt, 1, 1d-12, 2)
call test_single(integ_euler2, 'euler2.out')
print *, ''

call orbit_sympl_init(integ_gauss2, f, z0, dt, 1, 1d-12, 4)
call test_single(integ_gauss2, 'gauss2.out')
call orbit_sympl_init(integ_midpoint, f, z0, dt, 1, 1d-12, 3)
call test_single(integ_midpoint, 'midpoint.out')
call orbit_sympl_init_verlet(verlet, f, z0, dt, 1, 1d-12)
call test_multi(verlet, 'verlet.out')
print *, ''

call orbit_sympl_init(integ_gauss4, f, z0, dt, 1, 1d-12, 5)
call test_single(integ_gauss4, 'gauss4.out')
call orbit_sympl_init_mclachlan4(mclachlan4, f, z0, dt, 1, 1d-12)
call test_multi(mclachlan4, 'mclachlan4.out')
call orbit_sympl_init_blanes4(blanes4, f, z0, dt, 1, 1d-12)
call test_multi(blanes4, 'blanes4.out')
call orbit_sympl_init(integ_lobatto4, f, z0, dt, 1, 1d-12, 15)
call test_single(integ_lobatto4, 'lobatto4.out')
call orbit_sympl_init(integ_gauss6, f, z0, dt, 1, 1d-12, 6)
call test_single(integ_gauss6, 'gauss6.out')
print *, ''

call orbit_sympl_init(integ_gauss4, f, z0, dt, 1, 1d-12, 5)
call test_quasi(integ_gauss4, 'gauss4q.out')
call orbit_sympl_init(integ_gauss6, f, z0, dt, 1, 1d-12, 6)
call test_quasi(integ_gauss6, 'gauss6q.out')
call orbit_sympl_init(integ_lobatto4, f, z0, dt, 1, 1d-12, 15)
call test_quasi(integ_lobatto4, 'lobatto4q.out')
print *, ''

call orbit_sympl_init(integ_gauss4, f, z0, dt, 1, 1d-12, 5)
call test_single(integ_gauss4, 'gauss4.out')
print *, ''

call orbit_sympl_init_kahan6(kahan6, f, z0, dt, 1, 1d-12)
call test_multi(kahan6, 'kahan6.out')
call orbit_sympl_init(integ_gauss6, f, z0, dt, 1, 1d-12, 6)
call test_single(integ_gauss6, 'gauss6.out')
print *, ''

call orbit_sympl_init_kahan8(kahan8, f, z0, dt, 1, 1d-12)
call test_multi(kahan8, 'kahan8.out')
call orbit_sympl_init(integ_gauss8, f, z0, dt, 1, 1d-12, 7)
call test_single(integ_gauss8, 'gauss8.out')
print *, ''

call orbit_sympl_init(integ_euler1, f, z0, dt, 1, 1d-12, 1)
call test_quasi(integ_euler1, 'euler1q.out')
call orbit_sympl_init(integ_euler2, f, z0, dt, 1, 1d-12, 2)
call test_quasi(integ_euler2, 'euler2q.out')
print *, ''

call orbit_sympl_init_verlet(verlet, f, z0, dt, 1, 1d-12)
call test_multi_quasi(verlet, 'verletq.out')
call orbit_sympl_init(integ_midpoint, f, z0, dt, 1, 1d-12, 3)
call test_quasi(integ_midpoint, 'midpointq.out')
call orbit_sympl_init(integ_gauss2, f, z0, dt, 1, 1d-12, 4)
call test_quasi(integ_gauss2, 'gauss2q.out')
print *, ''

call orbit_sympl_init(integ_lobatto4, f, z0, dt, 1, 1d-12, 15)
call test_quasi(integ_lobatto4, 'lobatto4q.out')
call orbit_sympl_init_mclachlan4(mclachlan4, f, z0, dt, 1, 1d-12)
call test_multi_quasi(mclachlan4, 'mclachlan4q.out')
call orbit_sympl_init_blanes4(blanes4, f, z0, dt, 1, 1d-12)
call test_multi_quasi(blanes4, 'blanes4q.out')
call orbit_sympl_init(integ_gauss4, f, z0, dt, 1, 1d-12, 5)
call test_quasi(integ_gauss4, 'gauss4q.out')
print *, ''

call orbit_sympl_init_kahan6(kahan6, f, z0, dt, 1, 1d-12)
call test_multi_quasi(kahan6, 'kahan6q.out')
call orbit_sympl_init(integ_gauss6, f, z0, dt, 1, 1d-12, 6)
call test_quasi(integ_gauss6, 'gauss6q.out')
print *, ''

call orbit_sympl_init_kahan8(kahan8, f, z0, dt, 1, 1d-12)
call test_multi_quasi(kahan8, 'kahan8q.out')
call orbit_sympl_init(integ_gauss8, f, z0, dt, 1, 1d-12, 7)
call test_quasi(integ_gauss8, 'gauss8q.out')
print *, ''

contains

subroutine test_single(si, outname)
    type(SymplecticIntegrator) :: si
    character(*) :: outname

    integer :: nt
    double precision :: starttime, endtime
    double precision, allocatable :: out(:, :)

    nt = nbounce*steps_per_bounce
    allocate(out(5,nt))

    icounter = 0
    out(:,1:)=0d0

    out(1:4,1) = z0
    out(5,1) = f%H
    starttime = omp_get_wtime()
    do kt = 2, nt
        ierr = 0
        call orbit_timestep_sympl(si, f, ierr)
        if (.not. ierr==0) then
            print *, si%z
            exit
        endif
        out(1:4,kt) = si%z
        out(5,kt) = f%H
    end do
    endtime = omp_get_wtime()
    print *, outname(1:10), endtime-starttime, icounter

    open(unit=20, file=outname, action='write', recl=4096)
    do kt = 1, nt
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)
end subroutine test_single

subroutine test_quasi(si, outname)
    type(SymplecticIntegrator) :: si
    character(*) :: outname

    integer :: nt
    double precision :: starttime, endtime
    double precision, allocatable :: out(:, :)

    nt = nbounce*steps_per_bounce
    allocate(out(5,nt))

    icounter = 0
    out(:,1:)=0d0

    out(1:4,1) = z0
    out(5,1) = f%H

    f_quasi = f
    si_quasi = si

    starttime = omp_get_wtime()
    do kt = 2, nt
        ierr = 0
        call orbit_timestep_quasi(ierr)
        if (.not. ierr==0) then
            print *, si%z
            exit
        endif
        out(1:4,kt) = si_quasi%z
        out(5,kt) = f_quasi%H
    end do
    endtime = omp_get_wtime()
    print *, outname(1:10), endtime-starttime, icounter

    open(unit=20, file=outname, action='write', recl=4096)
    do kt = 1, nt
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)

end subroutine test_quasi

subroutine test_multi(mi, outname)
    type(MultistageIntegrator) :: mi
    character(*) :: outname

    integer :: nt
    double precision :: starttime, endtime
    double precision, allocatable :: out(:, :)

    nt = nbounce*steps_per_bounce
    allocate(out(5,nt))

    icounter = 0
    out(:,1:nbounce*steps_per_bounce)=0d0

    out(1:4,1) = z0
    out(5,1) = f%H

    starttime = omp_get_wtime()
    do kt = 2, nbounce*steps_per_bounce
        ierr = 0
        call orbit_timestep_sympl_multi(mi, f, ierr)
        if (.not. ierr==0) stop 'Error'
        out(1:4,kt) = mi%stages(1)%z
        out(5,kt) = f%H
    end do
    endtime = omp_get_wtime()
    print *, outname(1:10), endtime-starttime, icounter

    open(unit=20, file=outname, action='write', recl=4096)
    do kt = 1, nbounce*steps_per_bounce
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)
end subroutine test_multi

subroutine test_multi_quasi(mi, outname)
    type(MultistageIntegrator) :: mi
    character(*) :: outname

    integer :: nt
    double precision :: starttime, endtime
    double precision, allocatable :: out(:, :)

    f_quasi = f

    nt = nbounce*steps_per_bounce
    allocate(out(5,nt))

    icounter = 0
    out(:,1:nbounce*steps_per_bounce)=0d0

    out(1:4,1) = z0
    out(5,1) = f_quasi%H

    starttime = omp_get_wtime()
    do kt = 2, nbounce*steps_per_bounce
        ierr = 0
        call orbit_timestep_multi_quasi(mi, ierr)
        if (.not. ierr==0) stop 'Error'
        out(1:4,kt) = mi%stages(1)%z
        out(5,kt) = f_quasi%H
    end do
    endtime = omp_get_wtime()
    print *, outname(1:10), endtime-starttime, icounter

    open(unit=20, file=outname, action='write', recl=4096)
    do kt = 1, nbounce*steps_per_bounce
       write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)
end subroutine test_multi_quasi

end program test_sympl
