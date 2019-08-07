program test_sympl

use orbit_symplectic
use field_can_mod
use neo_orb
use new_vmec_stuff_mod, only: rmajor
use diag_mod, only : icounter
use omp_lib

implicit none
save

integer :: ierr, kt

double precision :: z0(4), vpar0, dt

type(FieldCan) :: f
type(SymplecticIntegrator) :: euler1, euler2, midpoint, gauss4
type(MultistageIntegrator) :: verlet, order4, mclachlan4, blanes4
type(NeoOrb) :: norb

integer :: npoiper2
real(8) :: rbig, dtau, dtaumax

integer :: nt

call init_field(norb, 5, 5, 3, 2)

npoiper2 = 64
rbig = rmajor*1.0d2
dtaumax = twopi*rbig/npoiper2
dtau = dtaumax

call init_params(norb, 2, 4, 3.5d6, dtau, dtaumax, 1d-8)  ! fusion alphas)

! Initial conditions
z0(1) = 0.4d0  ! r
z0(2) = 0.7d0  ! theta
z0(3) = 0.1d0  ! phi
vpar0 = 0.8d0  ! parallel velocity
call eval_field(f, z0(1), z0(2), z0(3), 0)

f%mu = .5d0**2*(1.d0-vpar0**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
f%ro0 = ro0/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
f%vpar = vpar0*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

print *, f%ro0, f%mu
z0(4) = vpar0*f%hph + f%Aph/f%ro0  ! p_phi

nt = 10000


npoiper2 = 48
dt = twopi*rbig/npoiper2
call orbit_sympl_init(gauss4, f, z0, dt, 1, 1d-12, 4, 2)
call test_single(gauss4, 'gauss4.out')

npoiper2 = 48
dt = twopi*rbig/npoiper2
call orbit_sympl_init(midpoint, f, z0, dt, 1, 1d-12, 3, 2)
call test_single(midpoint, 'midpoint.out')

npoiper2 = 64
dt = twopi*rbig/npoiper2
call orbit_sympl_init(euler1, f, z0, dt, 1, 1d-12, 1, 2)
call test_single(euler1, 'euler1.out')

call orbit_sympl_init(euler2, f, z0, dt, 1, 1d-12, 2, 2)
call test_single(euler2, 'euler2.out')

npoiper2 = 48
dt = twopi*rbig/npoiper2
call orbit_sympl_init_verlet(verlet, f, z0, dt, 1, 1d-12)
call test_multi(verlet, 'verlet.out')

npoiper2 = 64
dt = twopi*rbig/npoiper2
call orbit_sympl_init_order4(order4, f, z0, dt, 1, 1d-12)
call test_multi(order4, 'order4.out')

npoiper2 = 64
dt = twopi*rbig/npoiper2
call orbit_sympl_init_mclachlan4(mclachlan4, f, z0, dt, 1, 1d-12)
call test_multi(mclachlan4, 'mclachlan4.out')

npoiper2 = 48
dt = twopi*rbig/npoiper2
call orbit_sympl_init_blanes4(blanes4, f, z0, dt, 1, 1d-12)
call test_multi(blanes4, 'blanes4.out')


contains

subroutine test_single(si, outname)
    type(SymplecticIntegrator) :: si
    character(len=*) :: outname

    double precision :: starttime, endtime
    double precision, allocatable :: out(:, :)
    
    allocate(out(5,nt))

    icounter = 0
    out=0d0

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
    print *, outname(1:10), ' Time: ', endtime-starttime, 's, Evaluations: ', icounter

    open(unit=20, file=outname, action='write')
    do kt = 1, nt
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)
end subroutine test_single

subroutine test_multi(mi, outname)
    type(MultistageIntegrator) :: mi
    character(len=*) :: outname

    double precision :: starttime, endtime
    double precision, allocatable :: out(:, :)
    
    allocate(out(5,nt))

    icounter = 0
    out=0d0

    out(1:4,1) = z0
    out(5,1) = f%H
    starttime = omp_get_wtime()
    do kt = 2, nt
        ierr = 0
        call orbit_timestep_sympl_multi(mi, f, ierr)
        if (.not. ierr==0) stop 'Error'
        out(1:4,kt) = mi%stages(1)%z
        out(5,kt) = f%H
    end do
    endtime = omp_get_wtime()
    print *, outname(1:10), ' Time: ', endtime-starttime, 's, Evaluations: ', icounter

    open(unit=20, file=outname, action='write')
    do kt = 1, nt
        write(20,*) out(:,kt)
    end do
    close(20)
    deallocate(out)
end subroutine test_multi

end program test_sympl
