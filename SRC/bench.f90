module simple_bench

use field_can_mod, only: FieldCan_init
use orbit_symplectic
use orbit_symplectic_quasi, only: f_quasi => f, si_quasi => si, &
    orbit_timestep_quasi, orbit_timestep_multi_quasi
use simple
use new_vmec_stuff_mod, only: rmajor
use diag_mod, only : icounter
use omp_lib

implicit none
save

integer :: ierr, kt

double precision :: z0(4), vpar0, dt

type(FieldCan) :: f
type(SymplecticIntegrator) :: si
type(MultistageIntegrator) :: mi
type(Tracer) :: norb

integer :: npoiper2
double precision :: rbig, dtau, dtaumax

integer :: nt
double precision, allocatable :: out(:, :)

double precision :: starttime, endtime


logical :: multi  ! use multi-stage scheme
logical :: quasi  ! use quasi-Newton (numerical Jacobian)
logical :: tok    ! tokamak (.True.) or stellarator (.False.)
integer :: integ_mode
integer :: nlag
integer :: nplagr_invar
integer :: ncut
character(*), parameter :: infile = 'wout.nc'
character(*), parameter :: outfile = '/tmp/out.txt'

integer, parameter :: n_tip_vars = 7
double precision, allocatable :: var_cut(:, :)  ! r, th, ph, pph, H, Jpar

double precision :: taub
double precision :: rtol

contains

subroutine init_bench()
    print *, multi, quasi, tok, integ_mode, npoiper2, nlag, nplagr_invar, &
      ncut
    if (tok) then
        ! Initial conditions
        z0(1) = 0.1d0  ! r
        z0(2) = 1.5d0  ! theta
        z0(3) = 0.0d0  ! phi
        vpar0 = 0.0d0  ! parallel velocity

        call FieldCan_init(f, 1d-5, 1d0, vpar0, -1)

        ! Compute toroidal momentum from initial conditions
        call eval_field(f, z0(1), z0(2), z0(3), 0)

        taub = 7800d0  ! estimated bounce time
    else
        call init_field(norb, infile, 5, 5, 3, 1)

        rbig = rmajor*1.0d2
        dtaumax = twopi*rbig/npoiper2
        dtau = dtaumax

        call init_params(norb, 2, 4, 3.5d6, npoiper2, 1, 1d-8)  ! fusion alphas)

        ! Initial conditions
        z0(1) = 0.5d0    ! r
        z0(2) = 0.0d0    ! theta
        z0(3) = 0.314d0  ! phi
        vpar0 = 0.0d0    ! parallel velocity

        call eval_field(f, z0(1), z0(2), z0(3), 0)

        f%mu = .5d0**2*(1.d0-vpar0**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
        f%ro0 = ro0/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
        f%vpar = vpar0*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules
    end if

    z0(4) = vpar0*f%hph + f%Aph/f%ro0  ! p_phi

    if (.not. allocated(out)) allocate(out(6,nt))
    if (.not. allocated(var_cut)) allocate(var_cut(ncut, n_tip_vars))
    out=0d0
    out(1:4,1) = z0
    out(5,1) = f%H
end subroutine init_bench

subroutine do_bench()
    if (tok) then
        dt = taub/npoiper2
    else
        dt = twopi*rbig/npoiper2
    endif

    if (multi) then
        select case (integ_mode)
            case (21)
                call orbit_sympl_init_verlet(mi, f, z0, dt, 1, rtol)
            case (22)
                call orbit_sympl_init_mclachlan4(mi, f, z0, dt, 1, rtol)
        end select
    else
        call orbit_sympl_init(si, f, z0, dt, 1, rtol, integ_mode, nlag)
    end if

    icounter = 0
    starttime = omp_get_wtime()
    if (ncut>0) then
        call test_cuts(nplagr_invar)
    elseif (quasi) then
        call test_quasi
    elseif (multi) then
        call test_multi
    else
        call test_orbit
    end if
    endtime = omp_get_wtime()
    print *, endtime-starttime, icounter

    open(unit=20, file=outfile, action='write', recl=4096)
    if (ncut>0) then
        do kt = 1, ncut
            write(20,*) var_cut(kt, :)
        end do
    else
        do kt = 1, nt
            write(20,*) out(:, kt)
        end do
    end if
    close(20)
end subroutine do_bench

subroutine cleanup_bench()
    if (allocated(var_cut)) deallocate(var_cut)
    if (allocated(out)) deallocate(out)
end subroutine cleanup_bench

subroutine test_cuts(nplagr)
    use plag_coeff_sub, only : plag_coeff

    integer, intent(in) :: nplagr
    integer :: itip
    double precision :: vpar_old, z(4)

    integer, parameter :: nder = 0
    double precision :: orb_sten(n_tip_vars,nplagr), coef(0:nder,nplagr)
    integer :: ipoi(nplagr)

    integer, parameter :: nstep_max = 1000000000
    integer :: i, kcut
    double precision :: par_inv
    double precision :: temptime

    si_quasi = si
    f_quasi = f
    call eval_field(f_quasi, z0(1), z0(2), z0(3), 0)
    call get_derivatives(f_quasi, z0(4))

    par_inv = 0d0
    itip = nplagr/2+1
    do i=1,nplagr
      ipoi(i)=i
    enddo
    vpar_old=vpar0

    do kcut=1, ncut
        do i=1, nstep_max
            if (multi) then
                call orbit_timestep_sympl_multi(mi, f, ierr)
                z = mi%stages(1)%z
            elseif (quasi) then
                call orbit_timestep_quasi(ierr)
                z = si_quasi%z
                ! evaluate field, but don't count
                temptime = omp_get_wtime()
                call eval_field(f, z(1), z(2), z(3), 0)
                call get_derivatives(f, z(4))
                icounter = icounter - 1
                starttime = starttime - (omp_get_wtime() - temptime)
            else
                call orbit_timestep_sympl(si, f, ierr)
                z = si%z
            end if
            if (.not. ierr==0) stop 'Error'


            par_inv = par_inv+f%vpar**2 ! parallel adiabatic invariant

            if(i.le.nplagr) then          !<=first nplagr points to initialize stencil
                orb_sten(1:4,i)=z
                orb_sten(5,i)=f%H
                orb_sten(6,i)=par_inv
                orb_sten(7,i)=f%vpar
            else                          !<=normal case, shift stencil
                orb_sten(1:4,ipoi(1))=z
                orb_sten(5,ipoi(1))=f%H
                orb_sten(6,ipoi(1))=par_inv
                orb_sten(7,ipoi(1))=f%vpar
                ipoi=cshift(ipoi,1)
            endif

            if(vpar_old.lt.0.d0.and.f%vpar.gt.0.d0) itip=0   !<=tip has been passed
            itip=itip+1
            vpar_old=f%vpar

            !<=use only initialized stencil
            if(i .gt. nplagr .and. itip .eq. nplagr/2) then
                !print *, orb_sten(6, ipoi)
                call plag_coeff(nplagr, nder, 0d0, orb_sten(7, ipoi), coef)
                var_cut(kcut, :) = matmul(orb_sten(:, ipoi), coef(0,:))
                var_cut(kcut, 2) = modulo(var_cut(kcut, 2), twopi)
                var_cut(kcut, 3) = modulo(var_cut(kcut, 3), twopi)
                par_inv = par_inv - var_cut(kcut, 6)
                exit
            end if
        end do
    end do
end subroutine test_cuts

subroutine test_orbit
    call test_single()
end subroutine test_orbit

subroutine test_single
    do kt = 2, nt
        ierr = 0
        call orbit_timestep_sympl(si, f, ierr)
        if (.not. ierr==0) stop 'Error'
        out(1:4,kt) = si%z
        out(5,kt) = f%H
        out(6,kt) = f%vpar
    end do
end subroutine test_single

subroutine test_quasi
    f_quasi = f
    si_quasi = si
    do kt = 2, nt
        ierr = 0
        call orbit_timestep_quasi(ierr)
        if (.not. ierr==0) stop 'Error'
        out(1:4,kt) = si_quasi%z
        out(5,kt) = f_quasi%H
        out(6,kt) = f_quasi%vpar
    end do
end subroutine test_quasi

subroutine test_multi
    do kt = 2, nt
        ierr = 0
        call orbit_timestep_sympl_multi(mi, f, ierr)
        if (.not. ierr==0) stop 'Error'
        out(1:4,kt) = mi%stages(1)%z
        out(5,kt) = f%H
        out(6,kt) = f%vpar
    end do
end subroutine test_multi

subroutine test_multi_quasi
    f_quasi = f
    si_quasi = si
    do kt = 2, nt
        ierr = 0
        call orbit_timestep_multi_quasi(mi, ierr)
        if (.not. ierr==0) stop 'Error'
        out(1:4,kt) = mi%stages(1)%z
        out(5,kt) = f_quasi%H
        out(6,kt) = f_quasi%vpar
    end do
end subroutine test_multi_quasi

subroutine minsqdist(za, zref, result)
    double precision, intent(in) :: za(:,:)
    double precision, intent(in) :: zref(:,:)
    double precision, intent(inout) :: result(:)

    integer :: k, l, ka, la
    double precision :: current

    la = size(za, 2)
    l = size(zref, 2)
    result = 1d30

!$omp parallel private(current)
!$omp do
    do ka = 1, la
        do k = 1, l
            current = sum( (za(:,ka) - zref(:,k))**2 )
            if (current < result(ka)) result(ka) = current
        end do
    end do
!$omp end do
!$omp end parallel
end subroutine minsqdist

end module simple_bench
