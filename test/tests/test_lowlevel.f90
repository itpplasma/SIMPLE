module test_lowlevel
    use util, only: pi
    use funit

    implicit none

double precision :: errtol
character(*), parameter :: filename = 'wout.nc'

contains

    @test
    subroutine test_vmec_allocate()
        use new_vmec_stuff_mod
        use vmec_alloc_sub
        print *,'test_vmec_allocate'

        call new_allocate_vmec_stuff
        call new_deallocate_vmec_stuff
    end subroutine test_vmec_allocate

    @test
    subroutine test_spline_vmec_data()
        use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp
        print *,'test_spline_vmec_data'

        netcdffile = filename
        ns_s = 5
        ns_tp = 5
        multharm = 7
        !self%integmode = aintegmode

        call spline_vmec_data

    end subroutine test_spline_vmec_data

    @test
    subroutine test_vmecin()
        use new_vmec_stuff_mod
        use vector_potentail_mod, only : ns,hs,torflux,sA_phi
        double precision, dimension(:,:), allocatable :: splcoe
        double precision, dimension(:,:), allocatable :: almnc_rho,rmnc_rho,zmnc_rho
        double precision, dimension(:,:), allocatable :: almns_rho,rmns_rho,zmns_rho

        print *,'test_vmecin'


        call new_allocate_vmec_stuff
        call vmecin(rmnc,zmns,almns,rmns,zmnc,almnc,aiota,phi,sps,axm,axn,s,    &
              nsurfm,nstrm,kpar,torflux)

    end subroutine test_vmecin

    @test
    subroutine test_stevvo()
        use new_vmec_stuff_mod
        use vector_potentail_mod
        use vmecin_sub, only : stevvo
        integer             :: L1i
        double precision    :: RT0, R0i, cbfi, bz0i, bf0, volume, B00

        print *,'test_stevvo'

        call new_deallocate_vmec_stuff
        call spline_vmec_data ! initialize splines for VMEC field
        call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius

    end subroutine test_stevvo


end module test_lowlevel
