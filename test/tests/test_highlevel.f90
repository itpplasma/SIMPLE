module test_highlevel
    use simple
    use util, only: pi
    use funit
    use params, only : Tracer
    use new_vmec_stuff_mod

    implicit none

double precision :: errtol
character(*), parameter :: filename = 'wout.nc'

contains

    @test
    subroutine test_initfield()
        integer :: ans_s, ans_tp, amultharm, aintegmode
        type(Tracer) :: self
        ans_s = 5
        ans_tp = 5
        amultharm = 7
        aintegmode = 0
        
        call init_field(self, filename, ans_s, ans_tp, amultharm, aintegmode)
        
    end subroutine test_initfield

    @test
    subroutine test_volume_and_B00()
        integer             :: L1i
        double precision    :: RT0, R0i, cbfi, bz0i, bf0, volume, B00
        type(Tracer) :: self
        
        print *, 'test_volume_and_B00'
       
        call spline_vmec_data ! initialize splines for VMEC field
        call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius
        self%fper = (2.0*pi)/dble(L1i)   !<= field period
        print *, 'R0 = ', RT0, ' cm, fper = ', self%fper
        call volume_and_B00(volume,B00)
        print *,'volume = ',volume,' cm^3,  B_00 = ',B00,' G'
    end subroutine test_volume_and_B00

end module test_highlevel
