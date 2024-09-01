module field_vmec

use, intrinsic :: iso_fortran_env, only: dp => real64
use new_vmec_stuff_mod, only: rmnc, zmns, lmns => almns, rmns, zmnc, lmnc => almnc, &
    iota => aiota, psitor => phi, s_indices => sps, xm => axm, xn => axn, s, nsurfm, &
    n_mn => nstrm, ns => kpar, netcdffile
use vector_potentail_mod, only: psitor_a => torflux

implicit none

contains

subroutine init_field_vmec(filename)
    use vmec_alloc_sub, only: new_allocate_vmec_stuff
    use vmecin_sub, only: vmecin

    character(*), intent(in) :: filename

    integer :: dummy

    netcdffile = filename
    call new_allocate_vmec_stuff
    call vmecin(rmnc, zmns, lmns, rmns, zmnc, lmnc, iota, psitor, s_indices, xm, xn, &
        s, dummy, n_mn, ns, psitor_a)
end subroutine init_field_vmec


subroutine deinit_field_vmec()
    use vmec_alloc_sub, only: new_deallocate_vmec_stuff
    call new_deallocate_vmec_stuff
end subroutine deinit_field_vmec


subroutine vmec_to_cyl(ks, th, ph)
    integer, intent(in) :: ks
    real(dp), intent(in) :: th, ph
    real(dp) :: R, P, Z

    real(dp) :: cosphase, sinphase
    integer :: k_mn, m, n

    R = 0d0
    P = 0d0
    Z = 0d0

    do k_mn = 1, n_mn
        m = int(xm(k_mn))
        n = int(xn(k_mn))
        cosphase = cos(m*th - n*ph)
        sinphase = sin(m*th - n*ph)
        R = R + rmnc(k_mn, ks)*cosphase + rmns(k_mn, ks)*sinphase
        P = ph
        Z = Z + zmnc(k_mn, ks)*cosphase + zmns(k_mn, ks)*sinphase
    end do
end subroutine vmec_to_cyl

end module field_vmec
