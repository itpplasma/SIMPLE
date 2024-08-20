module field_vmec

use, intrinsic :: iso_fortran_env, only: dp => real64
use new_vmec_stuff_mod, only: rmnc, zmns, lmns => almns, rmns, zmnc, lmnc => almnc, &
    iota => aiota, psitor => phi, s_indices => sps, xm => axm, xn => axn, s, nsurfm, &
    n_mn => nstrm, n_s => kpar
use vector_potentail_mod, only: psitor_a => torflux

implicit none

contains

subroutine init
    use vmec_alloc_sub, only: new_allocate_vmec_stuff
    use vmecin_sub, only: vmecin
    integer :: dummy
    call vmecin(rmnc, zmns, lmns, rmns, zmnc, lmnc, iota, psitor, s_indices, xm, xn, &
        s, dummy, n_mn, n_s, psitor_a)
end subroutine init

! def vmec_to_cyl_ref(ks, th, ph):
! R = 0.0
! P = 0.0
! Z = 0.0

! for kmn in range(v.nstrm):
!     m = int(v.axm[kmn])
!     n = int(v.axn[kmn])
!     cosphase = np.cos(m*th - n*ph)
!     sinphase = np.sin(m*th - n*ph)
!     R += v.rmnc[kmn,ks]*cosphase + v.rmns[kmn,ks]*sinphase
!     P = ph
!     Z += v.zmnc[kmn,ks]*cosphase + v.zmns[kmn,ks]*sinphase
! return R, P, Z
subroutine vmec_to_cyl(ks, th, ph)
    integer, intent(in) :: ks
    real(dp), intent(in) :: th, ph
    real(dp) :: R, P, Z

    ! TODO
    R = 0.0
    P = 0.0
    Z = 0.0
end subroutine vmec_to_cyl

end module field_vmec
