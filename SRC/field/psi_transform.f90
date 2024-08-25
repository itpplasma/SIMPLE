module psi_transform
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use binsrc_sub
  use plag_coeff_sub
  implicit none

  contains

  subroutine grid_r_to_psi( n_r, n_phi, n_th, rmin, rmax, psi_inner, psi_outer, &
      psi_of_x, Aph_of_x, hph_of_x, hth_of_x, Bmod_of_x, &
      r_of_xc, Aph_of_xc, hph_of_xc, hth_of_xc, Bmod_of_xc)

      integer, intent(in) :: n_r, n_phi, n_th
      real(dp), intent(in) :: rmin, rmax, psi_inner, psi_outer
      real(dp), dimension(n_r, n_phi, n_th), intent(in) :: &
        psi_of_x, Aph_of_x, hph_of_x, hth_of_x, Bmod_of_x
      real(dp), dimension(n_r, n_phi, n_th), intent(out) :: &
          r_of_xc, Aph_of_xc, hph_of_xc, hth_of_xc, Bmod_of_xc

! ***** This should be an input variable:
!        integer,  intent(in) :: n_psi
! in case n_psi is different from n_r, n_r should be replaced
! with n_psi in the above declaration of output arrays
      real(dp) :: h_psi
      integer  :: n_psi
!
      integer, parameter  :: nder = 0 !(no derivatives)
! ***** This should be the input variable:
!        integer, intent(in) :: nplagr
! ***** For the moment it is a local parameter:
      integer, parameter :: nplagr = 6 !(4 points - cubic polynomial, 6: quintic)
!
      logical :: reverse
      integer :: i_r, i_psi, ibeg, iend, nshift
      integer :: i_th, i_phi
      real(dp) :: psi_fix
      real(dp), dimension(:),   allocatable :: xp, psi_loc, r
      real(dp), dimension(:,:), allocatable :: coef
!
      allocate(xp(nplagr),coef(0:nder,nplagr),psi_loc(n_r),r(n_r))
      nshift = nplagr/2

      do i_r = 1, n_r
        r(i_r) = rmin + dble(i_r-1)*(rmax-rmin)/dble(n_r-1)
      enddo
!
! Equidistant grid in psi:
!
      n_psi = n_r   !This needs not to be the same*****
!
      i_phi = n_phi/2
      i_th = n_th/2
      if(psi_of_x(n_r,i_phi,i_th).gt.psi_of_x(1,i_phi,i_th)) then
        reverse = .false.
      else
        reverse = .true.
      endif

      h_psi=(psi_outer-psi_inner)/dble(n_psi-1)
!
      do i_phi = 1,n_phi
        do i_th = 1,n_th
          if(reverse) then
            psi_loc = psi_of_x(n_r:1:-1, i_phi, i_th)
          else
            psi_loc = psi_of_x(:, i_phi, i_th)
          endif
          do i_psi = 1,n_psi
            psi_fix = psi_inner+h_psi*dble(i_psi-1)
!
            call binsrc(psi_loc,1,n_r,psi_fix,i_r)
!
            ibeg = i_r - nshift
            iend = ibeg+nplagr-1
            if(ibeg.lt.0) then
              ibeg = 1
              iend = ibeg+nplagr-1
            elseif(iend.gt.n_r) then
              iend = n_r
              ibeg = iend-nplagr+1
            endif
!
            call plag_coeff(nplagr,nder,psi_fix,psi_loc(ibeg:iend),coef)
!
            if(reverse) then
              ibeg = n_r - ibeg +1
              iend = n_r - iend +1
              r_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*r(ibeg:iend:-1))
              Aph_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*Aph_of_x(ibeg:iend:-1, i_phi, i_th))
              hph_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*hph_of_x(ibeg:iend:-1, i_phi, i_th))
              hth_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*hth_of_x(ibeg:iend:-1, i_phi, i_th))
              Bmod_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*Bmod_of_x(ibeg:iend:-1, i_phi, i_th))
            else
              r_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*r(ibeg:iend))
              Aph_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*Aph_of_x(ibeg:iend, i_phi, i_th))
              hph_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*hph_of_x(ibeg:iend, i_phi, i_th))
              hth_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*hth_of_x(ibeg:iend, i_phi, i_th))
              Bmod_of_xc(i_psi, i_phi, i_th) = sum(coef(0,:)*Bmod_of_x(ibeg:iend, i_phi, i_th))
            endif
          enddo
        enddo
      enddo
!
      deallocate(xp,coef,psi_loc)
  end subroutine grid_r_to_psi

end module psi_transform
