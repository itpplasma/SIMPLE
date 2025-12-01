module magfie_coils_sub

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use field_coils, only : CoilsField, create_coils_field
  use interpolate, only : evaluate_batch_splines_3d_der

  implicit none

  class(CoilsField), allocatable :: coils_field_gc

contains

  subroutine init_magfie_coils_from_file(filename)
    character(*), intent(in) :: filename

    if (allocated(coils_field_gc)) then
      deallocate(coils_field_gc)
    end if

    call create_coils_field(filename, coils_field_gc)
  end subroutine init_magfie_coils_from_file


  subroutine magfie_coils(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

    real(dp) :: Ar, Ath, Aphi, hr, hth, hphi
    real(dp), dimension(3) :: dAr, dAth, dAphi, dhr, dhth, dhphi
    real(dp) :: y_batch(7), dy_batch(3, 7)
    real(dp) :: sqrtg_bmod

    if (.not. allocated(coils_field_gc)) then
      print *, 'magfie_coils: coils field not initialized'
      error stop
    end if

    call evaluate_batch_splines_3d_der(coils_field_gc%spl_coils_batch, x, y_batch, dy_batch)

    Ar = y_batch(1);    dAr = dy_batch(:, 1)
    Ath = y_batch(2);   dAth = dy_batch(:, 2)
    Aphi = y_batch(3);  dAphi = dy_batch(:, 3)
    hr = y_batch(4);    dhr = dy_batch(:, 4)
    hth = y_batch(5);   dhth = dy_batch(:, 5)
    hphi = y_batch(6);  dhphi = dy_batch(:, 6)
    bmod = y_batch(7);  bder = dy_batch(:, 7)
    bder = bder / bmod

    sqrtg_bmod = dAth(1)*hphi - dAphi(1)*hth + &
                 dAphi(2)*hr  - dAr(2)*hphi  + &
                 dAr(3)*hth   - dAth(3)*hr

    sqrtg = sqrtg_bmod / bmod

    hcovar(1) = hr
    hcovar(2) = hth
    hcovar(3) = hphi

    hctrvr(1) = (dAphi(2) - dAth(3)) / sqrtg_bmod
    hctrvr(2) = (dAr(3)   - dAphi(1)) / sqrtg_bmod
    hctrvr(3) = (dAth(1)  - dAr(2))   / sqrtg_bmod

    hcurl(1) = (dhphi(2) - dhth(3)) / sqrtg
    hcurl(2) = (dhr(3)   - dhphi(1)) / sqrtg
    hcurl(3) = (dhth(1)  - dhr(2))   / sqrtg
  end subroutine magfie_coils

end module magfie_coils_sub
