module field_can_boozer
use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_base, only: field_can_t, n_field_evaluations, twopi

implicit none

contains

subroutine evaluate_boozer(f, r, th_c, ph_c, mode_secders)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    call eval_field_booz(f, r, th_c, ph_c, mode_secders)

    n_field_evaluations = n_field_evaluations + 1
end subroutine evaluate_boozer


subroutine integ_to_ref_boozer(xinteg, xref)
    use boozer_sub, only: boozer_to_vmec
    real(dp), intent(in) :: xinteg(3)
    real(dp), intent(out) :: xref(3)

    xref(1) = xinteg(1)
    call boozer_to_vmec( &
        xinteg(1), mod(xinteg(2), twopi), mod(xinteg(3), twopi), xref(2), xref(3))
    xref(2) = mod(xref(2), twopi)
    xref(3) = mod(xref(3), twopi)
end subroutine integ_to_ref_boozer


subroutine ref_to_integ_boozer(xref, xinteg)
    use boozer_sub, only: vmec_to_boozer
    real(dp), intent(in) :: xref(3)
    real(dp), intent(out) :: xinteg(3)

    xinteg(1) = xref(1)
    call vmec_to_boozer( &
        xref(1), mod(xref(2), twopi), mod(xref(3), twopi), xinteg(2), xinteg(3))
    xinteg(2) = mod(xinteg(2), twopi)
    xinteg(3) = mod(xinteg(3), twopi)
end subroutine ref_to_integ_boozer


subroutine eval_field_booz_many(npts, r, th_c, ph_c, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod)
    use boozer_sub, only: splint_boozer_coord_many

    integer, intent(in) :: npts
    real(dp), intent(in) :: r(npts), th_c(npts), ph_c(npts)
    real(dp), intent(out) :: Ath(npts), Aph(npts)
    real(dp), intent(out) :: dAth_dr(npts), dAph_dr(npts), d2Aph_dr2(npts)
    real(dp), intent(out) :: hth(npts), hph(npts)
    real(dp), intent(out) :: dhth(3, npts), dhph(3, npts)
    real(dp), intent(out) :: d2hth_dr2(npts), d2hph_dr2(npts)
    real(dp), intent(out) :: Bmod(npts), dBmod(3, npts), d2Bmod(6, npts)

    real(dp), allocatable :: r_local(:), d3Aph_dr3(:)
    real(dp), allocatable :: Bth(:), Bph(:), dBth(:), dBph(:), d2Bth(:), d2Bph(:)
    real(dp) :: bmod2
    integer :: i

    allocate(r_local(npts), d3Aph_dr3(npts))
    allocate(Bth(npts), Bph(npts), dBth(npts), dBph(npts), d2Bth(npts), d2Bph(npts))
    r_local = r

    call splint_boozer_coord_many(npts, r_local, th_c, ph_c, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, d3Aph_dr3, &
        Bth, dBth, d2Bth, Bph, dBph, d2Bph, &
        Bmod, dBmod, d2Bmod)

    do i = 1, npts
        bmod2 = Bmod(i)**2

        hth(i) = Bth(i) / Bmod(i)
        hph(i) = Bph(i) / Bmod(i)

        dhth(1, i) = dBth(i) / Bmod(i) - Bth(i) * dBmod(1, i) / bmod2
        dhph(1, i) = dBph(i) / Bmod(i) - Bph(i) * dBmod(1, i) / bmod2
        dhth(2, i) = -Bth(i) * dBmod(2, i) / bmod2
        dhth(3, i) = -Bth(i) * dBmod(3, i) / bmod2
        dhph(2, i) = -Bph(i) * dBmod(2, i) / bmod2
        dhph(3, i) = -Bph(i) * dBmod(3, i) / bmod2

        d2hth_dr2(i) = d2Bth(i) / Bmod(i) - 2.0d0 * dBth(i) * dBmod(1, i) / bmod2 &
            + Bth(i) / bmod2 * (2.0d0 * dBmod(1, i)**2 / Bmod(i) - d2Bmod(1, i))
        d2hph_dr2(i) = d2Bph(i) / Bmod(i) - 2.0d0 * dBph(i) * dBmod(1, i) / bmod2 &
            + Bph(i) / bmod2 * (2.0d0 * dBmod(1, i)**2 / Bmod(i) - d2Bmod(1, i))
    end do
end subroutine eval_field_booz_many


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine eval_field_booz(f, r, th_c, ph_c, mode_secders)
    use boozer_sub, only: splint_boozer_coord
    use vector_potentail_mod, only : torflux
    !
    ! Evaluates magnetic field in Boozer canonical coordinates (r, th_c, ph_c)
    ! and stores results in variable f
    ! Works for A_th linear in r (toroidal flux as radial variable)
    !
    ! mode_secders = 0: no second derivatives
    ! mode_secders = 1: second derivatives only in d/dr^2
    ! mode_secders = 2: all second derivatives, including mixed
    !
    ! tested in test_magfie.f90, 2018-10-23, C. Albert <albert@alumni.tugraz.at>
    !

    type(field_can_t), intent(inout) :: f
    double precision, intent(in) :: r, th_c, ph_c
    integer, intent(in) :: mode_secders

    double precision :: Bctr_vartheta, Bctr_varphi, bmod2, sqg, &
        d3Aphdr3, dummy, dummy3(3), dummy6(6), &
        Bth, Bph, dBth, dBph, d2Bth, d2Bph


    ! Count evaluations for profiling

    ! initialize to zero - no angular derivatives will be set due to straight field line Ath(r) Aph(r)
    f%dAth = 0d0
    f%dAph = 0d0

    ! initialize all 2nd derivatives to zero, as the mode decides, which one to use
    f%d2Ath = 0d0
    f%d2Aph = 0d0
    f%d2hth = 0d0
    f%d2hph = 0d0
    f%d2Bmod = 0d0

    call splint_boozer_coord(r, th_c, ph_c, &
        f%Ath, f%Aph, f%dAth(1), f%dAph(1), f%d2Aph(1), d3Aphdr3, &
        Bth, dBth, d2Bth, &
        Bph, dBph, d2Bph, &
        f%Bmod, f%dBmod, f%d2Bmod, dummy, dummy3, dummy6)

    bmod2 = f%Bmod**2
    sqg=(-f%dAph(1)/f%dAth(1)*Bth+Bph)/bmod2*torflux
    Bctr_vartheta = -f%dAph(1)/sqg
    Bctr_varphi = f%dAth(1)/sqg

    f%hth = Bth/f%Bmod
    f%hph = Bph/f%Bmod
    f%dhth(1) = dBth/f%Bmod - Bth*f%dBmod(1)/bmod2
    f%dhph(1) = dBph/f%Bmod - Bph*f%dBmod(1)/bmod2
    f%dhth(2:3) = -Bth*f%dBmod(2:3)/bmod2
    f%dhph(2:3) = -Bph*f%dBmod(2:3)/bmod2

    if(mode_secders > 0) then
        f%d2hth(1) = d2Bth/f%Bmod - 2d0*dBth*f%dBmod(1)/bmod2 + Bth/bmod2*(2d0*f%dBmod(1)**2/f%Bmod - f%d2Bmod(1))
        f%d2hph(1) = d2Bph/f%Bmod - 2d0*dBph*f%dBmod(1)/bmod2 + Bph/bmod2*(2d0*f%dBmod(1)**2/f%Bmod - f%d2Bmod(1))
    endif

    if(mode_secders.eq.2) then
        f%d2hth((/4,6/)) = Bth/bmod2*(2d0*f%dBmod((/2,3/))**2/f%Bmod - f%d2Bmod((/4,6/)))
        f%d2hph((/4,6/)) = Bph/bmod2*(2d0*f%dBmod((/2,3/))**2/f%Bmod - f%d2Bmod((/4,6/)))

        f%d2hth(2) =  -dBth*f%dBmod(2)/bmod2 &
        + Bth/bmod2*(2d0*f%dBmod(1)*f%dBmod(2)/f%Bmod - f%d2Bmod(2))
        f%d2hph(2) =  -dBph*f%dBmod(2)/bmod2 &
        + Bph/bmod2*(2d0*f%dBmod(1)*f%dBmod(2)/f%Bmod - f%d2Bmod(2))

        f%d2hth(3) =  -dBth*f%dBmod(3)/bmod2 &
        + Bth/bmod2*(2d0*f%dBmod(1)*f%dBmod(3)/f%Bmod - f%d2Bmod(3))
        f%d2hph(3) =  -dBph*f%dBmod(3)/bmod2 &
        + Bph/bmod2*(2d0*f%dBmod(1)*f%dBmod(3)/f%Bmod - f%d2Bmod(3))

        f%d2hth(5) = Bth/bmod2*(2d0*f%dBmod(2)*f%dBmod(3)/f%Bmod - f%d2Bmod(5))
        f%d2hph(5) = Bph/bmod2*(2d0*f%dBmod(2)*f%dBmod(3)/f%Bmod - f%d2Bmod(5))
    endif
end subroutine eval_field_booz

end module field_can_boozer
