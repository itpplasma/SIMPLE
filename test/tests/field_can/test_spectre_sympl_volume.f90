program test_spectre_sympl_volume
    !> Direct witnesses for the per-volume Meiss canonical coordinates of SPECTRE
    !> (#439). Two modes selected by the first argument:
    !>
    !>   construction <spectre.h5>  Build the per-volume charts and check that the
    !>       transformed covariant radial 1-form components vanish. A'_r (the
    !>       rho_par-independent part) and B'_r (the rho_par-linear part) are
    !>       exactly zero at the construction nodes; away from them the residual is
    !>       the spline-interpolation error, required below 5e-5 relative to the
    !>       poloidal components. This mirrors the ca/05 identities.
    !>
    !>   pphi <tok2vol.h5>  Trace passing orbits in an axisymmetric two-volume
    !>       tokamak with the canonical symplectic integrator. The toroidal
    !>       canonical momentum p_phi = m v_par h_phi + (e/c) A_phi is the fourth
    !>       symplectic coordinate; for an exactly axisymmetric field it is
    !>       conserved to machine precision, the witness that the canonical
    !>       structure (lambda AND chi gauge) is correct.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: read_config, netcdffile, ns_s, ns_tp, multharm
    use field_can_mod, only: field_can_t, field_can_init, eval_field => evaluate
    use field_can_spectre, only: spectre_volumes, spectre_mvol
    use field_spectre, only: spectre_field_t, create_spectre_field
    use interpolate, only: evaluate_batch_splines_3d_der
    use orbit_symplectic, only: symplectic_integrator_t, orbit_sympl_init, &
                                orbit_timestep_sympl, MIDPOINT
    use new_vmec_stuff_mod, only: nper
    use util, only: twopi

    implicit none

    character(len=32) :: mode
    character(len=1024) :: h5file
    character(len=256) :: config_file
    type(tracer_t) :: norb

    if (command_argument_count() < 2) then
        print *, 'usage: test_spectre_sympl_volume <construction|pphi> <file.h5>'
        error stop 1
    end if
    call get_command_argument(1, mode)
    call get_command_argument(2, h5file)

    call write_input(trim(h5file))
    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, MIDPOINT)

    select case (trim(mode))
    case ('construction')
        call check_construction(trim(h5file))
    case ('pphi')
        call check_pphi(norb)
    case default
        print *, 'unknown mode ', trim(mode)
        error stop 1
    end select

contains

    subroutine write_input(h5)
        character(*), intent(in) :: h5
        integer :: unit

        open (newunit=unit, file='simple.in', status='replace', action='write')
        write (unit, '(A)') '&config'
        write (unit, '(A)') "  field_input = '"//h5//"'"
        write (unit, '(A)') '  integ_coords = 6'
        write (unit, '(A)') '  integmode = 3'
        write (unit, '(A)') '  ntestpart = 1'
        write (unit, '(A)') '  trace_time = 1.0d-4'
        write (unit, '(A)') '  sbeg = 1.5d0'
        write (unit, '(A)') '/'
        close (unit)
    end subroutine write_input


    subroutine check_construction(h5)
        !> A'_r = A_s + A_zeta dlam/dr - dchi/dr, with A_s = 0 in the SPECTRE gauge,
        !> so the chi gauge must cancel A_zeta dlam/dr. B'_r = h_s + h_zeta dlam/dr,
        !> which the lambda transform must cancel. Both are zero at the nodes; the
        !> residual on the verification grid is the spline error.
        character(*), intent(in) :: h5

        integer, parameter :: NR = 6, NTH = 8, NPH = 8
        real(dp), parameter :: TOL = 5.0d-5
        real(dp), parameter :: R_LO = 0.1_dp, R_HI = 0.9_dp

        type(spectre_field_t) :: src
        integer :: ierr, lvol, ir, ith, iph
        real(dp) :: rmin, rmax, period, r, th, ph, frac
        real(dp) :: x(3), y(2), dy(3, 2), lam, dlam(3), dchi(3)
        real(dp) :: xref(3), Acov(3), hcov(3), Bmod
        real(dp) :: Ar_res, Br_res, Ath, Aph, Bth, Bph
        real(dp) :: maxA, maxB, max_endpoint_B, max_endpoint_h, ratioA, ratioB
        real(dp) :: hscale
        real(dp) :: field_values(5)
        real(dp) :: field_derivatives(3, 5)

        call create_spectre_field(src, h5, ierr)
        if (ierr /= 0) error stop 'check_construction: create_spectre_field failed'

        period = twopi/real(nper, dp)
        maxA = 0.0_dp
        maxB = 0.0_dp
        max_endpoint_B = 0.0_dp
        max_endpoint_h = 0.0_dp

        ! The transformed radial 1-form components A'_r and B'_r vanish at the
        ! construction nodes; on the interior verification grid the residual is the
        ! spline-derivative error. It is measured against the retained angular scale
        ! sqrt(theta^2 + phi^2) of each 1-form. Deviation from the issue's "relative
        ! to |A_theta|, |B_theta|": the covariant poloidal component reverses sign
        ! (passes through zero) inside a volume, so a pointwise |theta| denominator
        ! is a numerically singular normaliser, not a construction defect; the
        ! toroidal component keeps the retained-1-form scale finite.
        do lvol = 1, spectre_mvol
            rmin = spectre_volumes(lvol)%rmin
            rmax = spectre_volumes(lvol)%rmax
            do ir = 1, NR
                frac = R_LO + (R_HI - R_LO)*(real(ir, dp) - 0.5_dp)/real(NR, dp)
                if (ir == NR) frac = 0.95_dp
                r = rmin + (rmax - rmin)*frac
                do ith = 1, NTH
                    th = twopi*(real(ith, dp) - 0.5_dp)/real(NTH, dp)
                    do iph = 1, NPH
                        ph = period*(real(iph, dp) - 0.5_dp)/real(NPH, dp)

                        x = [r, th, ph]
                        call evaluate_batch_splines_3d_der( &
                            spectre_volumes(lvol)%spl_transform, x, y, dy)
                        lam = y(1)
                        dlam = dy(:, 1)
                        dchi = dy(:, 2)

                        xref = [r, th, ph + lam]
                        call src%evaluate(xref, Acov, hcov, Bmod)

                        Ar_res = Acov(1) + Acov(3)*dlam(1) - dchi(1)
                        Br_res = hcov(1) + hcov(3)*dlam(1)
                        Ath = Acov(2) + Acov(3)*dlam(2) - dchi(2)
                        Aph = Acov(3)*(1.0_dp + dlam(3)) - dchi(3)
                        Bth = hcov(2) + hcov(3)*dlam(2)
                        Bph = hcov(3)*(1.0_dp + dlam(3))

                        ratioA = abs(Ar_res)/hypot(Ath, Aph)
                        ratioB = abs(Br_res)/hypot(Bth, Bph)
                        maxA = max(maxA, ratioA)
                        maxB = max(maxB, ratioB)
                    end do
                end do
            end do
        end do

        do lvol = 1, spectre_mvol - 1
            r = spectre_volumes(lvol)%rmax
            do ith = 1, NTH
                th = twopi*(real(ith, dp) - 0.5_dp)/real(NTH, dp)
                ph = 0.0_dp
                x = [r, th, ph]
                call evaluate_batch_splines_3d_der( &
                    spectre_volumes(lvol)%spl_transform, x, y, dy)
                xref = [nearest(r, -1.0_dp), th, ph + y(1)]
                call src%evaluate(xref, Acov, hcov, Bmod)
                call evaluate_batch_splines_3d_der( &
                    spectre_volumes(lvol)%spl_field, x, field_values, &
                    field_derivatives)
                max_endpoint_B = max(max_endpoint_B, &
                    abs(field_values(5) - Bmod)/Bmod)
                Bth = hcov(2) + hcov(3)*dy(2, 1)
                Bph = hcov(3)*(1.0_dp + dy(3, 1))
                hscale = hypot(Bth, Bph)
                max_endpoint_h = max(max_endpoint_h, &
                    hypot(field_values(3) - Bth, field_values(4) - Bph)/hscale)
            end do
        end do

        print '(A,I0,A)', 'construction: ', spectre_mvol, ' volumes checked'
        print '(A,ES12.4)', 'construction: max |A''_r|/|A_ang| = ', maxA
        print '(A,ES12.4)', 'construction: max |B''_r|/|B_ang| = ', maxB
        print '(A,ES12.4)', 'construction: max lower endpoint |dB|/B = ', &
            max_endpoint_B
        print '(A,ES12.4)', 'construction: max lower endpoint |dh|/|h| = ', &
            max_endpoint_h
        if (maxA >= TOL .or. maxB >= TOL .or. max_endpoint_h >= TOL .or. &
            max_endpoint_B >= TOL) then
            print '(A,ES10.2)', 'FAIL: construction identity above tolerance ', TOL
            error stop 1
        end if
        print *, 'construction identity PASS'
    end subroutine check_construction


    subroutine check_pphi(self)
        !> Trace passing orbits in volume 2 (axisymmetric tok2vol); p_phi is the
        !> fourth symplectic coordinate and must stay constant to ~machine roundoff.
        type(tracer_t), intent(in) :: self

        integer, parameter :: NORBIT = 4, NSTEP = 100000
        real(dp), parameter :: MU = 1.0d-6, VPAR0 = 1.0d0
        real(dp), parameter :: DT = 1.0d0, TOL = 1.0d-10

        type(field_can_t) :: f
        type(symplectic_integrator_t) :: integ
        integer :: iorb, k, ierr
        real(dp) :: z0(4), pphi0, drift, maxdrift, th0, ph0, ro0
        real(dp) :: rmin2, rmax2, r0

        rmin2 = spectre_volumes(2)%rmin
        rmax2 = spectre_volumes(2)%rmax
        r0 = 0.5_dp*(rmin2 + rmax2)
        maxdrift = 0.0_dp

        ! ro0 = m c v0 / e (a length x field, cm*Gauss). Unit magnitude keeps the
        ! Larmor radius ro0/Bmod tiny so the passing orbit stays inside volume 2 and
        ! the Newton solve is well conditioned; p_phi = z(4) is conserved by the
        ! symplectic map for the exactly axisymmetric field regardless of ro0.
        ro0 = 1.0_dp
        call eval_field(f, r0, 0.0_dp, 0.0_dp, 0)
        print '(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)', 'pphi: Ath=', f%Ath, &
            ' Aph=', f%Aph, ' hph=', f%hph, ' Bmod=', f%Bmod

        do iorb = 1, NORBIT
            th0 = twopi*real(iorb - 1, dp)/real(NORBIT, dp)
            ph0 = 0.1_dp*real(iorb - 1, dp)

            call field_can_init(f, MU, ro0, VPAR0)
            z0 = [r0, th0, ph0, 0.0_dp]
            call eval_field(f, z0(1), z0(2), z0(3), 0)
            z0(4) = f%vpar*f%hph + f%Aph/f%ro0
            pphi0 = z0(4)

            call orbit_sympl_init(integ, f, z0, DT, 1, 1.0d-13, MIDPOINT)
            do k = 1, NSTEP
                ierr = 0
                call orbit_timestep_sympl(integ, f, ierr)
                if (ierr /= 0) then
                    print '(A,I0,A,I0,A,ES12.4)', 'pphi: step error ', ierr, &
                        ' orbit ', iorb, ' at r=', integ%z(1)
                    error stop 1
                end if
            end do

            drift = abs(integ%z(4) - pphi0)/max(abs(pphi0), 1.0d-30)
            maxdrift = max(maxdrift, drift)
            print '(A,I0,A,ES12.4,A,ES12.4)', 'pphi orbit ', iorb, &
                ': p_phi = ', pphi0, ' rel drift = ', drift
        end do

        print '(A,ES12.4)', 'pphi: max relative drift = ', maxdrift
        if (maxdrift >= TOL) then
            print '(A,ES10.2)', 'FAIL: p_phi drift above tolerance ', TOL
            error stop 1
        end if
        print *, 'p_phi conservation PASS'
    end subroutine check_pphi

end program test_spectre_sympl_volume
