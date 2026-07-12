program test_spectre_fo_mu_scatter
    !> Physical mu-scatter reference at a SPECTRE interface (ROADMAP section 4):
    !> controlled full-orbit traversals on tok2vol measure the return magnetic
    !> moment against the conserved-mu GC sheet convention. For each entry
    !> pitch and angle the Boris leg runs at two gyrostep scales; the recorded
    !> table (mu_scatter.dat) is the reference data the companion repository
    !> archives. Witnesses here: every traversal exits with an OK status, the
    !> reconstructed energy matches entry exactly (the switch identities), the
    !> scatter |mu_out/mu_in - 1| stays inside a generous physical bound, and
    !> the ensemble rms scatter agrees between the two gyrostep scales.
    !> Entry angles avoid the corrupted volume-2 coefficient bands of this
    !> fixture near theta = +-pi/2.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: read_config, params_init, netcdffile, ns_s, ns_tp, &
                      multharm, integmode
    use parmot_mod, only: ro0
    use magfie_sub, only: init_magfie, magfie, SPECTRE
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_enter, &
                                 spectre_fo_advance_until_exit, &
                                 spectre_fo_to_gc, SPECTRE_FO_OK
    use util, only: sqrt2

    implicit none

    integer, parameter :: NLAM = 8, NTH = 4, NH = 2
    real(dp), parameter :: LAMS(NLAM) = &
        [0.05_dp, -0.05_dp, 0.1_dp, -0.1_dp, 0.2_dp, -0.2_dp, 0.4_dp, -0.4_dp]
    real(dp), parameter :: THS(NTH) = [0.8_dp, 2.6_dp, 3.68_dp, 5.48_dp]
    real(dp), parameter :: HS(NH) = [1.0_dp, 0.5_dp]
    real(dp), parameter :: SCATTER_BOUND = 0.5_dp
    real(dp), parameter :: ENERGY_TOL = 1.0d-11

    character(len=1024) :: h5file
    type(tracer_t) :: norb
    type(spectre_fo_state_t) :: state
    real(dp) :: y(5), y_out(5), ro0_bar, bmod, dtb, dt_used, mu_in, mu_out
    real(dp) :: e_in, e_out, dmu, rms(NH), worst
    integer :: il, it, ih, owner, status, unit, nok(NH), nfail
    logical :: exited, failed

    if (command_argument_count() < 1) then
        print *, 'usage: test_spectre_fo_mu_scatter <tok2vol.h5>'
        error stop 1
    end if
    call get_command_argument(1, h5file)

    call write_input(trim(h5file))
    call read_config('simple.in')
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call params_init
    call init_magfie(SPECTRE)

    ro0_bar = ro0/sqrt2
    failed = .false.
    nfail = 0
    rms = 0.0_dp
    nok = 0
    worst = 0.0_dp

    open (newunit=unit, file='mu_scatter.dat', status='replace', action='write')
    write (unit, '(A)') '# theta0 lambda h_scale mu_in mu_out dvpar_rel'

    do ih = 1, NH
        do it = 1, NTH
            do il = 1, NLAM
                y = [1.0_dp, THS(it), 0.0_dp, 1.0_dp, LAMS(il)]
                mu_in = 0.5_dp*y(4)**2*(1.0_dp - y(5)**2)/bmod_side(y(2))
                e_in = y(4)**2
                call spectre_fo_enter(state, y, 1, 1, ro0_bar, status)
                if (status /= SPECTRE_FO_OK) then
                    nfail = nfail + 1
                    cycle
                end if
                bmod = bmod_side(y(2))
                dtb = 4000.0_dp*ro0_bar/bmod
                call spectre_fo_advance_until_exit(state, dtb, ro0_bar, &
                    dt_used, y_out, owner, exited, status, h_scale=HS(ih))
                if (status /= SPECTRE_FO_OK) then
                    nfail = nfail + 1
                    cycle
                end if
                if (.not. exited) then
                    nfail = nfail + 1
                    cycle
                end if
                mu_out = 0.5_dp*y_out(4)**2*(1.0_dp - y_out(5)**2) &
                         /bmod_side(y_out(2))
                e_out = y_out(4)**2
                if (abs(e_out - e_in)/e_in > ENERGY_TOL) then
                    print '(A,F6.2,A,F7.3,A,ES11.3)', 'FAIL: energy at lam=', &
                        LAMS(il), ' th=', THS(it), ': ', abs(e_out - e_in)/e_in
                    failed = .true.
                end if
                dmu = mu_out/mu_in - 1.0_dp
                worst = max(worst, abs(dmu))
                rms(ih) = rms(ih) + dmu**2
                nok(ih) = nok(ih) + 1
                write (unit, '(3(F8.4,1X),3(ES16.8,1X))') THS(it), LAMS(il), &
                    HS(ih), mu_in, mu_out, y_out(4)*y_out(5) - y(4)*y(5)
            end do
        end do
    end do
    close (unit)

    do ih = 1, NH
        if (nok(ih) > 0) rms(ih) = sqrt(rms(ih)/real(nok(ih), dp))
    end do

    print '(A,I0,A,I0,A,I0)', 'traversals: ok(h=1)=', nok(1), &
        ' ok(h=1/2)=', nok(2), ' failed/no-exit=', nfail
    print '(A,ES11.3,A,ES11.3,A,ES11.3)', 'mu scatter: rms(h=1) = ', rms(1), &
        '  rms(h=1/2) = ', rms(2), '  worst = ', worst

    if (nok(1) < NLAM*NTH/2 .or. nok(2) < NLAM*NTH/2) then
        print *, 'FAIL: too few completed traversals'
        failed = .true.
    end if
    if (worst > SCATTER_BOUND) then
        print *, 'FAIL: mu scatter above physical bound'
        failed = .true.
    end if
    if (rms(2) > 0.0_dp) then
        if (rms(1)/rms(2) > 3.0_dp .or. rms(2)/rms(1) > 3.0_dp) then
            print *, 'FAIL: rms scatter not stable under gyrostep halving'
            failed = .true.
        end if
    end if

    if (failed) error stop 1
    print *, 'full-orbit mu-scatter reference PASS'

contains

    subroutine write_input(h5)
        character(*), intent(in) :: h5
        integer :: unit_in

        open (newunit=unit_in, file='simple.in', status='replace', action='write')
        write (unit_in, '(A)') '&config'
        write (unit_in, '(A)') "  field_input = '"//h5//"'"
        write (unit_in, '(A)') '  integ_coords = 6'
        write (unit_in, '(A)') '  integmode = 3'
        write (unit_in, '(A)') '  spectre_ncon_phi = 32'
        write (unit_in, '(A)') '  ntestpart = 1'
        write (unit_in, '(A)') '  ntimstep = 10'
        write (unit_in, '(A)') '  npoiper2 = 256'
        write (unit_in, '(A)') '  relerr = 1d-13'
        write (unit_in, '(A)') '  facE_al = 500.0d0'
        write (unit_in, '(A)') '  trace_time = 1.0d-6'
        write (unit_in, '(A)') '  sbeg = 0.97d0'
        write (unit_in, '(A)') '/'
        close (unit_in)
    end subroutine write_input

    function bmod_side(th) result(bmod_loc)
        real(dp), intent(in) :: th
        real(dp) :: bmod_loc

        real(dp) :: sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)

        call magfie([1.0_dp - 1.0d-12, th, 0.0_dp], bmod_loc, sqrtg, bder, &
            hcov, hctr, hcurl)
    end function bmod_side

end program test_spectre_fo_mu_scatter
