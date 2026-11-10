!> Guiding-centre collision operator of a test particle a in a Maxwellian
!> field species b, derived from the Rosenbluth potentials, and identity proofs
!> for the coefficients that SIMPLE (collis_alphas.f90) and ASCOT5
!> (mccc_coefs.h) implement.
!>
!> Symbols: v test speed, vb = sqrt(2 T_b/m_b), x = v/vb, ma, mb masses,
!> cab = c_1 L = n_b q_a^2 q_b^2 L/(4 pi eps0^2) (SI; cgs: 4 pi n_b e_a^2 e_b^2 L),
!> L the Coulomb logarithm and Lp = dL/dv. SIMPLE works with p = v/v0 and
!> tau = v0 t, so its rates convert as dpp*v0**3 = D_par, dhh*v0 = nu/2 and
!> fpeff*v0**2 = K.
!>
!> Usage: derive_collision_operator [results.json]
program derive_collision_operator
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, exp, erf, sqrt, pi_expr, &
                            operator(+), operator(-), operator(*), operator(/), &
                            operator(**)
    use fortsym_diff, only: diff
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(suite_t) :: suite
    type(expr_t) :: v, vb, ma, mb, cab, nb, L, Lp, v0, x, gam, ex, sqpi
    type(expr_t) :: h, g, friction, dpar, dperp, k_ito, nu
    type(expr_t) :: mu0, mu1, dmu0, a_dpar, a_dperp, a_dd, a_q, a_k, a_nu
    type(expr_t) :: tb, zero_flux, dpar_total_dv
    type(expr_t) :: p, velrat, c_ref, efcolf, dp, dh, dpd, enrat, f, df
    type(expr_t) :: s_dpp, s_dhh, s_fpeff
    character(4096) :: out

    out = 'collision_operator_identities.json'
    if (command_argument_count() > 0) call get_command_argument(1, out)
    call arena%init()
    engine = make_symengine_engine(arena)
    call suite_begin(suite, 'guiding-centre collision operator, SIMPLE and ASCOT5')

    v = sym(arena, 'v'); vb = sym(arena, 'v_b'); ma = sym(arena, 'm_a')
    mb = sym(arena, 'm_b'); nb = sym(arena, 'n_b')
    L = sym(arena, 'L'); Lp = sym(arena, 'L_p'); v0 = sym(arena, 'v_0')
    cab = sym(arena, 'c_1')*L   ! cab is proportional to L
    x = v/vb
    sqpi = sqrt(pi_expr(arena))
    ex = exp(-(x**2))

    ! Rosenbluth potentials of a Maxwellian (h carries the mass factor) and
    ! Gamma = q_a^2 q_b^2 L/(4 pi eps0^2 m_a^2) = cab/(n_b m_a^2).
    gam = cab/(nb*ma**2)
    h = (1 + ma/mb)*nb*erf(x)/v
    g = nb*vb*((x + 1/(2*x))*erf(x) + ex/sqpi)
    friction = gam*diff(h, v)              ! <dv_par>/dt
    dpar = gam*diff(diff(g, v), v)/2       ! <dv_par^2>/(2 dt)
    dperp = gam*diff(g, v)/(2*v)           ! <dv_perp,k^2>/(2 dt), one direction
    ! Ito drift of |v|: friction plus the two perpendicular half-rates over v.
    k_ito = friction + 2*dperp/v
    nu = 2*dperp/v**2                      ! d xi = -nu xi dt + sqrt((1-xi^2) nu) dW

    ! ASCOT5 mccc_coefs.h (mufun: mu0, mu1 and mu0').
    mu0 = (erf(x) - 2*x*ex/sqpi)/x**2
    mu1 = erf(x) - mu0/2
    dmu0 = 4*ex/sqpi - 2*mu0/x
    a_dpar = cab*mu0/(2*ma**2*v)
    a_dperp = cab*mu1/(2*ma**2*v)
    a_dd = cab*(dmu0/vb - mu0/v)/(2*ma**2*v)
    a_q = -cab*mu0/(ma*mb*vb**2)
    a_k = a_q + a_dd + 2*a_dpar/v
    a_nu = 2*a_dperp/v**2

    call check_identity(suite, engine, 'ASCOT5 Dpara = Rosenbluth D_par', a_dpar - dpar)
 call check_identity(suite, engine, 'ASCOT5 Dperp = Rosenbluth D_perp', a_dperp - dperp)
    call check_identity(suite, engine, 'ASCOT5 dDpara = d D_par/dv at fixed L', &
                        a_dd - diff(dpar, v))
call check_identity(suite, engine, 'ASCOT5 K = Ito speed drift at fixed L', a_k - k_ito)
    call check_identity(suite, engine, 'ASCOT5 nu = pitch rate', a_nu - nu)

    ! Zero flux of a Maxwellian of the test species at T_b: with v-dependent L
    ! the drift needs dL/dv, the term ASCOT5 omits.
    tb = mb*vb**2/2
    dpar_total_dv = diff(dpar, v) + diff(dpar, L)*Lp
    zero_flux = dpar_total_dv + dpar*(2/v - ma*v/tb)
    call check_identity(suite, engine, 'fixed-L Ito drift = zero-flux drift at Lp=0', &
                        k_ito - (zero_flux - diff(dpar, L)*Lp))
    call check_identity(suite, engine, 'zero-flux drift - ASCOT5 K = D_par Lp/L', &
                        zero_flux - a_k - dpar*Lp/L)

    ! SIMPLE collis_alphas.f90, central branch of onseff, p = v/v0,
    ! velrat = v0/vb, enrat = E_a/T_b. efcolf*f carries L(v) via
    ! c_ref = cab/(2 m_a^2 v0^4) at the local L; df = (dL/dp)/L_ref so that
    ! efcolf*df = c_ref*velrat*Lp*v0/L.
    p = v/v0
    velrat = v0/vb
    c_ref = cab/(2*ma**2*v0**4)
    efcolf = c_ref*velrat
    f = num(arena, 1)
    df = Lp*v0/L
    dp = erf(x)/x**3 - 2*ex/(sqpi*x**2)
    dh = erf(x)*(1 - 1/(2*x**2))/x + ex/x**2/sqpi
    dpd = 4*ex/sqpi - dp
    enrat = (ma*v0**2/2)/tb
    s_dpp = f*dp*efcolf
    s_dhh = f*dh*efcolf/p**2
    s_fpeff = (f*(dpd/p - 2*dp*p*enrat) + df*dp)*efcolf

    call check_identity(suite, engine, 'SIMPLE dpp*v0^3 = D_par', s_dpp*v0**3 - dpar)
    call check_identity(suite, engine, 'SIMPLE dhh*v0 = nu/2', s_dhh*v0 - nu/2)
  call check_identity(suite, engine, 'SIMPLE fpeff*v0^2 = zero-flux drift with dL/dv', &
                        s_fpeff*v0**2 - zero_flux)
    call check_identity(suite, engine, 'SIMPLE fpeff*v0^2 - ASCOT5 K = D_par Lp/L', &
                        s_fpeff*v0**2 - a_k - dpar*Lp/L)

    ! itpplasma/ascot5 adds Dpara*dclog/clog to dDpara (fork PR #23), which
    ! makes its K the zero-flux drift and equal to SIMPLE's.
    call check_identity(suite, engine, 'ASCOT5 fork K = SIMPLE fpeff*v0^2', &
                        (a_q + a_dd + a_dpar*Lp/L + 2*a_dpar/v) - s_fpeff*v0**2)

    call suite_end(suite, trim(out))
end program derive_collision_operator
