# E: second-order GC production prerequisites (stopped early)

The controller stopped this item after a scope change: second-order GC is now
concept-only. **No generator, kernel or machine check exists in this directory.**
What follows is a desk analysis from reading the sources. None of it has been
checked in fortsym.

Sources read: the prototype at `SIMPLE-gcho/concept/gc_higher_order`, commit
598f041 (README, gen/derive_canonical.py, gen/derive_gc.f90 lines 100-200,
gen/gc_second_order.f90); SIMPLE `src/field/field_can_base.f90`,
`src/field_can.f90` (get_derivatives, get_derivatives2),
`src/field/field_can_meiss.f90`, DOC/coordinates-and-fields.md section 6;
libneo `spl_three_to_five.f90`.

## Unit mapping (READ FROM CODE, not proved)
SIMPLE code units match the prototype's units: eps = f%ro0 (= ro0/sqrt2),
mu = f%mu, u = vpar, unit mass. These follow from `pth = vpar*hth + Ath/ro0`,
`H = vpar^2/2 + mu*Bmod` and `mu = z4^2(1-lambda^2)/B`. The following
substitutions therefore apply:
  A_eff,a = A_a - ro0^2*mu*C_a,  H = vpar^2/2 + mu*B + ro0^2*H2(x,vpar,mu).
The positive-charge and orientation conventions have not been checked.

## Analysis status per claim
| Claim | Status |
|---|---|
| The H2 formula is identical in derive_canonical.py and derive_gc.f90 (alpha_sq, beta_perp, beta_par, vd, Pi1.vd); H2 = mu^2 c1(x) + mu u^2 c2(x) + u^4 c3(x) | READ FROM CODE / NOT MACHINE-CHECKED |
| Chart form with b^i = g^ij h_j and kappa_i = b^j(d_j h_i - d_i h_j) (exact form of kappa = -b x curl b): H2 uses Bmod and h to order 2 and the metric only to order 1, via Christoffel symbols, div and curl. b.curl R depends only on first derivatives, through F_jk = -b.(nabla_j b x nabla_k b) | ANALYSIS / OPEN (planned fortsym proof not run) |
| Consequence: grad H2 needs (B,h) to order 3 and g to order 2; Hess H2 needs (B,h) to order 4 and g to order 3 | ANALYSIS / OPEN |
| Representation-independent necessity of 4th derivatives: slab field B=(1+c x^4/24) z_hat gives H2 = mu^2/2[(lnB)''/2-((lnB)')^2], so d_xx H2(0) = mu^2 c/4 != 0, although the 3-jet is uniform. This formula reproduces the prototype's -3g^2mu^2/(4B^2) for linear B | HAND DERIVATION / NOT MACHINE-CHECKED |
| d_r C_a = d_r R*_a - d_a R*_r (local); d_a K = integral of d_a R*_r dr (nonlocal); C_a needs d and d2 for SIMPLE's pth Newton, hence a spline at least C^2 | ANALYSIS / OPEN |
| Whether metric order-3 jets are needed independently of the constraints (flatness, abs(b)=1, div B=0) | OPEN |
| Cross-checks against the prototype's Python/Fortran kernels | NOT STARTED |
| Generated fortsym kernel (C_a, K, H2, grad, Hess) and FD check | NOT STARTED |

## Production requirements (ANALYSIS)
- New tables:
  - six metric components g_ij (or the embedding), splined with derivatives
    to order 3;
  - C_theta and C_phi, mu-independent, with derivatives to order 2;
  - optionally K, which is needed only for phase reconstruction.
  - Jets to order 4 for Bmod, h_theta and h_phi, or precomputed c1..c3 with
    derivatives to order 2.
- Spline degree: SIMPLE's canonical tables already default to order 5.
  libneo spl_five_reg/per uses the 13+-sqrt(105) recurrence of quintic
  interpolation (C^4). Its continuity was not verified numerically. libneo
  evaluates batch splines only up to der2, so der3 and der4 evaluators are
  missing.
- Plug-in point: get_derivatives/get_derivatives2 need A -> A_eff and an
  added ro0^2*H2 term, with the chain rule through vpar. Also export dH/dmu
  and dpth/dmu at fixed (r,th,ph,pphi) for zeta_dot.
- Remaining work: the magnetic axis, where the radial frame and the
  K reference r_ref are singular; chart joins and gyrogauge patches;
  interfaces; and the D2 reconstruction, which needs third derivatives.
