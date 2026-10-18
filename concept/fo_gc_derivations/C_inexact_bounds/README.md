# C: a-posteriori symplecticity bound for the inexact-Newton Gauss2 step

Stopped early after the controller changed scope: Gauss2 FO is now concept-only.
This directory has only a written derivation. No code, Lean file or numerical
result was produced. Nothing outside this directory was edited. A
`make CONFIG=Fast` in the worktree was stopped part-way through. The partial,
gitignored `build/` it left can be deleted or reused.

## Solver facts (checked against the source)

| Stepper | Form | Newton | Stop rule |
|---|---|---|---|
| `concept/fo_symplectic_proof/generated_stepper.f90` (strict) | stages Y | full Newton, Jacobian rebuilt every iteration, LAPACK `dgesv` | `maxval(abs(r)) <= tol` (absolute), at most 20 iterations |
| `SIMPLE-fo-gauss2` 3ddc30b `gauss2_step` (shipped) | increments dy=Y-z | full Newton, fresh Jacobian | absolute max-abs residual, tol 1e-12 by default, at most 20 iterations |
| `reference_gauss2.f90` (frozen 89babb30) | stages Y | full Newton, fresh Jacobian | scaled backward-error floor or update < 3e-12·scale |

All three steppers evaluate the output z+h·b·F(Y_k) at the final iterate Y_k.
They all return the tangent T = I + h Σ b_i DF_i(Y_k) [A_k^{-1}(1⊗I)]_i,
where A_j = I - h(A_butcher⊗I)DF(Y_j). None of them uses simplified Newton.
In exact arithmetic the dy form gives the same iterates and the same derivative
formula as the Y form, because R_z = (A-I)(1⊗I) cancels. Only rounding and the
stopping index can differ. (Status: PROVED by hand, not machine-checked.)

## Derivation (PROVED on paper, not machine-checked)

Assume an open set of z on which the iteration count k is constant and F is C².
Write δ_j = Y_{j+1}-Y_j = -A_j^{-1}R(Y_j,z) and R_z = -(1⊗I). Differentiating
the Newton step N(Y,z) = Y - A(Y)^{-1}R(Y,z) gives the exact recursion

    DY_{j+1} = A_j^{-1}[(1⊗I) - dA_j[DY_j] δ_j],
    dA_j[w]δ_j = -h(A_butcher⊗I) blockdiag(D²F(Y_{j,i})[w_i, δ_{j,i}]).

The controller's lead is correct: N_Y = A^{-1}dA[·]A^{-1}R = O(|R(Y_{k-1})|).
Comparing with the computed-stage tangent DY^tan = A_k^{-1}(1⊗I) is sharper than
comparing with DY* on the exact branch. DY^tan is the quantity that is exactly
symplectic (Lean `Gauss2.symplectic_tangent`/`gauss2_symplectic` applied at Y_k).
Its error e_k := DY_k - DY^tan satisfies the exact identity

    e_k = A_{k-1}^{-1}[(A_k - A_{k-1}) DY^tan - C_{k-1} DY_{k-1}],
    C_{k-1} w := dA_{k-1}[w] δ_{k-1}.

Every term except C can be computed from the solver's own matrices. Because D²F
is symmetric, C_{k-1} = (A_k - A_{k-1}) + O(|δ_{k-1}|²), so to leading order

    e_k ≈ -A_{k-1}^{-1}(A_k - A_{k-1}) e_{k-1}.

The C¹ error therefore contracts by about γ‖ΔA‖ per iteration. It is a product
of the previous step size and the previous C¹ error, not O(tol). The bound is

    ‖e_k‖ ≤ γ_{k-1}( ‖A_k-A_{k-1}‖ ‖DY^tan‖ + ‖C_{k-1}‖ ‖DY_{k-1}‖ )
          ≤ γ_{k-1}² |h| ‖A_butcher‖ L_F √12 r_{k-1} (2‖DY^tan‖ + ‖e_{k-1}‖),

where γ = ‖A^{-1}‖ and L_F bounds ‖D²F‖ on the segment from Y_{k-1} to Y_k.
The factor √12 converts the max-norm residual r_{k-1} into a 2-norm.

For the one-step map, let η = ‖DΦ - T‖ ≤ ‖h[b_1DF_1, b_2DF_2]‖ ‖e_k‖. The
perturbation inequality, or the README quadratic-defect identity with
Dr = A_k e_k, gives

    ‖DΦᵀΩDΦ - Ω‖ ≤ d_T + 2‖T‖η + η²
    ‖DΦᵀΩDΦ - Ω‖ ≤ 2|h| Σ|b_i| ‖DF_i DY_i‖ ‖Dr_i‖ + (floating-point floor).

Here d_T is the floating-point defect of T itself, measured at about 1e-16 per step.

**Iteration-count boundary (PROVED on paper).** The terminated map is piecewise
smooth. Where the iteration count changes, the output jumps by
h Σ b_i DF_i δ_{k,i} + O(δ²). Since |δ_k| ≤ γ_k r_k and r_k ≈ tol at the switch,
the jump is at most about |h|‖b⊗DF‖ γ tol. A 4-point central-difference stencil
of spacing ε that straddles such a jump errs by up to 1.5·jump/ε. After N steps
this error is multiplied by the growth factor of the map.

**N steps (PROVED on paper).** Let P_n = M_n⋯M_1. The exact telescoping identity is

    P_NᵀΩP_N - Ω = Σ_n P_{n-1}ᵀ (M_nᵀΩM_n - Ω) P_{n-1},
    ‖·‖ ≤ Σ_n ‖P_{n-1}‖² d_n ≤ N G² max_n d_n,       G = max_n ‖P_n‖.

This corrects the brief in two ways. The weights are prefix products, not
suffix products, and they enter squared. The general form
‖P_{n-1}‖ ≤ Π_{m<n}‖M_m‖ gives a Gronwall-type exponential bound. No separate
recursion is needed, because the defects add exactly.

## Status per claim

| Claim | Status |
|---|---|
| Newton variant, residual norm and stop rule of all three steppers | Verified by reading the source |
| Exact derivative recursion and e_k identity; leading-order contraction | PROVED (paper only) |
| Explicit e_k, η and one-step defect bounds | PROVED (paper only), conditional on L_F; no machine checks |
| Jump size at an iteration-count switch; finite-difference noise model a·ε⁴ + b/ε | PROVED (paper) / model only |
| N-step telescoping identity and bound | PROVED (paper only) |
| Lean `InexactBound.lean`: (a) perturbation, (b) telescoping, (c) inexact-stage defect | NOT STARTED (file not created) |
| Numerical probe: logs of γ, r_{k-1}, r_k and k; spacing sweep; Richardson extrapolation; fit; bound/measured ratio | NOT STARTED (the build was interrupted, no driver written) |
| Whether the NCSX 100x 1.614e-7 defect is a finite-difference artefact | OPEN. Hint: the same value 1.6158e-7 appears for the frozen reference stepper, and the spacing disagreement of 1.18e-6 exceeds it, which is consistent with a finite-difference artefact. Not established. |

## Suggested next step if resumed

Write a driver that mirrors `generated_step` bit for bit. It would propagate
DY_j through the recursion above, getting D²F[·,δ] from a fourth-order
directional difference of DF along δ. That yields the derivative of the
terminated map with no map-level finite differences. Comparing the defect of
that derivative's product with the finite-difference sweeps would settle the
1.6e-7 question directly.
