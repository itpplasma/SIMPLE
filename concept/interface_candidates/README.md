# Charged-particle interface candidates

The executable methods are now **FO**, **C3** (carried phase, orders 0/1), and
**GL2** (resolved axisymmetric smooth-layer control). C1L0/C1L1, C2fixed/C2pphi,
C4 and C5 were removed with their settings and implementations. The exact planar
helix/cycle oracle remains; its successful fixed-model tests do not validate the
removed curved-field patch refresh. See the
[rejection ledger](../../DOC/spectre-rejected-interface-models.md#implementation-removal-ledger).

The tables and model descriptions below are the **historical comparison** from
`7ebcf4c`, including failed methods; they are not a list of current CLI options.
Their data are retained as evidence. Fresh `run_study.py` output goes to
`build/study/`, leaving those records intact. The historical report tools read
`results/`; they do not run or restore retired methods.

Standalone experiment on branch `feat/spectre-gc-fo-continuous`; no SIMPLE
source or build changes. Units are positive charge = mass = 1, no electric
field. This study compares **continuous formulations** using the same adaptive
Dormand–Prince RK45 except for the explicitly identified helix, Boris and
Gauss2 controls. A good loss count alone does not validate a formulation.

The preferred tested hybrid is **C3/T1: carried phase with one energy-matched
lift and its Newton inverse**. Its local errors have the expected quadratic
scaling, and its excursion distributions agree closely with full orbits.
At switch distances of 3 and 4 local gyroradii it spends 94–97% of the grid in full
orbits; measured wall costs range from 0.76 to 2.04 times FO on this shared
host, so no consistent speed gain is established. Full orbit remains the
bounded-time baseline. Resolved Boris and Gauss2 layers are secondary controls.

Long strong T-B has converged loss counts but a pointwise reversal limit
even for the tight full-orbit reference. Two C4 trajectories fail near the
first-order reconstruction fails for almost parallel velocity at very small mu.
Order 2 remains
unavailable. These limits are included in the measured tables and verdicts.

## Reproduce

Run from the repository root; gfortran, Python 3, NumPy and SymPy are required.

```sh
make -f concept/interface_candidates/Makefile test
OMP_NUM_THREADS=16 make -f concept/interface_candidates/Makefile accuracy
OMP_NUM_THREADS=16 make -f concept/interface_candidates/Makefile study
make -f concept/interface_candidates/Makefile summarize
```

The matrix takes tens of minutes on this shared host. `run_study.py --resume`
resumes complete configurations of the **same source revision**, checking their
parameters. Rebuild and restart after a physics change; there is no automatic
source-revision cache validation. Thread counts above 16 are rejected by the
ensemble executable. All computation is local. The copied older hybrid driver has been removed; the maintained C3 driver is
the executable hybrid comparison.

The original experiment copied `sheet_field.f90`, `rk_reference.f90`,
`fo_sheet.f90`, `gc_hybrid.f90`, `study_hybrid.f90` and `sheet_setup.f90` from
[the original prototype](../symplectic_sheet/README.md), repository commit
`a7788ee`. `candidate_field`, `transforms`, `candidate_orbits` and `compare`
implement this study's 3-D formulations. `proofs/generate_torsion.py` generates
`torsion_kernel.inc`; the emitted gradient has an independent finite-difference
oracle. Results and their complete parameters are in `results/main_manifest.json`
and `results/long_manifest.json`; marker CSVs are retained losslessly as gzip.

## Analytic fields and ensemble

Write cylindrical position as `(R,varphi,Z)`, `s=R-R0`, `r²=s²+Z²`, with
`R0=10`, sheet `r0=0.5`, `F=10`, `a_-=0.25`. In T-A,

```
A_R=0, A_varphi=psi/R, A_Z=-F log R,
psi_-=a_- r², psi_+=a_+ r²+(a_--a_+) r0²,
B=(-2a Z/R, F/R, 2a s/R).
```

`a_+=0.30,0.75,2.00` gives weak, moderate and strong sheets. At the outer
midplane their direction jumps are approximately 0.0050, 0.0499, 0.1724 rad;
the magnitude ratios are 1.00014, 1.00250, 1.01949. These are analytic sheet
models, not force-balanced MRxMHD equilibria or numerical SPECTRE fixtures.
The sheet has both a direction and a magnitude jump. `p_varphi=R v_varphi+psi`
is an exact invariant across it because A is continuous and axisymmetric.

T-B adds `delta A=delta H d(r²/2)`,
`H=s cos(n varphi)+Z sin(n varphi)`, `delta=0.1`, `n=3`. Setting
`K=-s sin(n varphi)+Z cos(n varphi)` gives

```
delta B=(delta n Z K/R, -delta K, -delta n s K/R).
```

The perturbation is smooth for `R>0`, its potential has the same trace on
both sides, `div B=0`, and `B.grad(r)=0` on **every** minor-radius surface.
The poloidal dependence and toroidal helicity remove the cyclic-coordinate
invariant: T-B's changing p_varphi is not counted as an error. SymPy checks
curl, divergence, tangency and sheet continuity exactly. The construction
is an original analytic test, not an approximation imported from SPECTRE.

There are 256 equal-weight markers. A fixed Park–Miller seed 20261014 generates
uniform poloidal angle, toroidal angle, velocity gyrophase and isotropic pitch.
All methods use the same physical initial `(x,v)` on `|v|=v` at
`r=r0-0.35 v`, with `v=0.01,0.02,0.04`. Define nominal `rho=v/Bref=v`,
`L=r0=0.5`, hence `epsilon=0.02,0.04,0.08`; the local smooth-field ordering
also depends on direction gradients and the toroidal major radius. Switching
uses **local** `v/B`, including parallel velocity in its conservative radius
bound, rather than the nominal radius or just perpendicular gyroradius.

The main horizon is `T=80/v`, i.e. 8000, 4000, 2000 (about 1273, 637,
318 gyroperiods at Bref). Loss means reaching the **fixed** radius 0.55 by T.
Tracing continues after a loss to collect invariant and excursion errors.
Radius is observed every 0.5 time units; refinements to 0.25 and 0.125 test
sampling bias. Radial excursion is `max_t r-min_t r`, including the start.
Additional return-orbit runs use T=16000 at v=.02 in moderate T-A and strong
T-B. This loss surface is a diagnostic absorbing surface, not a material wall.

Per-marker CSV data include maximum and final relative energy error, maximum
and final absolute p_varphi error, radial range, loss, event/switch counts,
failures and full-orbit fraction. Tables normalize momentum errors by `R0 v`,
position/reversal errors by nominal rho, velocity errors by v. Thirty-two
markers undergo the entire forward integration followed by **negative-time**
integration with the same B and retained solver mode. Reversing v alone in a
fixed B is not time reversal; that experiment would also require B -> -B.
For Gauss2 only the reduced `(R,Z,p_R,p_Z)` return is compared.

## Equations and switches

The physical reference is `xdot=v`, `vdot=v cross B(x)`; x and v remain
continuous at a bounded sheet. Each RK stage uses one smooth branch. A crossing
is landed by a scalar root of the continued branch trajectory and restarted
on the other branch, with no artificial impulse. Endpoint, midpoint and radial
extremum brackets resolve crossing-return pairs, including a dedicated shallow
pair regression. Landing time brackets are below 1e-11. This is numerical
event resolution, not an exact-arithmetic assertion or a uniform grazing proof.
The reference is tested at tolerances 1e-7/1e-9 everywhere and 1e-11 on the
strong sheet. At long strong T-B, loss=87/256 at all three tolerances, but
velocity return RMS is 0.219 at 1e-9 and 0.0231 at 1e-11. This limits
pointwise accuracy claims; converged counts do not certify individual paths. A 3-D trajectory oracle agrees with the copied reduced RK solver.

For smooth-volume GC dynamics define `b=B/B`, `kappa=(b.grad)b`,
`tau=b.curl b`, `R_i=e2.dot(partial_i e1)`, `e1=(b cross ez)/|b cross ez|`,
`e2=b cross e1`. Our phase is the angle of velocity
`w=cos(phi)e1+sin(phi)e2`, so uniform `phi_dot=-B`.
For T0 use `Rstar=R`; for T1 use Littlejohn's parallel symplectic correction
`Rstar=R+tau b/2` (perpendicular symplectic correction set to zero). Then

```
Gamma=(A+u b-mu Rstar).dX - mu dphi, Hgc=u²/2+mu B,
Bstar=B+u curl b-mu curl Rstar, Bpar=b.Bstar,
Xdot=(u Bstar+mu b cross grad B)/Bpar,
udot=-mu Bstar.grad B/Bpar, mu_dot=0,
phi_dot=-B-Rstar.Xdot.
```

T0 drops `mu curl Rstar` in Bstar, an order-two completion; this is an
order-zero comparison, not an exact six-dimensional Hamiltonian claim.
T1 includes it. The energy identity follows by dotting the equations into
`(mu grad B,u,B,0)`; the two cross terms cancel. Bpar must be nonzero and
the GC ordering applies only inside each smooth volume, away from the sheet.

* **C1L0/C1L1, negative controls.** GC everywhere; mu stays fixed. Level 0
  sets `u_+=sign(u_-) sqrt(2(H_--mu B_+))`. A forbidden radicand mirrors u
  on the home side. If that reset still points out of its assigned volume,
  the reset/flow contract fails and the marker is flagged, not counted as a
  numerical loss. Level 1 solves energy along the two-side averaged GC
  interface-generator direction, with midpoint refresh and a bounded
  tangential kick; unsuccessful solves fall back to Level 0. These are
  Cartesian analogues of the old rules, not bitwise ports of SIMPLE's chart
  code. The later special equal-B relocation option is not included.
* **C2fixed/C2pphi.** FO near the sheet and GC outside. At exit, use the
  prototype's radius shift and an energy-matched mu; discard phase. In T-A
  project u to the GC momentum root, as SIMPLE level 3 does; in T-B use the
  parallel projection. At entry reconstruct phase zero, or the nearest
  p_varphi root in T-A, with phase-zero fallback recorded on root failure.
  GC phase advanced for diagnostic reconstruction is reset at each entry:
  it carries no particle phase memory. This reproduces the mechanism, not
  SIMPLE's coordinates, interpolation or complete implementation.
* **C3.** Carry phase throughout. Use **one** T_n and its Newton inverse
  `P_n=T_n^{-1}`; switches preserve kinetic energy to the Newton tolerance.
  Enter FO at `|d(X)|=d_in v/B(X)`; leave when the inverse GC has
  `|d(X)|>d_out v/B(X)`, with `d_out=d_in+1` and scans 2,3,5,8.
  FO-to-GC exit is checked after accepted FO steps; entry is root-located.
  Hysteresis avoids chatter but does not make a switch reversible.
* **C4.** C3/T1 outside, analytic helices inside, in patches with frozen
  tangent plane and both continued fields at the patch start. Scan maximum
  `Omega dt=1,0.5,0.25,0.125`; region changes are analytic helix-plane roots,
  and repeated crossings within one patch are retained. Complete uniform-sheet
  cycles are compressed exactly, with explicit propagation of the remainder. Refreshing patches
  is an approximation to curved fields, not the exact planar map.
* **C5.** Symmetric drift–Boris-rotation–drift on the copied compact C2
  potential layer. Scan `w/rho=0.2,0.1,0.05` and `dt=0.025,0.0125,0.00625`.
  Here `v dt/w=dt/(w/rho)` is the resolution parameter. Exact speed does not
  certify p_varphi or the sharp-field limit; both width and dt must converge.

## The carried-phase transform

T0 is `x=X+rho0`, `v=u b+sqrt(2mu B) w`, `rho0=b cross v/B`.
For T1 set `rhat=b cross w`, `a1=-(rhat w+w rhat)/2`,
`a2=(w w-rhat rhat)/4`, `alpha_j=a_j:grad b`. The first generators are

```
G1u=-u rho0.kappa+mu(tau+alpha1),
G1mu=rho0.(mu grad B+u² kappa)/B-u mu(tau+alpha1)/B,
G1zeta=-rho0.R+u alpha2/B
        +sqrt(2mu/B) w.(grad B+u² kappa/(2mu))/B.
```

`zeta=-phi-pi/2`. T1 uses `x=X+rho0+rho1`, with
`rho1=-G2x-(G1.d rho0)/2`; `G2x` is the spatial generator with
`Pi1=-mu tau b/2`. Its explicit contractions are in `transforms.f90`.
Local particle variables are `u_p=u-G1u`, `mu_p=mu-G1mu`,
`phi_p=phi+G1zeta`; reconstruct v using B and the frame at x. The identity
`u G1u+B G1mu=mu rho0.grad B` cancels the first-order energy difference.
An explicit order-two completion rescales v to exactly Hgc. Newton in
`(X,u,log mu,phi)` inverts this **completed** lift, rather than mixing an
unrelated projection with an independently truncated inverse. Failed Newton
inversions keep FO active and are recorded; no clipping of a negative mu.
The measured maximum numerical switch-energy residual in the short grid
is 5.46e-11 for T0 and 4.76e-11 for T1, set by the Newton tolerance; the
completed lift itself matches energy algebraically.
A GC state whose local particle mu becomes nonpositive in T1 is a failed
reconstruction, not a valid full-orbit exit. Two long C4/T-B markers encounter
this near-parallel chart limit; their incomplete distributions are censored.

These formulas are the Littlejohn first-order symplectic representation;
[Tronko–Brizard (2015), equations 29,41–45,55](https://arxiv.org/html/1511.01005)
provide the generator convention. The momentum/phase expressions also appear
in [Hirvijoki et al. (2015), Appendix A](https://arxiv.org/html/1412.1966).
The energy completion is this experiment's choice, not a quoted theorem.

## What is proved, measured, or unresolved

**Exact algebra checked by SymPy:** T-B's field identities; first-order energy
cancellation; helix Lorentz equations/speed/plane-root identity; canonical
saltation `S^T J S=J` for a continuous-A planar event; and the T1 two-form and
energy pullback through first order in a nonuniform linear slab. This slab
witness is not an exhaustive symbolic proof for arbitrary fields. In the
uniform planar sheet, exact piecewise helices composed at true event times
are exactly symplectic at a fixed physical endpoint, away from grazing.
The boundary-to-boundary map preserves the induced section form, not an
unqualified six-dimensional canonical matrix at variable event times.

**Measured checks:** inverse/lift round trip; switch energy; 3-D/reduced FO
agreement; generated curl Rstar against finite differences; planar repeated
crossings and inverse; asymptotic two-form scaling in T-B. The scaled pullback
defects for v=.01/.02/.04 are T0 `1.255e-3/2.510e-3/5.021e-3` and T1
`5.365e-6/2.139e-5/8.495e-5`, consistent with O(epsilon)/O(epsilon²).

A separate smooth-segment study uses 256 phases/pitches at fixed **gyrotime**
T=50, identical physical initial conditions, and no sheet. Normalized position
exponents are T-A `0.996/2.002` and T-B `0.998/2.006` for T0/T1; velocity
exponents are `0.999/2.039` and `1.002/1.987`. See `results/accuracy.csv`.
The near-sheet/fixed-phase experiments are in `results/sheet_accuracy.csv`;
its T-B p_varphi change is physical variation, not an invariant error.
A power fit is descriptive over these three scales, not a proof or a uniform
long-time estimate. Main runs instead hold a transit-scale horizon fixed:
missing order-two gyrofrequency can accumulate phase error over O(1/epsilon)
gyroperiods. Local transform order does not require a chaotic loss statistic,
or its difference from a 256-marker sample, to obey the same exponent.

**Analytic estimate, conditional:** a locally frozen field makes relative
velocity errors O(rho/Lsmooth) over a gyrotime; curvature also displaces the
plane by O(rho²/Lsheet). Constants are not uniform at grazing, long residence,
or many repeated crossings. A tangent uniform field has no curvature-driven
normal GC migration, so an entire frozen planar “passage” may never exit.
No general full-passage or long-time C4 error bound is asserted. State-dependent
patch refresh destroys the exact planar symplectic guarantee. Complete repeated uniform-sheet cycles are compressed after two excursions;
the 128-root guard applies to unresolved explicit searches, not to the number
of cycles represented. A 2000-time planar test agrees with 8000 short direct
event steps and has more than 128 crossings.

The completed parallel FLR note motivated this exact **cycle compression**:
in a stationary uniform tangent sheet the tangential canonical momenta fix
the same tangential velocity at every sheet crossing, while normal velocity
alternates sign. A complete cycle therefore repeats velocity and adds a
fixed tangential displacement; N cycles are skipped analytically in the implemented solver. This
is an exact acceleration of the planar oracle, not a new curvature model.
A sheet-cycle slow Hamiltonian or a canonical smooth-residual correction
remains a research candidate; this study does not label it viable on T-A/T-B.

**Step 2:** the FO Hamiltonian has a natural canonical Gauss2 discretization;
Poincare Gauss2 is tested independently on a resolved smooth layer in T-A.
T1's smooth GC one-form supports a noncanonical variational method or a local
Darboux chart with Gauss2. Its finite-order switches have a residual form defect:
energy matching and symplectic segment solvers do not make the whole hybrid
exactly symplectic. C1's resets have the physical fixed-mu obstruction and no
repair from a finer RK tolerance. Boris is volume/speed preserving, generally
not canonical symplectic. See the
[rejected-model proofs](../../DOC/spectre-rejected-interface-models.md)
for the phase-dependent finite mu jump and event-composition criteria.

**Order 2 pending:** Task D committed `c9f23b1` in SIMPLE-gcho; its API refuses
n=2 and its n=1 velocity-order acceptance failed at fixed transit time.
It uses a different symplectic completion and time/epsilon scaling. That
result cannot be silently substituted for this study's tested lift/RHS.
`transforms.f90` exposes serial registration of matched T2/RHS/energy callbacks;
order 2 fails explicitly until all compatible kernels exist. Registration
requires the stated q,m, frame, energy and phase conventions. No order-2
results or arbitrary-order claim are supplied.

## Statistical interpretation

Loss intervals are 95% Wilson binomial intervals for the specified pseudorandom
launch distribution. Paired loss-difference intervals combine exact 97.5% Clopper–Pearson
intervals for the two discordance probabilities by Bonferroni, giving at least
95% simultaneous coverage. This remains nonzero when no discordance is seen.
Paired 4000-resample bootstrap differences and quantile intervals are also
retained; a degenerate bootstrap interval alone does not certify zero bias. No multiple-testing
correction is claimed. Failures bound loss by `[observed/N,(observed+failed)/N]`;
a Wilson interval alone is not sufficient for those rows. Distribution metrics
are omitted for failed trajectories, rather than reporting a censored ensemble
as complete. `W1` is the mean absolute difference of sorted excursions; the
summary also retains Kolmogorov distance, paired excursion errors and quantiles.
Confidence intervals measure sampling uncertainty, not model bias or a proof
of equivalence. A 256-marker experiment cannot exclude a few-percent loss bias.

Energy columns are RMS of each marker's maximum relative error; momentum
columns are RMS maximum absolute errors divided by R0 v. Switch energy is
separately recorded. Forward wall cost is divided by marker count and simulated
T and compared with FO at the **same** 1e-9 tolerance; reversing is separately
timed. The host is shared and no dedicated timing/scaling certification is
claimed. Failed rows divide cost by the intended horizon and are censored
work estimates, not a speed advantage over complete trajectories.
GL2 reduction fixes p_varphi by construction and does not integrate
the ignorable azimuth; its reduced return norm is not the Cartesian return norm.

<!-- MEASURED_RESULTS -->
## Measured comparison

Completed 373 configurations × 256 markers. All tables below are generated by `render_report.py` from the marker records.

The complete grid, confidence intervals and distribution metrics are in [main_summary.csv](results/main_summary.csv) and [main_summary.json](results/main_summary.json).

### Full-orbit statistical convergence

Counts are at tolerances 1e-7 / 1e-9 / 1e-11 (the last on strong sheets). The interval and excursion quantiles refer to the tightest reference. This grid uses observation dt=0.5; the finer observer below changes strong T-A from 11 to 12 losses.

| Test | a+ | v | Lost counts | Loss, 95% CI | Excursion q10/q50/q90 |
|---|---:|---:|---|---|---|
| T-A | 0.3 | 0.01 | 8 / 8 | 0.031 [0.016,0.060] | 0.0335/0.0459/0.0577 |
| T-A | 0.3 | 0.02 | 94 / 94 | 0.367 [0.311,0.428] | 0.0648/0.0913/0.1149 |
| T-A | 0.3 | 0.04 | 156 / 156 | 0.609 [0.548,0.667] | 0.1322/0.1826/0.2298 |
| T-A | 0.75 | 0.01 | 1 / 1 | 0.004 [0.001,0.022] | 0.0282/0.0409/0.0544 |
| T-A | 0.75 | 0.02 | 41 / 41 | 0.160 [0.120,0.210] | 0.0563/0.0815/0.1076 |
| T-A | 0.75 | 0.04 | 142 / 142 | 0.555 [0.493,0.614] | 0.1135/0.1618/0.2132 |
| T-A | 2 | 0.01 | 0 / 0 / 0 | 0.000 [0.000,0.015] | 0.0220/0.0345/0.0528 |
| T-A | 2 | 0.02 | 11 / 11 / 11 | 0.043 [0.024,0.075] | 0.0439/0.0696/0.1043 |
| T-A | 2 | 0.04 | 91 / 91 / 91 | 0.355 [0.299,0.416] | 0.0868/0.1362/0.2066 |
| T-B | 0.3 | 0.01 | 14 / 14 | 0.055 [0.033,0.090] | 0.0352/0.0461/0.0595 |
| T-B | 0.3 | 0.02 | 100 / 100 | 0.391 [0.333,0.452] | 0.0706/0.0916/0.1180 |
| T-B | 0.3 | 0.04 | 161 / 161 | 0.629 [0.568,0.686] | 0.1365/0.1827/0.2355 |
| T-B | 0.75 | 0.01 | 10 / 10 | 0.039 [0.021,0.070] | 0.0300/0.0433/0.0573 |
| T-B | 0.75 | 0.02 | 55 / 55 | 0.215 [0.169,0.269] | 0.0594/0.0861/0.1149 |
| T-B | 0.75 | 0.04 | 145 / 145 | 0.566 [0.505,0.626] | 0.1153/0.1681/0.2286 |
| T-B | 2 | 0.01 | 3 / 3 / 3 | 0.012 [0.004,0.034] | 0.0254/0.0391/0.0552 |
| T-B | 2 | 0.02 | 39 / 39 / 39 | 0.152 [0.113,0.201] | 0.0506/0.0779/0.1099 |
| T-B | 2 | 0.04 | 109 / 109 / 109 | 0.426 [0.367,0.487] | 0.0973/0.1517/0.2155 |

### Candidate comparison: moderate sheet, v=0.02, T=4000

E = RMS marker maximum relative energy error; P = RMS maximum p_varphi error/(R0 v); W1/rho compares the excursion distributions. Rev = RMS velocity return error/v over valid reversals. Cost excludes reversals; FO share identifies bands retaining FO.

| Test/method | Lost/256 [95% CI] | Failed | E | P | W1/rho | Rev | Cost/FO | FO share |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| A: FO 1e-09 | 41 [0.120,0.210] | 0 | 1.39e-06 | 1.33e-09 | 0.00e+00 | 2.78e-05 | 1.00 | 1.000 |
| A: C1L0 | 35 [0.100,0.184] | 9 | 4.07e-12 | 3.02e-02 | — | 1.20e-11 | 0.25 | 0.000 |
| A: C1L1 | 36 [0.103,0.189] | 9 | 3.52e-11 | 2.85e-02 | — | 5.04e-02 | 0.25 | 0.000 |
| A: C2fixed | 41 [0.120,0.210] | 0 | 1.36e-06 | 4.15e-04 | 5.99e-04 | 2.98e-01 | 1.01 | 0.967 |
| A: C2pphi | 41 [0.120,0.210] | 0 | 1.36e-06 | 4.15e-04 | 5.99e-04 | 3.13e-01 | 0.93 | 0.967 |
| A: C3/T0 d=3 | 41 [0.120,0.210] | 0 | 1.36e-06 | 5.24e-04 | 5.20e-04 | 6.61e-02 | 0.90 | 0.967 |
| A: C3/T1 d=3 | 41 [0.120,0.210] | 0 | 1.36e-06 | 5.39e-07 | 3.86e-06 | 1.57e-04 | 1.15 | 0.967 |
| A: C4/T1 d=3 h=0.25 | 40 [0.117,0.206] | 0 | 1.75e-12 | 6.03e-01 | 9.64e-02 | 1.32e+00 | 0.11 | 0.965 |
| A: C5 w/rho=0.2 h=0.025 | 41 [0.120,0.210] | 0 | 4.04e-14 | 6.78e-05 | 3.79e-03 | 6.85e-11 | 0.50 | 1.000 |
| B: FO 1e-09 | 55 [0.169,0.269] | 0 | 1.39e-06 | — | 0.00e+00 | 9.27e-05 | 1.00 | 1.000 |
| B: C1L0 | 50 [0.151,0.248] | 7 | 2.42e-12 | — | — | 5.84e-12 | 0.24 | 0.000 |
| B: C1L1 | 51 [0.155,0.252] | 7 | 2.69e-11 | — | — | 2.76e-01 | 0.24 | 0.000 |
| B: C2fixed | 55 [0.169,0.269] | 0 | 1.33e-06 | — | 9.68e-04 | 5.45e-01 | 1.20 | 0.951 |
| B: C3/T0 d=3 | 55 [0.169,0.269] | 0 | 1.33e-06 | — | 6.13e-04 | 1.20e-01 | 1.16 | 0.951 |
| B: C3/T1 d=3 | 55 [0.169,0.269] | 0 | 1.33e-06 | — | 4.89e-05 | 2.56e-03 | 2.04 | 0.951 |
| B: C4/T1 d=3 h=0.25 | 53 [0.162,0.261] | 0 | 5.79e-12 | — | 1.26e-01 | 1.38e+00 | 0.20 | 0.939 |
| B: C5 w/rho=0.2 h=0.025 | 55 [0.169,0.269] | 0 | 4.28e-14 | — | 3.16e-03 | 1.44e-10 | 0.69 | 1.000 |

### Candidate comparison: strong sheet, v=0.02, T=4000

E = RMS marker maximum relative energy error; P = RMS maximum p_varphi error/(R0 v); W1/rho compares the excursion distributions. Rev = RMS velocity return error/v over valid reversals. Cost excludes reversals; FO share identifies bands retaining FO.

| Test/method | Lost/256 [95% CI] | Failed | E | P | W1/rho | Rev | Cost/FO | FO share |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| A: FO 1e-09 | 11 [0.024,0.075] | 0 | 1.41e-06 | 2.72e-09 | 2.33e-06 | 1.07e-04 | 1.00 | 1.000 |
| A: C1L0 | 10 [0.021,0.070] | 22 | 4.04e-12 | 1.06e-01 | — | 3.46e-11 | 0.16 | 0.000 |
| A: C1L1 | 10 [0.021,0.070] | 22 | 3.65e-11 | 1.02e-01 | — | 2.60e-01 | 0.18 | 0.000 |
| A: C2fixed | 11 [0.024,0.075] | 0 | 1.38e-06 | 3.97e-04 | 5.72e-04 | 3.14e-01 | 0.74 | 0.968 |
| A: C2pphi | 11 [0.024,0.075] | 0 | 1.38e-06 | 3.97e-04 | 5.72e-04 | 2.98e-01 | 0.81 | 0.968 |
| A: C3/T0 d=3 | 11 [0.024,0.075] | 0 | 1.38e-06 | 5.11e-04 | 4.32e-04 | 6.28e-02 | 0.75 | 0.968 |
| A: C3/T1 d=3 | 11 [0.024,0.075] | 0 | 1.38e-06 | 5.43e-07 | 5.90e-06 | 1.72e-04 | 0.76 | 0.968 |
| A: C4/T1 d=3 h=0.25 | 3 [0.004,0.034] | 0 | 1.64e-12 | 8.84e-01 | 1.41e-01 | 1.27e+00 | 0.12 | 0.966 |
| A: C5 w/rho=0.2 h=0.025 | 12 [0.027,0.080] | 0 | 4.32e-14 | 2.95e-04 | 3.32e-03 | 1.05e-09 | 0.30 | 1.000 |
| A: GL2 w/rho=0.05 h=0.025 | 11 [0.024,0.075] | 0 | 1.01e-07 | 1.81e-16 | 6.87e-04 | 7.69e-06 | 1.93 | 1.000 |
| B: FO 1e-09 | 39 [0.113,0.201] | 0 | 1.40e-06 | — | 3.52e-06 | 2.10e-04 | 1.00 | 1.000 |
| B: C1L0 | 36 [0.103,0.189] | 21 | 2.49e-12 | — | — | 3.80e-10 | 1.08 | 0.000 |
| B: C1L1 | 36 [0.103,0.189] | 20 | 3.45e-11 | — | — | 3.61e-01 | 0.23 | 0.000 |
| B: C2fixed | 39 [0.113,0.201] | 0 | 1.36e-06 | — | 6.61e-04 | 5.93e-01 | 1.04 | 0.957 |
| B: C3/T0 d=3 | 39 [0.113,0.201] | 0 | 1.36e-06 | — | 5.86e-04 | 2.12e-01 | 1.06 | 0.957 |
| B: C3/T1 d=3 | 39 [0.113,0.201] | 0 | 1.36e-06 | — | 3.00e-05 | 2.78e-03 | 1.33 | 0.957 |
| B: C4/T1 d=3 h=0.25 | 21 [0.054,0.122] | 0 | 5.55e-12 | — | 3.49e-01 | 1.24e+00 | 0.18 | 0.944 |
| B: C5 w/rho=0.2 h=0.025 | 39 [0.113,0.201] | 0 | 4.46e-14 | — | 5.06e-03 | 9.04e-09 | 0.59 | 1.000 |

### Longer return-orbit test: v=0.02, T=16000

The earlier grid's short transit horizon can hide phase resets because many GC exits have no subsequent FO entry. These runs test returns; mean switches and FO residence are included.

| Test/method | Lost/256 | Delta loss [paired 95% CI] | W1/rho | P | Rev(v) | Cost | Switches | FO share |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| A: FO 1e-09 | 54 | 0.0000 [-0.0170,0.0170] | 1.75e-05 | 5.01e-09 | 8.27e-04 | 1.00 | 0.00 | 1.000 |
| A: C2fixed | 54 | 0.0000 [-0.0170,0.0170] | 2.48e-03 | 6.57e-04 | 7.57e-01 | 0.99 | 0.73 | 0.843 |
| A: C2pphi | 54 | 0.0000 [-0.0170,0.0170] | 2.54e-03 | 6.57e-04 | 7.41e-01 | 1.00 | 0.73 | 0.843 |
| A: C3/T0 d=3 | 54 | 0.0000 [-0.0170,0.0170] | 1.91e-03 | 7.98e-04 | 4.21e-01 | 1.04 | 0.73 | 0.843 |
| A: C3/T1 d=3 | 54 | 0.0000 [-0.0170,0.0170] | 6.31e-05 | 9.42e-07 | 7.26e-04 | 1.20 | 0.73 | 0.843 |
| A: C4/T1 d=3 h=0.125 | 70 | 0.0625 [-0.0387,0.1617] | 3.20e-01 | 8.51e-01 | 1.33e+00 | 0.44 | 0.91 | 0.846 |
| A: C5 w/rho=0.05 h=0.00625 | 54 | 0.0000 [-0.0170,0.0170] | 7.86e-04 | 3.04e-05 | 6.93e-09 | 2.18 | 0.00 | 1.000 |
| B: FO 1e-09 | 87 | 0.0000 [-0.0170,0.0170] | 4.87e-04 | — | 2.19e-01 | 1.00 | 0.00 | 1.000 |
| B: C2fixed | 87 | 0.0000 [-0.0170,0.0170] | 4.08e-03 | — | 7.41e-01 | 0.97 | 0.86 | 0.817 |
| B: C3/T0 d=3 | 87 | 0.0000 [-0.0170,0.0170] | 4.07e-03 | — | 5.27e-01 | 0.93 | 0.86 | 0.817 |
| B: C3/T1 d=3 | 87 | 0.0000 [-0.0170,0.0170] | 8.72e-04 | — | 2.20e-01 | 1.38 | 0.85 | 0.817 |
| B: C4/T1 d=3 h=0.125 | 86 (+2 failed) | — | — | — | 1.40e+00 | 0.45 | 1.04 | 0.832 |
| B: C5 w/rho=0.05 h=0.00625 | 87 | 0.0000 [-0.0170,0.0170] | 1.50e-02 | — | 1.72e-04 | 2.46 | 0.00 | 1.000 |

### Excursion uncertainty in the longer return tests

Medians and paired mean differences include 95% bootstrap intervals. Full q10/q50/q90 intervals and Kolmogorov distances are in the JSON summary.

| Test/method | Median excursion [95% CI] | Mean excursion difference from FO [95% CI] |
|---|---|---|
| A: FO 1e-09 | 0.09709 [0.09254,0.10379] | -2.92e-07 [-3.86e-07,-2.16e-07] |
| A: C2fixed | 0.09709 [0.09254,0.10379] | 1.36e-05 [-5.34e-06,3.29e-05] |
| A: C2pphi | 0.09709 [0.09254,0.10379] | 1.29e-05 [-6.06e-06,3.23e-05] |
| A: C3/T1 d=3 | 0.09709 [0.09254,0.10379] | -3.23e-07 [-7.52e-07,8.51e-08] |
| A: C4/T1 d=3 h=0.125 | 0.10526 [0.10128,0.11072] | 2.15e-03 [-6.77e-03,1.13e-02] |
| A: C5 w/rho=0.05 h=0.00625 | 0.09709 [0.09254,0.10379] | -1.88e-06 [-5.95e-06,2.37e-06] |
| B: FO 1e-09 | 0.12011 [0.11207,0.12799] | 6.50e-06 [-7.08e-07,1.96e-05] |
| B: C2fixed | 0.12011 [0.11207,0.12795] | 2.96e-05 [-4.47e-06,7.82e-05] |
| B: C3/T1 d=3 | 0.12011 [0.11207,0.12800] | 5.84e-06 [-3.35e-06,2.02e-05] |
| B: C5 w/rho=0.05 h=0.00625 | 0.12021 [0.11113,0.12837] | -1.14e-04 [-3.39e-04,8.78e-05] |

### Radius-observer convergence: strong sheet, v=0.02, T=4000

These rows test the loss observer separately from ODE tolerance. Excursion medians are physical minor-radius ranges.

| Test | Observation dt | FO lost | C3/T1 lost | FO median excursion | C3/T1 median excursion |
|---|---:|---:|---:|---:|---:|
| A | 0.5 | 11 | 11 | 0.069585 | 0.069585 |
| A | 0.25 | 11 | 11 | 0.069589 | 0.069589 |
| A | 0.125 | 12 | 12 | 0.069592 | 0.069592 |
| A | 0.0625 | 12 | 12 | 0.069602 | 0.069602 |
| A | 0.03125 | 12 | 12 | 0.069604 | 0.069604 |
| B | 0.5 | 39 | 39 | 0.077928 | 0.077928 |
| B | 0.25 | 39 | 39 | 0.077933 | 0.077933 |
| B | 0.125 | 39 | 39 | 0.078025 | 0.078025 |
| B | 0.0625 | 39 | 39 | 0.078031 | 0.078031 |
| B | 0.03125 | 39 | 39 | 0.078032 | 0.078032 |

### Negative-control tolerance checks: strong sheet, v=0.02, T=4000

Counts at RK tolerance 1e-7 / 1e-9 / 1e-11. A solver-converged count cannot repair missing sheet or phase physics.

| Test/method | Lost counts | Failed counts | Velocity return RMS at 1e-11 |
|---|---|---|---:|
| A: C1L0 | 10/10/10 | 22/22/22 | 3.46e-11 |
| A: C1L1 | 10/10/10 | 22/22/22 | 2.60e-01 |
| A: C2fixed | 11/11/11 | 0/0/0 | 3.14e-01 |
| A: C2pphi | 11/11/11 | 0/0/0 | 2.98e-01 |
| A: C3/T0 | 11/11/11 | 0/0/0 | 6.30e-02 |
| B: C1L0 | 36/36/36 | 21/21/21 | 3.80e-10 |
| B: C1L1 | 36/36/36 | 20/20/20 | 3.61e-01 |
| B: C2fixed | 39/39/39 | 0/0/0 | 5.97e-01 |
| B: C3/T0 | 39/39/39 | 0/0/0 | 2.16e-01 |

Rows with failures require the loss bounds in the JSON summary; their distributions are deliberately omitted. Statistical agreement at a wide band can simply mean almost all the trajectory used FO.

### Measured epsilon exponents

Fits use three radii at fixed gyrotime T=50. C1 slopes are conditional on its valid markers; C2 starts in the outer smooth volume to exercise discarding phase. C4/C5 start next to the strong sheet. Values near zero for C1/C2 describe finite normalized phase-space error, not a convergence order.

| Candidate | T-A p(position/rho), p(velocity/v) | T-B p(position/rho), p(velocity/v) |
|---|---:|---:|
| C3/T0 | 0.996, 0.999 | 0.998, 1.002 |
| C3/T1 | 2.002, 2.039 | 2.006, 1.987 |
| C1L0 | 0.021, 0.048 | 0.012, 0.056 |
| C1L1 | -0.060, 0.049 | -0.037, 0.058 |
| C2fixed | 0.020, 0.001 | 0.001, -0.000 |
| C2pphi | 0.020, 0.001 | — |
| C4 | 1.013, 1.061 | 1.014, 1.057 |
| C5 | -0.267, -0.179 | -0.403, -0.323 |

Negative C5 slopes at fixed w/rho and fixed gyro step do not establish a sharp-target order; its phase-space error is about 1% here and does not vanish with decreasing rho at these fixed layer parameters.

<!-- END_MEASURED_RESULTS -->

## Ranked verdict and limits

This ranking concerns the tested physical target and finite horizons, not a
production certification for SPECTRE. Loss equivalence is not inferred merely
from overlapping confidence intervals, and exact switch energy is not a
surrogate for correct sheet dynamics.

| Rank | Method | Verdict | Reason and next discretization step |
|---|---|---|---|
| 1 | Event-resolved FO/RK45 | Viable bounded-time baseline | Resolves the physical force jump; converged ensemble statistics and invariant/tolerance checks. Ordinary adaptive RK is not symplectic; Gauss2 is the independent route below. |
| 2 | C3/T1, carried phase and one inverse | Viable controlled hybrid on these tests | First-order lift has measured quadratic local errors and much smaller momentum/return errors than discarded phase. Scan band and tolerance; no uniform long-time phase theorem or demonstrated production speed gain. A consistent noncanonical variational discretization needs equally consistent transition forms. |
| 3 | Resolved layer + Poincare Gauss2 | Viable independent T-A control | Four-dimensional canonical reduction, exact momentum by symmetry, converging width/step statistics. The finite layer and sampled endpoint are explicit limitations; no T-B implementation is supplied by the copied prototype. |
| 4 | C5, resolved thin-layer Boris | Viable secondary control under convergence checks | Good speed/volume properties and cheap field samples; both width and step must converge to the sharp reference. Neither exact energy nor reversibility removes layer bias. |
| — | C3/T0 | Not a phase-accurate default | Correct uniform lift, but leading smooth-field phase/transform errors remain; widening the band often gives agreement by retaining FO. A low-order approximation remains useful only with a stated error budget. |
| — | C2 fixed/pphi reconstruction | Not viable for phase-faithful orbit reconstruction | Phase is discarded irreversibly; matching one symmetry invariant does not fix the full orbit. Loss counts agree on this finite ensemble, which is a restricted empirical result, not a general rejection of every loss-only use. |
| — | C4 refreshed frozen helices | Not competitive at the tested curved-sheet accuracy | Exact planar subflow is viable and symplectic; the explicit state-dependent curved-field refresh produces secular momentum/return errors and parameter-sensitive losses. A canonical smooth-residual or sheet-cycle extension is still unverified. |
| — | C1 old GC resets | Not viable through a sharp sheet | Finite phase-dependent mu jump is missing. Forbidden mirror/outgoing-flow failures and momentum errors survive solver refinement. |
| pending | C3/T2, higher-order phase theory | Unavailable | Task D supplies no accepted T2; callbacks explicitly refuse it. Higher-order coefficients or consensus cannot substitute for a verified transform and RHS. |

The tested loops do not eliminate full-orbit cost for persistent sheet residents.
An accurate higher-order phase model and a canonical curved-sheet/cycle method
remain separate research requirements. No SIMPLE option, production default,
or deployment is changed by this experiment.
