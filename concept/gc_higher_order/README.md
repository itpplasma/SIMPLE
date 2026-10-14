# Standalone guiding-centre Lie-transform prototype

**Task E passes local canonical and symplectic checks** for midpoint and
two-stage Gauss–Legendre, retaining D2's orbit orders. See the
[canonical formulation and measured results](#task-e-canonical-symplectic-integration).

**Task D2 passes the stated orbit-order checks.** Selector 2 reconstructs
position with at least order three and velocity with at least order two on
the test fields over t in [0,1]; use `GC_FIRST_PLUS` (11, label 1+) for the
phase-complete first-order selector. Selector 1 remains the original baseline.
SIMPLE's `src/` and build are untouched.
Physics and its symbolic derivation live in [gen/](gen/); generic Lie algebra
and code generation use the local fortsym branch `feat/lie-pullback` at
`1372fbd`. No arbitrary-order guiding-centre recurrence is claimed.

## Task E: canonical symplectic integration

The standalone prototype now has canonical implicit midpoint (time order 2)
and two-stage Gauss–Legendre (time order 4), with analytic Hamiltonian
derivatives, radial/parallel Newton inversion and analytic stage Newton.
The physics generator is [derive_canonical.py](gen/derive_canonical.py), its
committed analytic field kernels are [canonical_kernels.py](generated/canonical_kernels.py),
and the integrator is [gc_canonical.py](gen/gc_canonical.py). The existing
fortsym-generated Fortran D2 transform supplies particle reconstruction.
Python is the standalone driver; no production SIMPLE integration is claimed.

### Decision and local one-form proof

Choose the additional canonicalization route, retaining D2's representation.
The literature's **Hamiltonian representation does not eliminate the
gyrogauge connection**: BT 2016 (3), (35) retains `J(dtheta-R.dX)` even when
the symplectic momentum has no higher-order corrections. In fixed spatial
coordinates a scalar gauge cannot remove a connection with nonzero curl.
For example, magnetic shear `B=(0,k x,1)` gives
`[curl(R*)]_z(x=0)=k^2/2` in D2's frame. Changing the Hamiltonian alone cannot
change this two-form. A further local phase-space transform is permitted;
this is distinct from merely selecting a Lie-transform scalar gauge.

There is also a useful correction to the task's radial obstruction: a
**single radial component can be removed locally** while carrying phase.
Take any SIMPLE-style spatial chart `(r,q1,q2)` with `A_r=b_r=0`, on a
contractible radial segment. In the README units the physical gyroaction is
`J=epsilon*mu` because `Omega=B/epsilon`. Define

```
K(r,q) = integral_[r_ref,r] R*_r(s,q) ds
zeta = -phi - K;  C_a = R*_a - partial_a K
P_a = A_a/epsilon + u*b_a - J*C_a                 (a=1,2)
H_can(q,P,J) = [u^2/2 + (J/epsilon)*B + epsilon^2*H2]_(r,u=inverse(q,P,J))
```

Substitution, including variable J, gives the exact local identity for the
truncated D2 theory
`(A/epsilon+u b).dX-J(R*.dX+dphi)=P1 dq1+P2 dq2+J dzeta`.
The radial term cancels because `K_r=R*_r`; its angular curvature is retained
in `C_a`, and its effect on the Hamiltonian enters through the inverse map.
Gyrophase is integrated as `zetadot=partial_J H_can`, rather than discarded.
The spatial `A_r` gauge is chosen first; adding `d(chi/epsilon)` to the
one-form has no effect on the dynamics.

This is an exact extended Darboux completion, rather than a claim that every
new variable changes at order epsilon squared. At fixed `(q,P,J)`, let
`a_a=A_(a,r)+epsilon*u*b_(a,r)` and `D=a1*b2-a2*b1 != 0`. Relative to the
map without `-J C_a`, its leading inverse displacement is

```
delta_r = epsilon^2*mu*(C1*b2-C2*b1)/D
delta_u = epsilon*mu*(a1*C2-a2*C1)/D.
```

Thus the slow spatial correction is O(epsilon^2), while the parallel
momentum generally changes at O(epsilon); the carried phase offset K is
O(1). A requirement that all six variables change only at O(epsilon^2)
would be incompatible with a generic parallel connection. Keeping the
implicit map exact avoids introducing additional truncation error into H.
The transform requested here is
`T_2^can = C o T_2^(-1)` from particle state to canonical GC, and its
reconstruction is `T_2 o C^(-1)`. Here T2 is D2's GC-to-particle map, and C
is the map above. For a natural frame rotated by psi from D2's Cartesian
frame, use `phi_natural=phi_D2-psi` in C; reconstruction uses
`phi_D2=-zeta-K+psi`. Canonical angles, including zeta, are carried without
wrapping. The frame angle psi is recovered modulo 2pi for the periodic
particle reconstruction; a caller needing a continuous D2 phase must unwrap
psi across frame branches.

D2 uses `R*=R+curl(b)/2`, including the perpendicular polarization momentum.
The task's `R+(b.curl b)b/2` corresponds to a different perpendicular
symplectic representation. It must come with the corresponding change in
H2 and particle pullback. Task E preserves D2's representation throughout.

[prove_canonical.py](gen/prove_canonical.py) checks 29 exact local one-/two-form,
inverse-displacement, radial-primitive, obstruction and slab identities.
The sheet chart uses `(r,theta,varphi)` with `b_r=0`, radial natural frame,
`R=10+r cos(theta)`, `c=5/a`, `d=sqrt(100+c^2 cos(theta)^2)` and
`v=r/(sqrt(r^2+c^2)+c)`. Its radial primitive is
`K=c sin(theta)/d [atanh((10v-c cos(theta))/d)-atanh(-c cos(theta)/d)]`.
The radial gauge of A_theta uses 20 terms of its convergent geometric
series; omitted derivatives through order three are below 1e-15 on
`0<r<=1`. Each sheet branch is analytic and used only in its own smooth
volume. No crossing rule or sheet-interface approximation is introduced.

### Analytic derivatives and verification protocol

Let `w=(r,u,mu,q1,q2)`, `v=(q1,q2,P1,P2,J)`, `M=partial_w v`,
`W=M^(-1)` and `g_v=W^T grad_w H`. Then
`H_vv=W^T [H_ww-sum_i (g_v)_i (v_i)_ww] W`.
Symbolically differentiated P and H kernels supply both Newton systems;
there are no finite-difference derivatives in integration. Midpoint uses
`A=(1/2), b=(1)`. Gauss–Legendre uses diagonal `1/4`, off-diagonal
`1/4 +/- sqrt(3)/6`, and weights `(1/2,1/2)`. Both satisfy
`b_i A_ij+b_j A_ji=b_i b_j`. The six-dimensional discrete tangent is
obtained by differentiating the converged stage equations, and is also
checked independently by finite differences of the actual solved map.

Measured Task E results are retained in `results/canonical_*.json`.
Acceptance covers the original six fields and selectors 1, 1+ and 2,
with the unchanged D2 full-orbit protocol: five epsilon scales, 401 raw
samples on `[0,1]`, last-four log fits, matched physical initial states,
and no phase alignment, time rescaling or energy normalization. Symplectic
steps are refined until RMS errors agree with D2 RK45 within 1%; the finest
sheet scale is also checked by step halving. This deliberately resolves
time error separately from the asymptotic GC error. Selector 1 preserves
its baseline phase limitation; phase-complete first order uses selector 1+.

| Field | 1+ p(x),p(v), both methods | T2 midpoint p(x),p(v) | T2 Gauss–Legendre p(x),p(v) |
|---|---:|---:|---:|
| gradient | 2.986, 2.092 | 3.985, 2.986 | 3.985, 2.986 |
| curvature | 3.016, 3.001 | 3.542, 3.003 | 3.542, 3.003 |
| shear | 3.004, 1.988 | 3.009, 1.960 | 3.009, 1.960 |
| sheet inner | 3.010, 2.126 | 2.970, 1.970 | 2.971, 1.969 |
| sheet outer | 3.019–3.020, 2.218–2.219 | 3.036, 2.571 | 3.037, 2.571 |

Every nonuniform fit differs from its D2 RK45 exponent by less than 0.001.
The largest relative difference in RMS errors is 0.5831%; the largest finest
step-halving change is 0.5719%. Uniform-field motion is exact at solver
roundoff. Resolving the small outer-sheet T2 error required 45,600 midpoint
steps on `[0,1]` at epsilon=.005, with 91,200 for the precision check;
Gauss–Legendre used 400. These are accuracy settings for reconstruction,
not recommended production time steps.

Long-time checks use `epsilon=.04`, `t<=4000`, fixed steps 2 and 1, and
eight equal time windows. They check bounded energy envelopes, the expected
step powers on the sheets, exact cyclic momenta, and remaining in the smooth
volume. Uniform/slab/shear/circular-curvature slow motions in these charts
have constant momenta and energy at roundoff, so their energy test alone
does not demonstrate a nontrivial modified Hamiltonian. The sheet tests do.

For both GC orders the measured sheet energy envelopes at `t<=4000` are:

| Field | Midpoint h=2 / h=1 | Gauss–Legendre h=2 / h=1 |
|---|---:|---:|
| inner | 1.42e-5 / 3.55e-6 | 1.11e-9 / 6.93e-11 |
| outer | 1.01e-4 / 2.52e-5 | 1.80e-8 / 1.13e-9 |

The largest final/first energy-window ratio is 1.000032; halving h gives
the expected factors four and sixteen. J and the
axisymmetric canonical p_phi have **zero measured change**. The analytic
six-dimensional map defect `maxabs(T^T Omega T-Omega)` is at most 2.22e-16;
independent finite differences of the solved map give at most 1.03e-9.
These are finite numerical checks, not all-time error bounds.

Cost was measured on `[0,20]`, epsilon=.04, one thread, excluding initial
inversion and particle reconstruction. Representative selector-2 measurements
are seconds / RMS slow-position error against a tighter canonical RK45 run:

| Method | inner | outer |
|---|---:|---:|
| midpoint, h=.1 | 0.737 / 1.54e-6 | 0.960 / 1.44e-5 |
| Gauss–Legendre, h=.1 | 1.44 / 4.46e-13 | 1.46 / 4.38e-12 |
| RK45, same analytic Python kernels | 0.118 / 6.93e-10 | 0.175 / 6.54e-10 |
| RK45, existing Fortran D2 RHS | 0.0814 / 1.31e-10 | 0.0902 / 4.85e-10 |

The RK45 comparison uses rtol=1e-9, atol=1e-11. Carried-phase RMS errors are
1.52e-5 / 1.02e-4 for midpoint, 3.36e-12 / 5.54e-11 for Gauss–Legendre,
and 2.17e-9 / 4.43e-9 for the Fortran RK45 on inner / outer volumes.
This prototype demonstrates structure and accuracy, **not a speedup**:
the settings have different errors, and Python plus repeated analytic
Hessian/inverse evaluation adds substantial overhead. The complete selector-1/2
cost, call counts and errors are in [canonical_cost.json](results/canonical_cost.json).

### What production SIMPLE would need

SIMPLE's canonical construction can retain `A_r=b_r=0` and its linear
parallel dependence. It must add the radial phase primitive K and angular
connection tables C_a, use `A_eff,a=A_a-epsilon^2*mu*C_a`, modify the
canonical momentum derivatives and inverse, add H2 consistently, and carry
the `(zeta,J)` pair for full-orbit reconstruction. The geometry tables K
and C_a can remain independent of mu; effective potentials and the inverse
are mu dependent and could instead be tabulated per mu. Keeping the phase
pair avoids needing a separate mu-dependent spatial canonicalization.
Axisymmetric tables must preserve exact toroidal independence.
Where `b_phi != 0`, the existing scalar radial Newton construction remains
available: `u=(P_phi-A_phi/epsilon+J*C_phi)/b_phi` is still linear in P_phi,
and substituting it into P_theta leaves one radial root. The prototype solves
the equivalent two-variable `(r,u)` system, which also supports the circular
curvature chart where the other angular component must eliminate u.

Third magnetic derivatives are needed for the H2 RHS and D2 transform;
**fourth derivatives generally enter an analytic Hessian of H2** for stage
Newton. Production would need sufficiently smooth higher-derivative splines
or another consistent derivative representation, not just existing first-
and second-derivative field tables. Phase patches must be joined with the
matching gyrogauge change; neither a single global phase chart nor a global
nondegenerate momentum inverse is guaranteed. The prototype certifies local
smooth charts, not the discontinuous sheet, magnetic-axis chart singularities,
or a global toroidal Darboux theorem. SIMPLE `src/` and its build are untouched.

Sources checked for this decision: [Littlejohn 1983](https://doi.org/10.1017/S002237780000060X)
(variational conservation laws); [Cary & Brizard 2009](https://doi.org/10.1103/RevModPhys.81.693)
(Hamiltonian and gyrogauge framework); [BT 2016, (3), (35), Appendix A](https://arxiv.org/html/1606.06534v1)
(the connection remains in the Hamiltonian representation);
[Brizard & Tronko 2011](https://arxiv.org/abs/1105.1145)
(gyrokinetic momentum conservation, not the second-order GC transform);
[Burby–Squire–Qin 2013](https://arxiv.org/abs/1312.3973)
(representation freedom); [Burby–Qin](https://arxiv.org/abs/1312.3975)
(global phase-chart obstruction); [Albert–Kasilov–Kernbichler, JCP 403,
109065 (2020), section 3](https://arxiv.org/html/1903.06885)
(implicit canonical momenta and radial solve); and the later
[Albert et al., general toroidal canonicalization](https://arxiv.org/html/2503.15278v2)
(spatial canonicalization beyond nested flux surfaces). These sources do
not establish the new numerical orbit bounds; those are measured here.

## Units, phase and field contract

The particle equations are `xdot=v`, `vdot=v cross B(x)/epsilon`, with positive
charge, unit mass, no electric field and no time dependence. Supplied B,
length and speed scales remain order one as epsilon decreases. The physical
cyclotron frequency is `|B|/epsilon`; epsilon is proportional to rho/L.

```
Bmag=|B|; b=B/Bmag; mu=vperp^2/(2 Bmag)
e1=normalized projection of Cartesian gauge_axis perpendicular to b
e2=b cross e1; a=cos(phi)*e1+sin(phi)*e2
v0=u*b+sqrt(2*mu*Bmag)*a; rho0=b cross (v0-u*b)/Bmag
R_j=e2 dot partial_j e1
```

Angles are radians, carried without wrapping; `phi_dot=-Bmag/epsilon` in a
uniform field. Paper theta is `-phi`. A smooth frame rotation psi gives
`phi'=phi-psi`, `R'=R+grad(psi)`, preserving `dphi+R.dX`. Use one nonsingular
chart per segment and change phase when changing charts. Nonpositive mu is
excluded by this polar chart.

[gc_api.f90](gc_api.f90) exports `dp=real64`, `field_derivs_t`, and

```fortran
type :: field_derivs_t
    real(dp) :: origin(3), b(3), db(3,3), d2b(3,3,3), d3b(3,3,3,3)
    real(dp) :: epsilon = 1.0_dp
    integer :: gauge_axis = 1
end type
call gc_transform(n, f, X,u,mu,phi, x_particle,v, ierr)
call gc_transform_jacobian(n, f, X,u,mu,phi, jac, ierr)
call gc_inverse(n, f, x_particle,v, X,u,mu,phi, ierr)
call gc_rhs(n, f, X,u,mu,phi, X_dot,u_dot,mu_dot,phi_dot, &
            bstar,bstar_parallel,drift, ierr)
```

Inputs have `intent(in)`, returned values `intent(out)`. Numerical arguments
are real(dp), selectors/status integers. Vectors have shape (3); `jac(6,6)`
has columns `(X1,X2,X3,u,mu,phi)` and rows `(x1,x2,x3,v1,v2,v3)`.
`db(i,j)=partial_j B_i`, with successive derivative indices in d2b/d3b.
The field must be solenoidal and mixed derivatives symmetric.

The type defines one **cubic Taylor field about origin**, reanchored internally
at X, including db/d2b. Third derivatives default to zero; supply actual third derivatives for
nonquadratic fields when using 1+ or 2. Newton inverts that
fixed field with the analytic Jacobian and a damped pivoted solve. Refresh
analytic jets at each external evaluation. For nonpolynomial inverse conversion,
refresh about the returned X and repeat; check against the actual field.
The inverse returns a principal phase branch, which a phase-carrying caller
may unwrap relative to its preceding value.

Status: `GC_OK=0`, `GC_ORDER=1` (unsupported selector), `GC_FIELD=2` (invalid
field, frame, Bstar denominator or phase solvability), `GC_NEWTON=3`,
`GC_POLAR=4`. Consumers must check ierr; availability is not an order proof.
The integer constant `GC_FIRST_PLUS=11` selects the label **1+**.

## Derivation and selectors

Start with `Gamma_p=(A/epsilon+v).dx-(|v|^2/2)dt`, `curl A=B`.
The generic finite scalar Lie exponential lives in fortsym. Consumer generators
use Tronko–Brizard (TB) 2015 (41)–(49), (55)–(56), and Brizard–Tronko (BT)
2016 (97), (100), (104), (112)–(114). The representation is
`Pi1=-mu*curl(b)/2`, `Pi2=0`, `Rstar=R+curl(b)/2`.

BT (114) determines `partial_theta sigma4=Q`; choose `<sigma4>=0`. This fixes
`G2_theta` and `G3_parallel_x` through BT (112)–(113). Q is a trigonometric
polynomial of degree at most six: G1 has degree <=2, D1(Pn) <=3,
D1^2(Pn) <=5, and F1a <=4. Sixteen equispaced phase samples recover the
finite Fourier polynomial exactly in exact arithmetic; its periodic primitive
and momentum derivatives supply the phase closure. A nonzero mean is rejected,
not removed. Generated derivatives also enter the analytic map Jacobian.

| Selector | Reconstruction | Reduced Hamiltonian |
|---|---|---|
| 0 | x=X+eps*rho0; v=v0 | H0; parallel/mirror motion, connection phase |
| 1 | x adds eps^2*rho1; v=v0-eps*G1.dv0 | H0; first-order symplectic form |
| 1+ (11) | same position as 1; second-order velocity pullback | H0+eps^2*H2 |
| 2 | x also adds eps^3*rho2; second-order velocity pullback | H0+eps^2*H2 |

Here `rho1=-G2_x-G1.d(rho0)/2` and
`rho2=-G3_x-G2.d(rho0)+G1.d(G1.d(rho0))/6`.
The second-order velocity is `v0-eps*G1.dv0+eps^2*(-G2.dv0+G1.d(G1.dv0)/2)`.
The third spatial displacement contributes to velocity at second order because
`phi_dot=O(1/eps)`. Merely retaining G1 twice is insufficient.

For 1+ and 2, the complete Hamiltonian derivatives are used:

```
H=H0+eps^2*H2; H0=u^2/2+mu*Bmag
Bstar=B+eps*u*curl(b)-eps^2*mu*curl(Rstar); Bpar=b.Bstar
Xdot=(H_u*Bstar+eps*b cross grad(H))/Bpar
udot=-Bstar.grad(H)/Bpar; mudot=0
phidot=-H_mu/eps-Rstar.Xdot; drift=Xdot-u*b
```

The missing eps^2 H2 in selector 1 produces an O(eps) frequency error and
therefore O(eps) velocity error over t=O(1). Selector 1+ includes the
`-eps*partial_mu H2` phase correction **and second-order velocity/inverse
closure**: phase-only correction cannot generally prevent a second-order
initial slow-coordinate bias from being amplified by B(X)/eps.
It is consequently more than a frequency-only patch to selector 1.

With a consistent T2 and H through eps^2, omitted H3 gives an O(eps^2)
phase error over fixed slow time. Perpendicular velocity has order-one
amplitude, so expected velocity exponent is **2**; the gyroradius is O(eps),
so its phase contribution to position has exponent **3**. Remaining slow
transport errors must also be O(eps^3). A third-order velocity claim would
require further phase and slow-dynamics analysis. Selector 1+ targets full
phase-space exponent at least 2; selector 1 retains its measured baseline.

`gc_diagnostics(...,coefficients(14),ierr)` retains rho0 (1:3), G1_u (4),
G1_mu (5), G1_theta (6), rho1 (7:9), H2 (10), averaged G2_mu (11), tau (12),
alpha1 (13), alpha2 (14). Averaged G2_mu alone is not the full transformation.

## Verification and measured orbit orders

| Field | n=0 p(x),p(v) | n=1 p(x),p(v) | 1+ p(x),p(v) | n=2 p(x),p(v) |
|---|---:|---:|---:|---:|
| uniform | exact; solver floor | exact; solver floor | exact; solver floor | exact; solver floor |
| constant gradient | 0.992, 0.988 | 2.002, 1.001 | 2.986, 2.092 | 3.985, 2.986 |
| curvature on fieldline | 0.989, 1.015 | 2.020, 1.006 | 3.016, 3.001 | 3.542, 3.003 |
| magnetic shear | 1.011, -0.004 | 1.977, 0.935 | 3.004, 1.988 | 3.009, 1.960 |
| sheet inner | 1.044, 0.021 | 1.991, 1.015 | 3.010, 2.126 | 2.971, 1.969 |
| sheet outer | 1.026, 0.014 | 1.957, 1.007 | 3.019, 2.219 | 3.037, 2.571 |

The unchanged fitting protocol uses epsilon `{0.08,0.04,0.02,0.01,0.005}`,
401 samples on t in [0,1], matched physical initial conditions, refreshed-jet
Newton initialization, and reconstruction at every sample. Fits use the last
four scales of RMS Euclidean position/velocity errors. There is no averaging,
phase alignment, time rescaling or energy renormalization. RK45 uses
rtol=2e-12, atol=2e-14
(tighter than Task D); a ten-times tighter run checks the finest scale. The
largest nonuniform error change is 0.514%. `check-orders` accepts the requested
1+ bounds (2,2) and T2 bounds (3,2), with the original 0.3 slope allowance and
a 1% precision gate. It also requires the complete four-selector study.
Selectors 0/1 retain their measured phase limitations; their labels do not
certify full phase-space orders. Higher observed powers in special fields
do not strengthen the general claim. Uniform-field slopes measure the
solver floor; closed-form Fortran tests verify exactness independently.

The curvature field `B=e_theta` has constant magnitude and constant curvature
magnitude on each circular fieldline. Sheet inner/outer tests use smooth
branches of [sheet_field.f90](../symplectic_sheet/sheet_field.f90), R0=10,
F=10, a=1/1.3, sheet radius 0.5; full orbits must stay in their own volumes.
Raw measurements and dependency/check evidence are in [results/](results/).

Independent exact slab checks use `B=(0,0,B0+g*x)` and canonical momentum:
`H2=-3*g^2*mu^2/(4*B^2)` and
`J=Eperp/B-eps*g*Eperp*v_y/B^3+eps^2*g^2*(3*Eperp*v_y^2/2+3*Eperp^2/4)/B^5`.
They check reconstructed energy through eps^2, canonical momentum through
eps^3, magnetic action through eps^2, sigma4 solvability and its primitive.
The existing action study observed exponents 0.997, 1.982, 2.978.
General EL/energy identities and the BT/BSQ beta-coefficient equivalence are
checked separately. BSQ's Pi1_perp=0 differs from the present polarization
representation by the Pi1_perp.vgc term in H2; their gauges are not mixed.
These identities and finite-field studies do not prove uniform error bounds.
Fortran also checks cubic-jet Jacobians/inverses, chart covariance, and local
position/velocity and Lorentz reconstruction residuals of orders three/two.
The 16 phase samples run independently under OpenMP (at most 16 threads);
one- and sixteen-thread results agree bit-for-bit on 48 field/phase states.
The retained Lean sandwich theorem only proves an abstract pairing statement
under stated hypotheses; it proves no GC field identity or accuracy bound.

## Build and regenerate

Run from the SIMPLE repository root. Normal compilation uses only gfortran
and committed generated kernels; the standalone Makefile never includes
SIMPLE's build. Generated CSE leaves use packed inputs and array temporaries,
compiled with `GEN_FFLAGS=-O0` to bound GCC register-allocation time; runtime/API
use `-O2`, with `-fPIC` throughout. Regeneration uses fortsym
`build/include` and `build/lib`;
override `FORTSYM=/path/to/fortsym` when needed.

```sh
make -f concept/gc_higher_order/Makefile test
make -f concept/gc_higher_order/Makefile study
make -f concept/gc_higher_order/Makefile check-orders
make -f concept/gc_higher_order/Makefile regen
make -f concept/gc_higher_order/Makefile check-generated
make -f concept/gc_higher_order/Makefile lean
make -f concept/gc_higher_order/Makefile canonical-test
make -f concept/gc_higher_order/Makefile canonical-check-generated
make -f concept/gc_higher_order/Makefile canonical-study
make -f concept/gc_higher_order/Makefile canonical-check-results
```

## Sources checked on 2026-10-14

* Tronko & Brizard (2015), PoP 22, 112507,
  [full preprint](https://arxiv.org/html/1511.01005v1).
* Brizard & Tronko (2016), *Equivalent Higher-order Guiding-center Hamiltonian
  Theories*, [full fourth-order derivation](https://arxiv.org/html/1606.06534v1).
* Burby, Squire & Qin (2013), *Automation of the guiding center expansion*,
  PoP 20, 072105, [author preprint](https://arxiv.org/html/1312.3973).
* Brizard & Tronko (2011), PoP 18, 082307,
  [author preprint](https://arxiv.org/abs/1105.1145): gyrokinetic momentum
  conservation; this is not the general second-order GC reconstruction paper.
