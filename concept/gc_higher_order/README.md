# Standalone guiding-centre Lie-transform prototype

**Task D2 passes the stated orbit-order checks.** Selector 2 reconstructs
position with at least order three and velocity with at least order two on
the test fields over t in [0,1]; use `GC_FIRST_PLUS` (11, label 1+) for the
phase-complete first-order selector. Selector 1 remains the original baseline.
SIMPLE's `src/` and build are untouched.
Physics and its symbolic derivation live in [gen/](gen/); generic Lie algebra
and code generation use the local fortsym branch `feat/lie-pullback` at
`1372fbd`. No arbitrary-order guiding-centre recurrence is claimed.

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
