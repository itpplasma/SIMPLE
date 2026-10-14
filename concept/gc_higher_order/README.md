# Standalone guiding-centre Lie-transform research prototype

**Task D is incomplete.** This directory supplies a reproducible phase-carrying
baseline, generated analytic Jacobians, exact algebra checks, a second-order
slab magnetic-action derivation, and measured failures of the requested orbit
order. It does **not** supply a general second-order GC transformation or an
arbitrary-order GC solver. Every `n=2` API call returns `GC_ORDER=1`.
Task C must check `ierr` and must not interpret these selectors as certified
accuracy orders. SIMPLE's source and build are untouched.

## Units, phase and local field contract

The full-orbit equations are `xdot=v`, `vdot=v cross B(x)/epsilon`, with
positive charge and mass normalized to one. The supplied `B` has order-one
magnitude; the physical cyclotron frequency is `|B|/epsilon`. Length and speed
scales remain fixed while epsilon decreases. Thus epsilon is proportional to
rho/L; the orbit-dependent proportionality factor is not fitted away.
There is no electric field or time dependence.

Write `Bmag=|B|`, `b=B/Bmag`, `mu=vperp²/(2 Bmag)` and

```
e1 = normalized projection of Cartesian gauge_axis perpendicular to b
e2 = b cross e1
a  = cos(phi)*e1 + sin(phi)*e2
v0 = u*b + sqrt(2*mu*Bmag)*a
rho0 = b cross (v0-u*b)/Bmag
R_j = e2 dot partial_j e1
```

Angles are radians and never wrapped during evolution. `phi_dot=-Bmag/epsilon`
in a uniform field; the angle theta in Tronko–Brizard is `theta=-phi`.
A smooth frame rotation by psi gives `phi'=phi-psi`, `R'=R+grad(psi)`;
`dphi+R.dX` is invariant. Use a fixed chart over an integration segment.
Changing charts requires changing the phase; a chart singularity is refused.
The polar chart excludes `mu<=0`.

`gc_api` exports `dp=real64` and

```fortran
type :: field_derivs_t
    real(dp) :: origin(3), b(3), db(3,3), d2b(3,3,3)
    real(dp) :: epsilon = 1.0_dp
    integer :: gauge_axis = 1
end type
call gc_transform(n, field_derivs, X,u,mu,phi, x_particle,v, ierr)
call gc_transform_jacobian(n, field_derivs, X,u,mu,phi, jac, ierr)
call gc_inverse(n, field_derivs, x_particle,v, X,u,mu,phi, ierr)
call gc_rhs(n, field_derivs, X,u,mu,phi, X_dot,u_dot,mu_dot,phi_dot, &
            bstar,bstar_parallel,drift, ierr)
```

All inputs have `intent(in)` and all returned variables have `intent(out)`;
`n,ierr` are integers, all numerical arguments are `real(dp)`. Positions,
velocities, `bstar` and `drift` have shape `(3)`; `jac` has shape `(6,6)`,
with columns `(X1,X2,X3,u,mu,phi)` and rows `(x1,x2,x3,v1,v2,v3)`.
See [gc_api.f90](gc_api.f90) for the complete exact declarations.
`db(i,j)=partial_j B_i`, `d2b(i,j,k)=partial_j partial_k B_i`.
Mixed derivatives must be symmetric and the physical field solenoidal.

The type represents a **quadratic Taylor jet about origin**, not a callback.
The API evaluates/reanchors this jet at X. Refresh it for every external
RHS or transform evaluation. Newton inverts the fixed jet with the generated
Jacobian and a damped pivoted solve. For a nonquadratic field, refresh the jet
about the returned X and repeat; the orbit study verifies the residual against
the actual analytic field. The inverse returns a principal phase branch, which
a caller carrying phase can unwrap relative to its previous value.

| Status | Meaning |
|---|---|
| `GC_OK=0` | Supported calculation completed; no accuracy-order certification |
| `GC_ORDER=1` | Unsupported order, including 2 |
| `GC_FIELD=2` | Invalid field/epsilon, singular frame or Bstar denominator |
| `GC_NEWTON=3` | Inverse failed or its linear system was singular |
| `GC_POLAR=4` | Nonpositive mu in this polar chart |

## Method and the implemented truncation

The particle phase-space form is
`Gamma_p=(A/epsilon+v).dx-(|v|²/2)dt`, with `curl A=B`.
Varying v imposes `xdot=v`; varying x gives the Lorentz equation. For scalar
functions the Lie pullback is the finite sum
`exp(epsilon L_G)f = f+epsilon G.df+epsilon² G.d(G.df)/2+...`.
The generic finite exponential lives in fortsym. The physical generator,
jet dependencies, frame and representation choices live in [gen/](gen/).

The first generating field is TB (41)–(44), in coordinates `(x,u,mu,theta)`:

```
G1_x = -rho0
G1_u = -u*rho0.kappa + mu*(tau+alpha1)
G1_mu = rho0.(mu*grad(log Bmag)+u²*kappa/Bmag)
        -mu*u*(tau+alpha1)/Bmag
G1_theta = (vperp/Bmag).(grad(log Bmag)+u²*kappa/(2*mu*Bmag))
           -rho0.R + u*alpha2/Bmag
```

Here `kappa=(b.grad)b`, `tau=b.curl(b)`,
`alpha1=-(rhohat*a+a*rhohat):grad(b)/2`,
`alpha2=(a*a-rhohat*rhohat):grad(b)/4` and `rhohat=b cross a`.
The implementation uses `G1_phi=-G1_theta` and the inverse pullback of v0.
It includes the spatial correction `rho1=-G2_x-G1.d(rho0)/2` from TB (45),(55),
with `Pi1=-mu*curl(b)/2`.

The selectors are:

* `n=0`: `x=X+epsilon*rho0`, `v=v0`; parallel/mirror RHS, connection phase.
* `n=1`: `x=X+epsilon*rho0+epsilon²*rho1`,
  `v=v0-epsilon*G1.dv0`; the reduced Lagrangian below.
* `n=2`: refused in transform, Jacobian, inverse and RHS.

For selector 1 the implemented reduced form is

```
Gamma_gc = (A/epsilon+u*b-epsilon*mu*Rstar).dX
           -epsilon*mu*dphi - H0*dt
Rstar = R+curl(b)/2; H0=u²/2+mu*Bmag
Bstar = B+epsilon*u*curl(b)-epsilon²*mu*curl(Rstar)
Bpar  = b.Bstar
Xdot  = (u*Bstar+epsilon*mu*b cross grad(Bmag))/Bpar
udot  = -mu*Bstar.grad(Bmag)/Bpar
mudot = 0; phidot=-Bmag/epsilon-Rstar.Xdot
```

Exact inversion of this truncated form explains its exact reduced-energy
identity. Retaining its rational denominators does not make it a second-order
asymptotic theory. In particular, the missing `epsilon² H2` contributes an
`O(epsilon)` gyrofrequency correction and produces measured velocity errors.
`drift=Xdot-u*b`; its first-order coefficients give the curvature and grad-B
drifts. The shear-line phase test verifies the Baños term `-u*tau/2` in this
representation; no extra parallel drift is inserted by hand.

`gc_diagnostics(...,coefficients(14),ierr)` returns rho0 (1:3), G1_u (4),
G1_mu (5), G1_theta (6), rho1 (7:9), H2 (10), averaged G2_mu (11), tau (12),
alpha1 (13), alpha2 (14). H2 and averaged G2_mu use TB (33)–(38), with
Pi2_parallel=0 and Pi1=-mu*curl(b)/2. They are diagnostic transcriptions;
H2 is **not** applied in `gc_rhs`. Averaged G2_mu is not a full second-order
particle-space magnetic moment. Its scalar Hamiltonian constraint is checked
exactly; that check is not an independent general derivation.

## Independent second-order result on the linear slab

For `B=(0,0,B0+g*x)`, conserve canonical momentum to define X by
`A_y(X)=A_y(x)+epsilon*v_y`. Set `x=X+epsilon*q`. The perpendicular Hamiltonian
has harmonic part `(p_x²+B(X)²*q²)/2`, cubic perturbation
`epsilon*B(X)*g*q³/2`, and quartic perturbation `epsilon²*g²*q⁴/8`.
Solving the canonical angle homological equation and averaging gives

```
H(J)=B(X)*J - epsilon²*3*g²*J²/(4*B(X)²) + O(epsilon³)
J(x,v)=Eperp/B - epsilon*g*Eperp*v_y/B³
       +epsilon²*g²*(3*Eperp*v_y²/2+3*Eperp²/4)/B⁵ + O(epsilon³)
Eperp=(v_x²+v_y²)/2; B=B(x)
```

The derivation script proves the homological equation, the averaged Hamiltonian
and each cancellation through order one in the time derivative of J. The
Fortran diagnostic H2 agrees with this independent normal-form result. The
full-orbit magnetic-action study observes exponents **0.997, 1.982, 2.978**
for truncations 0, 1, 2. This is a special-field invariant result, not T2.

## Verification and measured orbit orders

| Field | n=0 p(x), p(v) | n=1 p(x), p(v) | n=2 |
|---|---:|---:|---|
| uniform | exact; solver floor | exact; solver floor | unsupported |
| constant gradient | 0.992, 0.988 | 2.002, 1.001 | unsupported |
| constant curvature on fieldline | 0.989, 1.015 | 2.020, 1.006 | unsupported |
| magnetic shear | 1.011, -0.004 | 1.977, 0.935 | unsupported |
| sheet inner | 1.044, 0.021 | 1.991, 1.015 | unsupported |
| sheet outer | 1.026, 0.014 | 1.957, 1.007 | unsupported |

The study uses analytic Cartesian derivatives, adaptive RK45 full orbits,
identical physical initial conditions for both selectors, refreshed-jet Newton
initial conversion, and T_n reconstruction at all 401 samples on `t in [0,1]`.
It fits the last four of epsilon `{0.08,0.04,0.02,0.01,0.005}` using RMS
Euclidean position and velocity errors. There is no gyroaveraging, phase
alignment, time rescaling, or energy renormalization of the comparison.
A ten-times tighter RK tolerance is checked on the finest scale. The largest
relative change in nonuniform errors is 6.16e-5. Selector 1 fails the requested
phase-space exponent two in every nonuniform test; selector 0 also fails
exponent one in shear and both torus branches. `make ... check-orders` reports
these failures with a nonzero exit status using the retained measurements.
Uniform-field fitted slopes measure the numerical solver floor; the independent
closed-form Fortran test verifies exactness to roundoff instead.

The curvature field `B=e_theta` has constant magnitude and curvature `1/R`
along each circular field line, not a globally constant curvature vector.
The sheet tests use each smooth branch of
[the existing field](../symplectic_sheet/sheet_field.f90), with R0=10, F=10,
a=1 or 1.3, and sheet radius 0.5. Full-orbit trajectories are checked to stay
inside their selected volumes. No sheet jump or nonsmooth Taylor jet is used.
Raw samples and precision checks are in [results/](results/).

| Evidence | Exact scope |
|---|---|
| Fortran `test_gc` | Uniform exactness; nonuniform inverse/Jacobian; chart covariance; slab H2; shear phase; explicit refusals |
| SymPy `prove_identities.py` | EL inversion and reduced energy; first drifts; alpha and Baños averages; gyrogauge; scalar H2 constraint; slab derivation; polynomial Lie composition |
| Lean `sandwich_iff` | Abstract invariant-pairing sandwich iff reset preserves the pairing, assuming both surrounding maps preserve it and the before map is surjective |
| RK45 studies | Numerical asymptotics on the stated samples, not proofs of uniform error bounds |

Lean uses no added axioms or `sorry`. It does not prove differentiability,
nondegeneracy, a GC error bound, or any physical field identity. To interpret
its maps geometrically one must supply appropriate tangent lifts and pairings.
The unrestricted proposed necessity of Hamiltonian preservation is false:
canonical translation `(q,p)->(q+1,p)` is symplectic but changes `H=q`.
Event-dependent flight times and energy-surface resets require extra hypotheses.

## Obstruction and remaining work

The finite Lie exponential is implemented with order as a parameter (bounded
at 12 for resource control); the GC homological recurrence is not. TB's
G2_theta and G3_parallel_x are determined at fourth Lagrangian order. Computing
H2 or only averaged G2_mu does not supply them. A full second-order velocity
reconstruction also needs the corresponding third spatial displacement through
the Lagrangian constraint. Those generators, their gyrophase-dependent gauges,
and a verified representation-consistent inverse are missing here.

Adding general H2 forces to the RHS needs third Cartesian B derivatives through
`grad(H2)`; the present public field type deliberately supports only second
jets. Extending it and the complete reconstruction is necessary before enabling
selector 2. This is unfinished physics/implementation work, not a theorem that
arbitrary order is impossible, and not a claim that existing literature lacks
an algorithm. The requested general n=1 phase-space and n=2 orbit-order goals
have not been established.

## Build and regenerate

Run from the SIMPLE repository root; this Makefile never includes SIMPLE's.
Normal compilation needs only gfortran and committed generated code.

```sh
make -f concept/gc_higher_order/Makefile test
make -f concept/gc_higher_order/Makefile study
make -f concept/gc_higher_order/Makefile check-orders  # expected to FAIL
make -f concept/gc_higher_order/Makefile lean
make -f concept/gc_higher_order/Makefile regen
make -f concept/gc_higher_order/Makefile check-generated
```

Python verification uses NumPy, SciPy and SymPy; Lean uses 4.14.0 without
Mathlib. Regeneration links the Fortran consumer generator against fortsym's
`build/include` and `build/lib`. Override `FORTSYM=/path/to/fortsym` as needed.
The local fortsym branch `feat/lie-pullback` contains the new generic finite
Lie scalar pullback and a blocking CSE guard repair. The latter avoids an
exponential tree-cost traversal in full CSE; its regression compiles and runs
a deep shared-DAG recurrence against an independent numerical loop.
No GC physics was added to fortsym. The affected upstream native CTest gates
passed. The resident Fo lane later returned a state-provider error for status
and stop; a complete Fo gate and verified lane cleanup are not claimed.
Both repositories are committed locally;
no push is part of this task. Exact dependency revision/check evidence is in
[results/provenance.json](results/provenance.json).

## References checked on 2026-10-14

* Burby, Squire & Qin, *Automation of the guiding center expansion*,
  PoP **20**, 072105 (2013), [DOI](https://doi.org/10.1063/1.4813247),
  [author preprint](https://arxiv.org/abs/1312.3973). Bibliography verified;
  its arbitrary-order algorithm has not been implemented here.
* Littlejohn, *Variational principles of guiding centre motion*, JPP **29**,
  111–125 (1983), [DOI](https://doi.org/10.1017/S002237780000060X).
  Bibliography verified; a line-by-line formula comparison was not performed.
* Brizard & Tronko, PoP **18**, 082307 (2011),
  [author preprint](https://arxiv.org/abs/1105.1145): exact gyrokinetic momentum
  conservation, not the general second-order GC reconstruction.
* Tronko & Brizard, *Lagrangian and Hamiltonian constraints for guiding-center
  Hamiltonian theories*, PoP **22**, 112507 (2015),
  [full preprint used here](https://arxiv.org/html/1511.01005v1),
  [DOI](https://doi.org/10.1063/1.4935925). Equations cited above were checked.
* Parra & Calvo, *Phase-space Lagrangian derivation of electrostatic
  gyrokinetics in general geometry*, PPCF **53**, 045001 (2011),
  [institutional record](https://sites.fusion.ciemat.es/wp-publications/nr184/),
  [DOI](https://doi.org/10.1088/0741-3335/53/4/045001). Bibliography verified;
  equivalence of its full transformation is unverified in this implementation.
