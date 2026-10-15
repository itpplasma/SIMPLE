# Rejected SPECTRE interface models

This is the mathematical and numerical record for removing the options below.
It supersedes the crossing draft's claims of a regular GC thin-layer limit and
automatic symplecticity from midpoint evaluation. Removal is a SIMPLE policy
for sharp SPECTRE interfaces, not an impossibility theorem for all integrators.
Several requested stronger claims are false; their counterexamples are retained.

## Setting and verification scope

Use SI conventions: charge $q\ne0$, mass $m>0$, parallel velocity $u$, magnetic
moment $\mu=m v_\perp^2/(2B)$, and $B=|\mathbf B|>0$. In Gaussian units replace
$q$ in the GC one-form and gyroradii by $q/c$. The stationary interface
$\Sigma$ has $\mathbf B_\pm\cdot n=0$, with a finite tangential jump. Assume a
continuous, locally Lipschitz $A$, piecewise smooth up to $\Sigma$, with
$\nabla\times A=\mathbf B_\pm$. The supplied physical premise reports
gauge-matching residuals of $10^{-18}$ on tok2vol and NCSX; this task does not
rerun that calculation. The README's raw per-volume gauge jumps are compatible
with this premise. Derivatives of $A$, not its matched trace, jump.

All crossing theorems below concern isolated transversal events, finitely many
on a bounded time interval, with no grazing or simultaneous boundary events.
Algebra is checked exactly by
[`verify_rejected_models.py`](../concept/symplectic_sheet/proofs/verify_rejected_models.py).
Analytic estimates, physical ordering assumptions, recorded measurements, and
cited backward-error theorems are justified separately; SymPy does not prove them.

## T1. GC resets and the sharp-layer obstruction

**T1a (event composition criterion).** On the fixed-$\mu$ four-dimensional GC
phase spaces, let $\Gamma_\pm=(qA+m u b_\pm)\cdot dx$,
$\omega_\pm=-d\Gamma_\pm$, and $H_\pm=m u^2/2+\mu B_\pm$.
Assume $\omega_\pm$ nondegenerate, exact per-volume flows, transversal incoming
and outgoing flow, and a smooth invertible reset $R:S_-\to S_+$ between event
hypersurfaces. The fixed-physical-time flow-reset-flow map is symplectic iff

$R^*(\omega_+|_{S_+})=\omega_-|_{S_-}$ and
$d(H_+\circ R-H_-)|_{TS_-}=0$.

Thus, **within the class of energy-preserving resets**, it is symplectic iff
the reset is presymplectic on the event hypersurface. This preserves the
two-form derived from the one-form, not necessarily the one-form itself;
adding an exact gauge differential changes neither condition.

*Proof.* Write an incoming variation uniquely as $v=t+aX_-$, $t\in TS_-$,
using transversality; the event-time derivative (saltation) sends it to
$Sv=DRt+aX_+$. For $w=s+bX_-$ and $\Delta H=H_+\circ R-H_-$,

$\omega_+(Sv,Sw)-\omega_-(v,w)=(R^*\omega_+-\omega_-)(t,s)+a\,d\Delta H(s)-b\,d\Delta H(t)$.

Here $\iota_{X_\pm}\omega_\pm=dH_\pm$. Varying $t,s,a,b$ independently proves
necessity and sufficiency. The surrounding exact fixed-time flows preserve
their forms. Energy *equality* is not necessary for symplecticity: replacing
$H_+$ by $H_++C$ leaves its vector field unchanged. This disproves the
unqualified requested “iff presymplectic and energy-preserving” statement.

**T1b (point-fixed reset obstruction).** Suppose $R$ keeps the same physical
surface point and $\mu$, changing only $u$ to $f(x,u)$ on an open set. If it
is presymplectic, then $f_u b_+=b_-$, so $b_+=\pm b_-$ and $f_u=\pm1$.
For $\mu>0$, energy preservation also forces $B_+=B_-$ there.

*Proof.* Evaluate the pulled-back two-form on $(\partial_u,t)$ for arbitrary
surface tangent $t$. The matched $A$ has no $u$ dependence; the coefficient
is $-m f_u b_+\cdot t$ on the target and $-m b_-\cdot t$ on the source.
Both unit vectors are tangent, so equality for every $t$ gives the vector
identity, then the unit norms give $f_u^2=1$. Differentiating energy equality
gives $f f_u=u$; hence $f^2=u^2$ and $\mu(B_+-B_-)=0$.
Level 0's root $f^2=u^2-2\mu(B_+-B_-)/m$ conserves energy but fails this
necessary condition at a genuine direction jump (and generally a magnitude jump).

This obstruction cannot be extended to **arbitrary tangential resets**.
For a planar sheet $z=0$, take equal constant magnitudes,
$b_-=(1,0,0)$, $b_+=(\cos\phi,\sin\phi,0)$, and
$A_-=B(0,-z,0)$, $A_+=B(z\sin\phi,-z\cos\phi,0)$.
These potentials match at the sheet. With $Q_\phi$ the tangential rotation,
$R(x,y,u)=(Q_\phi(x,y),u)$ satisfies
$b_+\cdot d(Q_\phi x)=b_-\cdot dx$ and preserves $H$.
It is presymplectic for every $\phi$, including $b_+\ne\pm b_-$.
This is a counterexample to the requested universal reset prohibition,
not a physically justified instantaneous displacement of a particle.

Level 1 permits tangential displacement, so T1b alone does not rule out every
one of its crossings. Its actual scalar energy solve supplies no general
two-form identity: in the constant-field slab above the side curls vanish,
$X_u=0$, and its energy residual is zero at $\lambda=0$ for every surface
state. `solve_crossing` therefore returns the identity reset, whose two-form
defect is $m(b_+-b_-)\cdot dx\wedge du\ne0$. This is an algebraic witness
against general presymplecticity of that reset rule; this uniform slab itself
has no GC normal drift. Its Level-0 fallback also inherits T1b. Midpoint
evaluation of a varying GC bracket does not establish the pullback identity.

**T1c (loss of uniform GC ordering).** In a smooth planar layer take
$b(z)=(\cos\alpha(z),\sin\alpha(z),0)$ and constant $B$. Then

$b\cdot\nabla\times b=-\alpha'(z)$,
$B^*_\parallel=B+m u\,b\cdot\nabla\times b/q=B(1-\rho_\parallel\alpha')$, $\rho_\parallel=m u/(qB)$,
and $\det\omega=(m q B^*_\parallel)^2$ in $(x,y,z,u)$.

For $\alpha=\alpha_-+\Delta\alpha\,s(z/L)$ on $0\le z\le L$, with a monotone
$C^1$ profile, $s(0)=0$, $s(1)=1$, and $s'(0)=s'(1)=0$, a zero of
$B^*_\parallel$ exists when $\rho_\parallel\Delta\alpha>0$ and
$L<\rho_\parallel\Delta\alpha$. At any fixed scaled interior position with
$s'>0$, $B^*_\parallel/B$ has no finite limit as $L\to0$ at fixed nonzero
$\rho_\parallel\Delta\alpha$.

*Proof.* $\int_0^1s'\,d\xi=1$ implies $\max s'\ge1$. Under the stated signed
condition, $1-\rho_\parallel\Delta\alpha\,s'/L$ is negative somewhere and
equals one at both endpoints; the intermediate value theorem gives zeros.
The determinant formula shows degeneracy, so a Hamiltonian GC vector field
is not defined for general gradients there. Special removable cancellations
do not make the two-form invertible. Divergence at fixed $\xi$ follows directly
from the $1/L$ factor. For opposite sign, the ratio is $1+|\rho_\parallel\Delta\alpha|s'/L\ge1$:
there is **no sign change**. Thus the requested unsigned
criterion $L<\rho_\parallel|\Delta\alpha|$ and a universal “no GC limit”
theorem are false. What fails generally is a *regular, uniformly ordered*
GC approximation, not the existence of every specially regularized limit.

Away from interfaces, the necessary slow-field conditions include
$\rho|\nabla b|\ll1$ and $\rho|\nabla\log B|\ll1$, where $\rho$ bounds the
perpendicular and absolute parallel gyroradii. For a profile with order-one
shape derivatives these require $\rho/L_{\rm jump}\ll1$, with
$L_{\rm jump}=L/\max(|\Delta\alpha|,|\log(B_+/B_-)|)$.
Other GC curvature/drift orderings must also hold. A sharp SPECTRE interface
has $L=0$, so every nonzero direction jump violates this condition for finite
nonzero $\rho$, however small the angle. Smoothing within GC does not restore
a controlled zero-width approximation at fixed gyroradius.

**T1d (physical witness against conserved-$\mu$ crossing).** For a bounded
magnetic step, full-orbit $x,v$ are continuous: the impulse over a vanishing
passage time vanishes. In the equal-$B$ slab of T1b, let $v=(v_t,0,v_n)$,
$v_n\ne0$. Direct projection gives

$\mu_+-\mu_-=m v_t^2\sin^2\phi/(2B)$.

This finite jump is independent of any integration step. A fixed-$\mu$ reset
therefore cannot reproduce general full-orbit crossings even if its two-form
condition were repaired. The held-back `feat/spectre-symplectic-crossing` GC
reset falls under T1c–d; formal GC symplecticity cannot repair the physical
ordering. No blanket impossibility theorem about its tangential reset is asserted.

## T2. A Boris smoothing layer changes the target problem

**Theorem (transversal thin-layer error and fixed-step cost).** Consider a
bounded slab field $\mathbf B_w(z)$ transitioning over width $w$, compared
with the sharp field $\mathbf B_0$. Assume fixed jump/profile shape, matched
entry $x,v$, $|v|=v_0$, and $|v_n|\ge\eta v_0>0$ throughout both passages.
Put $\Omega=|q|B_{\rm ref}/m$, $\rho=v_0/\Omega$, and $\epsilon=w/\rho$.
Constants below depend on $\eta$ and field ratios; they are not uniform at grazing.
The normalized exit-velocity/impulse discrepancy is $O(\epsilon)$ as
$\epsilon\to0$, and can be leading order $\epsilon$. At $w\sim\rho$ this
bound supplies no small error; an order-one discrepancy is possible, not mandatory.

*Proof.* Use $z$ as independent variable:
$dv/dz=(q/m)(v\times\mathbf B(z))/v_n$. Its state derivative is bounded
on the assumed transversal domain independently of $w$. An integral
inequality (or Gronwall on this interval of length $w$) bounds the difference
by a constant times $\Omega w/\eta$, divided by $v_0$ giving $O(\epsilon)$.
To first order, with entry velocity $v_0^{\rm vec}$ and entry normal speed
$v_{n0}$, the impulse difference is

$\Delta I/m=(q/(m v_{n0}))v_0^{\rm vec}\times\int(\mathbf B_w-\mathbf B_0)\,dz+O(v_0\epsilon^2)$.

For a layer on $[0,w]$ with
$\mathbf B_w=\mathbf B_-+\Delta\mathbf B\,s(z/w)$,
$s(\xi)=6\xi^5-15\xi^4+10\xi^3$, and the sharp jump at zero,
$\int(\mathbf B_w-\mathbf B_0)dz=-w\Delta\mathbf B/2$.
This verifies the nonzero leading coefficient for velocities not parallel to
the jump. A symmetric profile can cancel this coefficient: $O(\epsilon)$ is
an upper bound, not a universal lower bound. There is no delta-function
magnetic kick at the sharp sheet; both impulses are integrated between the
same entry and exit sections.

An exact finite-width witness uses $\mathbf B=B_y(z)e_y$: then
$dv_x/dz=-qB_y/m$. Take $q>0$, $B_-=B_{\rm ref}$, $B_+=2B_{\rm ref}$,
entry $v=(0,0,v_0)$, and the preceding profile. At $w=\rho/3$, the sharp
exit has $v_x=-2v_0/3$ and the layer exit $v_x=-v_0/2$, a difference $v_0/6$;
both retain $v_n\ge\sqrt5\,v_0/3$. Thus an order-one normalized difference
at fixed $w/\rho$ is possible while transversality remains uniform.

The requested **$O(1)$ Hamiltonian change inside the layer is false** in
matched gauges. Integrating the bounded field difference gives
$\delta A=O(B_{\rm ref}w)$, and with $\pi=p-qA_0$,

$H_w-H_0=-q\pi\cdot\delta A/m+q^2|\delta A|^2/(2m)$.

For $|\pi|=O(mv_0)$ this is $O(\epsilon)$ relative to kinetic energy, although
the *force* changes by order one relative to the sheet force jump. An arbitrary
gauge comparison at the same numerical $p$ would not measure a physical error.

Let $h$ denote time, $\delta=\Omega h$, and $\kappa=v_0h/w=\delta/\epsilon$.
Resolved second-order Boris stepping needs $\kappa\ll1$ (spatially
$v_0h\ll w$, not the dimensionally incorrect $h\ll w$) and $\delta\ll1$.
For a $C^2$ bounded-shape profile, force time derivatives in the layer bound
the accumulated crossing discretization error by $O(\epsilon\kappa^2)$:
sum $O(h^3\Omega v_0(v_0/w)^2)$ over $O(w/(v_0h))$ steps.
Thus the per-crossing velocity error bound relative to the sharp target is
$O(\epsilon+\delta^2/\epsilon)$, with smooth-region error $O(\delta^2)$ over
a fixed gyro-time. A global fixed step costs $O(1/\delta)$ steps per gyro-time;
at $\delta=\kappa\epsilon$ the bound is $O(\epsilon(1+\kappa^2))$ and cost
$O(1/(\kappa\epsilon))$. Local refinement changes this cost model.
For general crossings both width and step must converge; fixed nonzero width
converges to a different field. This excludes an uncontrolled smoothing knob,
not controlled thin-layer full-orbit approximations.

## T3. An unresolved single-branch step loses its smooth order

**Theorem.** Let a time step $h$ have one transversal crossing at fraction
$\theta\in(0,1)$. Suppose the method uses one side's field for the entire
step and is at least first-order consistent on that smooth side. With error
defined as numerical minus exact, and $[[f]]=f_+-f_-$, its leading defect is
$-(1-\theta)h[[f]]$ for the incoming branch, or $\theta h[[f]]$ for the
outgoing branch, plus $O(h^2)$. Smooth-side order does not remove this term.

*Proof for uniform magnetic fields.* Define $G_\pm v=(q/m)v\times\mathbf B_\pm$.
The exact velocity step is
$e^{(1-\theta)hG_+}e^{\theta hG_-}v_0$
$=[I+h(\theta G_-+(1-\theta)G_+)]v_0+O(h^2)$.
Boris uses the Cayley rotation
$C_h(G)=(I-hG/2)^{-1}(I+hG/2)=I+hG+h^2G^2/2+O(h^3)$.
Subtracting gives the stated defect with
$[[f_v]]=(q/m)v_0\times[[\mathbf B]]$ (or $[[F]]/m$).
Its position defect is $O(h^2)$; the velocity defect already fixes the leading
phase-space order. $C_h^TC_h=I$ proves that exact speed conservation does not
remove this error. The rejected Boris prototype samples the midpoint branch: the
coefficient is $\theta$ for $\theta<1/2$, $-(1-\theta)$ otherwise, provided
$\theta$ stays separated from $1/2$ as $h\to0$. In the shrinking interval
$|\theta-1/2|=O(h)$, the drift midpoint can choose either branch and the
remainder of this coefficient formula is not uniformly $O(h^2)$.
The leading uniform-phase moments below are unaffected by that interval.

For uniform random $\theta$, the incoming-branch coefficient has mean $-1/2$;
the midpoint coefficient has mean zero and variance $1/12$. **If** projected
kicks have equal magnitude $D=|[[f]]|$ and are independent with no transport
amplification, their accumulated rms is $\sqrt{N}hD/\sqrt{12}$.
Generally the variance includes covariance terms and tangent-flow weights;
equidistribution alone does not imply a random walk. Coherent errors can be
$O(NhD)$; cancellations are possible. These are conditional accumulation
models, not universal theorems for deterministic orbits.

The requested statement about **any fixed-step method** is false: exact
event splitting has no missed-branch term. Even unresolved single-branch
Boris can converge at fixed finite time as $h\to0$, with constants depending
on crossing count and stability. Its long-orbit requirement in the independent
kick model is $hD\sqrt N\lesssim$ the tolerated error. The rejected implementation
offers neither exact event splitting nor a width/step convergence contract;
the evidence below justifies its removal, not mathematical nonconvergence of Boris.

## T4. RK45 GC cannot repair the crossing model

**Theorem (conditional on the selected interface options).** For `integmode = 0`
using GC Level 0 or 1 resets, decreasing RK45 tolerance cannot remove their
model error. The limiting smooth-volume solver still composes the same reset
with the GC flows. T1c's ordering failure and T1d's fixed-$\mu$ counterexample
survive exact per-volume integration. Thus this path has no generally valid
sharp-interface crossing model. This does not exclude GC away from the sheet
joined to a separately validated full-orbit transition and higher-order lift.

## T5. Valid full-orbit alternatives and their limits

**T5a (exact full-orbit flow).** Under the stated regularity and transversal
event hypotheses, $H=|p-qA|^2/(2m)$ has a symplectic, energy-preserving
fixed-time flow across the sheet. If a coordinate is cyclic in a
symmetry-compatible matched gauge, its canonical momentum is conserved.

*Proof.* In each volume, $\dot Z=J\nabla H$ with
$J=\left(\begin{smallmatrix}0&I\\-I&0\end{smallmatrix}\right)$ gives
$\dot H=\nabla H^TJ\nabla H=0$ and
$(JH'')^TJ+J(JH'')=0$. Continuity of $A$ makes $x,p,H$ continuous at an event,
so there is no reset impulse. The common trace of $H$ has equal tangential
derivatives; its gradient jump is normal to the event surface. Consequently
the saltation matrix has the form
$S=\left(\begin{smallmatrix}I&0\\c\,nn^T&I\end{smallmatrix}\right)$,
which satisfies $S^TJS=J$. Composition proves the result; in a cyclic
coordinate $\dot p_j=-\partial_jH=0$. At grazing, a classical differentiable
event map need not exist; this theorem does not assert one.

**T5b (resolved smooth Poincare Gauss2).** For smooth $K=g(x)(H-E)$,
$0<g_{\min}\le g\le g_{\max}$, exact initial shell $H=E$, and a fixed fictitious
step $h_\tau$, two-stage Gauss-Legendre is symplectic and order four. It has
Butcher weights $b_1=b_2=1/2$ and coefficients
$a_{11}=a_{22}=1/4$, $a_{12}=1/4-\sqrt3/6$, $a_{21}=1/4+\sqrt3/6$;
$b_i a_{ij}+b_j a_{ji}=b_i b_j$ and all order-four conditions are checked.
Symplecticity here concerns the map at fixed fictitious time; sampling at a
state-dependent physical endpoint requires its own event/extended-phase construction.

The standard backward-error result gives $|K_n-K_0|\le C h_\tau^4$ over
exponentially long fictitious times for analytic $K$, a sufficiently small
step, a globally defined modified Hamiltonian, and numerical orbits confined
to a compact domain of uniform analytic bounds. Sufficient finite differentiability
gives only the corresponding finite-order, polynomial-time result. See
Hairer, Lubich and Wanner, *Geometric Numerical Integration*, 2nd ed. (2006),
chapters VI and IX ([publisher](https://doi.org/10.1007/3-540-30666-8));
the compactness and smoothness conditions are also explicit in
[Hairer's lectures, section 3](https://www.unige.ch/~hairer/poly_geoint/week3.pdf).

The product rule gives $X_K=gX_H+(H-E)X_g$, so on $K=0$, $dt/d\tau=g$
and the exact trajectory is physical. On a residual shell $K=k$,
$H-E=k/g$ and $X_K/g=X_{H-k/g}$: the effective extra potential has sign
**$-k/g$**. Therefore $|H_n-E|\le C h_\tau^4/g_{\min}$ under the preceding
bound, not exact energy conservation. Hairer's *Variable time step integration
with symplectic methods*, *Applied Numerical Mathematics* **25** (1997),
219–227, is verified against the
[author's abstract](https://www.unige.ch/~hairer/preprints/varsymp.html) and
[publisher](https://doi.org/10.1016/S0168-9274(97)00061-5).
No uncertain theorem number or unverified sharp-sheet extension is attributed to it.

The step must resolve **every** scale of $K$, including the transition of $g$:
locally $h_\tau\ll L/(g|v_n|)$. Constants need not be uniform as layer width
or $g_{\min}$ tends to zero. Positive $g$ alone does not smooth a discontinuous
gradient of $H$; standard smooth backward-error analysis does not cover a
sharp-sheet Gauss step. A finite smoothing layer still needs T2's target-limit
check. A smooth-layer Gauss2 experiment must validate its interface construction
and nonlinear solve accuracy before these bounds can certify it. No production
Gauss2 full-orbit option is retained; this document does not certify the experiment.

**T5c (event-resolved full-orbit RK45).** On a fixed finite time interval with
finitely many uniformly transversal crossings, piecewise smooth Lipschitz
fields and stable smooth-side solves, RK45 with exact event landing and
branch restart converges as the maximum step/tolerance tends to zero.
Smooth-side local errors tend to zero, event times/states depend continuously
on incoming states by transversality, and a finite induction over events
proves convergence. Fifth-order global accuracy requires sufficiently smooth
side fields and event errors commensurate with that order; exact landing alone
does not guarantee fifth order at grazing or with inaccurate dense output.

RK45 has no general symplectic guarantee. For explicit RK, the diagonal
symplectic coefficient condition is $2b_i a_{ii}=b_i^2$; $a_{ii}=0$ would
force every $b_i=0$, contradicting consistency. Ordinary adaptive stepping
does not restore a general guarantee. Energy error **need not** grow linearly:
the free-particle Hamiltonian is integrated exactly. As a witness that drift
is possible, Dormand–Prince's fifth-order stability polynomial is
$R(z)=\sum_{j=0}^5z^j/j!+z^6/600$. For a harmonic oscillator,

$|R(i\eta)|^2=1-\eta^6/1800+\eta^8/1600-\eta^{10}/14400+\eta^{12}/360000$.

Its energy multiplier after $N$ equal steps is this factor to power $N$;
for small $N\eta^6$ the leading drift is linear in $N$ (hence time).
For stable event-resolved runs, summing bounded per-step energy defects gives
an upper bound proportional to elapsed time and error budget, not a universal
nonzero linear law. The full-orbit implementation proposed in PR #530 uses this method for
bounded spans with tolerance and invariant checks; very long-time bounded energy needs a validated
symplectic interface method. No claim that RK45 is always unsuitable at long
times follows merely from its lack of symplecticity.

## Recorded evidence (not symbolic proof)

The [prototype README](https://github.com/itpplasma/SIMPLE/blob/4ce6d88/concept/symplectic_sheet/README.md) records 64 NCSX
3.5 keV alpha markers near each interface, traced for 2 ms. Counting uncertainty
is approximately 2/4 markers at interfaces 1/2, and individual chaotic fates
are not comparable. The outer jump reaches 0.32 rad and -25 to +38% in magnitude.

| Model | Lost, interface 1 | Lost, interface 2 | Interpretation |
| --- | ---: | ---: | --- |
| Full orbit everywhere, RK45 reference | 6 | 39 | Sharp-field baseline |
| Level 3, RK45 | 7 | 41 | Agrees within quoted counting uncertainty |
| Level 1, symplectic GC | 7 | 29 | Outer-interface under-loss, consistent with T1 |
| Level 1, RK45 GC | 1 | 7 | 32 reference losses missed; T1/T4 |
| Level 3, Boris, layer 0.2 rho, rotation 0.0125 | 5 | 48 | Over-loss; consistent with T2/T3, not a causal proof |

The [prototype crossing record](https://github.com/itpplasma/SIMPLE/blob/4ce6d88/DOC/spectre-interface-crossing.md) gives tok2vol canonical
$p_\zeta$ relative spreads: sharp Boris **14%** (rotation 0.05), **3.2%**
(0.0125); layer 0.2 rho **4.0%** (0.05), **0.93%** (0.0125);
event-resolved full-orbit RK45 **$6\times10^{-8}$** at tolerance $10^{-9}$
(a tangential touch reaches $1.2\times10^{-6}$). Speed conservation does not
guarantee canonical-momentum conservation, as T3 shows.

The prototype also records RK45 relative energy errors $2\times10^{-4}$ and
$2\times10^{-3}$ at times $4\times10^4$ and $4\times10^5$ (tolerance $10^{-6}$),
consistent with secular drift. Resolved analytic-layer Poincare Gauss2 records
$5.4\times10^{-7}$ and $5.1\times10^{-7}$, consistent with T5b. These figures
do not prove universal drift or a uniform sharp-width Gauss2 bound.

## Do not reintroduce

Multi-volume SPECTRE orbit tracing now uses only `crossing_level = 3` with
`integmode > 0`: sharp-field RK45 full orbits with exact interface crossings
and matched-gauge momentum bookkeeping. The rejected GC and smoothed-layer
Boris prototypes remain historical evidence, not available configurations.
T5b–c state the verification conditions; see
[the interface status](spectre-interface-crossing.md).

The removed GC code last existed on main at
[`b113786`](https://github.com/itpplasma/SIMPLE/commit/b113786), the branch-point
revision before this removal. The removal commit's parent `61316ee` adds only
this proof record and also retains the old code. Its introducing PRs were
#443 (Level 0), #440
(Level 1), and #441 (symplectic flow/map composition); #459 refined symplectic
event landing and #477 added tangential relocation. The code remains in Git
history for inspection, not reuse.

| Rejected option/model | Reason to reject as a sharp-interface prescription | Main provenance / future treatment |
| --- | --- | --- |
| `crossing_level = 0` | T1b point-fixed obstruction; T1c–d physical ordering and mu failure | `b113786`, #443; removed without aliases; use level-3 full orbit |
| `crossing_level = 1` | T1c–d; energy solve alone fails the general form guarantee in T1a/b | `b113786`, #440, #441, #459, #477; removed without aliases; use level-3 full orbit |
| RK45 GC crossing (`integmode = 0`, GC reset) | T4 inherits T1 irrespective of tolerance | `b113786`, #443/#440; removed without aliases; use level-3 full orbit |
| Prototype `spectre_fo_integrator = 'boris'` | T2 finite-width model error; T3 unresolved crossing error and measured invariant spread | Never a main interface option; event-resolved full orbit |
| Prototype `spectre_fo_layer_larmor` | T2 uncontrolled fixed-width target change | Never a main option; a new convergence contract is required |
| Prototype `spectre_fo_rotation` | T3 gyromotion bound alone does not locate or resolve a sheet | Never a main option; event landing and tolerance control |

The held-back `feat/spectre-symplectic-crossing` must not be promoted as a
physical sharp-interface cure: GC symplecticity does not bypass T1c–d.
Reopening requires a new physical/asymptotic contract and independent
full-orbit convergence/invariant evidence, not a renamed knob or energy check.

## Implementation removal ledger

Updated 2026-10-14. This ledger distinguishes a demonstrated failure of a
specific formulation from an unfinished proof, insufficient resolution, or a
method rejected for a different target. A small energy error or matching loss
count does not repair the failures below.

The newly collected comparison was added to main's cleanup branch in
[`4a1f842`](https://github.com/itpplasma/SIMPLE/commit/4a1f842), preserving the
original GPT worker result `345bfa8`. The removal commit
[`3a6fd37`](https://github.com/itpplasma/SIMPLE/commit/3a6fd37) deletes its retired
executable paths and settings, plus the earlier sharp-sheet branch/landing code
already on main. Both commits are retained by a merge commit,
rather than squash-merging away the add/remove history. Historical tables,
compressed marker records and mathematical witnesses remain evidence, not
selectable methods. F1–F5 are formulation/domain failures or measured failures
of the named implementation. F6 is retirement from the requested symplectic
target; it is not a claim that every Boris calculation is wrong. F7 records a
numerical chart failure in the later performance candidate, not a rejection of
full-orbit physics. F8 records an inconsistent absorbing boundary in the
multi-volume hybrid.

| ID | Failed contract and decisive evidence | Removed implementation/settings |
| --- | --- | --- |
| F1: fixed-mu GC sheet reset | **Physical counterexample**, T1d: continuous particle velocity at a direction/magnitude jump generally changes mu. T1b also obstructs point-fixed presymplectic energy resets. No smaller RK tolerance repairs the model. Historical moderate T-A C1L0/L1 runs each also reported 9 failed trajectories. | Comparison `reset_old`, C1L0/C1L1 dispatch and `old_level`; production levels 0/1 and GC-only multi-volume RK45 were already removed/rejected by #531/#530. |
| F2: phase-discarding hybrid as a reversible/symplectic formulation | **Information loss and measured reversal failure.** Multiple phases are mapped to one state; choosing phase zero or a momentum root cannot invert that projection. Moderate T-A reverse velocity RMS/v was about 0.30/0.31 for C2fixed/C2pphi, versus 2.8e-5 for FO at the same tolerance. This rejects that formulation as a reversible/symplectic candidate; it does not prove every finite-time hybrid loss statistic wrong. | C2fixed/C2pphi dispatch, `discarded_projection`, `match_phase`, `phase_rule`, and their matrix entries. |
| F3: refreshed frozen-plane patches as a curved-field orbit solver | **Measured failure of this implementation.** Moderate T-A C4/T1 at patch step .25 had normalized toroidal-momentum RMS error about .60 versus 1.3e-9 for FO, and reverse velocity error about 1.32. The long strong T-B run also had two reconstruction failures. Refreshing a state-selected plane has no inherited exact planar symplectic guarantee. | C4 dispatch, `local_patch`, `patch_step` and its scans. The exact constant-field `helix`/`plane_step`/cycle oracle stays: it passes inverse, repeated-event and compression checks. No impossibility of curvature corrections is claimed. |
| F4: unresolved fixed-branch sharp-sheet stepping | **Local error counterexample**, T3: crossing error is generically O(h times force jump), irrespective of the smooth method's nominal order. Historical naive sharp stepping gave rms relative mu change about .85 against a .017 reference. Midpoint branch choice reduced bias but retained excess scattering. | `SCHEME_NAIVE`, `SCHEME_MIDPOINT_REGION`, `midpoint_step` sharp-branch dispatch, their CLI selectors and dedicated stationary experiment. Resolved smooth midpoint is not rejected. |
| F5: energy-kick landing as a general sheet integrator | **Domain and accuracy failure.** The normal kick fails near grazing when the energy-matching kick exceeds available normal velocity. Recorded rms mu changes were 1.9 and .24 at h=.1/.05 against the .017 reference; a smaller step helped but did not supply grazing coverage. Successful tangent checks on accepted transversal events were not a general physical-accuracy proof. | `SCHEME_LANDING`, `sheet_step`, landing/root/kick helpers, tangent/record machinery and the dedicated landing test. This rejects this implemented kick rule, not variational event integrators generally. |
| F6: Boris advertised as a canonical symplectic sheet option | **Wrong guarantee**, T3 and T5: speed/volume preservation does not establish canonical symplecticity; finite smoothing also changes the target unless width and step converge. Resolved Boris can approximate a smooth field accurately, so it is not universally an inaccurate method. It is retired here as a candidate for the requested symplectic sharp-sheet solver. | Comparison C5 dispatch, `boris`, copied `boris_step`, `boris_step` setting and width/step scans. SPECTRE production Boris/layer controls were already removed by #530; generic full-orbit Boris elsewhere is unaffected. |

F7: `2e243d4`/`6bbdef0` applied the polar FO splines at the magnetic axis.
In the 32-marker, 3.5 keV NCSX volume pilot (200 microseconds), marker 31
reached rho=6.19e-13. FO entry returned `SPECTRE_FO_LOCATE_FAIL`; the default and
refined grids each recorded 4,666 warning-held intervals. Its apparent
confinement is invalid even though aggregate loss counts match main. The
independent opposite-point raw-field limit gives relative field error 49.0 at
rho=1e-12 for the default polar spline. The field derivative error is amplified
by the vanishing polar Jacobian; refining the integrator cannot fix it.
The replacement removes the polar-field/axis-inverse route near the axis and
uses regular xi/eta splines for geometry and all potential components, joined
to the bulk potential before taking curl. It adds no compatibility setting or
raw spectral fallback. The new axis-limit error is 4.21e-9 on NCSX and the
pole-seeded inverse passes. A replay of marker 31 has one FO entry/exit and
zero holds/failures. The retained initial condition and failure counters are
in `concept/spectre_loss_study/results/axis-failure.json` on the study branch.
Larger ensemble reruns and the existing invariant/event checks still apply.

F8: the multi-volume hybrid excluded `rho_g = Mvol` from FO proximity entry
and otherwise absorbed its guiding centre there. The FO reference absorbs the
particle, which can reach that same surface while its centre remains inside.
The matched 256-marker, 2 ms G3V3L3Fi replay lost 48 hybrid versus 53 FO markers.
Entering FO at the outer edge with the existing entry/exit distances recovers
markers 62, 124, 126 and 247, giving 52 losses on both default and refined
grids, with zero numerical stops, FO failures or warning holds. The old
multi-volume GC-only absorber is removed, without a compatibility setting.
The single-volume GC model retains its separate GC-edge contract.
An independent direct-FO wall witness has an interior guiding centre at the
particle's terminal event; raw boundary geometry also verifies immediate
loss for exterior launch particles. Both regressions fail with the old outer
exclusion and pass after repair. Source history retains the previous path.

The remaining marker 129 is not accepted as a physical hybrid/FO difference.
FO launch phases span loss times 1.910–2.055 ms across the 2 ms cutoff.
Field-grid and construction refinement do not remove its discrepancy;
GC-step controls do not establish trajectory convergence. Keeping FO active
throughout the hybrid pipeline recovers the reference loss time to about
1e-16 s. This supports the existing F2 switch/phase limitation and does not
justify tuning the entry distance until counts agree. Retained evidence is
`concept/spectre_loss_study/results/loss-diagnosis/` on the research branch.

Current comparison executables accept only FO, C3 (phase-carrying concept), and
GL2 (resolved smooth-layer control). The retired candidate dispatches and their settings are absent; no compatibility
handling or new rejection checks/tests are added.
The copied older hybrid driver and the accuracy driver dedicated solely to the
retired candidates were removed with their build targets. No aliases or shims
map old names onto a remaining method. Fresh study data use `build/study/` so
historical measurements are not overwritten.

### What this does not establish

- The current production event-resolved RK45 GC/FO bridge remains the finite-time
  baseline. It is not claimed reversible or globally symplectic. Its zeroth-order
  switch must be replaced coherently; removing it without a replacement would
  remove working multi-volume tracing, not validate another formulation.
- The slow NCSX smooth-layer Gauss2 implementation is **not a production candidate**
  in its present form. Its dispatch/namelist option remains absent. Its retained
  fixed-parameter concept has successful resolved tests; the timing comparison
  and open whole-map proof do not establish that all Gauss integrators fail.
- Higher-order GC, the canonical Lie switch, spline-provider FO and a variational
  event construction are unfinished or limited concepts, not demonstrated dead
  ends. No production flags are added for them by this cleanup.
- Random phase, ring/Bessel averaging alone at a sheet, and indefinite frozen
  planar cycles are not justified replacements for phase-dependent curved-sheet
  dynamics. They are not implemented as production options. Reopening a retired
  path needs a materially changed formulation and independent full-orbit evidence
  for its stated domain, not a renamed setting or a weaker tolerance.

### Verification of the removal

The retained exact symbolic rejection verifier and comparison behavioural tests
run unchanged in strength. The latter still check the sharp FO reference,
forward/backward events, shallow crossing-return detection, lift/inverse and
energy identities, nonaxisymmetric two-form order, and exact planar compression.
The native SIMPLE suite and before/after survivor comparison are recorded
in the cleanup PR. No production equations or namelist defaults change.

## Appendix: exact symbolic verification output

Run from the repository root with Python 3 and SymPy:
`python3 concept/symplectic_sheet/proofs/verify_rejected_models.py`.
The checks include the counterexamples as well as the supported identities.

```text
PASS T1a saltation pairing and constant-energy-offset exception
PASS T1b point-fixed obstruction, continuous A, level-1 slab witness
PASS T1b energy-root derivative and magnitude-jump obstruction
PASS T1b counterexample: rotating surface reset preserves form and energy
PASS T1c curl, GC determinant, signed layer-zero criterion
PASS T1d physical magnetic-moment jump
PASS T2 Hamiltonian difference, C2 layer area, combined error/cost
PASS T3 Boris norm preservation and explicit uniform-field crossing defect
PASS T3 phase means, variance, and independent-kick random-walk scaling
PASS T5a Hamiltonian energy, canonical symplecticity, saltation, symmetry
PASS T5b Poincare identities, residual sign, Gauss2 symplecticity/order four
PASS T5c explicit-RK obstruction, possible secular drift, free-particle exception
PASS all symbolic checks (analysis hypotheses remain in the document)
```
