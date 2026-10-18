# Why guiding-centre interface maps fail at SPECTRE sheets

Multi-volume SPECTRE tracing requires `integmode > 0` and
`crossing_level = 3`: guiding-centre (GC) motion inside each volume, full
orbit (FO) near each current sheet (`DOC/spectre-interface-crossing.md`).
This note derives why the alternatives are not valid models of a sharp
sheet: fixed-$\mu$ GC resets across the sheet (the former
`crossing_level = 0` point-fixed reset and `crossing_level = 1` tangential
reset, also with RK45 GC), a smoothed Boris layer, and unresolved
single-branch steps. The algebra is checked exactly by
`concept/symplectic_sheet/proofs/verify_rejected_models.py` (SymPy); analytic
estimates and cited theorems are argued here.

## Setting

SI units: charge $q\neq0$, mass $m$, parallel velocity $u$,
$\mu = m v_\perp^2/(2B)$, $B = |\mathbf B| > 0$ (Gaussian: $q \to q/c$). The
stationary sheet $\Sigma$ has $\mathbf B_\pm\cdot n = 0$ and a finite
tangential jump. $A$ is continuous and piecewise smooth up to $\Sigma$ with
$\nabla\times A = \mathbf B_\pm$; only its derivatives jump. Crossings are
isolated and transversal, finitely many per bounded time, without grazing.

## T1. GC resets

**T1a (composition criterion).** On the fixed-$\mu$ GC phase spaces let
$\Gamma_\pm = (qA + m u b_\pm)\cdot dx$, $\omega_\pm = -d\Gamma_\pm$,
$H_\pm = m u^2/2 + \mu B_\pm$, with nondegenerate $\omega_\pm$, transversal
flow and a smooth invertible reset $R: S_-\to S_+$ between event
hypersurfaces. The flow–reset–flow map at fixed physical time is symplectic
iff

$$R^*(\omega_+|_{S_+}) = \omega_-|_{S_-} \quad\text{and}\quad
d(H_+\circ R - H_-)|_{TS_-} = 0 .$$

*Proof.* Split a variation $v = t + aX_-$, $t\in TS_-$. The saltation map
sends it to $Sv = DR\,t + aX_+$. With $w = s + bX_-$ and
$\Delta H = H_+\circ R - H_-$,

$$\omega_+(Sv,Sw) - \omega_-(v,w) = (R^*\omega_+ - \omega_-)(t,s)
  + a\,d\Delta H(s) - b\,d\Delta H(t),$$

using $\iota_{X_\pm}\omega_\pm = dH_\pm$; $t,s,a,b$ are independent. Within
energy-preserving resets, symplecticity is presymplecticity on the event
surface. Energy equality itself is not necessary (a constant shift of $H_+$
leaves its flow unchanged).

**T1b (point-fixed resets).** Let $R$ keep the surface point and $\mu$ and
map $u\mapsto f(x,u)$. Presymplecticity forces $f_u b_+ = b_-$, hence
$b_+ = \pm b_-$ and $f_u = \pm1$; with energy conservation and $\mu>0$ also
$B_+ = B_-$.

*Proof.* Evaluate the pulled-back two-form on $(\partial_u, t)$ for tangent
$t$: $-m f_u b_+\cdot t = -m b_-\cdot t$ for all $t$. Unit norms give
$f_u^2 = 1$; differentiating energy equality gives $f f_u = u$, so
$f^2 = u^2$ and $\mu(B_+ - B_-) = 0$. The energy root
$f^2 = u^2 - 2\mu(B_+ - B_-)/m$ therefore fails at any direction jump.

Tangential resets escape T1b: for the sheet $z=0$ with
$b_- = (1,0,0)$, $b_+ = (\cos\phi,\sin\phi,0)$, $A_- = B(0,-z,0)$,
$A_+ = B(z\sin\phi,-z\cos\phi,0)$, the in-plane rotation
$R(x,y,u) = (Q_\phi(x,y),u)$ is presymplectic and energy preserving. It is an
algebraic counterexample, not a physical displacement. The tangential reset
with a scalar energy solve still has no general two-form identity: in this
uniform slab the energy residual vanishes at zero displacement, the reset is
the identity, and its two-form defect is $m(b_+ - b_-)\cdot dx\wedge du \neq 0$.

**T1c (loss of GC ordering).** For $b(z) = (\cos\alpha(z),\sin\alpha(z),0)$
and constant $B$,

$$b\cdot\nabla\times b = -\alpha', \qquad
B^*_\parallel = B + \frac{m u}{q}\,b\cdot\nabla\times b = B(1 - \rho_\parallel\alpha'),
\qquad \rho_\parallel = \frac{m u}{qB},$$

and $\det\omega = (mqB^*_\parallel)^2$ in $(x,y,z,u)$. Let
$\alpha = \alpha_- + \Delta\alpha\,s(z/L)$ on $[0,L]$ with monotone $C^1$
$s$, $s(0)=0$, $s(1)=1$, $s'(0)=s'(1)=0$. If $\rho_\parallel\Delta\alpha > 0$
and $L < \rho_\parallel\Delta\alpha$, $B^*_\parallel$ has a zero (since
$\max s'\ge1$), where the GC two-form degenerates; at fixed scaled position
with $s'>0$, $B^*_\parallel/B$ diverges as $L\to0$. For
$\rho_\parallel\Delta\alpha < 0$ there is no zero, but the slow-variation
conditions $\rho|\nabla b|\ll1$, $\rho|\nabla\ln B|\ll1$ require
$\rho \ll L/\max(|\Delta\alpha|, |\ln(B_+/B_-)|)$. A sharp sheet has $L = 0$,
so any finite gyroradius violates the GC ordering at any nonzero jump.

**T1d (μ is not conserved).** $x$ and $v$ of the full orbit are continuous
through a bounded magnetic step. In the equal-$B$ slab of T1b with
$v = (v_t, 0, v_n)$, $v_n\neq0$,

$$\mu_+ - \mu_- = \frac{m v_t^2\sin^2\phi}{2B},$$

a finite jump independent of any step size. No fixed-$\mu$ reset reproduces
this, whatever its two-form properties.

## T2. A smoothed Boris layer changes the problem

Replace the jump by a layer of width $w$ and compare with the sharp field at
matched entry $x, v$, $|v| = v_0$, $|v_n|\ge\eta v_0$. With
$\Omega = |q|B_\text{ref}/m$, $\rho = v_0/\Omega$, $\epsilon = w/\rho$, the
exit-velocity discrepancy is $O(\epsilon)$. Using $z$ as independent
variable, $dv/dz = (q/m)(v\times\mathbf B)/v_n$, so to first order

$$\frac{\Delta I}{m} = \frac{q}{m v_{n0}}\,v_0\times\int(\mathbf B_w - \mathbf B_0)\,dz
  + O(v_0\epsilon^2).$$

For $\mathbf B_w = \mathbf B_- + \Delta\mathbf B\,s(z/w)$,
$s(\xi) = 6\xi^5 - 15\xi^4 + 10\xi^3$, the integral is $-w\Delta\mathbf B/2$:
the leading term does not vanish in general. Exact witness:
$\mathbf B = B_y(z)e_y$, $q>0$, $B_- = B_\text{ref}$, $B_+ = 2B_\text{ref}$,
entry $v = (0,0,v_0)$, $w = \rho/3$: the sharp exit has $v_x = -2v_0/3$, the
layer exit $v_x = -v_0/2$. The Hamiltonian change in matched gauges is only
$O(\epsilon)$ ($H_w - H_0 = -q\pi\cdot\delta A/m + q^2|\delta A|^2/(2m)$,
$\delta A = O(B_\text{ref}w)$), while the force changes at order one.

With step $h$, $\delta = \Omega h$ and $\kappa = v_0h/w$, resolved Boris
needs $\kappa\ll1$ and $\delta\ll1$; the per-crossing error relative to the
sharp target is $O(\epsilon + \delta^2/\epsilon)$ at cost $O(1/\delta)$ steps
per gyration. Width and step must both converge; a fixed width converges to
a different field.

## T3. A single-branch step across a sheet is first order

Let a step $h$ cross at fraction $\theta\in(0,1)$ and use one side's field
throughout. With $[[f]] = f_+ - f_-$ the leading local error is
$-(1-\theta)h[[f]]$ (incoming branch) or $\theta h[[f]]$ (outgoing branch).
For uniform fields with $G_\pm v = (q/m)v\times\mathbf B_\pm$, the exact map
is $e^{(1-\theta)hG_+}e^{\theta hG_-} = I + h(\theta G_- + (1-\theta)G_+)
+ O(h^2)$, while Boris applies the Cayley rotation
$C_h(G) = (I - hG/2)^{-1}(I + hG/2) = I + hG + O(h^2)$ of one side.
$C_h^TC_h = I$: exact speed conservation does not remove the error. A
midpoint-branch choice has coefficient $\theta$ or $-(1-\theta)$, with mean
zero and variance $1/12$ for uniform $\theta$. If such kicks of size
$D = |[[f]]|$ were independent, $N$ crossings would accumulate
$\sqrt N hD/\sqrt{12}$; coherent accumulation $O(NhD)$ is also possible.
Exact event splitting has no such term.

## T4. RK45 GC does not repair the reset

Decreasing the RK45 tolerance with a fixed-$\mu$ reset converges to the exact
GC flows composed with the same reset. T1c and T1d hold for exact
per-volume flows, so the model error remains.

## T5. Full-orbit alternatives

**T5a (exact FO flow).** $H = |p - qA|^2/(2m)$ has a symplectic,
energy-preserving flow across the sheet. In each volume $\dot Z = J\nabla H$.
Continuity of $A$ makes $x$, $p$, $H$ continuous at the event; the gradient
jump of $H$ is normal to the event surface, so the saltation matrix is
$S = \bigl(\begin{smallmatrix}I&0\\ c\,nn^T&I\end{smallmatrix}\bigr)$ with
$S^TJS = J$. A cyclic coordinate in a symmetry-compatible matched gauge keeps
its canonical momentum; this is why SIMPLE maps $p_\zeta$ across sheets by
the $A_\zeta$ jumps in axisymmetric fields.

**T5b (resolved Poincaré Gauss 2).** For smooth $K = g(x)(H - E)$ with
$0 < g_\text{min} \le g$, two-stage Gauss–Legendre
($b_i = 1/2$, $a_{11} = a_{22} = 1/4$, $a_{12} = 1/4 - \sqrt3/6$,
$a_{21} = 1/4 + \sqrt3/6$) is symplectic of order four in fictitious time,
and under the usual analyticity and compactness assumptions
$|K_n - K_0| \le Ch_\tau^4$ over exponentially long times (Hairer, Lubich and
Wanner 2006). On a shell $K = k$, $X_K/g = X_{H - k/g}$, so
$|H_n - E| \le Ch_\tau^4/g_\text{min}$. The step must resolve every scale of
$K$, $h_\tau \ll L/(g|v_n|)$; smooth backward-error analysis does not cover a
sharp sheet.

**T5c (event-resolved FO RK45).** With finitely many uniformly transversal
crossings, piecewise smooth Lipschitz fields and exact event landing, RK45
converges as step and tolerance go to zero (induction over events, using
continuous dependence of event times on the incoming state). It is not
symplectic: explicit RK cannot satisfy $2b_ia_{ii} = b_i^2$ with
$a_{ii} = 0$. For the harmonic oscillator, Dormand–Prince
($R(z) = \sum_{j\le5}z^j/j! + z^6/600$) gives
$|R(i\eta)|^2 = 1 - \eta^6/1800 + \eta^8/1600 - \eta^{10}/14400 + \eta^{12}/360000$,
an energy drift linear in the number of steps for small $N\eta^6$. SIMPLE
uses this method for the bounded FO phases near sheets.

## References

- E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*,
  2nd ed., Springer (2006), https://doi.org/10.1007/3-540-30666-8
- E. Hairer, Variable time step integration with symplectic methods,
  *Appl. Numer. Math.* 25, 219 (1997),
  https://doi.org/10.1016/S0168-9274(97)00061-5
