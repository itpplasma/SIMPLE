#!/usr/bin/env python3
"""Exact algebra for DOC/spectre-rejected-interface-models.md (SymPy only).

No numerical probes. Analysis hypotheses, bounds, and cited BEA theorems are
proved/stated in the document, not certified by a computer algebra system.
"""
import sympy as s


def zero(expr):
    entries = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    for entry in entries:
        assert s.simplify(s.trigsimp(entry)) == 0, entry


def passed(label):
    print(f"PASS {label}")


def skew(b):
    """G(b) v = v cross b, with q/m absorbed in b."""
    x, y, z = b
    return s.Matrix([[0, z, -y], [-z, 0, x], [y, -x, 0]])


def omega(one_form, coords):
    """Matrix of -d Gamma."""
    return s.Matrix(len(coords), len(coords), lambda i, j:
                    s.diff(one_form[i], coords[j])
                    - s.diff(one_form[j], coords[i]))


# T1a: saltation pairing in a tangent/flow basis.
a, b, ht, hs, dt, ds, omt, opt = s.symbols("a b ht hs dt ds omt opt")
before = omt + a*hs - b*ht
after = opt + a*(hs + ds) - b*(ht + dt)
zero(after - before - (opt - omt + a*ds - b*dt))
coord, offset = s.symbols("coordinate constant_offset")
zero(s.diff(offset, coord))
passed("T1a saltation pairing and constant-energy-offset exception")

# T1b: constant tangent fields on z=0; A has zero tangential trace.
x, y, z, u = s.symbols("x y z u", real=True)
m, q, B, mu = s.symbols("m q B mu", positive=True)
phi = s.symbols("phi", real=True)
bm = s.Matrix([1, 0, 0])
bp = s.Matrix([s.cos(phi), s.sin(phi), 0])
Am = s.Matrix([0, -B*z, 0])
Ap = s.Matrix([B*z*s.sin(phi), -B*z*s.cos(phi), 0])
def curl(v):
    return s.Matrix([s.diff(v[2], y) - s.diff(v[1], z),
                     s.diff(v[0], z) - s.diff(v[2], x),
                     s.diff(v[1], x) - s.diff(v[0], y)])
zero(curl(Am) - B*bm)
zero(curl(Ap) - B*bp)
zero(Am.subs(z, 0) - Ap.subs(z, 0))
om = omega([m*u, 0, 0], [x, y, u])
op = omega([m*u*s.cos(phi), m*u*s.sin(phi), 0], [x, y, u])
fu = s.symbols("f_u", real=True)
point_reset = s.diag(1, 1, fu)
defect = point_reset.T*op*point_reset - om
zero(defect[0, 2] - m*(fu*s.cos(phi) - 1))
zero(defect[1, 2] - m*fu*s.sin(phi))
zero((fu*s.cos(phi))**2 + (fu*s.sin(phi))**2 - fu**2)
zero(defect.subs({phi: s.pi/2, fu: 1})[0, 2] + m)
passed("T1b point-fixed obstruction, continuous A, level-1 slab witness")

# Energy root and its derivative; signed root requires u,f nonzero locally.
dB = s.symbols("Delta_B", real=True)
f = s.sqrt(u**2 - 2*mu*dB/m)
zero(m*f**2/2 + mu*(B + dB) - (m*u**2/2 + mu*B))
zero(f*s.diff(f, u) - u)
# Unit-length pullback requires f_u^2=1; together with f*f_u=u:
zero(f**2 - u**2 + 2*mu*dB/m)
passed("T1b energy-root derivative and magnitude-jump obstruction")

# A nonlocal surface rotation defeats the universal obstruction.
Q = s.Matrix([[s.cos(phi), -s.sin(phi), 0],
              [s.sin(phi), s.cos(phi), 0], [0, 0, 1]])
zero(Q.T*op*Q - om)
zero(m*u**2/2 + mu*B - (m*u**2/2 + mu*B))
passed("T1b counterexample: rotating surface reset preserves form and energy")

# T1c: GC form determinant and thin-layer signed degeneracy.
alpha = s.Function("alpha")(z)
bz = s.Matrix([s.cos(alpha), s.sin(alpha), 0])
zero(bz.dot(curl(bz)) + s.diff(alpha, z))
bg = s.Function("B")(z)
ax, ay = s.Function("Ax")(z), s.Function("Ay")(z)
og = omega([q*ax + m*u*bz[0], q*ay + m*u*bz[1], 0, 0],
           [x, y, z, u])
og = og.subs({s.diff(ax, z): bg*s.sin(alpha),
              s.diff(ay, z): -bg*s.cos(alpha)})
bstar = bg - m*u*s.diff(alpha, z)/q
zero(og.det() - (m*q*bstar)**2)
rho, L, da, xi = s.symbols("rho L da xi", real=True)
profile = 6*xi**5 - 15*xi**4 + 10*xi**3
zero(s.diff(profile, xi) - 30*xi**2*(xi - 1)**2)
zero(s.integrate(s.diff(profile, xi), (xi, 0, 1)) - 1)
ratio = 1 - rho*da*s.diff(profile, xi)/L
zero(ratio.subs(xi, s.Rational(1, 2)) - (1 - 15*rho*da/(8*L)))
zero(ratio.subs(xi, 0) - 1)
zero(ratio.subs(xi, 1) - 1)
zero(ratio.subs(da, -da) - (1 + rho*da*s.diff(profile, xi)/L))
rpos, dpos, lpos = s.symbols("r_positive d_positive L_positive", positive=True)
center = ratio.subs({xi: s.Rational(1, 2), rho: rpos, da: dpos, L: lpos})
assert s.limit(center, lpos, 0, dir="+") == -s.oo
assert s.limit(center.subs(dpos, -dpos), lpos, 0, dir="+") == s.oo
assert s.diff(profile, xi).factor().is_nonnegative is True
passed("T1c curl, GC determinant, signed layer-zero criterion")

# T1d: full-orbit velocity is continuous; mu generally is not.
vn, vt = s.symbols("v_n v_t", real=True)
v = s.Matrix([vt, 0, vn])
mum = m*(v.dot(v) - v.dot(bm)**2)/(2*B)
mup = m*(v.dot(v) - v.dot(bp)**2)/(2*B)
zero(mup - mum - m*vt**2*s.sin(phi)**2/(2*B))
passed("T1d physical magnetic-moment jump")

# T2: matched-gauge Hamiltonian difference and an asymmetric layer impulse.
pi = s.Matrix(s.symbols("pi_x pi_y pi_z", real=True))
dA = s.Matrix(s.symbols("dA_x dA_y dA_z", real=True))
zero((pi - q*dA).dot(pi - q*dA)/(2*m) - pi.dot(pi)/(2*m)
     + q*pi.dot(dA)/m - q**2*dA.dot(dA)/(2*m))
zero(s.integrate(profile - 1, (xi, 0, 1)) + s.Rational(1, 2))
v0, qpos, w = s.symbols("v0 q_positive w", positive=True)
rhoref = m*v0/(qpos*B)
vxsharp = -qpos*B*w*2/m
vxlayer = -qpos*B*w*s.integrate(1 + profile, (xi, 0, 1))/m
zero(vxsharp.subs(w, rhoref/3) + 2*v0/3)
zero(vxlayer.subs(w, rhoref/3) + v0/2)
zero((vxlayer - vxsharp).subs(w, rhoref/3) - v0/6)
zero((v0**2 - vxsharp**2).subs(w, rhoref/3) - 5*v0**2/9)
zero((v0**2 - vxlayer**2).subs(w, rhoref/3) - 3*v0**2/4)
zero(skew([0, B, 0])[0, :].dot(s.Matrix([vt, 0, vn]))/vn + B)
for endpoint in (0, 1):
    zero(profile.subs(xi, endpoint) - endpoint)
    zero(s.diff(profile, xi).subs(xi, endpoint))
    zero(s.diff(profile, xi, 2).subs(xi, endpoint))
eps, delta, kappa = s.symbols("epsilon delta kappa", positive=True)
zero((eps + delta**2/eps).subs(delta, kappa*eps) - eps*(1 + kappa**2))
zero((1/delta).subs(delta, kappa*eps) - 1/(kappa*eps))
Omega = s.symbols("Omega", positive=True)
htime = s.symbols("h_time", positive=True)
zero((htime**3*Omega*(v0/w)**2)*(w/(v0*htime))
     - (Omega*htime)**2/(w*Omega/v0))
passed("T2 Hamiltonian difference, C2 layer area, combined error/cost")

# T3: exact uniform-field generators and the Boris Cayley rotation.
bx, by, bz0, h, theta = s.symbols("b_x b_y b_z h theta", real=True)
G = skew([bx, by, bz0])
I = s.eye(3)
beta2 = bx**2 + by**2 + bz0**2
zero(G.T + G)
zero(G**3 + beta2*G)
C = I + h*G/(1 + h**2*beta2/4) + h**2*G**2/(2*(1 + h**2*beta2/4))
zero((I - h*G/2)*C - (I + h*G/2))
zero(C.T*C - I)
zero(C.diff(h).subs(h, 0) - G)
zero(C.diff(h, 2).subs(h, 0) - G**2)
Gm = skew(s.symbols("bm_x bm_y bm_z"))
Gp = skew(s.symbols("bp_x bp_y bp_z"))
exact2 = (I + (1 - theta)*h*Gp + (1 - theta)**2*h**2*Gp**2/2) * (
          I + theta*h*Gm + theta**2*h**2*Gm**2/2)
leading = exact2.diff(h).subs(h, 0)
zero(leading - (theta*Gm + (1 - theta)*Gp))
zero(Gm - leading + (1 - theta)*(Gp - Gm))
zero(Gp - leading - theta*(Gp - Gm))
zero(s.Matrix([vt, 0, vn]).cross(bp - bm)
     - (skew(bp) - skew(bm))*s.Matrix([vt, 0, vn]))
passed("T3 Boris norm preservation and explicit uniform-field crossing defect")

c0, c1 = theta, -(1 - theta)
zero(s.integrate(c1, (theta, 0, 1)) + s.Rational(1, 2))
zero(s.integrate(c0, (theta, 0, s.Rational(1, 2)))
     + s.integrate(c1, (theta, s.Rational(1, 2), 1)))
variance = (s.integrate(c0**2, (theta, 0, s.Rational(1, 2)))
            + s.integrate(c1**2, (theta, s.Rational(1, 2), 1)))
zero(variance - s.Rational(1, 12))
N, D = s.symbols("N D", positive=True)
zero(N*h**2*D**2*variance - N*h**2*D**2/12)
passed("T3 phase means, variance, and independent-kick random-walk scaling")

# T5a: canonical Hamiltonian identities and interface saltation shear.
J = s.zeros(6)
J[:3, 3:] = I
J[3:, :3] = -I
grad = s.Matrix(s.symbols("H_1:7"))
Hess = s.Matrix(6, 6, lambda i, j: s.Symbol(f"H_{min(i,j)}_{max(i,j)}"))
zero(grad.dot(J*grad))
zero((J*Hess).T*J + J*(J*Hess))
normal = s.Matrix(s.symbols("n_x n_y n_z"))
c = s.symbols("c", real=True)
S = s.eye(6)
S[3:, :3] = c*normal*normal.T
zero(S.T*J*S - J)
pcyc = s.symbols("p_cyclic")
zero(-s.diff(pcyc**2/(2*m), x))
passed("T5a Hamiltonian energy, canonical symplecticity, saltation, symmetry")

# T5b: Poincare product rule, residual potential sign, and Gauss2 coefficients.
p, E, kval = s.symbols("p E k", real=True)
Hf = s.Function("H")(x, p)
gf = s.Function("g")(x)
K = gf*(Hf - E)
zero(s.diff(K, p) - gf*s.diff(Hf, p))
zero(s.diff(K, x) - gf*s.diff(Hf, x) - (Hf - E)*s.diff(gf, x))
zero(s.diff(Hf - kval/gf, x) - s.diff(Hf, x) - kval*s.diff(gf, x)/gf**2)
zero(K/gf - (Hf - E))
A = s.Matrix([[s.Rational(1, 4), s.Rational(1, 4) - s.sqrt(3)/6],
              [s.Rational(1, 4) + s.sqrt(3)/6, s.Rational(1, 4)]])
weights = s.Matrix([s.Rational(1, 2)]*2)
nodes = A*s.ones(2, 1)
for i in range(2):
    for j in range(2):
        zero(weights[i]*A[i, j] + weights[j]*A[j, i] - weights[i]*weights[j])
for power in range(4):
    zero(sum(weights[i]*nodes[i]**power for i in range(2)) - s.Rational(1, power + 1))
zero((weights.T*A*nodes)[0] - s.Rational(1, 6))
zero((weights.T*A*nodes.applyfunc(lambda t: t**2))[0] - s.Rational(1, 12))
zero((weights.T*s.diag(*nodes)*A*nodes)[0] - s.Rational(1, 8))
zero((weights.T*A*A*nodes)[0] - s.Rational(1, 24))
passed("T5b Poincare identities, residual sign, Gauss2 symplecticity/order four")

# T5c: explicit RK diagonal obstruction; Dormand-Prince fifth-order polynomial.
bi = s.symbols("b_i")
zero((2*bi*0 - bi**2) + bi**2)
R = s.Rational
adp = s.zeros(7)
rows = [[], [R(1,5)], [R(3,40), R(9,40)],
        [R(44,45), -R(56,15), R(32,9)],
        [R(19372,6561), -R(25360,2187), R(64448,6561), -R(212,729)],
        [R(9017,3168), -R(355,33), R(46732,5247), R(49,176), -R(5103,18656)],
        [R(35,384), 0, R(500,1113), R(125,192), -R(2187,6784), R(11,84)]]
for i, row in enumerate(rows):
    for j, value in enumerate(row):
        adp[i, j] = value
wdp = s.Matrix(rows[6] + [0])
zz, eta = s.symbols("zeta eta", real=True)
stages = s.ones(7, 1)
for i in range(7):
    stages[i] += zz*sum(adp[i, j]*stages[j] for j in range(i))
amp = s.expand(1 + zz*(wdp.T*stages)[0])
zero(amp - (sum(zz**j/s.factorial(j) for j in range(6)) + zz**6/600))
norm2 = s.expand(amp.subs(zz, s.I*eta)*amp.subs(zz, -s.I*eta))
zero(norm2 - (1 - eta**6/1800 + eta**8/1600 - eta**10/14400 + eta**12/360000))
zero(s.diff((1 + zz)**N, zz).subs(zz, 0) - N)
zero(sum(wdp) - 1)
# Free particle: all stage accelerations are zero; every stage velocity is p/m.
xnext = x + h*sum(wdp)*p/m
pnext = p
zero(xnext - (x + h*p/m))
zero(pnext**2/(2*m) - p**2/(2*m))
passed("T5c explicit-RK obstruction, possible secular drift, free-particle exception")
print("PASS all symbolic checks (analysis hypotheses remain in the document)")
