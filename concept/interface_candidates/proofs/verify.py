"""Exact identities and limited canonical witness, not a general GC proof."""
import sympy as s

R, Z, phi, R0, a, F, delta, n = s.symbols(
    "R Z phi R0 a F delta n", real=True
)
r2 = (R-R0)**2 + Z**2
H = (R-R0)*s.cos(n*phi) + Z*s.sin(n*phi)
Ar, Ap, Az = delta*H*(R-R0), a*r2/R, -F*s.log(R)+delta*H*Z
B = s.Matrix([s.diff(Az, phi)/R-s.diff(Ap, Z),
              s.diff(Ar, Z)-s.diff(Az, R),
              (s.diff(R*Ap, R)-s.diff(Ar, phi))/R])
div = s.diff(R*B[0], R)/R+s.diff(B[1], phi)/R+s.diff(B[2], Z)
assert s.simplify(div) == 0
assert s.simplify(B[0]*(R-R0)+B[2]*Z) == 0
am, ap, r0 = s.symbols("am ap r0", real=True)
psi_minus, psi_plus = am*r2, ap*r2+(am-ap)*r0**2
assert s.expand((psi_plus-psi_minus).subs(Z**2, r0**2-(R-R0)**2)) == 0
print("PASS T-A/T-B curl, div B=0, B.n=0, continuous vector potential")

# Energy cancellation in the first-order generators is algebraic for any
# smooth field. rk, rg are rho.kappa, rho.grad B; al is tau+a1:grad b.
u, mu, b, rk, rg, al = s.symbols("u mu b rk rg al", real=True)
gu = -u*rk+mu*al
gm = (mu*rg+u*u*rk-u*mu*al)/b
assert s.expand(u*gu+b*gm-mu*rg) == 0
print("PASS first-order energy expansion: u G1u+B G1mu=mu rho.grad B")

# Constant B: the exact helix and its time inverse.
t, om = s.symbols("t om", positive=True)
vx, vy, vz = s.symbols("vx vy vz", real=True)
v = s.Matrix([vx*s.cos(om*t)+vy*s.sin(om*t),
              vy*s.cos(om*t)-vx*s.sin(om*t), vz])
x = s.Matrix([(vx*s.sin(om*t)+vy*(1-s.cos(om*t)))/om,
              (vy*s.sin(om*t)-vx*(1-s.cos(om*t)))/om, vz*t])
assert all(s.simplify(q) == 0 for q in x.diff(t)-v)
assert all(s.simplify(q) == 0 for q in v.diff(t)-v.cross(s.Matrix([0, 0, om])))
assert s.trigsimp(v.dot(v)-(vx*vx+vy*vy+vz*vz)) == 0
aa, bb, alpha, theta = s.symbols("aa bb alpha theta", real=True)
assert s.expand_trig(s.sin(theta+alpha)).subs(
    {s.cos(alpha): aa, s.sin(alpha): -bb}) == aa*s.sin(theta)-bb*s.cos(theta)
print("PASS helix Lorentz equations, speed, analytic plane root identity")

# Matched continuous A => exact planar event flow has symplectic saltation.
nn = s.Matrix(s.symbols("nx ny nz"))
c = s.symbols("c")
I, O = s.eye(3), s.zeros(3)
J = O.row_join(I).col_join((-I).row_join(O))
salt = I.row_join(O).col_join((c*nn*nn.T).row_join(I))
assert salt.T*J*salt == J
print("PASS exact planar fixed-time event saltation is canonical symplectic")

Bx, By, z, Px, Py, Pz = s.symbols("Bx By z Px Py Pz", real=True)
mag2 = Bx**2+By**2
aa = s.Matrix([By, -Bx])
pt = s.Matrix([Px, Py])
center = pt.dot(aa)/mag2
Umin = (pt.dot(pt)-pt.dot(aa)**2/mag2)/2
Hplane = ((Px-By*z)**2+(Py+Bx*z)**2+Pz**2)/2
assert s.factor(Hplane-(Pz**2/2+mag2*(z-center)**2/2+Umin)) == 0
assert s.simplify(Hplane.subs(z, 0)-(Px**2+Py**2+Pz**2)/2) == 0
# Translation invariance preserves Px,Py; energy then fixes |Pz| at the
# sheet. Two returning excursions repeat velocity and translate tangentially.
print("PASS planar normal oscillator, crossing momenta, complete-cycle compression")

# A fully explicit nonuniform-field canonical witness, B=(0,1+k X,0).
# Convention: velocity phase, Gamma_gc=(A+u b).dX-mu dphase.
# This tests the spatial rho1 term that is absent from the older prototype.
X, Y, ZZ, U, M, ph, k = s.symbols("X Y ZZ U M ph k", real=True)
amp = s.sqrt(2*M)
sn, co = s.sin(ph), s.cos(ph)
q = s.Matrix([X, Y, ZZ, U, M, ph])
xx = X-amp*sn+k*(amp*sn*X/2-M*(1+co**2))
zz = ZZ-amp*co+k*(amp*co*X/2+M*co*sn)
vxx = amp*co+k*(amp*co*X/2-2*M*co*sn)
vzz = -amp*sn+k*(-amp*sn*X/2-2*M*co**2)
T = s.Matrix([xx, Y, zz, vxx, U, vzz])
D = T.jacobian(q)
omega = s.zeros(6)
# d[(A+v).dx], with Az=-x-k x^2/2.
omega[0, 2], omega[2, 0] = -(1+k*xx), 1+k*xx
for i in range(3):
    omega[i, i+3], omega[i+3, i] = -1, 1
gc = s.zeros(6)
gc[0, 2], gc[2, 0] = -(1+k*X), 1+k*X
gc[1, 3], gc[3, 1] = -1, 1
gc[4, 5], gc[5, 4] = -1, 1
def trunc(expr):
    return s.expand(expr).coeff(k, 0)+k*s.expand(expr).coeff(k, 1)
defect = D.T*omega*D-gc
for i in range(6):
    for j in range(6):
        assert s.trigsimp(trunc(defect[i, j])) == 0, (i, j)
energy = (vxx**2+U**2+vzz**2)/2-(U**2/2+M*(1+k*X))
assert s.trigsimp(trunc(energy)) == 0
print("PASS nonuniform slab T1 pullback and energy to O(k^2)")

# Gyrophase dependence of the physically continuous velocity reset.
ang, vp, vn, B0 = s.symbols("ang vp vn B0", real=True)
mu_m = vn**2/(2*B0)
mu_p = (vn**2+vp**2*s.sin(ang)**2)/(2*B0)
assert s.simplify(mu_p-mu_m-vp**2*s.sin(ang)**2/(2*B0)) == 0
print("PASS finite phase-dependent magnetic-moment jump; fixed-mu map fails")
print("PASS all stated algebraic identities (general-field ordering is separate)")
