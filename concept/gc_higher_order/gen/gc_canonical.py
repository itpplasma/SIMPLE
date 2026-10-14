"""Canonical six-dimensional GC flow and implicit symplectic RK steps.

Canonical ordering is (q1,q2,zeta,P1,P2,J), J=epsilon*mu. Physics is
evaluated in (r,u,mu,q1,q2); the radial inverse and both RK stages use
analytic Newton derivatives. This module does not modify SIMPLE.
"""
import importlib.util
from pathlib import Path
import numpy as np

_path = Path(__file__).resolve().parents[1]/'generated/canonical_kernels.py'
_spec = importlib.util.spec_from_file_location('canonical_kernels', _path)
_kernels = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_kernels)
OMEGA = np.block([[np.zeros((3, 3)), np.eye(3)],
                  [-np.eye(3), np.zeros((3, 3))]])
ACTIVE = np.array([0, 1, 3, 4, 5])


class CanonicalGC:
    def __init__(self, field, epsilon, order):
        if order not in (1, 11, 2):
            raise ValueError('canonical GC requires order 1, 1+ or 2')
        self.field, self.epsilon, self.order = field, epsilon, order
        self.kernel = getattr(_kernels, field.name)
        self.guess = None
        self.calls = self.inverse_iterations = 0

    def values(self, w):
        self.calls += 1
        return np.asarray(self.kernel(*w, self.epsilon, int(self.order != 1)))

    def chart(self, X):
        x, y, z = X
        if self.field.name.startswith('sheet'):
            R = np.hypot(x, y)
            r = np.hypot(R-10, z)
            th, ph = np.arctan2(z, R-10), np.arctan2(y, x)
            e1 = np.array([np.cos(th)*np.cos(ph), np.cos(th)*np.sin(ph), np.sin(th)])
        elif self.field.name == 'constant_curvature_on_fieldline':
            r, th, ph = np.hypot(x, y), np.arctan2(y, x), z
            e1 = np.array([0., 0., 1.])
        else:
            r, th, ph = x, y, z
            e1 = np.array([1., 0., 0.])
        return np.array([r, th, ph]), e1

    def cartesian(self, r, th, ph):
        if self.field.name.startswith('sheet'):
            R = 10+r*np.cos(th)
            return np.array([R*np.cos(ph), R*np.sin(ph), r*np.sin(th)])
        if self.field.name == 'constant_curvature_on_fieldline':
            return np.array([r*np.cos(th), r*np.sin(th), ph])
        return np.array([r, th, ph])

    def frame_angle(self, X, natural):
        b = self.field.b(X)
        b /= np.linalg.norm(b)
        axis = np.eye(3)[self.field.axis-1]
        e1 = axis-b*np.dot(axis, b)
        e1 /= np.linalg.norm(e1)
        return np.arctan2(np.dot(natural, np.cross(b, e1)), np.dot(natural, e1))

    def to_canonical(self, gc):
        q, natural = self.chart(gc[:3])
        r, th, ph = q
        w = np.array([r, gc[3], gc[4], th, ph])
        a = self.values(w)
        psi = self.frame_angle(gc[:3], natural)
        self.guess = w[:2].copy()
        return np.array([th, ph, -gc[5]+psi-a[93], a[0], a[1], self.epsilon*gc[4]])

    def inverse(self, y):
        if y[5] <= 0:
            raise ValueError('nonpositive gyroaction')
        if self.guess is None:
            raise ValueError('initialize from a GC state before inversion')
        w = np.array([*self.guess, y[5]/self.epsilon, y[0], y[1]])
        for _ in range(12):
            a = self.values(w)
            residual = a[:2]-y[3:5]
            D = a[2:12].reshape(2, 5)[:, :2]
            delta = np.linalg.solve(D, residual)
            self.inverse_iterations += 1
            if np.max(np.abs(delta)) < 3e-14*(1+np.max(np.abs(w[:2]))):
                self.guess = w[:2].copy()
                return w, a
            w[:2] -= delta
            if self.field.name.startswith('sheet') and not 0 < w[0] <= 1:
                raise RuntimeError('outside certified smooth sheet chart 0<r<=1')
        raise RuntimeError(('radial Newton failed', self.field.name, residual))

    def to_gc(self, y):
        w, a = self.inverse(y)
        X = self.cartesian(w[0], w[3], w[4])
        _, natural = self.chart(X)
        psi = self.frame_angle(X, natural)
        return np.r_[X, w[1], w[2], -y[2]-a[93]+psi]

    def evaluate(self, y, hessian=True):
        w, a = self.inverse(y)
        M = np.zeros((5, 5))
        M[0, 3] = M[1, 4] = 1
        M[2:4] = a[2:12].reshape(2, 5)
        M[4, 2] = self.epsilon
        W = np.linalg.inv(M)
        gv = W.T@a[63:68]
        g = np.zeros(6)
        g[ACTIVE] = gv
        f = OMEGA@g
        if not hessian:
            return a[62], f
        h = a[68:93].reshape(5, 5) - np.einsum(
            'a,aij->ij', gv[2:4], a[12:62].reshape(2, 5, 5))
        hv = W.T@h@W
        hy = np.zeros((6, 6))
        hy[np.ix_(ACTIVE, ACTIVE)] = hv
        return a[62], f, OMEGA@hy


def coefficients(method):
    if method == 'midpoint':
        return np.array([[.5]]), np.ones(1)
    if method == 'gauss2':
        v = np.sqrt(3)/6
        return np.array([[.25, .25-v], [.25+v, .25]]), np.full(2, .5)
    raise ValueError(method)


def step(model, y, dt, method='gauss2', tangent=False):
    A, b = coefficients(method)
    stages = len(b)
    _, f0 = model.evaluate(y, hessian=False)
    Z = np.tile(y, (stages, 1))+dt*A.sum(axis=1)[:, None]*f0
    identity = np.eye(6*stages)
    for iteration in range(10):
        evaluations = [model.evaluate(z) for z in Z]
        F = np.array([a[1] for a in evaluations])
        D = np.array([a[2] for a in evaluations])
        residual = Z-y-dt*A@F
        matrix = identity.copy().reshape(stages, 6, stages, 6)
        for i in range(stages):
            for j in range(stages):
                matrix[i, :, j, :] -= dt*A[i, j]*D[j]
        matrix = matrix.reshape(6*stages, 6*stages)
        correction = np.linalg.solve(matrix, residual.ravel()).reshape(stages, 6)
        # Test corrections relative to coordinate roundoff, not to a
        # large magnetic canonical momentum or accumulated fast phase.
        scale = 1+np.abs(Z)
        if np.max(np.abs(correction)/scale) < 2e-14:
            break
        Z -= correction
    else:
        raise RuntimeError(('stage Newton failed', method, residual))
    out = y+dt*b@F
    if not tangent:
        return out
    stage_deriv = np.linalg.solve(matrix, np.tile(np.eye(6), (stages, 1)))
    T = np.eye(6)+dt*sum(b[i]*D[i]@stage_deriv[6*i:6*i+6] for i in range(stages))
    return out, T, iteration+1
