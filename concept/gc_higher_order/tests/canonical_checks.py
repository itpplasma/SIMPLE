#!/usr/bin/env python3
"""Behavioral canonical-flow, analytic derivative and discrete-map oracles."""
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'gen'))
from gc_canonical import CanonicalGC, step, OMEGA
from orbit_study import Field, fields


def difference(fun, y, h):
    return np.column_stack([(fun(y+h*e)-fun(y-h*e))/(2*h) for e in np.eye(6)])


def main():
    rows = []
    for name, entry in fields.items():
        field = Field(name, entry)
        for order in (1, 2):
            model = CanonicalGC(field, .04, order)
            # Off-symmetry initial point and phase; action is varied by the
            # derivative oracle, so checking only a fixed-mu leaf cannot pass.
            gc = np.r_[field.initial, .35, .25, .7]
            y = model.to_canonical(gc)
            H, f, D = model.evaluate(y)
            J = difference(model.to_gc, y, 2e-6)
            rhs_error = np.max(np.abs(J@f-field.evaluate(2, order, .04, gc)))
            hessian_error = np.max(np.abs(D-difference(lambda z: model.evaluate(z)[1], y, 2e-6)))
            assert rhs_error < 2e-8, (name, order, rhs_error)
            assert hessian_error < 3e-7, (name, order, hessian_error)
            assert abs(H-field.evaluate(3, order, .04, gc)[0]) < 3e-14
            assert np.max(np.abs(model.to_gc(y)-gc)) < 3e-13
            for method in ('midpoint', 'gauss2'):
                out, T, iterations = step(model, y, .13, method, tangent=True)
                # Independent finite differences of the actual converged map.
                Td = difference(lambda z: step(model, z, .13, method), y, 2e-6)
                tangent_error = np.max(np.abs(T-Td))
                defect = np.max(np.abs(T.T@OMEGA@T-OMEGA))
                fd_defect = np.max(np.abs(Td.T@OMEGA@Td-OMEGA))
                assert tangent_error < 3e-7, (name, method, tangent_error)
                assert defect < 2e-12, (name, method, defect)
                assert fd_defect < 3e-7, (name, method, fd_defect)
                assert out[5] == y[5]
                if name.startswith('sheet') or name == 'constant_curvature_on_fieldline':
                    index = 4 if name.startswith('sheet') else 3
                    assert out[index] == y[index], 'axisymmetric canonical momentum'
                # Time reversal is a separate check of the solved map.
                back = step(model, out, -.13, method)
                assert np.max(np.abs(back-y)) < 3e-12
                rows.append(dict(field=name, order=order, method=method,
                    rhs_error=float(rhs_error), hessian_error=float(hessian_error),
                    tangent_fd_error=float(tangent_error), symplectic_defect=float(defect),
                    fd_symplectic_defect=float(fd_defect), newton_iterations=iterations))
        print(name, 'canonical oracles PASS', flush=True)
    (ROOT/'results/canonical_checks.json').write_text(json.dumps(dict(status='PASS', rows=rows), indent=2)+'\n')


if __name__ == '__main__':
    main()
