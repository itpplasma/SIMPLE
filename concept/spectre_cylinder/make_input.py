"""Write the SPECTRE input of the two-volume cylinder-limit Beltrami torus.

Design: core volume r < r1 with mu1 and B_s(0) = c1; outer volume with mu2
whose field at r1+ has pitch PITCH_OUT and |B| = RATIO |B(r1-)|. The pressure
step balances [[p + B^2/2]] = 0 at r1. Usage: make_input.py TEMPLATE OUT [R0]
"""

import sys

import numpy as np
import tomlkit
from scipy.special import j0, j1, y0, y1

from lundquist import Lundquist, flux_rows

A, R1 = 1.0, 0.5
C1, MU1, MU2 = 5.0, 1.5, 0.5
PITCH_OUT, RATIO = 0.1, 1.02


def design(r0):
    length = 2*np.pi*r0
    core = Lundquist([R1, A], [MU1, MU2], [C1, 0.0], [0.0, 0.0])
    bt, bs = core.field_in(0, R1)
    b_in = np.hypot(bt, bs)
    b_out = RATIO*b_in
    q = MU2*R1
    cj, cy = np.linalg.solve([[j0(q), y0(q)], [j1(q), y1(q)]],
                             [b_out*np.cos(PITCH_OUT), b_out*np.sin(PITCH_OUT)])
    tor0, _ = flux_rows(MU1, 0.0, R1, length)
    tor1, pol1 = flux_rows(MU2, R1, A, length)
    psi_t = np.array([C1*tor0[0], tor1 @ [cj, cy]])
    psi_p = pol1 @ [cj, cy]
    phiedge = float(f"{psi_t.sum():.10g}")
    return dict(phiedge=phiedge,
                tflux=[float(f"{psi_t[0]/phiedge:.10g}"), 1.0],
                pflux=[0.0, float(f"{psi_p/phiedge:.10g}")],
                pressure=[float(f"{(b_out**2 - b_in**2)/2:.10g}"), 0.0])


def main():
    template, out = sys.argv[1], sys.argv[2]
    r0 = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0
    d = design(r0)
    with open(template) as f:
        doc = tomlkit.load(f)
    ph = doc["physics"]
    ph["phiedge"] = d["phiedge"]
    ph["tflux"] = d["tflux"]
    ph["pflux"] = d["pflux"]
    ph["pressure"] = d["pressure"]
    ph["mu"] = [MU1, MU2]
    ph["rbc"] = {"(0, 0)": r0, "(1, 0)": A}
    ph["zbs"] = {"(1, 0)": -A}
    guess = tomlkit.table()
    guess["interface_1"] = {"rbc": {"(0, 0)": r0, "(1, 0)": R1},
                            "zbs": {"(1, 0)": -R1}}
    ph["allrzrz"] = guess
    doc["numeric"]["linitialize"] = 0
    with open(out, "w") as f:
        tomlkit.dump(doc, f)
    print(d)


if __name__ == "__main__":
    main()
