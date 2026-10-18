#!/usr/bin/env python3
"""Cross-check the RK/PRK tableaux in the Fortran source against the Lean ones.

Run from the SIMPLE root. The exact coefficients are read from the Lean
definitions (GcSymplectic.lean), including the irrational Gauss entries, and
enclosed in rational intervals with rigorously rounded square roots. The
Fortran literals are read from coeff_rk_gauss/coeff_rk_lobatto and converted
to the binary64 values gfortran assigns (correctly rounded decimal parsing).
Reports, with exact rational arithmetic:
  * max |a_code - a_exact| (upper bound from the enclosure),
  * the exact symplecticity defect b_i ah_ij + b_j a_ji - b_i b_j of the
    binary64 coefficients actually used by the code.
Exit status is nonzero if a bound is exceeded. No floating-point probe enters
any bound; floats only appear in the printed summary.
"""
import json
import re
import sys
from fractions import Fraction
from math import isqrt

LEAN = "concept/gc_symplectic_proof/lean/GcSymplectic.lean"
FORTRAN = "src/orbit_symplectic_base.f90"
OUT = "concept/gc_symplectic_proof/tableaux.json"
PREC = 2 ** 200
COEFF_BOUND = Fraction(1, 10 ** 14)   # code literal vs exact value
DEFECT_BOUND = Fraction(1, 10 ** 15)  # binary64 Butcher defect


class I:
    """Closed rational interval."""

    def __init__(self, lo, hi=None):
        self.lo = Fraction(lo)
        self.hi = Fraction(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def of(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, o):
        o = I.of(o)
        return I(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-I.of(o))

    def __rsub__(self, o):
        return I.of(o) - self

    def __mul__(self, o):
        o = I.of(o)
        p = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return I(min(p), max(p))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = I.of(o)
        assert o.lo > 0 or o.hi < 0, "division by interval containing 0"
        return self * I(1 / o.hi, 1 / o.lo)

    def __rtruediv__(self, o):
        return I.of(o) / self

    def width(self):
        return self.hi - self.lo


def S(x):
    x = I.of(x)
    assert x.lo >= 0
    lo = Fraction(isqrt(int(x.lo * PREC * PREC)), PREC)
    hi_num = x.hi * PREC * PREC
    hi = Fraction(isqrt(int(hi_num) + 1) + 1, PREC)
    return I(lo, hi)


def lean_defs(path):
    text = open(path).read()
    blocks = re.findall(r"^def (\w+)[^\n]*?:=\s*(.*?)(?=^\S)", text,
                        flags=re.S | re.M)
    env = {"S": S}
    for name, expr in blocks:
        if not re.match(r"(euler|gauss|g4|lobatto)", name):
            continue
        e = re.sub(r"Real\.sqrt (\d+)", r"S(\1)", expr)
        e = e.replace("Real.sqrt (", "S(").replace("![", "[")
        e = re.sub(r"(?<![\w.])(\d+)(?![\w.])", r"I(\1)", e)
        env[name] = eval(e, {"I": I}, env)
    return env


def fortran_tableaux(path):
    text = open(path).read()
    out = {}
    for sub in ("coeff_rk_gauss", "coeff_rk_lobatto"):
        body = re.search(r"subroutine %s\b(.*?)end subroutine %s" % (sub, sub),
                         text, flags=re.S).group(1)
        n = None
        for line in body.splitlines():
            m = re.search(r"if \(n == (\d+)\)", line)
            if m:
                n = int(m.group(1))
                out.setdefault((sub, n), {})
                continue
            m = re.match(r"\s*(a|ahat|b|c)\(([\d, ]+)\)\s*=\s*([-+0-9.eEdD]+)\s*$",
                         line)
            if m and n is not None:
                idx = tuple(int(k) - 1 for k in m.group(2).split(","))
                lit = m.group(3).lower().replace("d", "e")
                out[(sub, n)][(m.group(1),) + idx] = Fraction(float(lit))
    return out


def matrix(code, key, s):
    return [[code[(key, i, j)] for j in range(s)] for i in range(s)]


def main():
    env = lean_defs(LEAN)
    code = fortran_tableaux(FORTRAN)
    cases = [
        ("GAUSS1", ("coeff_rk_gauss", 1), "gauss1_a", "gauss1_a", "gauss1_b"),
        ("GAUSS2", ("coeff_rk_gauss", 2), "gauss2_a", "gauss2_a", "gauss2_b"),
        ("GAUSS3", ("coeff_rk_gauss", 3), "gauss3_a", "gauss3_a", "gauss3_b"),
        ("GAUSS4", ("coeff_rk_gauss", 4), "gauss4_a", "gauss4_a", "gauss4_b"),
        ("LOBATTO3", ("coeff_rk_lobatto", 3), "lobatto3_a", "lobatto3_ah",
         "lobatto3_b"),
    ]
    report = {"lean": LEAN, "fortran": FORTRAN, "variants": {}}
    ok = True
    for label, key, la, lah, lb in cases:
        c = code[key]
        s = len(env[lb])
        a = matrix(c, "a", s)
        ah = a if la == lah else matrix(c, "ahat", s)
        b = [c[("b", i)] for i in range(s)]
        ea, eah, eb = env[la], env[lah], env[lb]
        err = Fraction(0)
        width = Fraction(0)
        for i in range(s):
            pairs = [(b[i], I.of(eb[i]))]
            pairs += [(a[i][j], I.of(ea[i][j])) for j in range(s)]
            pairs += [(ah[i][j], I.of(eah[i][j])) for j in range(s)]
            for val, enc in pairs:
                err = max(err, abs(val - enc.lo), abs(val - enc.hi))
                width = max(width, enc.width())
        defect = max(abs(b[i] * ah[i][j] + b[j] * a[j][i] - b[i] * b[j])
                     for i in range(s) for j in range(s))
        passed = err <= COEFF_BOUND and defect <= DEFECT_BOUND
        ok = ok and passed
        report["variants"][label] = {
            "stages": s,
            "max_abs_code_minus_exact_upper_bound": float(err),
            "max_enclosure_width": float(width),
            "max_abs_binary64_symplectic_defect": float(defect),
            "binary64_defect_is_exactly_zero": defect == 0,
            "pass": passed,
        }
        print(f"{label:9s} s={s} |code-exact|<={float(err):.3e} "
              f"(enclosure {float(width):.1e})  "
              f"binary64 defect={float(defect):.3e}  {'PASS' if passed else 'FAIL'}")
    report["bounds"] = {"coefficient": float(COEFF_BOUND), "defect": float(DEFECT_BOUND)}
    if len(sys.argv) > 1 and sys.argv[1] == "--write":
        with open(OUT, "w") as f:
            json.dump(report, f, indent=1)
            f.write("\n")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
