"""Bind the joint FortSym scalar DAG to compact embedding arrays; no algebra."""

from pathlib import Path
import re
import sys
import textwrap

root = Path(sys.argv[1])
binding = {}
for r in range(4):
    for t in range(4 - r):
        for z in range(4 - r - t):
            for c in (1, 2):
                binding[f"j{r}{t}{z}{c}"] = f"jets({r},{t},{z},{c})"
for prefix, count, target in [("a", 3, "aj"), ("g", 6, "gj")]:
    for i in range(1, count + 1):
        for j in range(1, 11):
            binding[f"{prefix}{i}_{j}"] = f"{target}({i},{j})"
for prefix, count in [("p", 3), ("gh", 6), ("f", 6), ("t", 10)]:
    for i in range(1, count + 1):
        binding[f"{prefix}{i}"] = f"{'tj' if prefix == 't' else prefix}({i})"
for prefix in ("hh", "df"):
    for i in range(1, 7):
        for j in range(1, 7):
            binding[f"{prefix}{i}_{j}"] = f"{prefix}({i},{j})"
source = """! Generated array bindings for gen/gen_fo_joint.f90. Do not edit.
module fo_joint_kernel
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none
contains
"""
for name in ("fo_joint_first", "fo_joint_full", "fo_joint_flow"):
    scalar = (root / f"{name}.tmp").read_text()
    body = scalar[
        scalar.index("    real(dp) :: work") : scalar.rindex("end subroutine")
    ]
    body = re.sub(r"\b\w+\b", lambda m: binding.get(m[0], m[0]), body)
    if name == "fo_joint_flow":
        head = f"""subroutine {name}(energy,e0,tj,gh,hh,kenergy,f,df)
real(dp),intent(in)::energy,e0,tj(10),gh(6),hh(6,6)
real(dp),intent(out)::kenergy,f(6),df(6,6)
"""
    else:
        head = f"""subroutine {name}(p,aj,jets,gj,energy,gh,hh)
real(dp),intent(in)::p(3),aj(3,10),jets(0:3,0:3,0:3,2)
real(dp),intent(out)::gj(6,10),energy,gh(6),hh(6,6)
"""
        if name.endswith("first"):
            split = body.index("\n    work1 =")
            body = body[:split] + "\ngj=0; hh=0\n" + body[split:]
    source += head + body + f"end subroutine {name}\n"
source += "end module fo_joint_kernel\n"
logical, pending = [], ""
for line in source.splitlines():
    line = line.strip()
    pending += (" " if pending else "") + line.removesuffix("&").lstrip("&").strip()
    if not line.endswith("&"):
        logical.append(pending)
        pending = ""
assert not pending
lines = []
for line in logical:
    if not line.startswith("!"):
        line = re.sub(r"(\*\*|[+*/=])", r" \1 ", line)
        line = line.replace(",", ", ")
    pieces = textwrap.wrap(
        line, width=68, break_long_words=False, break_on_hyphens=False
    ) or [""]
    lines.extend(
        piece + (" &" if i < len(pieces) - 1 else "") for i, piece in enumerate(pieces)
    )
(root / "fo_joint_kernel.f90").write_text("\n".join(lines) + "\n")

for name in ("fo_joint_first", "fo_joint_full", "fo_joint_flow"):
    (root / f"{name}.tmp").unlink()
