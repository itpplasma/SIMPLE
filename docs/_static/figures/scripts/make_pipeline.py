#!/usr/bin/env python3
"""Generate docs/_static/figures/pipeline.svg, the SIMPLE workflow schematic.

The SVG is written directly (no plotting library) so that all labels stay
editable text. It has a transparent background; box fills are light tints of
the Okabe-Ito palette with dark text, and text outside boxes switches colour
with ``prefers-color-scheme`` so it reads on light and dark pages.

Every box names a capability present on main; the comment next to each entry
points to its source of truth.

Usage::

    python docs/_static/figures/scripts/make_pipeline.py
"""

from __future__ import annotations

import argparse
from pathlib import Path
from xml.sax.saxutils import escape

import simple_style as st

W = 980
COL_W = 176
COL_X = [2, 202, 402, 602, 802]
TOP = 44

# Stage tints (light) and strokes (Okabe-Ito).
STAGES = {
    "source": ("#E1F0FA", st.OKABE_ITO["sky"]),
    "field": ("#D6E6F2", st.OKABE_ITO["blue"]),
    "canon": ("#D7EFE8", st.OKABE_ITO["green"]),
    "orbit": ("#FBEBD0", st.OKABE_ITO["orange"]),
    "out": ("#F8DFD0", st.OKABE_ITO["vermillion"]),
    "py": ("#F4E1EC", st.OKABE_ITO["purple"]),
}
INK = "#1F1F1F"

out: list[str] = []


def box(x, y, w, h, stage, lines, bold_first=True, size=12.5, rx=7, dash=None):
    fill, stroke = STAGES[stage]
    d = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="1.6"{d}/>'
    )
    lh = size * 1.28
    y0 = y + h / 2 - (len(lines) - 1) * lh / 2 + size * 0.36
    for i, text in enumerate(lines):
        weight = ' font-weight="700"' if (bold_first and i == 0) else ""
        out.append(
            f'<text x="{x + w / 2}" y="{y0 + i * lh:.1f}" text-anchor="middle" '
            f'font-size="{size if i == 0 else size - 1}" fill="{INK}"{weight}>'
            f"{escape(text)}</text>"
        )


def arrow(x1, y1, x2, y2, dashed=False, cls="arrow"):
    d = ' stroke-dasharray="5 4"' if dashed else ""
    out.append(
        f'<line class="{cls}" x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
        f'stroke-width="1.6"{d} marker-end="url(#head)"/>'
    )


def path(d, dashed=False):
    da = ' stroke-dasharray="5 4"' if dashed else ""
    out.append(
        f'<path class="arrow" d="{d}" fill="none" stroke-width="1.6"{da} '
        f'marker-end="url(#head)"/>'
    )


def label(x, y, text, size=12, anchor="middle", weight=None, cls="label"):
    w = f' font-weight="{weight}"' if weight else ""
    out.append(
        f'<text class="{cls}" x="{x}" y="{y}" text-anchor="{anchor}" '
        f'font-size="{size}"{w}>{escape(text)}</text>'
    )


def build() -> str:
    fs = 14.5  # box title size; body lines use fs - 1
    for x, text in zip(COL_X, ["1  Sources", "2  Field", "3  Canonical", "4  Orbits",
                               "5  Outputs"]):
        label(x + 3, 30, text, size=16, anchor="start", weight="700")

    # 1  Sources: field_input / coord_input (examples/simple_full.in).
    bh, gap = 56, 12
    ys = [TOP + i * (bh + gap) for i in range(5)]
    srcs = [
        ["VMEC", "wout NetCDF"],
        ["GVEC", "state file"],
        ["booz_xform", "boozmn NetCDF"],
        ["Coil set", "coils file"],
        ["SPECTRE", "multi-volume HDF5"],
    ]
    for y, lines in zip(ys, srcs):
        box(COL_X[0], y, COL_W, bh, "source", lines, size=fs)

    # 2  Field representation: src/field/*, tools/*_to_boozer_chartmap.py.
    x2 = COL_X[1]
    reps = {
        "vmec": (ys[0], ["Splined VMEC", "R, Z, λ on (s, θ, φ)"]),
        "chart": ((ys[1] + ys[2]) / 2, ["Boozer chartmap", "via tools/ converter"]),
        "coils": (ys[3], ["Biot–Savart", "on VMEC reference"]),
        "spectre": (ys[4], ["Per-volume fields", "interface sheets"]),
    }
    for y, lines in reps.values():
        box(x2, y, COL_W, bh, "field", lines, size=fs)
    xa, xb = COL_X[0] + COL_W, x2
    arrow(xa, ys[0] + bh / 2, xb - 2, ys[0] + bh / 2)
    cy = reps["chart"][0] + bh / 2
    path(f"M{xa},{ys[1] + bh / 2} C{xa + 12},{ys[1] + bh / 2} {xb - 14},{cy - 8} {xb - 2},{cy - 8}")
    path(f"M{xa},{ys[2] + bh / 2} C{xa + 12},{ys[2] + bh / 2} {xb - 14},{cy + 8} {xb - 2},{cy + 8}")
    arrow(xa, ys[3] + bh / 2, xb - 2, ys[3] + bh / 2)
    arrow(xa, ys[4] + bh / 2, xb - 2, ys[4] + bh / 2)

    # 3  Canonical coordinates: integ_coords 0, 2, 3, 4, 6.
    x3 = COL_X[2]
    ctop, cbot = ys[0], ys[3] + bh
    box(x3, ctop, COL_W, cbot - ctop, "canon",
        ["Boozer", "canonical flux", "Meiss", "Albert", "", "radial gauge",
         "A_r = h_r = 0"], size=fs + 0.5)
    box(x3, ys[4], COL_W, bh, "canon", ["Per-volume", "Meiss (SPECTRE)"], size=fs)
    for key in ("vmec", "chart", "coils"):
        y = reps[key][0] + bh / 2
        arrow(x2 + COL_W, y, x3 - 2, y)
    arrow(x2 + COL_W, ys[4] + bh / 2, x3 - 2, ys[4] + bh / 2)

    # 4  Orbit integration: integmode, orbit_model, crossing_level, gpu_*.
    x4 = COL_X[3]
    gc_h = 112
    box(x4, ys[0], COL_W, gc_h, "orbit",
        ["Guiding centre", "symplectic Euler,", "midpoint, Gauss,", "Lobatto; RK45;",
         "collisions optional"], size=fs)
    fo_y = ys[0] + gc_h + gap
    fo_h = 92
    box(x4, fo_y, COL_W, fo_h, "orbit",
        ["Full orbit", "Boris (Boozer chart);", "SPECTRE GC↔FO", "near interfaces"],
        size=fs)
    be_y = fo_y + fo_h + gap
    half = (COL_W - 8) / 2
    be_h = ys[4] + bh - be_y
    box(x4, be_y, half, be_h, "orbit", ["CPU", "OpenMP"], size=fs)
    # CUDA supports guiding-centre RK in Boozer coordinates, not full orbit.
    box(x4 + half + 8, be_y, half, be_h, "orbit",
        ["GPU", "CUDA", "RK", "Boozer"], size=fs)
    arrow(x3 + COL_W, (ctop + cbot) / 2, x4 - 2, (ctop + cbot) / 2)
    arrow(x3 + COL_W, ys[4] + bh / 2 - 12, x4 - 2, fo_y + fo_h / 2 + 18)

    # 5  Outputs: README main outputs, results.nc, classification, wall_input.
    x5 = COL_X[4]
    outs = [
        (["Confined fraction", "vs time"], bh),
        (["Loss times,", "exit codes, results.nc"], bh),
        (["Orbit classes", "trapped / passing,", "regular / chaotic"], 74),
        (["Wall hits", "STL wall (chartmap)"], bh),
    ]
    y = ys[0]
    for lines, h in outs:
        box(x5, y, COL_W, h, "out", lines, size=fs)
        y += h + gap
    out_bot = y - gap
    mid = (ys[0] + fo_y + fo_h) / 2
    arrow(x4 + COL_W, mid, x5 - 2, mid)

    # 6  pysimple (python/pysimple/*) and the external optimisation loop.
    py_y = ys[4] + bh + 52
    fx, fw = COL_X[1], COL_X[4] + COL_W - COL_X[1]
    fill, stroke = STAGES["py"]
    out.append(f'<rect x="{fx}" y="{py_y}" width="{fw}" height="118" rx="9" '
               f'fill="none" stroke="{stroke}" stroke-width="1.6"/>')
    label(fx + 10, py_y + 22, "pysimple  (Python)", size=15, anchor="start",
          weight="700")
    iy, ih = py_y + 34, 72
    box(COL_X[1], iy, COL_W, ih, "py",
        ["Python API", "same core and", "simple.in"], size=fs)
    box(COL_X[2], iy, COL_W, ih, "py",
        ["plotting", "loss statistics,", "energy vs J_perp"], size=fs)
    box(COL_X[3], iy, COL_X[4] + COL_W - COL_X[3], ih, "py",
        ["engineering", "wall heat-flux maps from loss positions,",
         "port placement"], size=fs)
    box(COL_X[0], iy, COL_W, ih, "py",
        ["External", "optimiser", "loss objective"], size=fs, dash="6 4")
    # Outputs feed pysimple; pysimple feeds an external optimiser.
    arrow(x5 + COL_W / 2, out_bot, x5 + COL_W / 2, py_y - 2)
    arrow(fx, iy + ih / 2, COL_X[0] + COL_W + 2, iy + ih / 2)
    # Loop back to the sources.
    ax_ = COL_X[0] + COL_W / 2
    path(f"M{ax_},{iy} L{ax_},{ys[4] + bh + 2}", dashed=True)
    label(ax_ + 8, ys[4] + bh + 26, "new equilibrium", size=12.5, anchor="start",
          cls="note")
    label(ax_ + 8, ys[4] + bh + 42, "or coils", size=12.5, anchor="start",
          cls="note")
    ly = py_y + 118

    css = (
        "text{font-family:'DejaVu Sans',Helvetica,Arial,sans-serif}"
        ".label{fill:#333333}.note{fill:#555555;font-style:italic}"
        ".arrow{stroke:#7F7F7F}#head path{fill:#7F7F7F}"
        "@media (prefers-color-scheme: dark){.label{fill:#E6E6E6}"
        ".note{fill:#BDBDBD}.arrow{stroke:#A6A6A6}#head path{fill:#A6A6A6}}"
    )
    total_h = ly + 4
    head = (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {total_h}" '
        f'width="{W}" height="{total_h}" role="img" '
        f'aria-labelledby="title desc">\n'
        "<title id=\"title\">SIMPLE workflow</title>\n"
        "<desc id=\"desc\">Field sources (VMEC, GVEC and booz_xform via Boozer "
        "chartmap converters, coils, SPECTRE) feed a splined field representation, "
        "canonical coordinates, guiding-centre integration on CPU or Boozer RK "
        "on GPU, full-orbit integration on CPU, and outputs "
        "(confined fraction, loss times, orbit classes, wall hits) "
        "that pysimple post-processes and exposes to external optimisers.</desc>\n"
        f"<style>{css}</style>\n"
        '<defs><marker id="head" viewBox="0 0 10 10" refX="8" refY="5" '
        'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
        '<path d="M0,0 L10,5 L0,10 z"/></marker></defs>\n'
    )
    return head + "\n".join(out) + "\n</svg>\n"


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=here.parent / "pipeline.svg")
    args = parser.parse_args()
    args.out.write_text(build())
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
