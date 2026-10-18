# Documentation figures

The README and Sphinx documentation use the text-based workflow schematic.
The scripts also produce scientific plots locally; generated PNG and PDF
outputs are not committed.

| Figure | Script | Input |
|---|---|---|
| `pipeline.svg` | `scripts/make_pipeline.py` | none (schematic) |
| local loss plot | `scripts/make_front.py` | real SIMPLE run, see below |

Run all commands from the repository root.

## Shared style

`scripts/simple_style.py` holds the colours and matplotlib settings shared by
all figure scripts. Categorical colours use the Okabe-Ito palette and are
paired with a second cue (line style, marker or direct label); continuous data
use `viridis` or `cividis`. New plots (verification, benchmarks, scaling)
should import it:

```python
import simple_style
simple_style.use()
...
simple_style.save(fig, "build/figures/my_figure.png")
```

`simple_style.save` writes PNGs on an opaque white background so dark text
stays legible on dark page themes, and keeps SVG text as text.

## `pipeline.svg`

```bash
python docs/_static/figures/scripts/make_pipeline.py
```

The SVG has a transparent background. Box text is dark on light tints; text
outside boxes switches colour with `prefers-color-scheme`. Each box names a
capability on `main`; the comments in the script give the source of truth
(input keys in `examples/simple_full.in`, modules in `src/` and
`python/pysimple/`, converters in `tools/`). Update the script when a
capability is added or removed.

## Local loss plot

```bash
make CONFIG=Fast
python docs/_static/figures/scripts/make_front.py               # run and plot
python docs/_static/figures/scripts/make_front.py --plot-only   # replot only
```

The script downloads the public VMEC equilibrium
`wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc` (precise
quasi-axisymmetric configuration of Landreman & Paul, Phys. Rev. Lett. 128,
035001 (2022), scaled to reactor size) from the simsopt test data at a pinned
commit, checks its SHA-256, writes `simple.in` and runs `build/simple.x` in
`build/figures/front/`. The run traces 3000 collisionless 3.5 MeV alpha
guiding-centre markers from `s = 0.5` for 50 ms in Boozer coordinates with
symplectic Euler; it takes a few minutes on a 32-core workstation.

The left panel shows where lost markers cross the last closed flux surface
(`zend` in `results.nc`, evaluated on the VMEC boundary and cross-checked
against `xend_cart`), marked by their class at the start (`trap_par`): discs
for trapped, diamonds for passing. These are guiding-centre crossing points of
the LCFS, not wall impact points. The right panel shows the passing, trapped
and total confined fractions of all markers. The script computes them at full
time resolution from `times_lost` and checks them against SIMPLE's
`confined_fraction.dat`. Strongly passing markers with trapping parameter at
or below `contr_pp = -1` (the default) count as confined without tracing.

The example plot was generated with GNU 16.1 on macOS, the Fortran
sources of SIMPLE `fe2de27`, and libneo
`f6090e8150c2f5b1c0afde09baadf79b50b1ffe7` in the `Fast` profile. This run
had 194 LCFS losses (157 initially trapped, 37 initially passing), 1260
completed markers, 1546 skipped markers, and no numerical-failure exits.
The original generation reported 193 losses. `deterministic = .true.` fixes
the random seed; it does not guarantee identical trajectories across
toolchains and dependency revisions.
