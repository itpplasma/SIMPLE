# Documentation figures

The documentation build renders every figure from text sources in the
repository. PNG and PDF outputs are gitignored and are published with the
site only.

```bash
python docs/_static/figures/scripts/build_doc_figures.py
sphinx-build -W -b html docs docs/_build/html
```

Run all commands from the repository root.

| Figure | Script | Source |
|---|---|---|
| `pipeline.svg` | `scripts/make_pipeline.py` | none (schematic) |
| `front.png` | `scripts/make_front.py` via `build_doc_figures.py` | `data/front/` and the pinned public VMEC file |
| `strong-scaling.png`, `marker-throughput.png` in `docs/performance/efss-hpc/` | `tools/perf/plot_efss_scaling.py`, `tools/perf/plot_efss_throughput.py` | CSV files in `docs/performance/efss-hpc/` |

## Shared style

`scripts/simple_style.py` holds the colours and matplotlib settings for all
figure scripts. Categorical colours use the Okabe–Ito palette paired with a
second cue (line style, marker or label); continuous data use `viridis` or
`cividis`. `simple_style.save` writes PNGs on an opaque white background so
they stay legible in dark themes, and keeps SVG text as text.

```python
import simple_style
simple_style.use()
...
simple_style.save(fig, "build/figures/my_figure.png")
```

## `pipeline.svg`

```bash
python docs/_static/figures/scripts/make_pipeline.py
```

Each box names a capability on `main`; comments in the script give the source
of truth (keys in `examples/simple_full.in`, modules in `src/` and
`python/pysimple/`, converters in `tools/`). Update the script when a
capability is added or removed. Text outside the boxes follows
`prefers-color-scheme`.

## Front figure (`front.png`)

The figure shows where lost markers cross the last closed flux surface and the
confined fractions over time for the public quasi-axisymmetric equilibrium
`wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc`
(Landreman & Paul, Phys. Rev. Lett. 128, 035001 (2022), scaled to reactor
size). The run traces 3000 collisionless 3.5 MeV alpha guiding-centre markers
from $s = 0.5$ for 50 ms in Boozer coordinates with symplectic Euler.

`data/front/` stores the result as text: positions, times and initial classes
(`trap_par`) of the 194 lost markers in `losses.csv`, 101 confinement samples
in `confinement.csv`, and population counts and provenance in
`metadata.json` (SIMPLE `ffd7f4a`, libneo `f6090e8`, GNU Fortran 16.1 on
macOS, `Fast` profile). The other markers enter only through their counts.
`build_doc_figures.py` downloads the VMEC file from the simsopt test data at a
pinned commit, checks its SHA-256, checks the event-based confined fractions
against the samples and replots without running SIMPLE.

The left panel marks LCFS crossings (`zend`, checked against `xend_cart`) of
initially trapped markers with discs and of passing markers with diamonds;
these are guiding-centre crossings, not wall impacts. The right panel shows
passing, trapped and total confined fractions from `times_lost`. Strongly
passing markers with trapping parameter at or below `contr_pp = -1` count as
confined without tracing.

To regenerate the data with a new run (a few minutes on 32 cores):

```bash
make CONFIG=Fast
python docs/_static/figures/scripts/make_front.py               # run and plot
python docs/_static/figures/scripts/make_front.py --plot-only   # replot only
python docs/_static/figures/scripts/export_front_data.py \
  build/figures/front docs/_static/figures/data/front \
  --source-revision SOURCE_COMMIT --libneo-revision LIBNEO_COMMIT
```

`deterministic = .true.` fixes the random seed; trajectories still differ
between toolchains and dependency revisions.
