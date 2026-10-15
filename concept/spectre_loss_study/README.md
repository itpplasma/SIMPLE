# Matched SPECTRE loss and throughput study

This research driver uses the production GC/FO microstep pipeline at 3.5 keV.
It does not introduce another orbit integrator. Physical edge losses (88),
numerical stops (89), recovered failures and warning holds remain distinct.
`fo_elapsed` measures accepted FO time; the same two-line probe and driver are
applied to the detached main baseline. GC output positions represent guiding
centres, FO positions represent particles; compare like modes. FO `lambda` is
the most recent reconstructed GC pitch, not instantaneous particle pitch.
The FO reference explicitly enters every marker before stepping and disables
return to GC with a large exit distance. Enlarging the interface-entry distance
alone is insufficient: the production proximity screen excludes points whose
nearest integer is the magnetic axis. Before the outer-loss repair it also
excluded the outer edge, giving the hybrid a GC-centre absorber there.

From a configured native checkout:

```sh
python concept/spectre_loss_study/run.py --field ncsx --sampling sheet \
  --mode fo --output build/study/spectre/ncsx_sheet_fo
```

Use `--root` for the baseline checkout and `--starts .../study_starts.csv` for
every subsequent matched case. Baseline configs omit `--grid`; the spline
candidate accepts, for example, `--grid 32 64 32`. `--field g3v3` uses libneo's
documented G3V3L3Fi fixture. Each case stores the input, source revision,
equilibrium hash, native log, starts, marker outcomes, sparse traces and timing.
`process_s` includes native build/check overhead and sampling; `setup_s` includes
field initialization only; `trace_s` measures the parallel marker loop including
sparse diagnostics. Do not confuse the short pilot with a measured 2 ms run.

Default smoke: `make test CONFIG=Fast TEST='^study_spectre_orbits$' VERBOSE=1`.
The runner invokes that same native authority, never the executable directly.
Study output belongs under ignored `build/study`; retain compact summaries and
plot sources here when a milestone is complete. Parameter rejection tests are
outside this study's scope.

The controller records all solver/provider source hashes; completed cases are
reused only for matching settings and sources. Unresolved FO failures or held
intervals fail the native study after its diagnostic files have been saved.
They are never accepted just because the physical-loss count matches.

`pilot.py --baseline PATH --variants raw spline fine --markers 32
--duration 2e-4 --threads 16 --output build/study/spectre-long` runs the longer
paired controls. `full.py` then measures 256-marker, 2 ms volume controls on both
fields; `full_report.py INPUT OUTPUT` retains data and plots, and also replots a
retained results directory. These use the current entry phase and boundary
conventions, not a completed phase-consistent hybrid model.

`boozer_cost.py` compares GC, hybrid and FO cost with Boozer GC on the
source VMEC equilibrium; see [results/boozer-cost](results/boozer-cost/README.md).

## Targeted loss controls

The outer-loss repair resolves particle-wall events with the existing FO
entry/exit distances. `results/full-gc-fixed/` retains the earlier GC-edge
snapshot; `results/full-wall-fixed/` retains the particle-edge replay.
The diagnosis and current plot links are in [RND.md](RND.md). The initial five
G3V3L3Fi disagreements are markers 62, 124, 126, 129 and 247. Subset starts
must be renumbered consecutively for the driver; retain their original IDs
in the report.

Research controls accept `--samples` for denser diagnostics, `--npoiper2` for
GC step refinement, and `--construction-tolerance` for the existing Meiss ODE
tolerance. `--duration .002125 --macrosteps 68` extends the original 2 ms /
64-macrostep case while keeping its microstep exactly unchanged. `--phase`
rotates the initial perpendicular velocity and Larmor offset together around
the launch field, preserving the prescribed GC start and speed. FO mode with
`--sampling phase` applies offsets `phase + 2*pi*(marker-1)/markers`; repeat
the same GC start across markers for a launch-phase sweep. These are research
driver controls, not new production namelist settings or a phase-preserving
hybrid. Do not use them to tune the hybrid until a loss count matches.
