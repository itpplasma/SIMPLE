# SPECTRE orbit cost relative to Boozer GC, 2026-10-15

This compares the measured cost of the SPECTRE models with SIMPLE's Boozer GC
on the VMEC equilibrium that generated the NCSX SPECTRE fixture,
`wout_li383_low_res_reference.nc` (simsopt test file, SHA-256
`975b13acbcef2a9ef2c34ec565fc9799ec989e3acf37686e5208d048f960f92f`).
The fixture's plasma boundary is that equilibrium's s = 0.7 surface, extended
by 0.15 m of vacuum; both fields have B of about 1.6 T and Nfp = 3.

All cases use 3.5 keV alphas (`facE_al=1000`), volume-sampled starts, 2 ms,
16 OpenMP threads and the same deterministic Fast build at `70b5e6a` (libneo
`c92b336`, GNU Fortran 16.1.0) on one Apple M3 Ultra. The three SPECTRE cases
share the hybrid's 256 starts. Cost is normalized as core-seconds per simulated
particle-second, which removes differing loss times and thread imbalance.
SPECTRE uses per-marker wall times. For Boozer it is trace time times 16
threads, an upper bound that favours SPECTRE.

| Model | Integrator, step | Core-s / particle-s | vs Boozer default | vs Boozer same step |
| --- | --- | ---: | ---: | ---: |
| Boozer GC, production default | Euler1, `npoiper2=256` | 19.7 | 1 | 0.14 |
| Boozer GC | midpoint, 256 | 37.2 | 1.9 | 0.26 |
| Boozer GC, SPECTRE study settings | midpoint, 1024 | 142.8 | 7.3 | 1 |
| SPECTRE GC (77/151 markers never in FO) | midpoint, 1024 | 385–400 | ~20 | 2.7–2.8 |
| SPECTRE hybrid, entry 0.01 gyroradii (2.5% FO time) | | 700 | 36 | 4.9 |
| SPECTRE hybrid, default entry 3 (30% FO time) | | 3,852 | 196 | 27 |
| SPECTRE FO everywhere | RK45, 1e-9 | 11,629 | 591 | 81 |

Boozer rows use 2,560 markers; the earlier 256-marker runs gave 14 and 138.
All SPECTRE health counters are zero. Wall-clock trace times for 256 SPECTRE
markers are 32.4 s (entry 0.01), 191.5 s (hybrid) and 349.8 s (FO).

## Findings

- A pure-GC SPECTRE step costs 8.0–8.4 us at a 20.9 ns step. At the same
  integrator and `npoiper2`, this is 2.7x the Boozer GC cost. The remaining factor
  to production Boozer comes from the study's 4x smaller step and midpoint
  instead of Euler1. Midpoint at 256 cuts the SPECTRE GC part 3.7x without a
  resolvable loss change; Euler1 makes the switch chatter (`../fo-residence/`).
- FO costs about 30x GC on the same field. The default hybrid spends 30% of
  particle time in FO, which accounts for about 90% of its compute; it is 3.0x
  cheaper than FO everywhere in core-seconds.
- Hybrid thread efficiency is 42–45% versus 70% for FO everywhere. This came
  from the study driver's static OpenMP schedule; production already uses
  dynamic scheduling, and the driver now matches it (wall time 191 s to 94 s,
  identical outcomes; see `../fo-residence/`). The earlier wall-clock
  decomposition that implied FO/GC ratios of only 3–7 was wrong for the same reason.
- Residual FO at sheets and the wall (entry 0.01) already costs 36x Boozer.
  Reducing FO residence is therefore the main cost lever; faster GC steps alone
  cannot approach Boozer cost.

## Limits

These are single runs on a shared host (load average 9–15 on 28 cores).
Boozer losses occur at the VMEC s = 1 boundary, whereas SPECTRE includes a
vacuum region: 20% versus 38–41% losses. The normalization compares cost only,
not physical equivalence. The entry-0.01 case is a cost probe, not a validated
model. The Boozer step equals the SPECTRE step only up to the two fields'
major-radius normalization.

## Reproduction

From a configured checkout (on macOS, GCC with OpenMP and Homebrew NetCDF via
`CC=gcc-16 CXX=g++-16 CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib`):

```sh
O=build/study/spectre-mac
python concept/spectre_loss_study/run.py --field ncsx --mode hybrid \
  --markers 256 --duration 2e-3 --threads 16 --output $O/hybrid
for c in "hybrid 0.01 gcprobe" "fo 3 fo"; do set -- $c
  python concept/spectre_loss_study/run.py --field ncsx --mode $1 --entry $2 \
    --markers 256 --duration 2e-3 --threads 16 --output $O/$3 \
    --starts $O/hybrid/study_starts.csv
done
mkdir -p build/study/boozer
curl -L -o build/study/boozer/wout_li383_low_res_reference.nc <wout_url in summary.json>
B="python concept/spectre_loss_study/boozer_cost.py run"
$B setup_only 1 256 256 --duration 1e-7
$B x10_default 1 256 2560; $B x10_mid256 3 256 2560; $B x10_matched 3 1024 2560
python concept/spectre_loss_study/boozer_cost.py report \
  build/study/spectre-mac build/study/boozer concept/spectre_loss_study/results/boozer-cost
```

`summary.json` retains every case's counts, times, FO fraction, thread
efficiency and ratios. Run `make test` first so that `build/simple.x` is the
deterministic build the SPECTRE study uses.
