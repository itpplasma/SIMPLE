# Marker-level controls for the 48-versus-53 discrepancy

`diagnosis.json` records the original five marker IDs, all control summaries,
case-specific configurations and source hashes, phase offsets and paired status
statistics. `traces.npz` retains every start/outcome and plotting traces sampled
every 16 records plus radial extrema, endpoints and both sides of switches.
Subset marker IDs 1–5 map to originals 62, 124, 126, 129, 247. The 16 launch-phase
markers all repeat original marker 129. A one-marker hybrid control uses ID 1.

`verification.json` preserves the wall-exclusion RED witnesses, the initially
overbroad single-volume failure, its unchanged-test rerun and the isolated PR
checks. No failed outcome is silently replaced by the final passing case.

The full scientific interpretation, plotting commands and expiring slopbox
links are in `../../RND.md`. Reproduction arguments are in each case's `run`
metadata; extend 2 ms/64 macrosteps to 2.125 ms/68 to preserve the base step.
