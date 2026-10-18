# Multi-volume particle-wall replay

The two hybrid cases replay the same starts after the GC extrapolation and
outer-boundary repairs. Four unchanged FO references are reused from `../full/`.
Per-case run paths, source hashes and revisions remain separate in `summary.json`.
The measured multi-volume equations are unchanged by the subsequent correction
that preserves the separate single-volume GC boundary convention.

NCSX hybrid loses 97/256; G3V3L3Fi hybrid loses 52/256. All health counters are
zero. A separate G3V3L3Fi refined-hybrid replay also loses 52 with identical
final statuses; its outcomes/provenance are retained in `../loss-diagnosis/`.
See `../../RND.md` for the marker-level diagnosis, uncertainty and verification.
These shared-host accuracy timings overlap other local checks and are not a
new controlled throughput comparison. `full_report.py` preserves the explicit
particle-boundary scope when replotting this retained directory.
