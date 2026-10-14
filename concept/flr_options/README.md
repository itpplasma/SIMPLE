# FLR and sharp-sheet research note recovered from the GPT worker

`flr_options.tex` and `flr/` contain the original note and its algebra,
local-map and citation-metadata scripts, reconstructed from the worker's
2026-10-14 tool log after the scratch directory was deleted. Only the scripts'
result paths were changed to work from the SIMPLE root.

```sh
python3 concept/flr_options/flr/check_algebra.py
python3 concept/flr_options/flr/check_local_maps.py
latexmk -pdf -outdir=/tmp/simple-flr-note concept/flr_options/flr_options.tex
```

The original delivered PDF was 23 pages at
<https://box.sloppy.at/5bcb4.pdf> (ephemeral hosting). The editable source is
retained here so that delivery does not depend on the URL. Generated PDFs and
rendering products are not committed.

The note distinguishes exact planar constant-field sheet maps and cycle
compression from unproved curved-sheet and long-time extensions. Phase-blind
Bessel/ring averaging alone does not reproduce the phase-dependent sheet reset.
KIM/KiLCA source observations and citation verification are historical research,
not newly repeated audits. The Neishtadt 1986 source limitation is retained.
`verify_citations.py` is an optional network tool; it was not rerun during this
recovery. See `PLAN.md` for how these ideas relate to W3–W6.
