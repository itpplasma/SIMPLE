# Recovered FO and GC derivations

These are the stopped derivation workers A–E, collected with their original
sources and available results. `PLAN.md` distinguishes historical reports from
checks repeated in the consolidated checkout.

- A: fixed-parameter composition and selection obstruction. The checking
  declarations are in `../fo_symplectic_proof/lean/EntryFrozen.lean`; failed
  declarations remain explicitly separated in `EntryFrozenDraft.lean.txt`.
- B: provider-jet identities and a mutation-control draft. The baseline can
  run against the collected concept source and this checkout's libneo sources
  using `bash concept/fo_gc_derivations/B_provider_jets/run_provider_jets.sh`.
  The mutation harness is unfinished: its final mutation did not apply, and
  mutation runs overwrite the baseline JSON. Rerun the baseline afterward.
- C: written inexact-Newton derivative bounds; no validated implementation or
  machine-checked theorem.
- D: first-order particle/GC transform and Newton inverse, generated kernel,
  SPECTRE field-jet adapter and retained seven-case study. The standalone build
  script requires `build/neo/include` and `build/neo/lib` containing its libneo
  dependency snapshot. These ignored build products are not distributed.
  Completing the portable build and symbolic checks is remaining work.
- E: written second-order GC production prerequisites, without machine checks.

`../fo_symplectic_proof/lean/LieTruncation.lean` states generic truncation
results. Its hypotheses still have to be established for the emitted switch
and the production field; it does not itself certify that implementation.
