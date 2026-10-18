# Hybrid replay after the GC extrapolation repair

The two hybrids were rerun through native Make/CTest after integrating the
main-based fix `6970eff` as `0993478`. Their source hashes include the GC
integrator and canonical field provider. The four FO controls are retained
from `../full/`: explicit FO throughout does not call the repaired GC stepper,
and their field, starts and FO equations are unchanged. Case revisions and
original run paths remain separate in `summary.json`; this is not a claim
that every case was rerun at the same revision.

NCSX hybrid loses 98/256 (previously 99), versus 99 refined-FO losses, with
seven paired status disagreements. G3V3L3Fi remains 48/256 versus 53, with
five paired disagreements. Both replays have zero stops, FO failures or
warning-held intervals. Timings are shared-host accuracy replays; use the
original controlled paired runs for speedup claims. The existing production
entry phase, switch and GC-centre versus particle edge conventions persist.
