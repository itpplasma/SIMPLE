# SPECTRE orbit interface status

Multi-volume SPECTRE orbit tracing (`Mvol > 1`) is disabled on main until the
full-orbit interface treatment of PR #530 lands. Configuration fails with:

> no valid interface crossing on this version; see DOC/spectre-rejected-interface-models.md

The conserved-magnetic-moment guiding-centre crossing maps were removed because
a sharp sheet changes the physical magnetic moment by a gyrophase-dependent
amount and the maps produce spurious confinement. The
[rejection record](spectre-rejected-interface-models.md) retains the mathematical
arguments, symbolic verifier, recorded measurements, and code provenance.

Single-volume SPECTRE tracing, multi-volume field evaluation, canonical field
construction, geometry, field-line diagnostics, and start-point generation
remain supported. Sampling a multi-volume field requires
`generate_start_only = .true.`; it does not trace orbits.
