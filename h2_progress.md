# PR #530 rebase regression repair

- Starting branch: feat/spectre-fo-interfaces, 9f5f395; working tree clean.
- Read CLAUDE.md and relevant Fortran guidance; using make only, no delegation.
- review_h.md is absent; supplied findings checked against production source.
- Confirmed simple.x and simple_api share init_spectre_field.
- Existing direct volume test overrides sympl_rmax and masks the init regression.
- Added production-init volume-2 GC witness without overriding sympl_rmax.
- Fail before (9f5f395 sources): bound 1; first step enters FO; angular advance 0.
- Pass after: bound 3; 10000 steps; FO entries/stops/warning holds all zero;
  angular advance 47.010 rad. Outward edge marker loses at rho_g=2,
  status 88, terminal fraction 3.1826e-4; one landing, no FO or stops.
- Restored Mvol+1 trial margin; reset default bound on field initialization.
- Updated bound comments and coordinates documentation. Interface doc has no
  bound claim and needs no edit.
- Next: dead stop-record chain cleanup and full make test.
