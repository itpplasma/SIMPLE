# Lundquist cylinder limit of the SPECTRE interface crossing

Validation rung for SPECTRE support: a two-volume stepped-pressure Beltrami
torus at `R0/a = 100` (`test/tests/data/cyl2vol.h5`, see its
`PROVENANCE.md`) against the analytic two-volume Lundquist field of the
straight periodic cylinder with the same fluxes and `mu`. The sheet at
`r = 0.5 m` turns the field by 0.28 rad and raises `|B|` by 2 %.

- `lundquist.py`: analytic field from per-volume fluxes and SPECTRE's flux
  conventions.
- `make_input.py`: designs the fixture input (`R0` as optional argument).
- `test_spectre_cylinder_field` (ctest, `test/tests/field_can`): SIMPLE's
  full-orbit spline field and raw SPECTRE field against Lundquist, bound
  `2 a/R0`; writes `cylinder_field.dat`.
- `cylinder_markers.f90` (`cylinder_markers.x`, built with SIMPLE, not a
  test): production GC/FO microstep loop (`integmode = 3`,
  `crossing_level = 3`) for 3.5 MeV alphas whose guiding centres start inside
  the FO entry band, `|d| < 1.4` radial gyroradii from the sheet, uniform
  pitch `|lambda| < 0.9`. Microstep 1/32 gyro-period, 3 gyro-periods.
- `reference_fo.py`: independent full orbit in the analytic cylinder, DOP853
  (`rtol 1e-11`), the sheet as an exact terminal event switching the field.
  *Independent* orbits start from the same guiding-centre distance to the
  sheet, side and pitch with 4 random gyrophases and a plain Larmor lift;
  *paired* orbits start from SIMPLE's entry particle, moved so that its
  distance to the sheet is kept (the torus interface is Shafranov-shifted by
  6.6 mm).
- `analyze.py`: statistics and figure.

Observables per marker: whether the particle reaches the other volume, and
at the first microstep end there `mu_after/mu_before` (`mu = v_perp^2/(2B)`
with the local one-sided field, `mu_before` from the guiding-centre state)
and the guiding centre's distance from the sheet. In the cylinder `p_theta`
and `p_z` are conserved, so a particle that touches the sheet stays on it;
"one passage" is the first far-side arc, and the hybrid never returns such
markers to GC.

After `make CONFIG=Fast`, in an empty work directory (`$S` = SIMPLE root):

```sh
$S/build/test/tests/field_can/cylinder_markers.x $S/test/tests/data/cyl2vol.h5 800 markers.dat
python $S/concept/spectre_cylinder/reference_fo.py markers.dat $S/test/tests/data/cyl2vol.h5 ref.npz
```

The `R0 = 50` comparison uses the same `make_input.py ... 50` input, run
through SPECTRE outside the repository.

## Findings

Field (16 poloidal angles, both volumes): `|B|` within 0.96 % and pitch
within 5.2 mrad at `R0 = 100`; 1.97 % and 10.4 mrad at `R0 = 50`. The
deviation is the `cos theta` toroidal modulation and scales with `a/R0`.

Crossing (800 markers per `R0`, gyroradius 2.5-5.8 % of `a`):

| quantity | SIMPLE, R0 = 100 | reference | SIMPLE, R0 = 50 |
|---|---|---|---|
| crossing fraction | 0.743 +- 0.015 | 0.742 +- 0.008 | 0.744 +- 0.015 |
| mean mu_after/mu_before | 1.016 +- 0.008 | 1.012 +- 0.005 | 1.022 +- 0.009 |
| mean GC distance from sheet after passage | 0.5 +- 1.2 mm | 0.5 +- 0.6 mm | 0.2 +- 1.2 mm |

Paired orbits (same entry particle): the crossing flag agrees for 99.4 %
of markers; on common crossings the median `|d mu|/mu` is 1.9e-3 at
`R0 = 100` and 3.5e-3 at `R0 = 50`, the median GC offset 0.6 % and 1.0 % of
`rho_L`, i.e. linear in `a/R0`.

The `mu_after/mu_before` distributions differ in shape (KS p = 0.001, SIMPLE
narrower, means equal) while SIMPLE and the paired reference agree (KS
p = 1.0). The difference therefore comes from the entry, not the orbit:
SIMPLE's entry matches `p_zeta` through `v_par`, which moves `v_par` by
rms 0.09 `v` for guiding centres within 2 cm of the sheet (0.01 beyond
4 cm), against the plain Larmor lift of the reference. This is the documented
`O(rho*)` entry correction on a current sheet (SIMPLE#528).
