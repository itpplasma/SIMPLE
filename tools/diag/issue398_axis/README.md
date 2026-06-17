# Near-axis symplectic loss evidence (issue #398)

Harness for the figures in itpplasma/SIMPLE#398: the symplectic step's implicit
Newton Jacobian consumes the Boozer field's second radial derivatives, which
diverge at the axis in `s = rho^2` coordinates; RK never uses them.

## Regenerate

Two SIMPLE diagnostics produce the inputs (build them: `diag_field_deriv.x`,
`diag_traj.x`):

1. Field and radial derivatives at fixed Boozer angle, internal transform:
   ```
   diag_field_deriv.x <internal_config>.in field_internal.dat
   ```
2. One axis-grazing orbit (e.g. particle 12), three integrators on the internal
   field:
   ```
   for im in 0 1 3; do
     sed -i "s/^integmode = .*/integmode = $im/" simple.in
     diag_traj.x simple.in 12 5      # writes traj_p12_im$im.dat
   done
   ```
3. Plot:
   ```
   python plot_axis_evidence.py field_internal.dat \
       boozer_chartmap_w7x.nc chartmap_simple_axispow.nc <orbit_dir>
   ```

`boozer_chartmap_w7x.nc` (booz_xform) and `chartmap_simple_axispow.nc` (SIMPLE
internal transform) are the W7-X high-mirror chartmaps from the benchmark case.

## Figures

- `fig1_d2_divergence.png` - `d2Bmod/ds2` and `d2hth/ds2` vs s. Both fields
  diverge ~ `s^(-3/2)`; in the orbit grazing band (`s < 1e-4`) the internal
  field reaches ~4e11 while booz_xform saturates near ~2e9. The internal curve
  is the actual integrator input (`diag_field_deriv`); the booz_xform curve is a
  cubic-spline estimate from the chartmap values, since the current loader
  requires the new `A_phi` format and cannot read that file directly.
- `fig2_orbit_energy.png` - radius `s(t)` and energy `|dH/H|(t)` for one
  grazing orbit. RK conserves energy and is bounded; symplectic Euler1 and
  midpoint random-walk the energy in steps at each near-axis pass. The RK trace
  is short because this build's RK aborts at the deepest axis pass on the
  axispow field; the production old-clamp RK confines this orbit for the full
  slowing-down time (run 101).
- `fig3_values_match.png` - `iota(s)` and `Bmax/Bmin(s)` agree between the two
  fields to <0.2%. It is not the field values; it is their second derivatives.
