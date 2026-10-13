# NCSX free-boundary test fixture

This HDF5 file contains two plasma volumes and one external vacuum volume
(`Nvol=2`, `Mvol=3`) with `Nfp=3`. It was generated from SPECTRE's committed
`tests/test_files/wout_li383_low_res_reference.nc` at SPECTRE commit
`1b26166ec7b24d6b445a49490d36f4b242a84849` using simsopt 1.10.6.

The fixed-boundary VMEC conversion used normalized toroidal-flux surfaces 0.1
and 0.7, `mpol=4`, `ntor=2`, and `interface_guess="wout"`. The computational
boundary was extended outward by 0.15 m along the surface normal. SPECTRE's
`fix2free` converter used `nvns_mpol=4`, `nvns_ntor=4`, and `vacuum=True`.
The resulting input, `ncsx_freeboundary.toml`, preserves the NCSX `Nfp`,
toroidal-flux sign, and coil-linked current sign. Commands, with
`make_ncsx_freeboundary.py` from `plasma/proj/spectre-orbits` and its first
flux surface set to 0.1:

    make_ncsx_freeboundary.py wout_li383_low_res_reference.nc ncsx_freeboundary.toml
    run_case.sh ncsx_freeboundary.toml field 3

The three-rank SPECTRE field solve completed without overlapping interfaces.
This compact file is a field-evaluation and particle-transition fixture, not a
force-minimized free-boundary equilibrium. At this resolution the Beltrami
residual is large (largest component average 0.56, 2.8 and 12 in the three
volumes); `mpol=8`, `ntor=4` lowers it to 0.19 and 0.20 in the plasma volumes
at five times the file size and test time.

The interface blend of every volume keeps a single Jacobian sign. An earlier
version placed the first surface at 0.2, where the NCSX bean indentation is
already strong: the axis-volume blend `sbar**m` then folded near theta = pi for
`sbar > 0.95`, and libneo rejects such a chart.

SHA-256:
`225976b7c5e2f5cd9480b2e2d233fda2a68a590028f36f24cd63f17aa29b24d9`.
