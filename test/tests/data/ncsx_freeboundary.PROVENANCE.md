# NCSX free-boundary test fixture

This HDF5 file contains two plasma volumes and one external vacuum volume
(`Nvol=2`, `Mvol=3`) with `Nfp=3`. It was generated from SPECTRE's committed
`tests/test_files/wout_li383_low_res_reference.nc` at SPECTRE commit
`1b26166ec7b24d6b445a49490d36f4b242a84849` using simsopt 1.10.6.

The fixed-boundary VMEC conversion used normalized toroidal-flux surfaces 0.2
and 0.7, `mpol=4`, `ntor=2`, and `interface_guess="wout"`. The computational
boundary was extended outward by 0.15 m along the surface normal. SPECTRE's
`fix2free` converter used `nvns_mpol=4`, `nvns_ntor=4`, and `vacuum=True`.
The resulting input preserves the NCSX `Nfp`, toroidal-flux sign, and
coil-linked current sign.

The three-rank SPECTRE field solve completed without overlapping interfaces.
This compact file is a field-evaluation and particle-transition fixture, not a
force-minimized free-boundary equilibrium.

SHA-256:
`03e64793ca915199fad7823d8260448270d4c538786ff98d47fc7580bb25312a`.
