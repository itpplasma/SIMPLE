This file contains a (possibly non-exhaustive) list of important changes.

Changes marked with "Y" in the Differing Output column may yield different output when run on the same input
compared to earlier versions.

| Differing output | commit hash / tag and timestamp | Summary |
| ---------------- | ------------------------------- | ------- |
| Y | 8e19626 2022-06-29 | Removed bmod_ref that rescaled to match B_00 on flux surface |
| N | b5fc202 2022-06-29 | MPI working again, misses later columns in times_lost.dat |
| N | c796d90 2022-06-28 | Merge collisions branch back to master, cross-checked results |
| Y | 7ca8fad 2022-06-23 | Default to old VMEC near-axis healing again |
| N | 86c3a5d 2022-06-23 | simple.in is now in namelist format |
| Y | e264ef7 2022-05-31 | New healing procedure for VMEC near axis |
| Y | 468e317 2022-02-14 | Compute actual momentum in z(4) every timestep |
| Y | a247433 2022-02-13 | Disable extrapolation entirely in init_sympl |
| N | 6981791 2021-12-20 | start.dat is now always stored in VMEC coordinates. |
| Y | v1.1 2020-02-13 | Lambda is now evaluated on half mesh, improves B-field interpolation. |
