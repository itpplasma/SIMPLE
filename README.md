[![DOI](https://zenodo.org/badge/146476544.svg)](https://zenodo.org/badge/latestdoi/146476544)

# SIMPLE

**S**ymplectic **I**ntegration **M**ethods for **P**article **L**oss **E**stimation

SIMPLE traces energetic-particle guiding-center orbits in three-dimensional
magnetic fields to compute fusion alpha-particle losses in stellarators.
Symplectic integration and early orbit classification make it suitable for
configuration studies and stellarator optimization.

- VMEC equilibria, Boozer chartmaps from GVEC and booz_xform, coil fields,
  and SPECTRE equilibria.
- Loss statistics, orbit classification, collisions and STL wall intersections;
  available combinations are described in the [user guide](https://itpplasma.github.io/SIMPLE/latest/user-guide/index.html).
- OpenMP parallelism, an optional CUDA backend and the `pysimple` Python interface.

## Install and run

Install GNU Fortran, CMake, Ninja, NetCDF-C/Fortran and BLAS/LAPACK, then:

```bash
git clone https://github.com/itpplasma/SIMPLE.git
cd SIMPLE
make
mkdir -p run
cp examples/simple.in run/simple.in
curl -L https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -o run/wout.nc
(cd run && ../build/simple.x simple.in)
```

This example traces 128 alpha particles for 10 ms in a public
quasi-axisymmetric equilibrium. Start with the
[installation guide](https://itpplasma.github.io/SIMPLE/latest/getting-started/installation.html)
for platform packages and the Python interface, and the
[quick start](https://itpplasma.github.io/SIMPLE/latest/getting-started/quickstart.html)
for interpreting the output.

## Documentation

The [manual](https://itpplasma.github.io/SIMPLE/latest/) includes tutorials,
the physical model, input reference, verification and HPC instructions.
The source and release tags are also available on the
[EUROfusion GitLab mirror](https://gitlab.eufus.psnc.pl/itpplasma/simple).
For contributions, see the [development guide](https://itpplasma.github.io/SIMPLE/latest/development/index.html).

## Cite

Cite the methods used in your study:

- [Guiding-center symplectic integration (JCP, 2020)](https://doi.org/10.1016/j.jcp.2019.109065).
- [Accelerated alpha-loss computation (JPP, 2020)](https://doi.org/10.1017/S0022377820000203).
- [Orbit classification (JPP, 2023)](https://doi.org/10.1017/S0022377823000351).
- [Wall heat loads (Nuclear Materials and Energy, 2026)](https://doi.org/10.1016/j.nme.2026.102200).

The [publications page](https://itpplasma.github.io/SIMPLE/latest/reference/publications.html)
contains full references and related work; [CITATION.cff](CITATION.cff)
provides citation metadata.

## License

[MIT](LICENSE). The FIRM3D-compatible tracing path has separate
[provenance and license information](LICENSES/FIRM3D-MIT.txt).
