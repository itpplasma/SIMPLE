Tutorial 2: Wall heat loads
===========================

This tutorial turns lost alphas into a heat-load estimate in two ways.
**Part A** bins the exit points of lost guiding centres on the last closed
flux surface (LCFS) of the QA case of :doc:`loss-fraction` into a heat flux in
MW/m² (seconds of work after Tutorial 1). **Part B** sets up the chain for
losses on a triangulated first wall from an STL file (a second build of a few
minutes, about one minute of run time). Setup reference:
:doc:`/user-guide/wall-heat-loads`; estimator formulas:
:doc:`/model/heat-loads`.

Part A: heat flux on the last closed flux surface
-------------------------------------------------

``tutorial1/qa/results.nc`` holds for each marker the exit code and the final
state ``zend`` :math:`= (s, \vartheta, \varphi, v/v_0, \lambda)`. For a lost
marker this is where its guiding centre crossed the LCFS.

The LCFS area comes from a chartmap, a NetCDF file with Cartesian positions on
a :math:`(\rho, \vartheta, \zeta)` grid (:doc:`/boozer-chartmap-schema`).
libneo's ``vmec_to_chartmap.x``, built with SIMPLE, writes one from the VMEC
file; the last three arguments are the grid sizes in :math:`\rho`,
:math:`\vartheta`, :math:`\zeta`:

.. code-block:: bash

   source .venv/bin/activate
   mkdir -p tutorial2/lcfs
   cd tutorial2/lcfs
   ../../build/_deps/libneo-build/vmec_to_chartmap.x \
       ../../tutorial1/qa/wout.nc qa.chartmap.nc 17 33 33

The script bins the energy of the lost markers in 8 poloidal × 16 toroidal
bins. Each marker carries :math:`P_\alpha/N` at birth and the fraction
:math:`(v/v_0)^2` of it at the exit (1 without collisions):

.. literalinclude:: scripts/lcfs_heat_flux.py
   :language: python
   :caption: docs/tutorials/scripts/lcfs_heat_flux.py

.. code-block:: bash

   python ../../docs/tutorials/scripts/lcfs_heat_flux.py \
       ../../tutorial1/qa/results.nc qa.chartmap.nc

::

   lost markers: 31 of 256 (12.1%)
   lost power: 73 MW
   LCFS area: 887 m^2
   peak flux: 0.93 MW/m^2, mean flux: 0.088 MW/m^2

With 31 lost markers the map ``lcfs_heat_flux.png`` is Monte Carlo noise: one
marker in a bin already carries 2.3 MW, and a bin with :math:`k` markers has
relative error :math:`1/\sqrt{k}`. A converged pattern needs :math:`10^4` to
:math:`10^5` markers and finer bins (:doc:`/user-guide/performance`).

* The points are guiding-centre positions on the LCFS, not wall impact
  points; the gyroradius of a 3.5 MeV alpha at 5 to 6 T is a few centimetres.
* ``zend`` angles are VMEC angles, and ``vmec_to_chartmap.x`` uses the same
  angles, so bins and areas refer to the same surface patches.

Limitation: ``pysimple.engineering.WallHeatMap.from_netcdf`` selects losses by
``class_lost``, which only classification runs set, and indexes ``zend`` with
transposed axes; the script above reads ``orbit_exit_code`` and ``zend``
directly.

Part B: losses on an STL first wall
-----------------------------------

To follow lost markers to a wall outside the LCFS, SIMPLE needs a CGAL
build, a chartmap that continues beyond the LCFS (``rho_lcfs`` < 1), a field
defined there (the coil field), and the wall as a closed STL surface. This
part uses the public NCSX coils ``c09r00`` with a fixed-boundary vacuum VMEC
equilibrium of the same configuration, and a wall 11 to 12 cm outside the
LCFS.

.. warning::

   This example demonstrates the input chain; its loss numbers are not a
   physics result. The public coils give 1.3 to 1.9 T on the VMEC surfaces,
   the public equilibrium is scaled to 0.5 T, and the coil field has a normal
   component of up to 5 % on these surfaces. For a physics result use an
   equilibrium computed from the coil set.

Build with CGAL (header-only download; needs the Boost headers, ``boost`` on
Arch, ``libboost-dev`` on Ubuntu):

.. code-block:: bash

   cmake -S . -B build-cgal -G Ninja -DCMAKE_BUILD_TYPE=Release \
       -DSIMPLE_ENABLE_CGAL=ON
   cmake --build build-cgal --target simple.x

Fetch the inputs and make the libneo Python tools of the build importable:

.. code-block:: bash

   source .venv/bin/activate
   mkdir -p tutorial2/stl
   cd tutorial2/stl
   wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc \
       -O wout_ncsx.nc
   wget https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00
   export PYTHONPATH=$PWD/../../build/_deps/libneo-src/python

Convert the coils from STELLOPT/MAKEGRID to the SIMPLE format (points in m,
currents in A). The test-suite converter keeps the zero-current point that
closes each coil; ``libneo.convert_coils_to_simple`` drops it and joins
consecutive coils (:doc:`/user-guide/field-sources`):

.. code-block:: bash

   python ../../test/tests/convert_coils_stellopt_to_simple.py coils.c09r00 coils.simple

Write a chartmap that uses the VMEC coordinates inside the LCFS and extends
them by 15 cm with a cubic Hermite continuation:

.. code-block:: bash

   python -m libneo.chartmap_cli from-vmec-extended wout_ncsx.nc chartmap_wall.nc \
       --nrho 33 --ntheta 65 --nzeta 33 --boundary-offset 0.15

The file records ``rho_lcfs = 0.719``. Take the chartmap surface at
:math:`\rho = 0.94`, between LCFS and outer boundary, as the wall:

.. literalinclude:: scripts/make_wall_stl.py
   :language: python
   :caption: docs/tutorials/scripts/make_wall_stl.py

.. code-block:: bash

   python ../../docs/tutorials/scripts/make_wall_stl.py chartmap_wall.nc 0.94 wall.stl

In an application ``wall.stl`` is the first-wall CAD export; ``wall_units``
gives its length unit, because STL carries none. Write ``simple.in``:

.. code-block:: fortran

   &config
   netcdffile = 'wout_ncsx.nc'       ! VMEC equilibrium: start positions
   field_input = 'coils.simple'      ! Biot-Savart field of the coils
   coord_input = 'chartmap_wall.nc'  ! reference coordinates beyond the LCFS
   isw_field_type = 3                ! Meiss canonical coordinates
   integmode = 1                     ! symplectic Euler
   wall_input = 'wall.stl'           ! triangulated first wall
   wall_units = 'm'
   contr_pp = -1d10                  ! trace all markers
   facE_al = 100d0                   ! 35 keV alphas
   sbeg = 0.5d0
   ntestpart = 64
   trace_time = 1d-3
   deterministic = .true.
   output_results_netcdf = .true.
   /

.. code-block:: bash

   ../../build-cgal/simple.x simple.in > simple.log

Most of the run time is the Biot–Savart field on the spline grid and the
canonical transformation. Meiss coordinates (``isw_field_type = 3``) do not
assume flux surfaces and remain valid between LCFS and wall
(:doc:`/model/coordinates`). ``facE_al = 100`` lowers the energy to 35 keV:
with a minor radius of about 0.3 m and a 3.5 MeV gyroradius of about 18 cm at
1.5 T, full-energy alphas would be lost on their first orbit.

``results.nc`` contains per marker ``wall_hit``, the intersection point
``wall_hit_cart``, the facet normal and ``wall_hit_angle_rad`` (0 is normal
incidence); a wall hit has exit code 2:

.. literalinclude:: scripts/wall_hits.py
   :language: python
   :caption: docs/tutorials/scripts/wall_hits.py

.. code-block:: bash

   python ../../docs/tutorials/scripts/wall_hits.py results.nc

::

   exit codes: {0: 17, 2: 47}
   wall hits: 47 of 64
   hit points: R = 1.04..1.76 m, Z = -0.65..0.69 m
   median incidence angle: 67 deg
   median loss time: 1.65e-05 s

Chaotic orbits are not bitwise reproducible across compilers and dependency
versions, so expect counts within the Monte Carlo error of these. The hits,
the energy fraction :math:`(v/v_0)^2` (``zend`` column 4) and the facet areas
give the load per triangle with ``pysimple.wall_heat_load.TriangleHeatLoad``
(:doc:`/user-guide/wall-heat-loads`).

The same case with the field splined on the chartmap and adaptive RK45
(``isw_field_type = 5``, ``integmode = 0``) gives 44 hits of 64. Agreement
within counting error checks both input paths but does not bound
discretisation bias; a quantitative study also checks convergence in step size
and field resolution.
