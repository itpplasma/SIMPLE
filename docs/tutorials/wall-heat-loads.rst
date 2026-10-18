Tutorial 2: Wall heat loads
===========================

This tutorial turns lost alpha particles into a heat-load estimate. It has two
parts.

* **Part A** bins the exit points of lost guiding centres on the last closed
  flux surface (LCFS) of the QA configuration of :doc:`loss-fraction` and
  converts them into a heat flux in MW/m². This is the standard first
  estimate during configuration design.
* **Part B** sets up the full chain for losses on a triangulated first wall
  given as an STL file: a CGAL build, reference coordinates that extend beyond
  the LCFS, a field that is defined there, and the wall itself.

The physics of the heat-load estimate is described in
:doc:`/user-guide/wall-heat-loads` and :doc:`/model/wall-intersection`.

Runtime: Part A takes seconds after Tutorial 1. Part B needs a second build
(a few minutes) and about one minute of run time.

Part A: heat flux on the last closed flux surface
-------------------------------------------------

Tutorial 1 wrote ``tutorial1/qa/results.nc`` because of
``output_results_netcdf = .true.``. It contains, for each particle, the exit
code and the final phase-space position ``zend`` = :math:`(s, \vartheta,
\varphi, v/v_0, v_\parallel/v)`. For a lost particle this is the point where
its guiding centre crossed the LCFS.

The area of the LCFS comes from a chartmap file, a NetCDF file with the
Cartesian positions on a grid in :math:`(\rho, \vartheta, \zeta)`
(:doc:`/boozer-chartmap-schema`). libneo's ``vmec_to_chartmap.x``,
which is built together with SIMPLE, writes one from the VMEC file. The
arguments after the file names are the numbers of grid points in
:math:`\rho`, :math:`\vartheta` and :math:`\zeta`:

.. code-block:: bash

   source .venv/bin/activate
   mkdir -p tutorial2/lcfs
   cd tutorial2/lcfs
   ../../build/_deps/libneo-build/vmec_to_chartmap.x \
       ../../tutorial1/qa/wout.nc qa.chartmap.nc 17 33 33

The script ``docs/tutorials/scripts/lcfs_heat_flux.py`` bins the energy
carried by the lost markers in 8 poloidal by 16 toroidal bins. Each marker
represents the power :math:`P_\alpha / N`; it carries the fraction
:math:`(v/v_0)^2` of its birth energy, which is 1 without collisions:

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

The plot ``lcfs_heat_flux.png`` shows the flux over the toroidal and poloidal
angle. With 31 lost markers it is dominated by Monte Carlo noise: a bin with
one marker already carries 2.3 MW. The relative error of a bin with :math:`k`
markers is about :math:`1/\sqrt{k}`. A converged heat-load pattern needs
:math:`10^4` to :math:`10^5` markers and finer bins; use the GPU or a cluster
(:doc:`/user-guide/performance`).

Two properties of this estimate matter for its interpretation:

* The positions are guiding-centre positions on the LCFS, not impact points on
  a wall. For 3.5 MeV alphas at 5 to 6 T the gyroradius is a few centimetres.
* The angles in ``zend`` are VMEC angles, and the chartmap from
  ``vmec_to_chartmap.x`` uses the same angles, so the bins and the areas refer
  to the same surface patches.

.. note::

   ``pysimple.engineering.WallHeatMap.from_netcdf`` implements a similar
   binning, but at this revision it does not read the losses of a standard run
   correctly; see the open issues in :doc:`/user-guide/wall-heat-loads`. The
   script above reads ``orbit_exit_code`` and ``zend`` directly.

Part B: losses on an STL first wall
-----------------------------------

A first wall lies outside the LCFS. To follow lost particles to the wall,
SIMPLE needs

1. a build with CGAL, which provides the segment-triangle intersection,
2. reference coordinates that continue beyond the LCFS: a chartmap whose
   radial coordinate :math:`\rho` reaches the LCFS at ``rho_lcfs`` < 1,
3. a magnetic field that is defined between LCFS and wall, in practice the
   vacuum field of the coils,
4. the wall as a closed triangulated surface in an STL file.

This part uses the only public pair of coil set and VMEC equilibrium that is
also used by the SIMPLE test suite: the NCSX coils ``c09r00`` and a
fixed-boundary vacuum equilibrium of the same configuration. The wall is a
surface 11 to 12 cm outside the LCFS.

.. warning::

   This example demonstrates the input chain; its loss numbers are not a
   physics result. The public NCSX coil set produces 1.3 to 1.9 T on the
   VMEC flux surfaces, while the public VMEC equilibrium is scaled to 0.5 T,
   and the coil field has a normal component of up to 5 % of its magnitude on
   these surfaces (checked with an independent Biot-Savart evaluation). The
   VMEC surfaces are therefore not exact flux surfaces of the coil field. For
   a physical result use a coil set and an equilibrium computed from it, as in
   the published wall-load studies (:doc:`/verification/index`).

Build SIMPLE with CGAL
~~~~~~~~~~~~~~~~~~~~~~

CGAL is downloaded as a header-only library; the Boost headers must be
installed (``boost`` on Arch Linux, ``libboost-dev`` on Ubuntu):

.. code-block:: bash

   cmake -S . -B build-cgal -G Ninja -DCMAKE_BUILD_TYPE=Release \
       -DSIMPLE_ENABLE_CGAL=ON
   cmake --build build-cgal --target simple.x

Prepare the inputs
~~~~~~~~~~~~~~~~~~

From the repository root:

.. code-block:: bash

   source .venv/bin/activate
   mkdir -p tutorial2/stl
   cd tutorial2/stl
   wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc \
       -O wout_ncsx.nc
   wget https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00
   export PYTHONPATH=$PWD/../../build/_deps/libneo-src/python

Convert the coils from the STELLOPT/MAKEGRID format to the SIMPLE format
(points in metres and currents in ampere) with the converter of the SIMPLE
test suite:

.. code-block:: bash

   python ../../test/tests/convert_coils_stellopt_to_simple.py coils.c09r00 coils.simple

::

   Converted 18690 coil segments from coils.c09r00 to coils.simple

The converter keeps the zero-current point that closes each coil. Do not use
``python -m libneo.convert_coils_to_simple`` at this revision: it drops these
points, so that consecutive coils are joined by spurious current segments
(see :doc:`/user-guide/field-sources`).

The libneo Python tools used below come with the libneo source that the build
has downloaded; ``PYTHONPATH`` makes them importable.

Write a chartmap that uses the VMEC coordinates inside the LCFS and extends
them by 15 cm with a smooth (cubic Hermite) continuation:

.. code-block:: bash

   python -m libneo.chartmap_cli from-vmec-extended wout_ncsx.nc chartmap_wall.nc \
       --nrho 33 --ntheta 65 --nzeta 33 --boundary-offset 0.15

The file records ``rho_lcfs = 0.719``: the LCFS is at 72 % of the radial
coordinate. Write the wall as the chartmap surface at :math:`\rho = 0.94`,
which lies between the LCFS and the outer chartmap boundary:

.. literalinclude:: scripts/make_wall_stl.py
   :language: python
   :caption: docs/tutorials/scripts/make_wall_stl.py

.. code-block:: bash

   python ../../docs/tutorials/scripts/make_wall_stl.py chartmap_wall.nc 0.94 wall.stl

::

   rho = 0.9375: 12870 triangles written to wall.stl

In an application, ``wall.stl`` is the first-wall CAD geometry exported as STL.
``wall_units`` gives its length unit, because STL files carry none.

Run
~~~

Write ``simple.in``:

.. code-block:: fortran

   &config
   netcdffile = 'wout_ncsx.nc'       ! VMEC equilibrium: start positions
   field_input = 'coils.simple'      ! Biot-Savart field of the coils
   coord_input = 'chartmap_wall.nc'  ! reference coordinates beyond the LCFS
   isw_field_type = 3                ! Meiss canonical coordinates
   integmode = 1                     ! symplectic Euler
   wall_input = 'wall.stl'           ! triangulated first wall
   wall_units = 'm'
   contr_pp = -1d10                  ! trace all particles
   facE_al = 100d0                   ! 35 keV alphas
   sbeg = 0.5d0
   ntestpart = 64
   trace_time = 1d-3
   deterministic = .true.
   output_results_netcdf = .true.
   /

and run the CGAL build:

.. code-block:: bash

   ../../build-cgal/simple.x simple.in > simple.log

The run takes about 35 seconds, most of it for the Biot-Savart field on the
spline grid and the canonical transformation. ``isw_field_type = 3`` builds
canonical Meiss coordinates on the chartmap; they do not assume flux surfaces
and remain valid between LCFS and wall (:doc:`/model/coordinates`).
``facE_al = 100`` reduces the energy to 35 keV: NCSX has a minor radius of
about 0.3 m, and the gyroradius of a 3.5 MeV alpha at 1.5 T is about 18 cm,
so full-energy alphas are lost on their first orbit. ``contr_pp = -1d10``
traces all particles, including the strongly passing ones that the default
counts as confined without tracing.

Read the wall hits
~~~~~~~~~~~~~~~~~~

``results.nc`` contains for each particle the flag ``wall_hit``, the
intersection point ``wall_hit_cart``, the wall normal and the incidence angle
``wall_hit_angle_rad`` (0 is normal incidence). The exit code of a wall hit
is 2.

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

These are representative results. A fixed sampling seed does not make
chaotic orbit outcomes bitwise identical
across compilers and dependency versions.

The intersection points lie on the STL surface. Together with the energy
:math:`(v/v_0)^2` of each marker (column 4 of ``zend``) and the area of the
wall triangles they give the heat load on the wall geometry, as in Part A.
The incidence angle enters the load on tilted tiles.

As a check of the numerics, the same case with the field splined on the
chartmap and the adaptive Runge-Kutta integrator (``isw_field_type = 5``,
``integmode = 0``) gives 44 wall hits of 64 in the reference run. This checks
both input paths, but agreement within the Monte Carlo counting uncertainty
does not bound discretisation bias. A quantitative study must also check
step-size and field-resolution convergence.
