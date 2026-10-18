Inputs and field sources
========================

A run is configured by the ``&config`` namelist of ``simple.in`` and, for
collisions, the ``&profiles`` namelist. ``simple.x`` reads ``simple.in`` from
the working directory or the file given as the first argument. All keys have
defaults; ``examples/simple_full.in`` lists them
(:doc:`/reference/input-parameters`).

Three keys select the files:

``field_input``
   The magnetic field: a VMEC ``wout`` file, a Boozer chartmap, a coil file or
   a SPECTRE HDF5 file. The type is detected from the file.
``coord_input``
   The reference coordinates in which start points and results are given:
   a VMEC ``wout`` file or a chartmap. Defaults to ``netcdffile``, then to
   ``field_input``.
``netcdffile``
   A VMEC ``wout`` file. It is the fallback for the two keys above, and the
   default ``wout.nc`` is used when a VMEC file is required.

``integ_coords`` (or its older alias ``isw_field_type``) selects the
coordinates in which the orbits are integrated:

.. list-table::
   :header-rows: 1
   :widths: 8 16 28 48

   * - Value
     - Name
     - Coordinates
     - Use
   * - 0
     - CANFLUX
     - :math:`(s, \vartheta_c, \varphi_c)`
     - Canonical flux coordinates from VMEC (symplectic).
   * - 1
     - VMEC
     - :math:`(s, \vartheta, \varphi)`
     - VMEC coordinates, Runge-Kutta only.
   * - 2
     - BOOZER
     - :math:`(s, \vartheta_B, \varphi_B)`
     - Boozer coordinates from VMEC or a Boozer chartmap (symplectic).
       The standard choice for nested flux surfaces.
   * - 3
     - MEISS
     - :math:`(r, \vartheta_c, \varphi_c)`
     - Canonical coordinates without flux surfaces, for example for coil
       fields (symplectic).
   * - 4
     - ALBERT
     - :math:`(\psi, \vartheta_c, \varphi_c)`
     - Meiss coordinates with a radial coordinate proportional to
       :math:`A_\vartheta` (symplectic).
   * - 5
     - REFCOORDS
     - reference coordinates
     - ``field_input`` splined on ``coord_input``, Runge-Kutta only.
   * - 6
     - SPECTRE
     - :math:`(\rho_g, \vartheta, \zeta)`
     - Per-volume coordinates of SPECTRE/SPEC equilibria.
   * - -1
     - TEST
     - :math:`(r, \vartheta, \varphi)`
     - Analytic tokamak field for tests.

The coordinates and the canonical transformations are described in
:doc:`/model/coordinates`.

VMEC equilibria
---------------

A VMEC ``wout*.nc`` file in NetCDF format is the standard input:

.. code-block:: fortran

   &config
   netcdffile = 'wout.nc'
   isw_field_type = 2      ! Boozer
   /

SIMPLE splines :math:`R`, :math:`Z` and :math:`\lambda` of VMEC over
:math:`s` and the angles (orders ``ns_s``, ``ns_tp``; angular resolution
``multharm``), transforms to Boozer or canonical flux coordinates and splines
the canonical field components. Near the magnetic axis the VMEC harmonics are
continued smoothly (``axis_healing``; ``'polyfit'`` is recommended and will
become the default).

``vmec_B_scale`` and ``vmec_RZ_scale`` multiply the field and the size of the
device. They allow a quick scan of a configuration at another field or size
without a new equilibrium; a scaled equilibrium is still an MHD equilibrium
only if the pressure is scaled consistently.

Start positions, ``sbeg`` and ``start.dat`` refer to VMEC coordinates
:math:`(s, \vartheta, \varphi)` with :math:`s` the normalised toroidal flux.

Boozer chartmaps (GVEC, booz_xform)
-----------------------------------

A Boozer chartmap is a NetCDF file with the geometry and the Boozer field
components on a uniform grid. SIMPLE runs from it without a VMEC file:

.. literalinclude:: ../../examples/simple_chartmap.in
   :language: fortran
   :caption: examples/simple_chartmap.in

``field_input`` and ``coord_input`` both point to the chartmap. ``sbeg`` is
still the normalised toroidal flux :math:`s`; the chartmap radial coordinate is
:math:`\rho = \sqrt{s}`. The schema, sign conventions and scaling rules are in
:doc:`/boozer-chartmap-schema`.

Converters in ``tools/`` write chartmaps from other codes:

.. code-block:: bash

   # GVEC (needs the gvec Python package, installed by setup-venv.sh)
   python tools/gvec_to_boozer_chartmap.py parameter_final.ini State_final.dat \
       boozer_chartmap.nc

   # booz_xform output
   python tools/booz_xform_to_boozer_chartmap.py boozmn_case.nc boozer_chartmap.nc

Both accept ``--nrho``, ``--ntheta`` and a toroidal resolution option
(``--nphi`` or ``--nzeta``). The GVEC converter computes Boozer coordinates
with GVEC's own transformation and handles the orientation change to SIMPLE's
left-handed convention (option ``--flip``). SIMPLE does not link GVEC at run
time. The tested reference workflow is ``test/tests/run_gvec_qa_roundtrip.py``.

Coil fields
-----------

A coil file gives the field of filamentary coils by the Biot-Savart law. It is
needed when particles are followed outside the last closed flux surface, for
example to a wall. The SIMPLE coil format is a text file: the number of points
in the first line, then one line per point with :math:`x, y, z` in metres and
the current in ampere. Point :math:`i` and point :math:`i+1` are joined by a
straight segment that carries the current of point :math:`i`, so a current of
zero ends a coil. The converter of the test suite translates the
STELLOPT/MAKEGRID format and keeps these zero-current points:

.. code-block:: bash

   python test/tests/convert_coils_stellopt_to_simple.py coils.c09r00 coils.simple

At this revision, ``python -m libneo.convert_coils_to_simple`` drops the
zero-current points, so that consecutive coils are joined by spurious current
segments; do not use it for SIMPLE input.

Coil fields need canonical coordinates that do not assume flux surfaces
(``isw_field_type = 3`` or ``4``) or the Runge-Kutta path on reference
coordinates (``isw_field_type = 5``, ``integmode = 0``). A VMEC file
(``netcdffile``) is still used for the start positions and as reference
coordinates unless ``coord_input`` is a chartmap:

.. code-block:: fortran

   &config
   netcdffile = 'wout.nc'
   field_input = 'coils.simple'
   isw_field_type = 3      ! Meiss
   /

The coil set and the VMEC equilibrium must describe the same field; SIMPLE
does not check this; see the warning in :doc:`/tutorials/wall-heat-loads`.

Chartmaps beyond the LCFS
-------------------------

For wall losses, the reference coordinates must continue beyond the last
closed flux surface. libneo writes such chartmaps from a VMEC file:

.. code-block:: bash

   export PYTHONPATH=build/_deps/libneo-src/python
   python -m libneo.chartmap_cli from-vmec-extended wout.nc chartmap_wall.nc \
       --nrho 33 --ntheta 65 --nzeta 33 --boundary-offset 0.15

The global attribute ``rho_lcfs`` records the radial position of the LCFS.
``chart_boundary_kind`` defines what the outer boundary of a chartmap means: a
loss at the LCFS, a loss at a wall, or the end of the numerical domain
(:doc:`/reference/input-parameters`). Other subcommands of
``libneo.chartmap_cli`` build chartmaps from an STL boundary
(``from-stl``) or extend a VMEC chartmap to an STL wall
(``from-vmec-to-wall``); they need the optional Python packages
``map2disc`` and ``shapely``.

SPECTRE equilibria
------------------

SPECTRE/SPEC multi-region relaxed MHD equilibria (HDF5, ``.h5``) are detected
from ``field_input``. They require ``integ_coords = 6``. SIMPLE integrates
guiding-centre orbits per volume and switches to full orbits near the
interfaces between volumes, where the field jumps:

.. literalinclude:: ../../examples/simple_spectre.in
   :language: fortran
   :caption: examples/simple_spectre.in

``sbeg`` is the SPECTRE radial label :math:`\rho_g \in [0, M_\text{vol}]`.
The interface treatment and its controls (``crossing_level``,
``spectre_fo_entry_larmor``, ``spectre_fo_tolerance``, ``spectre_fo_grid``)
are described in :doc:`/reference/design-notes`.

Start positions
---------------

``startmode`` selects how start positions are generated:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Value
     - Start positions
   * - 1
     - Random points along one field line on the surface ``sbeg``, weighted
       by volume (default). ``grid_density`` places them on a grid instead.
   * - 2
     - Read ``start.dat``: :math:`(s, \vartheta, \varphi, v/v_0,
       v_\parallel/v)` per line in reference coordinates.
   * - 3
     - Read an ANTS-format start file.
   * - 4
     - A random batch of the lines of ``start.dat``.
   * - 5
     - Uniform random :math:`s` in a range set by ``num_surf`` and
       ``sbeg``, random angles.
   * - 6
     - Read ``start.dat`` with positions in integration coordinates.

Modes 1 and 5 write the generated positions to ``start.dat``. The pitch
:math:`v_\parallel/v` is uniform in :math:`[-1, 1]`, the speed is the birth
speed (:math:`v/v_0 = 1`). The particle species is set by ``n_e`` (charge
number), ``n_d`` (mass number) and ``facE_al`` (energy :math:`3.5\,\text{MeV}
/ \text{facE\_al}`). All angles are in radians and are not wrapped to
:math:`[0, 2\pi)`.
