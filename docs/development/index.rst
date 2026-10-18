Development
===========

SIMPLE is developed at https://github.com/itpplasma/SIMPLE. Start with
:doc:`/contributing` for the pull-request rules, :doc:`architecture` for the
code structure and :doc:`testing` for the tests a change needs.

.. toctree::
   :maxdepth: 1

   /contributing
   coding-standard
   architecture
   testing
   /collision-golden-contract
   design-notes
   changelog

Repository layout
-----------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Path
     - Content
   * - ``src/``
     - Fortran library and ``simple.x`` driver; ``src/field/`` field classes
       and canonical coordinates, ``src/coordinates/`` reference
       coordinates, ``src/wall/`` STL wall, ``*.cu`` CUDA kernels.
   * - ``app/``
     - Main program ``simple.x``, diagnostics, ``export_boozer_chartmap``,
       ``spectre_poincare``.
   * - ``python/pysimple/``
     - Python interface, plotting and heat-load post-processing.
   * - ``tools/``
     - Converters (GVEC, booz_xform, QUASR), input validator, reference
       generators, benchmark tools.
   * - ``examples/``
     - Example inputs and scripts; ``simple_full.in`` lists every key.
   * - ``test/``
     - Unit, integration, Python and golden-record tests (:doc:`testing`).
   * - ``DOC/``
     - Implementation notes shown under :doc:`design-notes`.
   * - ``docs/``
     - This documentation (Sphinx).
   * - ``RUN/``, ``examples/hpc/``
     - Cluster job scripts (:doc:`/hpc/index`).
   * - ``cmake/dependency_versions.json``
     - Pinned fortnum, fortplot and golden-record libneo revisions.

Build and test
--------------

From the repository root:

.. code-block:: bash

   make CONFIG=Fast          # build in build/
   ./setup-venv.sh           # once: .venv with the Python test dependencies
   source .venv/bin/activate
   make test-fast            # before every push
   make test TEST=name       # tests whose name matches a regex

:doc:`testing` lists all targets and labels. Dependencies and CMake options
are in :doc:`/getting-started/installation`. A libneo checkout or revision is
selected on the first configuration with ``make LIBNEO_PATH=/path/to/libneo``
or ``make LIBNEO_REF=<branch|tag|sha>``; run ``make clean`` before changing
it.

Build this documentation
------------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   python docs/_static/figures/scripts/build_doc_figures.py
   sphinx-build -W -b html docs docs/_build/html

The figure script downloads one pinned equilibrium on its first run; sources
and provenance of the figures are in ``docs/_static/figures/README.md``. The
API pages are generated from the ``pysimple`` docstrings with the compiled
extension mocked, so no build of SIMPLE is needed.
``.github/workflows/docs.yml`` runs the same commands on every push and pull
request and publishes ``main`` to GitHub Pages at ``latest/``.
