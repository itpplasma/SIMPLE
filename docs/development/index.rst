Development
===========

SIMPLE is developed on GitHub at https://github.com/itpplasma/SIMPLE. The
contribution guidelines are in ``CONTRIBUTING.md`` in the repository root
(https://github.com/itpplasma/SIMPLE/blob/main/CONTRIBUTING.md). Pull requests
follow the review policy below and fill ``.github/pull_request_template.md``.

.. toctree::
   :maxdepth: 1

   /review-policy
   /performance/non_soa_cpu_regression
   changelog

Repository layout
-----------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Path
     - Content
   * - ``src/``
     - Fortran sources of the library and of ``simple.x``; ``src/field/``
       field classes and canonical coordinates, ``src/wall/`` STL wall,
       ``*.cu`` CUDA kernels.
   * - ``app/``
     - Main programs (``simple.x``, diagnostics).
   * - ``python/pysimple/``
     - Python interface and post-processing (``plotting``, ``engineering``).
   * - ``tools/``
     - Converters (GVEC, booz_xform, QUASR) and benchmark tools.
   * - ``examples/``
     - Example inputs and scripts; ``simple_full.in`` lists all keys.
   * - ``test/``
     - Unit, integration, Python and golden-record tests.
   * - ``DOC/``
     - Technical design notes (:doc:`/reference/design-notes`).
   * - ``docs/``
     - This documentation (Sphinx).
   * - ``RUN/``
     - Cluster job scripts.

Build and test
--------------

Run all commands from the repository root:

.. code-block:: bash

   make                      # build
   make test-fast            # fast tests, deterministic build
   make test                 # all tests except regression tests
   make test TEST=name       # one test by name
   make test-golden-main     # compare with the main branch

``make test`` and ``make test-fast`` rebuild with
``SIMPLE_DETERMINISTIC_FP=ON``. The golden-record tests build a reference
version and compare the outputs of a set of cases.

Code that changes coordinates or fields must update
``DOC/coordinates-and-fields.md``.

Build this documentation
------------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   sphinx-build -W -b html docs docs/_build/html

``-W`` turns warnings into errors. The Python API pages are generated from the
docstrings of ``pysimple``; the compiled extension is not needed, because it
is mocked. Figures that are produced separately are referenced from
``docs/_static/figures/``; until a figure file exists, the build shows a note
in its place (``docs/conf.py``).

The workflow ``.github/workflows/docs.yml`` builds the documentation for every
push and pull request and publishes it on GitHub Pages: ``latest/`` for
``main``, ``branches/<name>/`` for branches and ``previews/pr-<n>/`` for pull
requests. The same two commands build the site on any other CI system, for
example a GitLab Pages job that copies ``docs/_build/html`` to ``public/``.
