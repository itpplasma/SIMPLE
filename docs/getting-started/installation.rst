Installation
============

This page builds ``build/simple.x``, the library ``libsimple`` and optionally
the Python package ``pysimple`` from source. Run all commands from the
repository root.

Prerequisites
-------------

CMake ≥ 3.22, Ninja, GNU Fortran, Git, NetCDF-C and NetCDF-Fortran, LAPACK
and BLAS (OpenBLAS for ``setup-venv.sh``), and Python ≥ 3.10 with ``venv`` for
``pysimple`` and the tools in ``tools/``.

.. code-block:: bash

   # Arch Linux
   sudo pacman -S gcc-fortran cmake ninja netcdf netcdf-fortran lapack openblas

   # Ubuntu or Debian
   sudo apt install gfortran cmake ninja-build libnetcdf-dev libnetcdff-dev \
       liblapack-dev libopenblas-dev python3-venv

   # macOS with Homebrew
   brew install gcc cmake ninja netcdf netcdf-fortran lapack openblas

On macOS, select the Homebrew GNU compilers before the first configuration;
Apple Clang has no OpenMP:

.. code-block:: bash

   export FC=gfortran
   export CC="gcc-$(gfortran -dumpversion | cut -d. -f1)"
   export CXX="g++-$(gfortran -dumpversion | cut -d. -f1)"
   export CPATH="$(brew --prefix)/include${CPATH:+:$CPATH}"
   export LIBRARY_PATH="$(brew --prefix)/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"

The first configuration downloads, with CMake ``FetchContent`` and therefore
with network access to GitHub, `libneo <https://github.com/itpplasma/libneo>`_
(fields, splines, coordinates, chartmaps), `fortnum
<https://github.com/lazy-fortran/fortnum>`_ (root finders, ODE integrators),
fortio (file I/O, through libneo) and `fortplot
<https://github.com/lazy-fortran/fortplot>`_. Pinned revisions are in
``cmake/dependency_versions.json``.

Build
-----

.. code-block:: bash

   git clone https://github.com/itpplasma/SIMPLE.git
   cd SIMPLE
   make                 # Release build in build/
   make CONFIG=Fast     # native CPU flags; use for production runs

``make`` runs ``cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release`` and
``cmake --build build``. Further options:

=================================  =============  ================================================
Option                             Default        Effect
=================================  =============  ================================================
``CMAKE_BUILD_TYPE``               ``Release``    ``Release``, ``Fast`` (native CPU flags),
                                                  ``Debug`` or ``Profile``
``ENABLE_OPENMP``                  ``ON``         OpenMP over particles
``ENABLE_PYTHON_INTERFACE``        ``ON``         ``pysimple`` bindings (needs numpy, f90wrap)
``SIMPLE_ENABLE_CGAL``             ``OFF``        STL walls with CGAL (downloads CGAL, needs
                                                  Boost headers); :doc:`/user-guide/wall-heat-loads`
``SIMPLE_ENABLE_CUDA``             ``OFF``        native CUDA tracing; :doc:`/user-guide/gpu`
``SIMPLE_DETERMINISTIC_FP``        ``OFF``        no fast-math, bitwise reproducible
                                                  (``make build-deterministic``)
``LIBNEO_REF``                     pinned         libneo branch, tag or commit
                                                  (``make LIBNEO_REF=...``)
``LIBNEO_SOURCE_DIR``              unset          local libneo checkout (``make LIBNEO_PATH=...``)
=================================  =============  ================================================

Builds with different options live in separate build directories, for
example with CGAL:

.. code-block:: bash

   cmake -S . -B build-cgal -G Ninja -DCMAKE_BUILD_TYPE=Release \
       -DSIMPLE_ENABLE_CGAL=ON
   cmake --build build-cgal --target simple.x

Python package
--------------

``setup-venv.sh`` creates ``.venv`` in the repository, installs the Python
requirements, builds SIMPLE and installs ``pysimple`` in editable mode:

.. code-block:: bash

   ./setup-venv.sh
   source .venv/bin/activate
   python -c "import pysimple"

In an existing environment, run ``python -m pip install -r requirements.txt``
and ``python -m pip install --no-build-isolation -e .``; pip builds the
Fortran code with scikit-build-core in ``_skbuild/`` (several minutes). ``python -m pip install '.[plotting]'``
adds Matplotlib and netCDF4 for :doc:`/user-guide/postprocessing`.

Check
-----

.. code-block:: bash

   make test-fast

builds with ``SIMPLE_DETERMINISTIC_FP=ON`` and runs all tests except those
labelled ``slow``, ``regression``, ``performance`` and ``scalability``. Next:
:doc:`quickstart`.
