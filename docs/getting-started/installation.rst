Installation
============

SIMPLE is built from source with CMake. The build produces the executable
``build/simple.x``, the Fortran library ``libsimple`` and, when the Python
prerequisites are present, the Python package ``pysimple``. All commands on
this page are run from the root of the repository.

Prerequisites
-------------

* CMake 3.22 or newer and Ninja
* GNU Fortran (gfortran) and Git
* NetCDF-C and NetCDF-Fortran
* LAPACK and BLAS (OpenBLAS for ``setup-venv.sh``)
* OpenMP (part of gfortran)
* Python 3.10 or newer with ``venv``, for ``pysimple`` and the tools

On common systems:

.. code-block:: bash

   # Arch Linux
   sudo pacman -S gcc-fortran cmake ninja netcdf netcdf-fortran lapack openblas

   # Ubuntu or Debian
   sudo apt install gfortran cmake ninja-build libnetcdf-dev libnetcdff-dev \
       liblapack-dev libopenblas-dev python3-venv

   # macOS with Homebrew
   brew install gcc cmake ninja netcdf netcdf-fortran lapack openblas

On macOS, select the Homebrew GNU toolchain before the first configuration
if other compilers are installed. This also avoids requiring OpenMP support
from Apple Clang. Expose the Homebrew headers and libraries to the compiler:

.. code-block:: bash

   export FC=gfortran
   export CC="gcc-$(gfortran -dumpversion | cut -d. -f1)"
   export CXX="g++-$(gfortran -dumpversion | cut -d. -f1)"
   export CPATH="$(brew --prefix)/include${CPATH:+:$CPATH}"
   export LIBRARY_PATH="$(brew --prefix)/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"

The build downloads further dependencies with CMake ``FetchContent``:
`libneo <https://github.com/itpplasma/libneo>`_ (field, spline and
coordinate routines), `fortnum <https://github.com/lazy-fortran/fortnum>`_
(root finders and ODE integrators), fortio (scientific file I/O, through
libneo) and `fortplot <https://github.com/lazy-fortran/fortplot>`_ (plots).
The fortnum and fortplot revisions are pinned in
``cmake/dependency_versions.json``. Unless ``LIBNEO_REF`` is set, libneo
follows the matching SIMPLE branch when it exists in libneo, otherwise
libneo's ``main`` branch. The build machine therefore needs network access
to GitHub on the first configuration.

Get the source
--------------

.. code-block:: bash

   git clone https://github.com/itpplasma/SIMPLE.git
   cd SIMPLE

Build with CMake
----------------

The ``Makefile`` wraps CMake and Ninja:

.. code-block:: bash

   make                 # Release build in build/
   make CONFIG=Fast     # native optimisation (recommended for production runs)

This is equivalent to

.. code-block:: bash

   cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
   cmake --build build

Useful CMake options:

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Option
     - Default
     - Purpose
   * - ``CMAKE_BUILD_TYPE``
     - ``Release``
     - ``Release``, ``Fast`` (native CPU flags), ``Debug`` or ``Profile``.
   * - ``ENABLE_OPENMP``
     - ``ON``
     - OpenMP parallelisation over particles.
   * - ``ENABLE_PYTHON_INTERFACE``
     - ``ON``
     - Build the ``pysimple`` bindings (needs numpy and f90wrap).
   * - ``SIMPLE_ENABLE_CGAL``
     - ``OFF``
     - STL wall intersection with CGAL (downloads CGAL headers, needs Boost
       headers). See :doc:`/user-guide/wall-heat-loads`.
   * - ``SIMPLE_ENABLE_CUDA``
     - ``OFF``
     - Native CUDA kernels for GPU tracing. See :doc:`/user-guide/gpu`.
   * - ``SIMPLE_DETERMINISTIC_FP``
     - ``OFF``
     - Disable fast-math for bitwise reproducible results
       (``make build-deterministic``).
   * - ``LIBNEO_REF`` / ``LIBNEO_SOURCE_DIR``
     - matching branch or main
     - Select a libneo branch or use a local checkout
       (``make LIBNEO_REF=...`` or ``make LIBNEO_PATH=...``).

Several builds can live side by side in different build directories, for
example a CGAL build:

.. code-block:: bash

   cmake -S . -B build-cgal -G Ninja -DCMAKE_BUILD_TYPE=Release \
       -DSIMPLE_ENABLE_CGAL=ON
   cmake --build build-cgal --target simple.x

Python package with pip
-----------------------

The script ``setup-venv.sh`` creates a virtual environment ``.venv`` in the
repository, installs the Python requirements, builds SIMPLE and installs
``pysimple`` into the environment in editable mode:

.. code-block:: bash

   ./setup-venv.sh
   source .venv/bin/activate
   python -c "import pysimple"

The script runs the following steps, which can also be used in an existing
environment:

.. code-block:: bash

   python -m pip install -r requirements.txt
   python -m pip install --no-build-isolation -e .

``pip`` builds the Fortran code with scikit-build-core in ``_skbuild/``. The
first build takes several minutes.

Check the installation
----------------------

.. code-block:: bash

   make test-fast

runs the fast part of the test suite with a deterministic floating-point
build. :doc:`quickstart` runs a first physics case.
