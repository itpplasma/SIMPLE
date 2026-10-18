Clusters
========

This page shows how to run SIMPLE under Slurm and which job scripts in the
repository fit which machine. One SIMPLE process traces all markers of a run
with OpenMP on one node, or on one GPU; there is no MPI decomposition. Scale
out by submitting independent jobs (equilibria, seeds, ``batch_size``
batches), each in its own run directory, because output file names are fixed.

Job pattern
-----------

.. code-block:: bash

   #SBATCH --nodes=1
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=<cores>
   export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
   export OMP_PLACES=cores OMP_PROC_BIND=spread OPENBLAS_NUM_THREADS=1
   srun --cpu-bind=cores /path/to/SIMPLE/build/simple.x simple.in

Submit from a run directory that contains ``simple.in`` and the files it
references. Build on the target architecture; ``CONFIG=Fast`` tunes for the
build host.

GPU jobs
--------

GPU runs use the native CUDA backend (``gpu_backend = 'cuda-native'``,
:doc:`/user-guide/gpu`), one process per GPU. Build with
``-DSIMPLE_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=<arch>``:

=========================== =================================== ============================
System                      GPU                                 ``CMAKE_CUDA_ARCHITECTURES``
=========================== =================================== ============================
aCluster                    NVIDIA Tesla T4                     ``75``
Perlmutter, CINECA Leonardo NVIDIA A100                         ``80``
CINECA Pitagora (Booster)   NVIDIA H100                         ``90``
sCluster                    NVIDIA RTX PRO 6000 Blackwell Max-Q ``120``
=========================== =================================== ============================

On Perlmutter the GNU environment suffices for the Fortran host code:

.. code-block:: bash

   module load PrgEnv-gnu cudatoolkit/12.9 craype-accel-nvidia80 cray-hdf5 cray-netcdf
   cmake -S . -B build-cuda -G Ninja -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn \
     -DSIMPLE_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=80
   cmake --build build-cuda -j

Request the GPU with ``--gpus-per-task=1`` (Perlmutter) or ``--gres=gpu:1``
(aCluster, sCluster).

Job scripts
-----------

* ``RUN/{Perlmutter,aCluster,sCluster}/submit_gpu.slurm``: one GPU; set
  ``SIMPLE_BUILD`` (default ``$HOME/code/SIMPLE/build-cuda``) and
  ``SIMPLE_INPUT``.
* ``examples/hpc/cineca/``, ``examples/hpc/pitagora/``: Leonardo smoke test
  and untested Pitagora templates (:doc:`cineca`).

.. toctree::
   :maxdepth: 1

   cineca
