Clusters
========

This page shows how to run SIMPLE under Slurm and which job scripts in the
repository fit which machine. One SIMPLE process traces its markers with
OpenMP on one node or on the GPUs it sees. A build with MPI spreads the
markers of one run over several processes and nodes (`MPI and several
GPUs`_). Independent runs (equilibria, seeds, ``batch_size`` batches) go into
separate jobs, each in its own run directory, because output file names are
fixed.

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

MPI and several GPUs
--------------------

Configure with ``-DSIMPLE_ENABLE_MPI=ON`` (needs an MPI Fortran compiler
wrapper; the option changes nothing when off). Marker ``i`` goes to rank
``mod(i - 1, nranks)``; every rank reads the input and sets up the field,
rank 0 broadcasts the start positions, collects the results and alone writes
the output files, which equal those of a run with one process byte for byte.

.. code-block:: bash

   #SBATCH --nodes=4
   #SBATCH --ntasks-per-node=1
   #SBATCH --cpus-per-task=<cores>
   export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
   mpirun -np $SLURM_NTASKS --map-by ppr:1:node:pe=$SLURM_CPUS_PER_TASK \
     -x OMP_NUM_THREADS /path/to/build-mpi/simple.x simple.in

Several ranks reject ``gpu_mode = 'benchmark'``, orbit output,
``class_plot``, ``restart`` and ``gpu_particle_profile``, and need at least
one marker per rank.

On the GPU (``gpu_mode = 'production'``) a process uses ``gpu_num_devices``
GPUs, all visible ones for ``0``, with one host thread per GPU and the
markers dealt to them in turn. Ranks on one node start at different GPUs, so
one rank per GPU (``gpu_num_devices = 1``) and one rank per node with all its
GPUs both work. The early stop (``gpu_maxloss``) is decided by each GPU and
rank for its own markers, so a split run equals the single-GPU run only with
``gpu_maxloss = 1``; with a smaller threshold a process traces on one GPU,
and several ranks stop at input validation.

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
