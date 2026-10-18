GPU tracing
===========

SIMPLE traces collisionless guiding-centre losses in Boozer coordinates on an
NVIDIA GPU with native CUDA kernels. The GPU run computes the segmented loss
objective of `arXiv:2606.19523 <https://arxiv.org/abs/2606.19523>`_, not the full set of CPU options.

Build and run
-------------

Needs the CUDA toolkit, a C++17 compiler and GNU or NVIDIA Fortran.

.. code-block:: bash

   cmake -S . -B build-cuda -G Ninja -DSIMPLE_ENABLE_CUDA=ON
   cmake --build build-cuda -j
   (cd run && ../build-cuda/simple.x simple.in)

.. code-block:: fortran

   &config
   isw_field_type = 2
   gpu_mode = 'production'
   gpu_backend = 'cuda-native'
   gpu_method = 'dopri'
   /

Requirements: ``integ_coords = 2``, no collisions, no wall, no orbit output,
no classification, guiding-centre orbits. Other combinations stop with an
error message.

Keys
----

==========================  =================  ==========================================
Key                         Default            Meaning
==========================  =================  ==========================================
``gpu_mode``                ``'off'``          ``'production'`` runs on the GPU
``gpu_backend``             ``'cuda-native'``  only ``'cuda-native'`` (= ``'auto'``) runs
``gpu_method``              ``'dopri'``        ``'dopri'`` (Dormand-Prince) or
                                               ``'cash-karp'``; tolerance ``relerr``
``gpu_landreman``           ``.true.``         segmented loss objective; must stay on
``gpu_t_block``             ``1d-3``           segment length in s, at most ``trace_time``
``gpu_loss_tau``            ``1d-1``           decay time :math:`\tau_L` (s) of the loss
                                               weight :math:`e^{-t/\tau_L}`
``gpu_maxloss``             ``2d-2``           stop threshold, see below
``gpu_min_timestep``        ``1d-10``          smallest Runge-Kutta step in s
``gpu_start_coordinates``   ``'reference'``    ``'boozer'``: ``start.dat`` is in Boozer
                                               coordinates
``gpu_compact_init``        ``.false.``        compact initializer; needs a Boozer chartmap
                                               and ``gpu_start_coordinates = 'boozer'``
``gpu_particle_profile``    ``''``             CSV with start/end state, field evaluations,
                                               loss step and loss time (s) per marker
==========================  =================  ==========================================

Early stop
----------

After each segment SIMPLE evaluates :math:`L = N^{-1}\sum_{\text{lost } i}
e^{-t_i/\tau_L}`. When :math:`L >` ``gpu_maxloss``, tracing stops and the
markers still in the plasma are written as confined to ``trace_time`` (exit
code 0). The log line ``observed until`` gives the time actually traced.
For loss fractions over the full ``trace_time``, set ``gpu_maxloss = 1``,
which :math:`L \le 1` never exceeds. The stop rule is in
``src/simple_cuda_native_rk54.cu``.

The run writes the usual :doc:`outputs` and uses one GPU per process. Tests
(with CUDA): ``test_cuda_native_rk54``, ``test_cuda_native_bridge``.
