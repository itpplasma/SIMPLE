"""Python interface to SIMPLE.

pysimple is a thin layer over the Fortran module ``simple_api``. A run is
configured by the ``&config`` namelist of ``simple.in``; every keyword argument
of :func:`init` is a ``&config`` name (see ``examples/simple_full.in``)::

    import pysimple

    pysimple.init(config="simple.in")             # exactly as simple.x
    pysimple.init("wout.nc", trace_time=1e-3)     # equilibrium plus settings
    pysimple.init("wout.nc", config="simple.in")  # file with overrides

    start = pysimple.sample()                     # as configured by startmode
    result = pysimple.trace_parallel(start)

States are ``(s, theta, phi, v/v0, v_par/v)`` in reference coordinates and
have shape ``(5,)`` for one particle or ``(5, n)`` for ``n`` particles.
Angles are radians. ``pysimple.params`` exposes the Fortran ``params`` module
for reading derived values such as ``params.v0``.
"""

from importlib import import_module

# Keep pure post-processing usable without loading the compiled tracer.
__all__ = [
    "init", "sample", "sample_surface", "sample_volume", "load_particles",
    "trace_parallel", "classify_parallel", "classify_fast", "classify_fractal",
    "write_output", "trace_orbit", "compute_canonical_frequencies",
    "invariants_from_state", "states_from_invariants",
    "states_from_potato_invariants", "trace_to_cut", "params",
    "RK45", "EXPL_IMPL_EULER", "IMPL_EXPL_EULER", "MIDPOINT", "GAUSS1",
    "GAUSS2", "GAUSS3", "GAUSS4", "LOBATTO3",
    "ORBIT_CLASS_UNKNOWN", "ORBIT_CLASS_TRAPPED", "ORBIT_CLASS_PASSING",
    "FREQ_SUCCESS", "FREQ_INVALID_INPUT", "FREQ_ORBIT_LOST",
    "FREQ_INTEGRATOR_ERROR", "FREQ_MAX_STEPS",
    "INVARIANT_SUCCESS", "INVARIANT_INVALID_INPUT", "INVARIANT_NOT_AXISYMMETRIC",
    "INVARIANT_NO_INTERSECTION", "INVARIANT_ROOT_FAILURE", "INVARIANT_CAPACITY",
    "SECTION_EXTREMUM_UNKNOWN", "SECTION_B_MINIMUM", "SECTION_B_MAXIMUM",
    "CUT_ANY", "CUT_TIP", "CUT_TOROIDAL",
]
__all__ += ["cylindrical_to_reference", "reference_to_cylindrical"]


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    tracing = import_module('._tracing', __name__)
    for export in __all__:
        globals()[export] = getattr(tracing, export)
    return globals()[name]


def __dir__():
    return sorted(set(globals()) | set(__all__))
