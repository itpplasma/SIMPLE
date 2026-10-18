Post-processing and visualization
=================================

``simple-plot`` reads existing SIMPLE output and never runs the tracer. It uses
Matplotlib, as do the Python plotting helpers. From a source checkout, install
it with the optional plotting dependencies::

   python -m pip install '.[plotting]'

Save figures with ``-o`` (the extension selects PNG, PDF, SVG, etc.), or omit it
to open an interactive plot. For a headless machine, set ``MPLBACKEND=Agg``::

   simple-plot losses run/ -o losses.png
   simple-plot losses run/ --confidence 0.99 -o losses-99.png
   simple-plot classes run/ -o classes.png
   simple-plot wall run/results.nc --color incidence -o wall.png
   simple-plot poincare run/ -o poincare.png

Losses and Monte Carlo uncertainty
----------------------------------

``losses`` reads ``confined_fraction.dat``: time in seconds, confined passing
fraction, confined trapped fraction, and resolved marker count. It plots the
confined total and both components, plus the lost fraction, on a logarithmic
time axis. The components use the same denominator and sum to the total; they
are not separately normalized to their initial populations. Zero-time rows are
omitted, and rows with no resolved markers appear as gaps.

The shaded band is a pointwise Wilson binomial confidence interval on the
lost fraction. At each time, ``n`` is column 4 (the resolved count), and ``k`` is
``n * (1 - passing - trapped)``, rounded to the nearest integer to recover the
marker count from the printed fractions. Legacy three-column files instead
use the number of rows in ``times_lost.dat`` for ``n``. With ``p = k/n`` and
standard-normal quantile ``z`` at ``(1 + confidence)/2``, the limits are::

   (p + z^2/(2n) +/- z*sqrt(p*(1-p)/n + z^2/(4n^2))) / (1 + z^2/n)

The default confidence is 95%; ``--confidence`` changes it. This assumes
independent, equally weighted markers. The interval describes sampling error
conditional on resolved outcomes, not numerical error or a simultaneous band
for the entire curve. Numerically unresolved markers are excluded rather than
counted as physical losses; inspect ``unresolved_fraction.dat`` separately.

Classes, wall hits and Poincare sections
----------------------------------------

``classes`` reads ``class_parts.dat`` (index, initial reference radial coordinate,
perpendicular invariant, three class codes), or falls back to ``iclass`` in
``results.nc``. Three panels
count the J_parallel, topology and fractal codes, including 0 (unclassified).
Class codes are shown as written by SIMPLE; a NetCDF run with classification
disabled can contain only zeros.

``wall`` requires ``results.nc``, produced with ``output_results_netcdf = .True.``.
It selects only markers with ``wall_hit = 1`` and uses the authoritative
``wall_hit_cart`` intersection, not the final orbit position or generic loss
flag. The cylindrical R--Z projection uses ``cart_units`` from the file and
combines toroidal positions. The default energy color is the dimensionless
``E/E0 = (v/v0)^2`` from ``zend`` (not energy in eV); ``--color incidence`` uses
``wall_hit_cos_incidence`` (0 = grazing, 1 = normal). Wall hits require a run
with STL wall support and ``wall_input`` configured; ordinary LCFS losses do
not supply a wall footprint.

``poincare`` discovers ``poincare_vol<N>.dat`` files from
``spectre_poincare.x`` (seed, section, volume-local s, theta, R, Z). It plots
R--Z in the equilibrium's native length units, colored by volume-local s on a
shared scale. Theta and SIMPLE's other angle outputs are in radians. These
are field-line sections, distinct from orbit-class statistics.

Each subcommand reports missing files or required variables with an error and
exit status 2. It reads no input configuration. For reference,
``times_lost.dat`` contains index, loss time, trapping parameter, initial
reference radial coordinate, perpendicular invariant and the five final
phase-space components; ``start.dat`` contains the five initial components
``(u1, theta, phi, v/v0, v_par/v)`` in reference coordinates. Here ``u1`` is
``s`` for VMEC and ``rho = sqrt(s)`` for Boozer chartmaps; other reference
charts retain their native radial coordinate. Neither file is required for
current four-column confinement output.
