Code traces a set of guiding center orbits starting from a given flux surface.
Starting values of coordinates and velocity are chosen randomly assuming
the isotropic mono-energetic distribution constant on a starting flux surface.
Orbits are traced until they leave the confinement volume or until tracing
time is exhausted.

Phase space variables are $(x,\varphi,z,\lambda,\hat p)$ where spatial coordinates
$x = N_p s \cos(\theta)$ and $z = N_p s \sin(\theta)$ are expressed via VMEC
coordinates $(s,\theta,\varphi)$ being the normalized toroidal flux, poloidal
angle and toroidal angle, and $N_p$ is the number of flux surfaces in VMEC file.
Velocity space variables are $\lambda = v_\parallel / v$ - pitch parameter
and $\hat p = v / v_0$ - velocity module normalized by the starting value $v_0$

Orbit parameter is the normalized time $\tau = v_0 t$ (has the dimension of length).

Boundary of the confinement volume is given by $s N_p = N_p - 1$.

All equations are in Gaussian units.

Random start from the flux surface is modelled by a start from a single, long enough
field line taking the specific volume (flux tube volume) into account.

For details of starting field line and other computation parameters see comments
on input and output files below.

Present code version has been compiled on Linux computers using gfortran.
Code does not use external libraries, the only external routine is "random_number"
which is intrinsic in fortran95 (present in gfortran, may be absent in other 
compilers).

All ODE are solved with an adaptive RK4-5 integrator using "odeint" set of routines
from Numerical Recipes. Relative error is set by parameter statements in the file
"sub_alpha_lifetime.f90" on lines 194,

      double precision, parameter :: relerr=1d-9

and 334,

      double precision, parameter :: relerr=1d-9

To compile the code on Linux machines execute in root directory

make -f Alpha_lifetime_vmec.mk

what produces and executable "alpha_lifetime.x".
To run the code change to subdirectory RUN and execute there

../alpha_lifetime.x

One of output files obtained for the present input, "confined_fraction.dat.test",
is in the root directory for the comparison (result may be different because of
computer architecture).

---------------------------------------------------------------------------
---------------------------------------------------------------------------

Input files:


1. "alpha_lifetime.inp" has comments inside. File is read by free format, 
[read(1,*) data]. In more details input variables are:

1) notrace_passing - integer, takes values 0 or 1. 
If notrace_passing=1 only trapped particles are traced, while passing particles 
are assumed to be ideally confined. If notrace_passing is different from 1 all
particles are traced.

2) nper - integer, number of toroidal field periods which are covered by a starting
field line representing the starting magnetic surface.

3) npoiper - integer, number of equidistant points per field period on a starting 
field line where data is stored. This data is used for random generation of initial
particle coordinates.

4) ntimstep - integer, number of "time steps" per slowing down time. After each
of this time steps number of particles remaining in the volume is recorded.

5) ntestpart - integer, initial number of test particles traced

6) bmod_ref - double precision, prescribed amplitude of $B_{00}$ harmonic of the 
magnetic field module in Boozer coordinates, Gauss. Magnetic field from VMEC file 
is rescaled to fulfill the prescription.

7) trace_time - double precision, maximum time test particles are followed, seconds.

8) sbeg - double precision, value of flux surface label (normalized toroidal flux)
for the starting magnetic surface.

9) phibeg - double precision, value of toroidal angle of VMEC coordinates at the 
initial point of the starting field line.

10) thetabeg - double precision, value of poloidal angle of VMEC coordinates at the 
initial point of the starting field line.

11) loopskip - integer, run number (starts from 0). Used to enable independent 
statistics if few runs are performed independently (e.g. in parallel computations
using CONDOR script). If loopskip > 0 random numbers are subsequently generated
for loopskip runs and ignored. Only after that random numbers are used for generation
of starting phase space coordinates of test particles (i.e. loopskip shifts the 
random sequence).

12) contr_pp - double precision, shifts trapped-passing boundary computed for zero
Larmor radius into the passing region so that all particles regarded as passing with
respect to a new boundary are indeed passing in case of finite Larmor radius.
Normally - B_min/(B_max-B_min) < contr_pp < 0. Set notrace_passing=0 to activate this.

13) facE_al - double precision, ratio of fusion alpha particle energy (3.5 MeV) 
to the desired test particle energy (specified in MeV).

14) npoiper2 - integer, additional split factor for maximum integration time step.
Maximum step in normalized time ($\tau=v_0 t$) is defined as $\Delta\tau=2\pi R/N$ 
where N = L1i*npoiper*npoiper2. Here, L1i - number of toroidal field periods and 
npoiper is described in 3) above.

15) n_e - integer, the same as Z - charge number of test particle (n_e=2 for alphas)

16) n_d - integer, the same as A - mass number of test particle (n_d=4 for alphas)

17) filed - character*24, ascii VMEC data file name

...................................................................................

2. "wout_26.scaled" - ascii VMEC data file, has comments inside. Numerical data is read
by routines from the file "vmecinm.f" using a free format [read(71,*) data,data,...].
If number of flux surfaces, number of Fourier harmonics or number of field
periods (respectively, these are variables nsurf, maximum value of nstr or nper 
defined in the list below) are different from a given example, one should modify 
variable values in the parameter statement in the file "parvmec.f" where:
nsurfm - number of flux surfaces (currently nsurfm=99)
nstrm  - number of Fourier harmonics (currently nstrm=288)
nper   - number of field periods (currently nper=5)

More detailed comments:

Global data:

1) nsurf - integer, number of flux surfaces.

2) nper  - integer, number of field periods.

3) flux  - double precision, toroidal flux at the last closed surface (SI).

4) Rmajor_p - double precision, major radius (SI).

5) Aminor_p - double precision, minor radius (SI).

Local data (per flux surface):
Scalars:

6) ns - integer, flux surface index (starts from 0).

7) iotaf - double precision, iota value at this flux surfaces.

8) phi - double precision, toroidal flux at this flux surfaces (SI).

Vectors:

9)  xm - double precision, poloidal Fourier index m

10) xn - double precision, toroidal Fourier index n. This index corresponds to the
expansion over toroidal angle changing between 0 and 2*pi (xn are multiples of nper).

11) rmnc - double precision, amplitudes of Fourier expansion of R over cos(m*theta-n*phi) 

12) zmns - double precision, amplitudes of Fourier expansion of Z over sin(m*theta-n*phi)

13) lmns - double precision, amplitudes of Fourier expansion of stream function lambda
over sin(m*theta-n*phi)

14) nstr - integer, (m,n) harmonic number in a general 1D list

--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------

Output files:

1. "alpha_lifetime.log" - contains input parameters from "alpha_lifetime.inp"
and some extra computation parameters. The later are:

1) Rbig - major radius (cm);
2) dphi - integration step in phi for the starting field line;
3) v0 - test particle velocity module (cm/s);
4) rlarm - Larmor radius of test particle in the reference magnetic field (cm);
5) dtau - step value in normalized time ($\tau=v_0 t$) for counting particles remaining
in the confinement volume (cm);
6) dtaumin - maximum normalized time step for orbit integration (cm).

......................................................................................

2. "confined_fraction.dat" has 4 columns:

1) time slice (seconds);
2) fraction of trapped particles remaining in the volume at this time slice (number
of trapped particles which have reached this time slice divided by ntestpart -total 
number of test particles to be traced);
3) fraction of passing particles remaining in the volume at this time slice (see above);
4) number of test particles traced so far. This is needed for recalculation of trapped 
and passing particle fractions during the run (equals ntestpart if the run is finished).

.......................................................................................

3. "start_p.dat" has 5 columns with the following orbit data for all passing test 
particles:

1) trapping parameter $\kappa = (1-\lambda^2)B_{max}/B-1)B_{min}/(B_{max}-B_{min})$
evaluated at the starting point where $B$ and $\lambda=v_\parallel/v$ are magnetic
field module and pitch parameter at this point, and $B_{min}$ and $B_{max}$ are
minimum and maximum values of $B$ at the starting surface. This parameter is
$0 < \kappa < 1$ for trapped particles and $-B_{min}/(B_{max}-B_{min}) < \kappa < 0$
for passing particles;
2) step number when particle has been lost (out of ntimstep steps described above);
3,4,5) - starting coordinates of test particle.

.......................................................................................

4. "start_t.dat" - the same as "start_p.dat" for trapped particles.

.......................................................................................

5. "bminmax.dat" - minimum $B$, maximum $B$ and Boozer $B_{00}$ harmonic of $B$
at the starting surface (Gauss).

.......................................................................................

6. "volstart.dat" - specific volume as function of field line length,
$\int_0^l \rd l / B$ normalized by its final value at $l=L$. Columns:

1) point index on a starting field line
2) normalized specific volume at this point

.......................................................................................

7. "starting_surface.dat" - Poincare plot of the starting field line at the poloidal
cross-section phi=phibeg (see phibeg above). Columns:

1) coordinate $x=N_s s \cos \theta$ where N_s$ - number of flux surfaces in VMEC file;
2) coordinate $z=N_s s \sin \theta$.

.......................................................................................

8. "chamb0.dat" (4 columns) - exit coordinates of lost particles. Columns:

1,2) coordinates x and y (see above);
3) index of lost particle;
4) toroidal coordinate phi.
