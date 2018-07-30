
Code of VER_1 does the following:

1. Reads magnetic field data from VMEC file in a standard NETCDF format

2. From $\iota$ profile given at VMEC flux surfaces (equidistant in $s$, firs point - axis, last point - LCMS)
builds 4-th order spline and integrates it giving 5-th order spline for vector potential component $A_\varphi$
(vector potential component $A_\vartheta$ is the toroidal flux - radial variable times flux at the edge surface).

2. Computes functions of VMEC coordinates $(s,\theta,\varphi)$, namely $R((s,\theta,\varphi))$, $Z(s,\theta,\varphi)$ 
and $\lambda(s,\theta,\varphi)$ on VMEC flux surfaces at nodes of the equidistant $(\theta,\varphi)$ grid.
Size of the grid is determined by maximum values of poloidal ("m_max") and toroidal ("n_max") mode numbers times
a given integer factor "multharm". Values tried so far are multharm=3 and multharm=7, for splining with computer accuracy
it is estimated to be around 20. The value of this parameter is set directly in the main (driver) routine 
"canonical_coordinates.f90":

  multharm=3 !7

3. Splines all three quantities in 3D using usual and periodic splines. Orders of splines can be 3 or 5. These are:
"ns_A" - spline order for vector potential (should be always 5 because second order spline needed for 3-rd order
splining is absent in the generic splining routine), "ns_s" - radial spline (over $s$) for $R$, $Z$ and $\lambda$ (can be 
between 3 and 5, fourth order has not been tried but is present), "ns_tp" - order of splines over the angles (3,...,5).
These parameters are also set in "canonical_coordinates.f90":

  ns_A=5
  ns_s=5
  ns_tp=5

4. Integrates the 1D (over radius $r=s$) ODE for all nodes of $(\vartheta_c,\varphi_c)$ grid and thus determines transformation
function $G_c(r,\vartheta_c,\varphi_c)$ on a 3D grid (currently grid size is set to be the same as of VMEC data grid).

5. Computes covariant canonical magnetic field components $B^c_\vartheta$ and $B^c_\varphi$ and metric determinant
of canonical coordinates $\sqqrt{g_c}$ (can be negative) using numerical derivatives of $G$ (7-th order).

6. Splines $B^c_\vartheta$, $B^c_\varphi$ and $\sqqrt{g_c}$ in 3D (splines also $G_c$ which is needed for transformation from
canonical to VMEC coordinates if a logical switch "fullset" is set to .true.; if one works only within canonical coordinates
this spline is not needed)

7. Computes the same guiding center orbit using different methods (two methods in VMEC coordinates and one in canonical 
coordinates). Uses for this the existing (since long) set of subroutines for guiding center orbit integration in arbitrary
curvilinear coordinates. This set interacts with magnetic field via a driver routine "magfie" which has a standard
input/output. There are 2 (or 3, depends on distribution, see below) such drivers in "magfie.f90": "magfie_vmec" - for the
field by 3D splines in VMEC coordinates and "magfie_can" - for the field in canonical coordinates. A particular routine
is selected by a control variable "isw_field_type" which is set directly in the main driver routine "canonical_coordinates.f90":

  isw_field_type=0
...
  isw_field_type=1

whee 0 is for canonical coordinates and 1 is for VMEC coordinates. See comments in "magfie" routines on input/output
variables.

8. There are two "distributions" - "VER_1_BAREBONES" and "VER_1_FULLTEST". Subroutines regarding the splining of VMEC 
coordinates and construction and splining of canonical coordinates are the same (the only difference is in one deallocation
command which is commented in "VER_1_FULLTEST"), but version "VER_1_FULLTEST" contains also the original routines by
Viktor Nemov where he integrates guiding center equations in VMEC coordinates differently (computes the field in cylindrical
coordinates and at the very end transforms the guiding center velocity to VMEC-associated coordinates; splining is used only
over $s$, over the angles Fourier series is used directly).
Distribution "VER_1_BAREBONES" uses less memory (does not keep 2 coordinate systems in the splined form at the same time).
It is recommended for the further development.

9. Usage: go to "VER_1_BAREBONES" or "VER_1_FULLTEST" and run

make -f Canonical_coordinates.mk 

change there to "RUN" subdirectory and run

../canonical_coordinates.x

There are two input files in "RUN" subdirectories - "alpha_lifetime_m.inp" - control file including various parameters 
for orbit integration and magnetic field location, "netcdf_file_for_test.nc" - VMEC output file.
Output files are "fort.3001" - orbit via splined VMEC coordinates, "fort.3002" - orbit by Viktor method (only in
"VER_1_FULLTEST") and "fort.3003" - orbit via canonical (splined) coordinates. See the main driver routine 
"canonical_coordinates.f90" for the meaning of columns.

