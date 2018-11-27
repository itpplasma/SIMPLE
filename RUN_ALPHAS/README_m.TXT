Version for VMEC file in netcdf format. 

Differs from version the original version for ascii input (see README.txt) by
the names of makefile, executable file and by input files.

Code uses netcdf libraries and include files. Paths in the makefile corresponds to DEBIAN
Linux:
NCINC = -I/usr/include
NCLIB = -L/usr/lib -lnetcdff -lnetcdf

To compile the code on Linux machines execute in root directory

make -f Alpha_lifetime_vmec_m.mk

what produces and executable "alpha_lifetime_m.x".
To run the code change to subdirectory RUN and execute there

../alpha_lifetime_m.x

---------------------------------------------------------------------------
---------------------------------------------------------------------------

Input files:


1. "alpha_lifetime_m.inp" differs from "alpha_lifetime.inp" (see README.txt) by the last 
line where the name of netcdf VMEC file is specified instead of ascii VMEC file name

2. VMEC file in netcdf format
