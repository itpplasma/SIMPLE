# Example for plotting classification over radius s and perpendicular invariant Jperp

FC = gfortran
FFLAGS = -g -fbacktrace

all: plot

plot: prompt1.dat regular1.dat stochastic1.dat prompt2.dat regular2.dat stochastic2.dat plot_classification.py
	python plot_classification.py

prompt1.dat regular1.dat stochastic1.dat prompt2.dat regular2.dat stochastic2.dat: class_parts.dat postprocess_class.x
	./postprocess_class.x

class_parts.dat: simple.in wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc
	../../build/simple.x

wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc:
	wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc

postprocess_class.x: postprocess_class.f90
	$(FC) $(FFLAGS) -o postprocess_class.x postprocess_class.f90

clean:
	rm -f *.dat fort.* *.pdf
