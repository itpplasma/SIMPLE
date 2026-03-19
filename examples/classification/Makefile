# Example for plotting classification over radius s and perpendicular invariant Jperp

all: plot

plot: class_parts.dat bminmax.dat plot_classification.py
	python plot_classification.py

class_parts.dat bminmax.dat: simple.in wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc
	../../build/simple.x

wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc:
	wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc

clean:
	rm -f *.dat fort.* *.pdf
