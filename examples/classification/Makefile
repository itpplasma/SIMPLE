# Classification example reproducing JPP 2020 paper results
# (Albert, Kasilov, Kernbichler, J. Plasma Phys. 86, 815860201)
#
# Targets:
#   make run    - run all three cases (s03, s06, volume) under /tmp
#   make plot   - generate all paper-style plots from results
#   make clean  - remove local output files
#
# The existing volume classification (plot_classification.py) is also retained.

SIMPLE_X = ../../build/simple.x
BASE_DIR = /tmp/simple_classification
WOUT = wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc
WOUT_URL = https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/$(WOUT)

all: plot

# Run all three cases in /tmp
run: $(SIMPLE_X) $(WOUT)
	bash run_all.sh

# Generate paper-style plots from /tmp results
plot: $(BASE_DIR)/s03/times_lost.dat
	python plot_paper_results.py $(BASE_DIR)

# Dependency: results must exist (run first if not)
$(BASE_DIR)/s03/times_lost.dat: $(SIMPLE_X) $(WOUT)
	bash run_all.sh

# Legacy volume classification (runs in-place)
plot-local: class_parts.dat bminmax.dat
	python plot_classification.py

class_parts.dat bminmax.dat: simple.in $(WOUT)
	$(SIMPLE_X)

$(WOUT):
	wget -q $(WOUT_URL)

clean:
	rm -f *.dat fort.* *.pdf
	rm -rf $(BASE_DIR)
