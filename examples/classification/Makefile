# Classification example reproducing JPP 2020 paper results
# (Albert, Kasilov, Kernbichler, J. Plasma Phys. 86, 815860201)
#
# Configurations: QI (Drevlak 2014), QH (Drevlak 2018), QA (Henneberg 2019)
#
# Targets:
#   make run    - run all 10 cases under /tmp/simple_classification
#   make plot   - generate all paper-style plots from results
#   make clean  - remove local output files and /tmp results

SIMPLE_X = ../../build/simple.x
BASE_DIR = /tmp/simple_classification

all: plot

run: $(SIMPLE_X)
	bash run_all.sh

plot:
	source ../../.venv/bin/activate && python plot_paper_results.py $(BASE_DIR)

# Legacy volume classification (runs in-place with simple.in)
plot-local: class_parts.dat bminmax.dat
	source ../../.venv/bin/activate && python plot_classification.py

class_parts.dat bminmax.dat: simple.in
	$(SIMPLE_X)

clean:
	rm -f *.dat fort.* *.pdf
	rm -rf $(BASE_DIR)
