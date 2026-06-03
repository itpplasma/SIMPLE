# Classification example reproducing JPP 2023 paper results
# (Albert, Buchholz, Kasilov, Kernbichler, Rath, J. Plasma Phys. 89, 955890301)
#
# Configurations: QI (Drevlak 2014), QH (Drevlak 2018), QA (Henneberg 2019)
#
# Targets:
#   make run    - run all 12 cases under /tmp/simple_classification
#   make plot   - generate all paper-style plots from results
#   make clean  - remove /tmp results

SIMPLE_X = ../../build/simple.x
BASE_DIR = /tmp/simple_classification

all: plot

run: $(SIMPLE_X)
	bash run_all.sh

plot:
	source ../../.venv/bin/activate && python plot_paper_results.py $(BASE_DIR)

clean:
	rm -rf $(BASE_DIR)
