# Independent of SIMPLE's CMake build. Invoke from the repository root.
FC = gfortran
FFLAGS ?= -O2 -march=native -fopenmp -Wall -Wextra -fimplicit-none -std=f2018
ROOT := concept/interface_candidates
BUILD := $(ROOT)/build
MODULES := sheet_field fo_sheet rk_reference sheet_setup gc_hybrid candidate_field transforms candidate_orbits
OBJECTS := $(addprefix $(BUILD)/,$(addsuffix .o,$(MODULES)))

all: $(BUILD)/compare.x $(BUILD)/test_candidates.x $(BUILD)/study_hybrid.x

$(BUILD):
	mkdir -p $@

$(BUILD)/%.o: $(ROOT)/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/fo_sheet.o $(BUILD)/rk_reference.o $(BUILD)/sheet_setup.o $(BUILD)/gc_hybrid.o: $(BUILD)/sheet_field.o
$(BUILD)/candidate_field.o: $(BUILD)/sheet_field.o
$(BUILD)/candidate_field.o: $(ROOT)/torsion_kernel.inc
$(BUILD)/transforms.o: $(BUILD)/candidate_field.o
$(BUILD)/candidate_orbits.o: $(BUILD)/transforms.o $(BUILD)/fo_sheet.o

$(BUILD)/%.x: $(ROOT)/%.f90 $(OBJECTS)
	$(FC) $(FFLAGS) -I$(BUILD) $< $(OBJECTS) -o $@

test: all
	python3 $(ROOT)/proofs/verify.py
	OMP_NUM_THREADS=4 $(BUILD)/test_candidates.x

study: all
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-16} python3 $(ROOT)/run_study.py
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-16} python3 $(ROOT)/run_study.py --long
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-16} python3 $(ROOT)/run_study.py --sampling

summarize:
	python3 $(ROOT)/pack_results.py
	python3 $(ROOT)/summarize.py
	python3 $(ROOT)/render_report.py

accuracy: $(BUILD)/accuracy.x $(BUILD)/sheet_accuracy.x
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-16} $(BUILD)/accuracy.x
	OMP_NUM_THREADS=$${OMP_NUM_THREADS:-16} $(BUILD)/sheet_accuracy.x

clean:
	rm -rf $(BUILD)

.PHONY: all test study summarize accuracy clean
