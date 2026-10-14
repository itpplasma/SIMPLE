# Invoke from the SIMPLE repository root. Entirely separate from its build.
ROOT := concept/gc_higher_order
BUILD := $(ROOT)/build
FORTSYM ?= /home/ert/code/fortsym
FC := gfortran
FFLAGS ?= -O2 -fopenmp -Wall -Wextra -fimplicit-none -std=f2018 -fPIC
# Large straight-line CSE leaves make GCC register allocation pathological.
# Their algebra is already optimized by fortsym; keep compilation bounded.
GEN_FFLAGS ?= -O0
export OMP_NUM_THREADS := 16

all: $(BUILD)/test_gc

$(BUILD):
	mkdir -p $@

$(BUILD)/derive_gc: $(ROOT)/gen/derive_gc.f90 | $(BUILD)
	$(FC) -O0 -g -fcheck=all -ffree-line-length-none -I$(FORTSYM)/build/include $< \
	    -L$(FORTSYM)/build/lib -Wl,-rpath,$(FORTSYM)/build/lib -lfortsym -o $@

regen: $(BUILD)/derive_gc
	$(BUILD)/derive_gc $(ROOT)/generated/gc_generated.f90

check-generated: $(BUILD)/derive_gc
	$(BUILD)/derive_gc $(BUILD)/regenerated.f90
	cmp $(ROOT)/generated/gc_generated.f90 $(BUILD)/regenerated.f90

$(BUILD)/gc_generated.o: $(ROOT)/generated/gc_generated.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(GEN_FFLAGS) -J$(BUILD) -c $< -o $@

$(BUILD)/gc_second_order.o: $(ROOT)/gen/gc_second_order.f90 $(BUILD)/gc_generated.o
	$(FC) $(FFLAGS) -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/gc_api.o: $(ROOT)/gc_api.f90 $(BUILD)/gc_second_order.o
	$(FC) $(FFLAGS) -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/c_bridge.o: $(ROOT)/tests/c_bridge.f90 $(BUILD)/gc_api.o
	$(FC) $(FFLAGS) -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/libgc.so: $(BUILD)/gc_generated.o $(BUILD)/gc_second_order.o \
                   $(BUILD)/gc_api.o $(BUILD)/c_bridge.o
	$(FC) $(FFLAGS) -shared $^ -o $@

$(BUILD)/test_gc: $(ROOT)/tests/test_gc.f90 $(BUILD)/gc_api.o $(BUILD)/gc_second_order.o $(BUILD)/gc_generated.o
	$(FC) $(FFLAGS) -fcheck=all -I$(BUILD) $^ -o $@

test: $(BUILD)/test_gc
	$(BUILD)/test_gc
	python3 $(ROOT)/gen/prove_identities.py
	python3 $(ROOT)/gen/prove_second_order.py

# Canonical kernels are committed, so these tests do not regenerate them.
canonical-test: $(BUILD)/libgc.so
	python3 $(ROOT)/gen/prove_canonical.py
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 $(ROOT)/tests/canonical_checks.py

canonical-study: $(BUILD)/libgc.so
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 $(ROOT)/tests/canonical_study.py energy
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 $(ROOT)/tests/canonical_study.py orders
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 $(ROOT)/tests/canonical_study.py cost

canonical-regen:
	python3 $(ROOT)/gen/derive_canonical.py

canonical-check-generated: | $(BUILD)
	python3 $(ROOT)/gen/derive_canonical.py $(BUILD)/canonical_regenerated.py
	cmp $(ROOT)/generated/canonical_kernels.py $(BUILD)/canonical_regenerated.py

canonical-check-results:
	python3 $(ROOT)/tests/check_canonical_results.py

study: $(BUILD)/libgc.so
	python3 $(ROOT)/tests/orbit_study.py
	python3 $(ROOT)/tests/slab_action.py

check-orders:
	python3 $(ROOT)/tests/check_orders.py

lean:
	cd $(ROOT)/lean && lake env lean Structure.lean

clean:
	rm -rf $(BUILD)

.PHONY: all regen check-generated test study check-orders lean clean \
        canonical-test canonical-study canonical-regen \
        canonical-check-generated canonical-check-results
