# Invoke from the SIMPLE repository root. Entirely separate from its build.
ROOT := concept/gc_higher_order
BUILD := $(ROOT)/build
# fortsym branch feat/lie-pullback at 1372fbd, built with CMake in build/.
FORTSYM ?= ../fortsym
FC := gfortran
FFLAGS ?= -O2 -fopenmp -Wall -Wextra -fimplicit-none -std=f2018 -fPIC
# Large straight-line CSE leaves make GCC register allocation pathological.
# Their algebra is already optimized by fortsym; keep compilation bounded.
GEN_FFLAGS ?= -O0
export OMP_NUM_THREADS := 16
# Generated kernels are build products; these hashes pin the reviewed output.
GC_GENERATED := $(BUILD)/gc_generated.f90
GC_GENERATED_SHA256 := 2067ae6dab198dae0a3e3ddf7d0d40d0c006899303729fc544d2c569c9a1c192
CANONICAL_KERNELS := $(BUILD)/canonical_kernels.py
CANONICAL_KERNELS_SHA256 := 8ae1ffa57d4dbfdb299797edcfc97181f5c35b70425207e7eb729f2b320cfe87

all: $(BUILD)/test_gc

$(BUILD):
	mkdir -p $@

$(BUILD)/derive_gc: $(ROOT)/gen/derive_gc.f90 | $(BUILD)
	@test -f $(FORTSYM)/build/lib/libfortsym.so || \
	    { echo 'Set FORTSYM to a built fortsym feat/lie-pullback checkout'; exit 1; }
	$(FC) -O0 -g -fcheck=all -ffree-line-length-none -I$(FORTSYM)/build/include $< \
	    -L$(FORTSYM)/build/lib -Wl,-rpath,$(abspath $(FORTSYM))/build/lib -lfortsym -o $@

$(GC_GENERATED): $(BUILD)/derive_gc
	$(BUILD)/derive_gc $@

regen: $(BUILD)/derive_gc
	$(BUILD)/derive_gc $(GC_GENERATED)

check-generated: $(GC_GENERATED)
	echo '$(GC_GENERATED_SHA256)  $<' | sha256sum -c -

$(BUILD)/gc_generated.o: $(GC_GENERATED)
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

$(CANONICAL_KERNELS): $(ROOT)/gen/derive_canonical.py | $(BUILD)
	python3 $< $@

canonical-test: $(BUILD)/libgc.so $(CANONICAL_KERNELS)
	python3 $(ROOT)/gen/prove_canonical.py
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 $(ROOT)/tests/canonical_checks.py

canonical-study: $(BUILD)/libgc.so $(CANONICAL_KERNELS)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 $(ROOT)/tests/canonical_study.py energy
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 $(ROOT)/tests/canonical_study.py orders
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 $(ROOT)/tests/canonical_study.py cost

canonical-check-generated: $(CANONICAL_KERNELS)
	echo '$(CANONICAL_KERNELS_SHA256)  $<' | sha256sum -c -

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
        canonical-test canonical-study \
        canonical-check-generated canonical-check-results
