CONFIG ?= Release
FLAGS ?=
BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

# Common ctest command with optional verbose and test name filtering
CTEST_CMD = cd $(BUILD_DIR) && ctest --test-dir test --output-on-failure $(if $(filter 1,$(VERBOSE)),-V) $(if $(TEST),-R $(TEST))

# NVIDIA HPC SDK paths for nvfortran builds
NVHPC_ROOT := /opt/nvidia/hpc_sdk/Linux_x86_64/25.11
NVHPC_HPCX := $(NVHPC_ROOT)/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi
NVHPC_BUILD_DIR := build_nvfortran

.PHONY: all configure reconfigure build build-deterministic build-deterministic-nopy test test-fast test-slow test-regression test-all test-golden-main test-golden-tag test-golden install clean nvfortran nvfortran-configure nvfortran-clean
all: build

$(BUILD_NINJA):
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON $(FLAGS)

configure: $(BUILD_NINJA)

reconfigure:
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

build: configure
	cmake --build $(BUILD_DIR) --config $(CONFIG)

# Test targets
# Usage: make [test-target] [TEST=test_name] [VERBOSE=1]
# Example: make test TEST=test_gvec VERBOSE=1
# Example: make test-golden-main TEST=classifier VERBOSE=1

# Run all tests except regression tests (default)
test: build
	$(CTEST_CMD) -LE "regression"

# Run only fast tests (exclude slow and regression tests)
test-fast: build
	$(CTEST_CMD) -LE "slow|regression"

# Run only slow tests
test-slow: build
	$(CTEST_CMD) -L "slow" -LE "regression"

# Run only regression tests (requires deterministic FP build without Python interface)
test-regression: build-deterministic-nopy
	$(CTEST_CMD) -L "regression"

# Build with deterministic floating-point for regression tests
# Keeps CODE intact so local dependency paths are respected
# Patches libneo to also use deterministic FP (no -ffast-math)
build-deterministic:
	@echo "Building with deterministic FP..."
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DSIMPLE_DETERMINISTIC_FP=ON -DCMAKE_COLOR_DIAGNOSTICS=ON $(FLAGS)
	@LIBNEO_CMAKE="$(BUILD_DIR)/_deps/libneo-src/CMakeLists.txt"; \
	if [ -f "$$LIBNEO_CMAKE" ] && grep -q "\-ffast-math" "$$LIBNEO_CMAKE" 2>/dev/null; then \
		echo "Patching libneo for deterministic FP..."; \
		sed 's/-ffast-math[[:space:]]*-ffp-contract=fast/-ffp-contract=off/g' "$$LIBNEO_CMAKE" > "$$LIBNEO_CMAKE.tmp" && mv "$$LIBNEO_CMAKE.tmp" "$$LIBNEO_CMAKE"; \
		sed 's/-ffast-math/-ffp-contract=off/g' "$$LIBNEO_CMAKE" > "$$LIBNEO_CMAKE.tmp" && mv "$$LIBNEO_CMAKE.tmp" "$$LIBNEO_CMAKE"; \
		echo "Reconfiguring after libneo patch..."; \
		cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DSIMPLE_DETERMINISTIC_FP=ON -DCMAKE_COLOR_DIAGNOSTICS=ON $(FLAGS); \
	fi
	cmake --build $(BUILD_DIR) --config $(CONFIG)

# Deterministic build matching golden record reference configuration (Python OFF)
build-deterministic-nopy:
	@echo "Building with deterministic FP and Python interface OFF (golden record reference)..."
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DSIMPLE_DETERMINISTIC_FP=ON -DENABLE_PYTHON_INTERFACE=OFF -DCMAKE_COLOR_DIAGNOSTICS=ON $(FLAGS)
	@LIBNEO_CMAKE="$(BUILD_DIR)/_deps/libneo-src/CMakeLists.txt"; \
	if [ -f "$$LIBNEO_CMAKE" ] && grep -q "\-ffast-math" "$$LIBNEO_CMAKE" 2>/dev/null; then \
		echo "Patching libneo for deterministic FP..."; \
		sed 's/-ffast-math[[:space:]]*-ffp-contract=fast/-ffp-contract=off/g' "$$LIBNEO_CMAKE" > "$$LIBNEO_CMAKE.tmp" && mv "$$LIBNEO_CMAKE.tmp" "$$LIBNEO_CMAKE"; \
		sed 's/-ffast-math/-ffp-contract=off/g' "$$LIBNEO_CMAKE" > "$$LIBNEO_CMAKE.tmp" && mv "$$LIBNEO_CMAKE.tmp" "$$LIBNEO_CMAKE"; \
		echo "Reconfiguring after libneo patch..."; \
		cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DSIMPLE_DETERMINISTIC_FP=ON -DENABLE_PYTHON_INTERFACE=OFF -DCMAKE_COLOR_DIAGNOSTICS=ON $(FLAGS); \
	fi
	cmake --build $(BUILD_DIR) --config $(CONFIG)

# Run all tests including regression tests
test-all: build
	$(CTEST_CMD)

# Golden record test target
# Compares current branch against main to enforce strict numerical reproducibility
# Any differences must be manually reviewed before merging
test-golden: build
	cd $(BUILD_DIR) && ctest --output-on-failure $(if $(filter 1,$(VERBOSE)),-V) -L "golden_record" $(if $(TEST),-R $(TEST))

doc: configure
	cmake --build --preset default --target doc

fpm:
	fpm build

clean:
	rm -rf $(BUILD_DIR)

# NVIDIA nvfortran build targets
# Uses NVIDIA HPC SDK with HPC-X MPI and proper CUDA setup
nvfortran-configure:
	NVHPC_CUDA_HOME=/opt/cuda \
	cmake -S . -B$(NVHPC_BUILD_DIR) -GNinja \
		-DCMAKE_BUILD_TYPE=$(CONFIG) \
		-DCMAKE_Fortran_COMPILER=nvfortran \
		-DCMAKE_C_COMPILER=nvc \
		-DMPI_HOME=$(NVHPC_HPCX) \
		-DMPI_C_COMPILER=$(NVHPC_HPCX)/bin/mpicc \
		-DMPI_Fortran_COMPILER=$(NVHPC_HPCX)/bin/mpifort \
		-DCMAKE_COLOR_DIAGNOSTICS=ON \
		$(FLAGS)

nvfortran: nvfortran-configure
	NVHPC_CUDA_HOME=/opt/cuda cmake --build $(NVHPC_BUILD_DIR) --config $(CONFIG)

nvfortran-clean:
	rm -rf $(NVHPC_BUILD_DIR)
