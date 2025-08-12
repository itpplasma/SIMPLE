CONFIG ?= Release
FLAGS ?=
BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

# Common ctest command with optional verbose and test name filtering
CTEST_CMD = cd $(BUILD_DIR) && ctest --test-dir test --output-on-failure $(if $(filter 1,$(VERBOSE)),-V) $(if $(TEST),-R $(TEST))

.PHONY: all configure reconfigure build test test-fast test-slow test-regression test-all install clean coverage coverage-build coverage-clean
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

# Run all tests except regression tests (default)
test: build
	$(CTEST_CMD) -LE "regression"

# Run only fast tests (exclude slow and regression tests)
test-fast: build
	$(CTEST_CMD) -LE "slow|regression"

# Run only slow tests
test-slow: build
	$(CTEST_CMD) -L "slow" -LE "regression"

# Run only regression tests
test-regression: build
	$(CTEST_CMD) -L "regression"

# Run all tests including regression tests
test-all: build
	$(CTEST_CMD)

doc: configure
	cmake --build --preset default --target doc

fpm:
	fpm build

# Coverage targets
# Build with coverage instrumentation using Profile build type  
coverage-build:
	$(MAKE) CONFIG=Profile FLAGS="-DENABLE_COVERAGE=ON"

# Generate coverage report after running tests
coverage: coverage-build
	$(MAKE) CONFIG=Profile test-all
	cmake --build $(BUILD_DIR) --target coverage
	@echo "✅ Coverage data generated in $(BUILD_DIR)/coverage_filtered.info"

# Clean coverage data
coverage-clean:
	find $(BUILD_DIR) -name "*.gcda" -delete
	find $(BUILD_DIR) -name "*.gcno" -delete
	rm -f $(BUILD_DIR)/coverage*.info

clean:
	rm -rf $(BUILD_DIR)
