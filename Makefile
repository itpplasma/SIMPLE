CONFIG ?= Release
BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all configure reconfigure build test test-fast test-regression install clean
all: build

$(BUILD_NINJA):
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

configure: $(BUILD_NINJA)

reconfigure:
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

build: configure
	cmake --build $(BUILD_DIR) --config $(CONFIG)

# Run tests including slow ones but excluding regression tests
# Usage: make test [TEST=test_name]
# Example: make test TEST=test_gvec
# Note: Verbose output is now default. Use VERBOSE=0 to disable.
test: build
	cd $(BUILD_DIR) && ctest --test-dir test --output-on-failure $(if $(filter 0,$(VERBOSE)),,-V) $(if $(TEST),-R $(TEST)) -LE "regression"

# Run only fast tests (exclude slow and regression tests)
test-fast: build
	cd $(BUILD_DIR) && ctest --test-dir test --output-on-failure $(if $(filter 0,$(VERBOSE)),,-V) $(if $(TEST),-R $(TEST)) -LE "slow|regression"

# Run all tests including regression tests
test-regression: build
	cd $(BUILD_DIR) && ctest --test-dir test --output-on-failure $(if $(filter 0,$(VERBOSE)),,-V) $(if $(TEST),-R $(TEST))

doc: configure
	cmake --build --preset default --target doc

fpm:
	fpm build

clean:
	rm -rf $(BUILD_DIR)
