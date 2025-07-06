CONFIG ?= Release
BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all configure reconfigure build test install clean
all: build

$(BUILD_NINJA):
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

configure: $(BUILD_NINJA)

reconfigure:
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

build: configure
	cmake --build $(BUILD_DIR) --config $(CONFIG)

# Run tests with optional verbose output and single test selection
# Usage: make test [VERBOSE=1] [TEST=test_name] [INCLUDE_SLOW=1]
# Example: make test VERBOSE=1 TEST=test_gvec
# By default, tests labeled "slow" are excluded. To include them: make test INCLUDE_SLOW=1
test: build
	cd $(BUILD_DIR) && ctest --test-dir test --output-on-failure $(if $(VERBOSE),-V) $(if $(TEST),-R $(TEST)) $(if $(INCLUDE_SLOW),,-L "^((?!slow).)*$$")

# Run all tests including slow ones
test-all: build
	cd $(BUILD_DIR) && ctest --test-dir test --output-on-failure $(if $(VERBOSE),-V) $(if $(TEST),-R $(TEST))

doc: configure
	cmake --build --preset default --target doc

fpm:
	fpm build

clean:
	rm -rf $(BUILD_DIR)
