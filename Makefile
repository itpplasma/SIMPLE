CONFIG ?= Release
BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all configure reconfigure build test install clean
all: build

$(BUILD_NINJA):
	cmake -S . -Bbuild -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

configure: $(BUILD_NINJA)

reconfigure:
	cmake -S . -Bbuild -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

build: configure
	cmake --build build --config $(CONFIG)

# Run tests with optional verbose output and single test selection
# Usage: make test [VERBOSE=1] [TEST=test_name]
# Example: make test VERBOSE=1 TEST=test_gvec
test: build
	cd $(BUILD_DIR) && ctest --test-dir test --output-on-failure $(if $(VERBOSE),-V) $(if $(TEST),-R $(TEST))

doc: configure
	cmake --build --preset default --target doc

fpm:
	fpm build

clean:
	rm -rf $(BUILD_DIR)
