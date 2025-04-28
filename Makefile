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

test: build
	cd $(BUILD_DIR) && ctest

doc: configure
	cmake --build --preset default --target doc

fpm:
	fpm build

clean:
	rm -rf $(BUILD_DIR)
