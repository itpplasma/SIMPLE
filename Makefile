CONFIG ?= Release
BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all ninja test install clean
all: ninja

$(BUILD_NINJA):
	cmake -S . -Bbuild -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

ninja: $(BUILD_NINJA)
	cmake --build build --config $(CONFIG)

test: ninja
	cd $(BUILD_DIR) && ctest

doc: $(BUILD_NINJA)
	cmake --build --preset default --target doc

clean:
	rm -rf $(BUILD_DIR)
