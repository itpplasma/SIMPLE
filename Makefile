CONFIG ?= Release
BUILD_DIR := build
BUILD_MAKE := $(BUILD_DIR)/Makefile

.PHONY: all configure build clean
all: build

$(BUILD_MAKE):
	cmake -S . -Bbuild -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

configure: $(BUILD_MAKE)

build: configure
	cmake --build build --config $(CONFIG) -j$(shell nproc)

clean:
	rm -rf $(BUILD_DIR)
