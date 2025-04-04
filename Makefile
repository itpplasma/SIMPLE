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

test: build
	cd $(BUILD_DIR) && ctest

doc: configure
	cmake --build --preset default --target doc

fpm:
	fpm build --flag "-I$HDF5_INCLUDE -I$NETCDF_FORTRAN_INCLUDE -I$OPENBLAS_INCLUDE -L$OPENBLAS_LIB"

clean:
	rm -rf $(BUILD_DIR)
