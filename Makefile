BUILD_DIR := build
CMAKE_CACHE := $(BUILD_DIR)/CMakeCache.txt
NUM_PROC := $(shell python -c 'import os; print(os.cpu_count())')

.PHONY: all compile install test fasttest clean
all: compile

install: compile
	cd $(BUILD_DIR) && make install

test: compile
	ctest --output-on-failure --test-dir $(BUILD_DIR)

compile: $(CMAKE_CACHE)
	cmake --build $(BUILD_DIR) -j$(NUM_PROC)

$(CMAKE_CACHE):
	cmake -B $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
