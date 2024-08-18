BUILD_DIR := build
CMAKE_CACHE := $(BUILD_DIR)/CMakeCache.txt
NUM_PROC := $(shell python -c 'import os; print(os.cpu_count())')

.PHONY: all compile install test clean
all: compile

install: compile
	cd $(BUILD_DIR) && make install

test: compile
	cd build/test/tests && ./utility_test.x

compile: $(CMAKE_CACHE)
	cd $(BUILD_DIR) && make -j$(NUM_PROC)

$(CMAKE_CACHE):
	cmake -B $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
