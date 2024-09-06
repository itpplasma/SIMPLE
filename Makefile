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

$(CMAKE_CACHE): libneo
	cmake -B $(BUILD_DIR)

# Clone libneo if it doesn't exist, or update it if it does
LIBNEO_PATH := SRC/libneo
LIBNEO_REPO := https://github.com/itpplasma/libneo.git
libneo:
	@if [ ! -d $(LIBNEO_PATH) ]; then \
		echo "Cloning libneo repository..."; \
		git clone $(LIBNEO_REPO) $(LIBNEO_PATH); \
	elif [ -L $(LIBNEO_PATH) ]; then \
		echo "libneo is a symlink. Skipping clone..."; \
	elif [ -d $(LIBNEO_PATH)/.git ]; then \
		echo "libneo exists, pulling latest changes..."; \
		cd $(LIBNEO_PATH) && git pull origin main; \
	else \
		echo "libneo already exists as a non-git directory. Skipping update..."; \
	fi

clean:
	rm -rf $(BUILD_DIR)
