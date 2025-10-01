# Transition libneo integration away from `find_or_fetch`

Goal: replace the in-tree `find_or_fetch(libneo)` call with an isolated build/install flow that keeps the `$CODE` override, prevents dependency tests from polluting SIMPLE's CTest, and reliably tracks upstream branch updates.

## Pre-flight checks
- [ ] Verify libneo currently exports modules into `${CMAKE_Fortran_MODULE_DIRECTORY}` and confirm location after `cmake --install`.
- [ ] Capture minimum libneo revision we depend on (branch or tag) and record in documentation.
- [ ] Ensure all developers know to keep `$CODE/libneo` up to date if they rely on the local override.

## Implementation steps
1. **Introduce dependency module**
   - Create `cmake/deps/libneo.cmake` (new directory if needed).
   - Detect source path:
     - If `$ENV{CODE}/libneo` exists and contains a `.git` directory, use it directly.
     - Otherwise prepare `${CMAKE_BINARY_DIR}/_deps/libneo-src` for fetched sources.
   - Expose cache entries (`LIBNEO_SOURCE`, `LIBNEO_BRANCH`) so users can override without editing CMake.
2. **Fetch sources when override absent**
   - Use `ExternalProject_Add` (preferred) or `FetchContent`+manual clone *without* `add_subdirectory`.
   - Pass `GIT_REPOSITORY` (https://github.com/itpplasma/libneo.git) and `GIT_TAG` (default `main`, override via cache/env).
   - Keep `UPDATE_DISCONNECTED` default so updates pull new commits.
   - Enable shallow clones for speed (`GIT_SHALLOW TRUE`) and document the behaviour.
3. **Configure isolated build**
   - Set `BINARY_DIR = ${CMAKE_BINARY_DIR}/_deps/libneo-build` and `INSTALL_DIR = ${CMAKE_BINARY_DIR}/_deps/libneo-stage`.
   - Supply `CMAKE_ARGS`:
     - `-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}`
     - `-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}`
     - Disable options we do not need (`-DENABLE_PYTHON_INTERFACE=OFF`, `-DENABLE_OPENMP=ON|OFF` follow SIMPLE’s choice).
     - Suppress dependency tests (`-DBUILD_TESTING=OFF` and, once available, `-DLIBNEO_BUILD_TESTS=OFF`).
   - Ensure environment mirrors `$ENV{CODE}` usage (export `CODE` when building staged copy if libneo scripts expect it).
4. **Expose imported targets**
   - After defining the external project, call `ExternalProject_Get_Property` to grab `BINARY_DIR` / `INSTALL_DIR`.
   - Define imported libraries, e.g.:
     ```cmake
     add_library(libneo::neo STATIC IMPORTED)
     set_target_properties(libneo::neo PROPERTIES
         IMPORTED_LOCATION "${INSTALL_DIR}/lib/libneo.a"
         INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include"
     )
     add_dependencies(libneo::neo libneo_external)
     ```
   - If install step does not copy `.mod` files, fall back to `${BINARY_DIR}/include` and add TODO to upstream.
   - Update all SIMPLE targets to link against `libneo::neo` instead of directory-relative paths.
5. **Remove legacy helper**
   - Delete `find_or_fetch` from `cmake/Util.cmake`.
   - Replace call sites (search for `find_or_fetch(libneo)` and import new module).
   - Guard the transition with CMake policy checks so consumers get a clear error if both systems are invoked.
6. **CTest hygiene**
   - Confirm `ctest -N` only lists SIMPLE tests after change.
   - Add regression guard by introducing a CTest label (e.g. `third_party`) if any dependency tests accidentally slip through.
7. **Documentation & Developer ergonomics**
   - Update README/CONTRIBUTING with new workflow (`cmake -DLIBNEO_SOURCE=...`).
   - Note that `$CODE/libneo` must be built automatically; no manual `cmake --build` calls required.
   - Add troubleshooting section: rebuild external project (`cmake --build build --target libneo_external`), cleaning staged dirs, handling branch mismatch.

## Testing checklist
- [ ] Configure from clean tree (`rm -rf build && cmake -S . -B build`).
- [ ] Build libneo via external project (`cmake --build build --target libneo_external`).
- [ ] Build SIMPLE (`cmake --build build`).
- [ ] Run focused tests (`ctest --test-dir build -R '^golden_record$'`).
- [ ] Repeat with `$CODE/libneo` override pointing at a different commit.
- [ ] Verify updating libneo branch (`git -C build/_deps/libneo-src pull`) triggers rebuild next run.

## Follow-up tasks
- [ ] Upstream feature request: add `LIBNEO_BUILD_TESTS` option that gates `enable_testing()`.
- [ ] Investigate providing a packaged libneo (fpm or CMake config) to simplify long-term consumption.
- [ ] Evaluate caching/stashing staged builds in CI to cut minutes off clean jobs.
