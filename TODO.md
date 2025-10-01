# TODO: Meiss Canonical Coordinates on EQDSK/Geoflux — Full Implementation & Robust Testing

_Last updated: 2025-09-30_

This is an exhaustive execution guide. Every step includes file names, code snippets, commands, and validation criteria. The emphasis is on **precise, fast, non-tautological tests** and **system-level regressions**.

---
## 0. Baseline Safety Net (Tests Only)
1. **Ensure clean trees** (SIMPLE + libneo):
   ```bash
   git status
   cd ../libneo && git status
   cd ../SIMPLE
   ```

2. **VMEC Meiss guard** (`test/tests/test_field_can_meiss_vmec.f90`)
   - Implement as described previously **but add analytic invariants**:
     - After `evaluate_meiss`, compute energy `H = 0.5*vpar^2 + mu*Bmod` and ensure matches reference value from prior run (hard-code expected value with tolerance 1e-10).
     - Validate orthogonality: `dot_product(hcov, tracer%f%dhth)` should be near zero (tolerance 1e-12).
   - Register in CMake and run: `ctest -R test_field_can_meiss_vmec`.

3. **EQDSK guard (WILL_FAIL initially)** (`test/tests/test_field_can_meiss_eqdsk.f90`)
   - Implement same checks as VMEC guard; currently expect failure.
   - Mark test `WILL_FAIL TRUE` in CMake.
   - Run: `ctest -R test_field_can_meiss_eqdsk` and confirm it fails but suite continues.

4. **Geoflux field smoke** (`ctest -R test_magfie_geoflux`) – must pass.

5. **Golden record awareness**
   - Inspect `test/golden_record/canonical/simple.in` – note it uses `isw_field_type=0` (canonical flux).
   - Plan to add new golden case for Meiss EQDSK later (see §4).

---
## 1. libneo Coordinate System (follow ../libneo/TODO.md)
Complete libneo work first. Do not proceed until `../libneo/build/ctest -R test_coordinate_systems` passes and exposes `make_vmec_coordinate_system` / `make_geoflux_coordinate_system`.

---
## 2. Integrate coordinate_system_t (no wrappers)

1. **Imports**
   - Add `use libneo_coordinates` wherever geometry is needed: `src/field/field_can_meiss.f90`, `src/field_can.f90`, `app/simple_diag_meiss.f90`, `src/diag/diag_meiss.f90`, and any helper modules use the geometry.

2. **Global storage** (`field_can_meiss`):
   - Replace `field_noncan` with:
     ```fortran
     class(coordinate_system_t), allocatable :: meiss_coords
     class(MagneticField), allocatable :: source_field
     ```

3. **Initialise coordinate system** (`init_meiss`)
   - After storing incoming magnetic field, allocate coordinate system directly:
     ```fortran
     if (allocated(meiss_coords)) deallocate(meiss_coords)
     if (geoflux_ready()) then
       call make_geoflux_coordinate_system(meiss_coords)
     else
       call make_vmec_coordinate_system(meiss_coords)
     end if
     ```
   - Note: `geoflux_ready()` should be a new public function from `field_geoflux`. Add it if missing (e.g., module procedure returning logical flag).

4. **Use coordinate system** across the module:
   - Replace all VMEC-specific data access (e.g., global arrays `lam_phi`, `h_cov`, etc.) with method calls:
     - `meiss_coords%evaluate_point(u, position)`.
     - `meiss_coords%covariant_basis(u, basis)` to compute derivatives.
     - `meiss_coords%metric_tensor(u, g, ginv, sqrtg)`.
   - Compute derivatives/invariants using the returned basis/metric—*no* direct `splint_vmec_data` or `geoflux_to_cyl` inside SIMPLE.

5. **Magnetic field sampling**
   - Continue using the existing `MagneticField` polymorphic object (`source_field`) for `evaluate` calls. Geometry now depends solely on `meiss_coords`.

6. **Entry points update** (`src/field_can.f90`)
   - Ensure `init_field_can` passes the magnetic field so `init_meiss` can allocate the right coordinate system.
   - No new SIMPLE types are created; everything uses libneo’s `coordinate_system_t` instance.

7. **Diagnostics**
   - In `app/simple_diag_meiss.f90`, after reading config, call `field_from_file` and then `init_meiss(field)`; remove direct `init_vmec` + manual cache calls.
   - `src/diag/diag_meiss.f90` simply relies on `rh_can`; confirm that `rh_can` internally uses `meiss_coords`.

8. **Upgrade EQDSK test**
   - Remove `WILL_FAIL TRUE` from EQDSK test once refactor is complete.
   - Enhance checks in both VMEC and EQDSK tests to compare outputs against reference values stored in a small JSON/Fortran data file (e.g., `test/reference/meiss_expected.dat`). Populate this file by running the old VMEC pipeline and recording `Bmod`, `hth`, `hph` for a few canonical points; these become the expected values with tolerance.

9. **Run focused unit tests**
   ```bash
   cmake --build build --target test_field_can_meiss_vmec.x test_field_can_meiss_eqdsk.x
   ctest -R test_field_can_meiss_vmec
   ctest -R test_field_can_meiss_eqdsk
   ctest -R test_magfie_geoflux
   ```
   All must pass (<1s each).

---
## 3. Regression Tests & Golden Record
1. **Add Meiss EQDSK golden test**
   - Create `test/golden_record/meiss_eqdsk/simple.in` with EQDSK-specific configuration (e.g., `netcdffile='EQDSK_I.geqdsk'`, `isw_field_type=3` to request Meiss, set deterministic seeds).
   - Update `test/golden_record/run_golden_tests.sh` so it copies the EQDSK file into each test case directory (download once to `$TEST_DATA_DIR/EQDSK_I.geqdsk`).
   - Ensure symlink/hard copy placed in each run directory (similar to `wout.nc`).
   - Add validation at the end of each run: check `diag_meiss_*` outputs exist if expected.

2. **Golden record reference update**
   - Generate reference outputs by running `test/golden_record/golden_record.sh main` (or chosen reference commit) after new code is in place. Keep `RUN_DIR_REF` results for diff.
   - Ensure `compare_golden_results.sh` accounts for new files (update script to include Meiss-specific outputs).

3. **System regression command**
   - After every change, run:
     ```bash
     ./test/golden_record/golden_record.sh main
     ```
     and ensure the current run matches reference (no diffs aside from known improvements).

---
## 4. Diagnostics Automation
1. **CTest driver** `test/diag/run_diag_meiss.cmake` (see prior template). Extend to check *both* VMEC and EQDSK outputs and to verify key numeric summaries (e.g., parse `diag_meiss_*.pdf` metadata using a small Python script or check associated data files if generated).
2. **Register test** `diag_meiss_runs` with `TIMEOUT 120`.
3. **Run** `ctest -R diag_meiss_runs`.

---
## 5. Examples & Documentation
1. **Makefile** `examples/tokamak/Makefile` – ensure targets `all`, `run`, `diag`, `clean` as previously outlined.
2. **README** `examples/tokamak/README.md` must include:
   - Download step (`make`).
   - Execution (`make run`), expected runtime, note about OpenMP threads.
   - Diagnostics (`make diag`), location of generated plots.
   - Troubleshooting tips (e.g., missing libneo branch).
3. **Manual run**: follow README to completion; confirm outputs align with documentation.

---
## 6. Final Regression Sweep
1. **Unit tests**: `ctest --output-on-failure` in SIMPLE repo.
2. **System tests**: `./test/golden_record/golden_record.sh main` (ensure 0 exit code).
3. **libneo tests**: `ninja -C ../libneo/build && ctest`.
4. **Clean status**: `git status` (SIMPLE + libneo).
5. **Document** results (changelog, PR summary).

Completion criteria: unit & system tests all green, new golden record case passes, diagnostics produce EQDSK Meiss plots, documentation updated.
