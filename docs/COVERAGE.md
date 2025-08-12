# Code Coverage Integration

This document describes the code coverage analysis setup for SIMPLE, integrated similarly to the Fortfront project.

## Overview

SIMPLE now includes automated code coverage analysis integrated into the CMake build system and GitHub Actions workflow. Coverage analysis helps identify which parts of the codebase are exercised by tests and which may need additional testing.

## Local Usage

### Basic Coverage Analysis

To generate coverage reports locally:

```bash
# Build with coverage instrumentation and run all tests
make coverage

# Or step by step:
make coverage-build      # Build with coverage flags
make CONFIG=Profile test-all  # Run tests with coverage collection
cmake --build build --target coverage  # Generate coverage report
```

### Cleaning Coverage Data

```bash
make coverage-clean      # Remove generated coverage files
```

### Manual Coverage Commands

```bash
# Build with coverage enabled
make CONFIG=Profile FLAGS="-DENABLE_COVERAGE=ON"

# Run tests to generate coverage data
make CONFIG=Profile test-fast

# Generate LCOV report
cmake --build build --target coverage
```

## CI/CD Integration

Coverage analysis is automatically performed on pull requests through GitHub Actions:

### What Gets Analyzed

- **Project Coverage**: Overall line coverage percentage across all source files
- **Patch Coverage**: Coverage of lines changed in the pull request (simplified implementation)
- **Branch Coverage**: Enabled via LCOV's `branch_coverage=1` flag

### Coverage Workflow

1. **Build with Coverage**: Uses Profile build type with coverage instrumentation
2. **Run Tests**: Executes all test suites to collect coverage data
3. **Generate Reports**: Creates LCOV and Cobertura XML reports
4. **Coverage Checks**: Creates GitHub status checks for project/patch coverage
5. **PR Comments**: Posts detailed coverage report as pull request comment

### Coverage Thresholds

- **Project Coverage**: 70% minimum threshold
- **Patch Coverage**: 70% minimum threshold (simplified - uses project coverage)

## Technical Details

### CMake Configuration

Coverage analysis is controlled by the `ENABLE_COVERAGE` CMake option:

- **Available Build Types**: Debug, Profile (coverage requires one of these)
- **Compiler Flags**: `--coverage -fprofile-arcs -ftest-coverage`
- **Linker Flags**: `--coverage -lgcov`

### File Exclusions

Coverage analysis excludes:
- Third-party dependencies (`build/dependencies/*`)
- Test files (`test/*`)
- System headers (`/usr/*`)

### Generated Files

Coverage analysis generates several files (automatically git-ignored):
- `*.gcda`, `*.gcno` - GCC coverage data files
- `coverage.info` - Raw LCOV coverage data
- `coverage_filtered.info` - Filtered LCOV data
- `cobertura.xml` - Cobertura XML format for GitHub Actions
- `coverage_html/` - HTML coverage reports (when generated)

### Dependencies

The system uses:
- **LCOV** - Coverage data collection and processing
- **gcov** - GCC's coverage analysis tool
- **lcov_cobertura** - Python package for XML conversion
- **coverage-action** - GitHub Action for PR reporting

## Makefile Targets

| Target | Description |
|--------|-------------|
| `coverage-build` | Build with coverage instrumentation |
| `coverage` | Full coverage analysis (build + test + report) |
| `coverage-clean` | Clean coverage data files |

## Troubleshooting

### Coverage Data Not Generated

Ensure you're using a supported build type:
```bash
make CONFIG=Profile FLAGS="-DENABLE_COVERAGE=ON"
```

### Missing LCOV

Install LCOV if not available:
```bash
# Ubuntu/Debian
sudo apt-get install lcov

# macOS
brew install lcov
```

### Low Coverage Numbers

Coverage analysis only measures code executed during test runs. To improve coverage:
1. Add more comprehensive tests
2. Ensure tests exercise edge cases
3. Check for unreachable or dead code

## Integration with fortfront

This implementation mirrors the coverage setup used in the fortfront project, including:
- Similar LCOV configuration and filtering
- Cobertura XML generation for GitHub Actions
- Coverage threshold enforcement
- Automated PR commenting
- Same directory exclusion patterns