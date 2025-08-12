#!/bin/bash

# Function to modernize type declarations in a Fortran module
modernize_module() {
    local module_file="$1"
    local issue_num="$2"
    local module_name=$(basename "$module_file" .f90)
    module_name=$(basename "$module_name" .F90)
    
    echo "Processing $module_file (Issue #$issue_num)..."
    
    # Create branch
    git checkout main
    git pull
    git checkout -b "type-modernization-${module_name}-${issue_num}"
    
    # Check if module has double precision declarations
    if ! grep -q "double precision" "$module_file"; then
        echo "  No double precision declarations found in $module_file"
        git checkout main
        return 0
    fi
    
    # Add the dp parameter and modernize types
    # First check if the module already uses iso_fortran_env
    if grep -q "iso_fortran_env" "$module_file"; then
        # Module already uses iso_fortran_env, need to be careful
        echo "  Module uses iso_fortran_env, using dp_kind alias"
        # Add dp_kind => real64 to existing use statement or add new one
        sed -i '1,/implicit none/{
            /use.*iso_fortran_env/! {
                /implicit none/i\
use, intrinsic :: iso_fortran_env, only: dp_kind => real64
            }
        }' "$module_file"
        # Replace double precision with real(dp_kind)
        sed -i 's/double precision/real(dp_kind)/g' "$module_file"
    else
        # Module doesn't use iso_fortran_env, we can add dp parameter safely
        echo "  Adding dp parameter and modernizing types"
        # Add dp parameter after implicit none
        sed -i '/implicit none/a\\ninteger, parameter :: dp = kind(1.0d0)' "$module_file"
        # Replace double precision with real(dp)
        sed -i 's/double precision/real(dp)/g' "$module_file"
    fi
    
    # Test build
    if make test-fast > /dev/null 2>&1; then
        echo "  Build successful!"
        
        # Commit changes
        git add "$module_file"
        git commit -m "refactor: Modernize type declarations in $module_name (#$issue_num Part 2)"
        
        # Push branch
        git push -u origin "type-modernization-${module_name}-${issue_num}"
        
        # Create PR
        gh pr create \
            --title "refactor: Modernize type declarations in $module_name (#$issue_num Part 2)" \
            --body "Part 2 of refactoring for issue #$issue_num

## Changes
- Replace all \`double precision\` with \`real(dp)\` or \`real(dp_kind)\`
- Add appropriate parameter definition for double precision kind
- Modernize type declarations to follow current Fortran standards
- No functional changes, purely type modernization

Part of the systematic refactoring effort."
        
        echo "  PR created successfully!"
    else
        echo "  Build failed, skipping..."
        git checkout main
        git branch -D "type-modernization-${module_name}-${issue_num}"
    fi
    
    return 0
}

# Process all modules
modernize_module "src/orbit_symplectic.f90" "99"
modernize_module "src/get_canonical_coordinates.F90" "100"
modernize_module "src/field_can.f90" "101"
modernize_module "src/collis_alphas.f90" "102"
modernize_module "src/sub_alpha_lifetime_can.f90" "103"
modernize_module "src/orbit_symplectic_quasi.f90" "104"
modernize_module "src/params.f90" "105"
modernize_module "src/simple_main.f90" "106"
modernize_module "src/samplers.f90" "107"
modernize_module "src/find_bminmax.f90" "108"
modernize_module "src/classification.f90" "109"
modernize_module "src/orbit_symplectic_base.f90" "110"
modernize_module "src/simple.f90" "111"
modernize_module "src/cut_detector.f90" "112"

echo "All modules processed!"
git checkout main