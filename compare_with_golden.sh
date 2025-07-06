#!/bin/bash
# Script to compare current implementation with golden reference

echo "=== Refactoring Comparison Tool ==="
echo

# Check if golden reference tag exists
if ! git tag | grep -q "golden-reference-v1"; then
    echo "ERROR: Golden reference tag 'golden-reference-v1' not found"
    echo "Please create the golden reference first"
    exit 1
fi

echo "1. Capturing current output..."
./build/simple.x | grep "DEBUG:" > current_debug.txt 2>&1

echo "2. Checking out golden reference..."
git stash push -m "Temporary stash for comparison"
git checkout golden-reference-v1 -q

echo "3. Building golden reference..."
make clean > /dev/null 2>&1
make > /dev/null 2>&1

echo "4. Capturing golden reference output..."
./build/simple.x | grep "DEBUG:" > golden_debug.txt 2>&1

echo "5. Comparing outputs..."
if diff -q current_debug.txt golden_debug.txt > /dev/null; then
    echo "✅ SUCCESS: Outputs are identical!"
    echo "   Refactoring preserved exact behavior"
else
    echo "❌ DIFFERENCE DETECTED:"
    echo "   Running detailed comparison..."
    echo
    echo "=== DIFF OUTPUT ==="
    diff current_debug.txt golden_debug.txt
    echo "==================="
fi

echo "6. Restoring current branch..."
git checkout - -q
git stash pop -q > /dev/null 2>&1

echo
echo "Files created for manual inspection:"
echo "- current_debug.txt  (your current implementation)"
echo "- golden_debug.txt   (golden reference output)"
echo
echo "Done."