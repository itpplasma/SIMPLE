#!/usr/bin/env python3
"""Compare multiple numerical files for bit-identical results.

With deterministic floating-point builds (-ffp-contract=off, no -ffast-math),
golden record tests should produce bit-identical results. By default, this
script requires exact equality. Use --rtol and --atol flags to allow tolerances
if needed for specific cases.
"""
import sys
import os
import numpy as np
import argparse
import json

def compare_numerical_files(old_file, new_file, rtol=0.0, atol=0.0):
    """Compare two numerical files. Default is bit-identical (rtol=0, atol=0).

    Lines starting with '=' or '#' are treated as comments and skipped.
    """
    try:
        # Skip lines starting with '=' (section headers) or '#' (comments)
        old_data = np.loadtxt(old_file, comments=['#', '='])
        new_data = np.loadtxt(new_file, comments=['#', '='])
        
        # Handle scalar values
        if old_data.ndim == 0:
            old_data = old_data.reshape(1, 1)
        elif old_data.ndim == 1:
            old_data = old_data.reshape(-1, 1)
            
        if new_data.ndim == 0:
            new_data = new_data.reshape(1, 1)
        elif new_data.ndim == 1:
            new_data = new_data.reshape(-1, 1)
        
        # Check shape
        if old_data.shape != new_data.shape:
            return False, f"Shape mismatch: {old_data.shape} vs {new_data.shape}"
        
        # Compare values
        non_match_count = 0
        for old_row, new_row in zip(old_data, new_data):
            for old_val, new_val in zip(old_row, new_row):
                if not np.isclose(old_val, new_val, rtol=rtol, atol=atol):
                    non_match_count += 1
        
        if non_match_count > 0:
            return False, f"Found {non_match_count} non-matching entries"
        
        return True, "Files match"
        
    except Exception as e:
        return False, f"Error comparing files: {str(e)}"

def compare_binary_files(old_file, new_file):
    """Compare two binary files byte by byte"""
    try:
        with open(old_file, 'rb') as f1, open(new_file, 'rb') as f2:
            return f1.read() == f2.read(), "Binary comparison"
    except Exception as e:
        return False, f"Error comparing binary files: {str(e)}"

def get_file_comparison_method(filename):
    """Determine comparison method based on filename"""
    # Text/numerical files that should be compared numerically
    numerical_extensions = ['.dat', '.txt']
    numerical_prefixes = ['fort.', 'times_lost', 'confined_fraction', 'class_parts', 
                         'avg_inverse_t_lost', 'healaxis', 'start']
    
    # Binary files
    binary_extensions = ['.nc']
    
    # Check if it's a binary file
    for ext in binary_extensions:
        if filename.endswith(ext):
            return 'binary'
    
    # Check if it's a numerical file
    for ext in numerical_extensions:
        if filename.endswith(ext):
            return 'numerical'
    
    for prefix in numerical_prefixes:
        if filename.startswith(prefix):
            return 'numerical'
    
    # Default to binary for unknown files
    return 'binary'

def compare_file_lists(ref_dir, cur_dir, file_list, skip_files=None, rtol=0.0, atol=0.0):
    """Compare a list of files between two directories.

    Default is bit-identical comparison (rtol=0, atol=0) for deterministic FP builds.
    """
    if skip_files is None:
        skip_files = ['simple.in', 'wout.nc']  # Default files to skip

    results = {}
    summary = {
        'total': 0,
        'passed': 0,
        'failed': 0,
        'missing': 0,
        'skipped': 0
    }

    for filename in file_list:
        if filename in skip_files:
            results[filename] = {'status': 'skipped', 'message': 'File in skip list'}
            summary['skipped'] += 1
            continue

        summary['total'] += 1

        ref_file = os.path.join(ref_dir, filename)
        cur_file = os.path.join(cur_dir, filename)

        # Check if files exist
        if not os.path.exists(ref_file):
            results[filename] = {'status': 'missing', 'message': f'Reference file missing: {ref_file}'}
            summary['missing'] += 1
            continue

        if not os.path.exists(cur_file):
            results[filename] = {'status': 'missing', 'message': f'Current file missing: {cur_file}'}
            summary['missing'] += 1
            continue

        # Determine comparison method
        method = get_file_comparison_method(filename)

        # Compare files
        if method == 'numerical':
            match, message = compare_numerical_files(ref_file, cur_file, rtol=rtol, atol=atol)
        else:
            match, message = compare_binary_files(ref_file, cur_file)
        
        if match:
            results[filename] = {'status': 'passed', 'message': message}
            summary['passed'] += 1
        else:
            results[filename] = {'status': 'failed', 'message': message}
            summary['failed'] += 1
    
    return results, summary

def main():
    parser = argparse.ArgumentParser(description='Compare multiple files between directories')
    parser.add_argument('ref_dir', help='Reference directory')
    parser.add_argument('cur_dir', help='Current directory to compare')
    parser.add_argument('--files', nargs='+', help='List of files to compare')
    parser.add_argument('--skip', nargs='+', default=['simple.in', 'wout.nc'], 
                       help='Files to skip (default: simple.in wout.nc)')
    parser.add_argument('--json', action='store_true', help='Output results as JSON')
    parser.add_argument('--rtol', type=float, default=0.0,
                       help='Relative tolerance (default: 0.0 for bit-identical)')
    parser.add_argument('--atol', type=float, default=0.0,
                       help='Absolute tolerance (default: 0.0 for bit-identical)')
    
    args = parser.parse_args()
    
    # If no files specified, use all files in current directory
    if args.files is None:
        if os.path.exists(args.cur_dir):
            args.files = [f for f in os.listdir(args.cur_dir) if os.path.isfile(os.path.join(args.cur_dir, f))]
        else:
            print(f"Error: Current directory not found: {args.cur_dir}")
            sys.exit(1)
    
    # Compare files (default: bit-identical for deterministic FP builds)
    results, summary = compare_file_lists(
        args.ref_dir, args.cur_dir, args.files, args.skip,
        rtol=args.rtol, atol=args.atol
    )
    
    # Output results
    if args.json:
        output = {
            'results': results,
            'summary': summary
        }
        print(json.dumps(output, indent=2))
    else:
        # Human-readable output
        print(f"Comparing files between:")
        print(f"  Reference: {args.ref_dir}")
        print(f"  Current: {args.cur_dir}")
        print()
        
        # Show individual results
        for filename, result in sorted(results.items()):
            status = result['status']
            message = result['message']
            
            if status == 'passed':
                symbol = '✓'
            elif status == 'failed':
                symbol = '✗'
            elif status == 'missing':
                symbol = '⚠'
            else:  # skipped
                symbol = '-'
            
            print(f"{symbol} {filename}: {message}")
        
        print()
        print("=" * 50)
        print("Summary:")
        print(f"  Total files checked: {summary['total']}")
        print(f"  Passed: {summary['passed']}")
        print(f"  Failed: {summary['failed']}")
        print(f"  Missing: {summary['missing']}")
        print(f"  Skipped: {summary['skipped']}")
        print("=" * 50)
    
    # Exit with appropriate code
    if summary['failed'] > 0 or summary['missing'] > 0:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == '__main__':
    main()