#!/usr/bin/env python3
"""
Test the improved BOOZXFORM field evaluation implementation.
This script tests that magfie_boozxform now reads real data and returns 
meaningful field values instead of dummy placeholders.
"""
import numpy as np
import sys
import os

# Add parent directory to path for pysimple
script_dir = os.path.dirname(os.path.abspath(__file__))
build_dir = os.path.join(script_dir, '../../build')
if os.path.exists(build_dir):
    sys.path.insert(0, build_dir)

try:
    import pysimple
except ImportError as e:
    print(f"Error: pysimple not found: {e}")
    sys.exit(1)


def test_boozxform_field_evaluation():
    """Test BOOZXFORM field evaluation using the magfie interface"""
    print("Testing improved BOOZXFORM field evaluation...")
    
    # Test coordinates (s=0.5, theta=0, zeta=0)
    x = np.array([0.5, 0.0, 0.0])
    
    # Temporary arrays for field evaluation
    bder = np.empty(3)
    hcovar = np.empty(3)
    hctrvr = np.empty(3)
    hcurl = np.empty(3)
    
    print(f"Evaluating field at coordinates: s={x[0]}, θ={x[1]}, ζ={x[2]}")
    
    try:
        # Call magfie_boozxform directly to bypass init_magfie issue
        modB, sqrtg = pysimple.magfie_sub.magfie_boozxform(x, bder, hcovar, hctrvr, hcurl)
        
        print(f"\nField evaluation results:")
        print(f"  |B| = {modB:.6f}")
        print(f"  √g = {sqrtg:.6e}")
        print(f"  ∂B/∂s = {bder[0]:.6e}")
        print(f"  ∂B/∂θ = {bder[1]:.6e}")
        print(f"  ∂B/∂ζ = {bder[2]:.6e}")
        print(f"  B_s = {hcovar[0]:.6e}")
        print(f"  B_θ = {hcovar[1]:.6e}")
        print(f"  B_ζ = {hcovar[2]:.6e}")
        print(f"  B^s = {hctrvr[0]:.6e}")
        print(f"  B^θ = {hctrvr[1]:.6e}")
        print(f"  B^ζ = {hctrvr[2]:.6e}")
        
        # Check if we're getting realistic values
        realistic = True
        issues = []
        
        if modB == 1.0:
            issues.append("Field strength still placeholder (exactly 1.0)")
            realistic = False
            
        if sqrtg == 1.0:
            issues.append("Jacobian still placeholder (exactly 1.0)")
            realistic = False
            
        if np.allclose(bder, 0.0):
            issues.append("All field derivatives are zero")
            
        if np.allclose(hcurl, 0.0):
            issues.append("All curl components are zero")
        
        if sqrtg <= 0:
            issues.append("Non-positive Jacobian")
            realistic = False
            
        if realistic:
            print(f"\n✅ SUCCESS: Field evaluation appears to use real BOOZXFORM data")
            print(f"   - Non-trivial field strength: {modB}")
            print(f"   - Non-trivial Jacobian: {sqrtg:.2e}")
            if not np.allclose(hcovar[1:], 0.0):
                print(f"   - Non-zero field components: B_θ={hcovar[1]:.3f}, B_ζ={hcovar[2]:.3f}")
        else:
            print(f"\n⚠️  PARTIAL: Field evaluation working but still has placeholder elements:")
            for issue in issues:
                print(f"   - {issue}")
                
        return modB, sqrtg, realistic
        
    except Exception as e:
        print(f"\n❌ ERROR: Field evaluation failed: {e}")
        return None, None, False


def test_coordinate_variations():
    """Test field evaluation at multiple coordinates"""
    print(f"\n" + "="*50)
    print("Testing field evaluation at multiple coordinates...")
    
    # Test at several points
    test_coords = [
        (0.1, 0.0, 0.0),    # Near axis
        (0.5, 0.0, 0.0),    # Mid-radius
        (0.9, 0.0, 0.0),    # Near edge
        (0.5, np.pi/2, 0.0), # Different theta
        (0.5, np.pi, 0.0),   # Different theta
    ]
    
    results = []
    
    for i, (s, theta, zeta) in enumerate(test_coords):
        x = np.array([s, theta, zeta])
        bder = np.empty(3)
        hcovar = np.empty(3)
        hctrvr = np.empty(3)
        hcurl = np.empty(3)
        
        try:
            modB, sqrtg = pysimple.magfie_sub.magfie_boozxform(x, bder, hcovar, hctrvr, hcurl)
            results.append((s, theta, zeta, modB, sqrtg, hcovar[1], hcovar[2]))
            print(f"  Point {i+1}: s={s:.1f}, θ={theta:.2f} → |B|={modB:.4f}, √g={sqrtg:.2e}")
            
        except Exception as e:
            print(f"  Point {i+1}: ERROR - {e}")
            results.append((s, theta, zeta, None, None, None, None))
    
    # Analyze results
    valid_results = [r for r in results if r[3] is not None]
    
    if len(valid_results) > 1:
        modB_values = [r[3] for r in valid_results]
        sqrtg_values = [r[4] for r in valid_results]
        
        modB_variation = np.std(modB_values) / np.mean(modB_values) if np.mean(modB_values) > 0 else 0
        sqrtg_variation = np.std(sqrtg_values) / np.mean(sqrtg_values) if np.mean(sqrtg_values) > 0 else 0
        
        print(f"\nVariation analysis:")
        print(f"  |B| coefficient of variation: {modB_variation:.1%}")
        print(f"  √g coefficient of variation: {sqrtg_variation:.1%}")
        
        if modB_variation > 0.01:  # >1% variation
            print(f"  ✅ Field strength shows spatial variation (good)")
        else:
            print(f"  ⚠️  Field strength shows little variation (may still be placeholder)")
            
        if sqrtg_variation > 0.01:  # >1% variation  
            print(f"  ✅ Jacobian shows spatial variation (good)")
        else:
            print(f"  ⚠️  Jacobian shows little variation (may still be placeholder)")
    
    return valid_results


def main():
    print("="*60)
    print("SIMPLE BOOZXFORM Field Evaluation Test")
    print("="*60)
    
    # Check that test data exists
    if not os.path.exists('boozmn_LandremanPaul2021_QA_lowres.nc'):
        print("❌ ERROR: BOOZXFORM test file not found")
        print("Run setup_booz_xform_test.sh first")
        return
    
    # Test single point evaluation
    modB, sqrtg, realistic = test_boozxform_field_evaluation()
    
    if modB is not None:
        # Test coordinate variations
        results = test_coordinate_variations()
        
        print(f"\n" + "="*60)
        print("SUMMARY")
        print("="*60)
        
        if realistic:
            print("✅ BOOZXFORM field evaluation is working with real data")
        else:
            print("⚠️  BOOZXFORM field evaluation partially working but needs improvement")
            
        print(f"Evaluated {len(results)} coordinate points successfully")
        
        if len(results) > 0:
            modB_range = [r[3] for r in results if r[3] is not None]
            if modB_range:
                print(f"|B| range: [{min(modB_range):.4f}, {max(modB_range):.4f}]")
        
    else:
        print("❌ BOOZXFORM field evaluation failed completely")
    
    print("="*60)


if __name__ == '__main__':
    main()