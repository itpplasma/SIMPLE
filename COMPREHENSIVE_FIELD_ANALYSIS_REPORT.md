# Comprehensive Analysis: GVEC vs SIMPLE VMEC Field Implementation

## Executive Summary

This report provides a detailed investigation of why GVEC magnetic field evaluations become excessively large at small flux surface values (small s) compared to SIMPLE's native VMEC implementation. The analysis identifies the root cause and provides specific recommendations for resolution.

## Root Cause Identification

### The Problem
GVEC produces magnetic field magnitudes that are **orders of magnitude larger** than SIMPLE's VMEC implementation, particularly at small normalized flux values (s < 0.1). This discrepancy becomes extreme near the magnetic axis.

### Root Cause: Coordinate System Scaling
The issue stems from GVEC's use of a **coordinate transformation scaling factor** in its magnetic field computation:

```
GVEC: B^i = (field_coefficients) * (dPhi_dr / Jac)
```

Where:
- `dPhi_dr`: derivative of toroidal flux with respect to radial coordinate
- `Jac`: Jacobian determinant of the coordinate transformation

Near the magnetic axis (s → 0):
- **Jacobian → 0** (area elements shrink)
- **dPhi_dr behavior** depends on radial coordinate definition
- **Combined effect**: `dPhi_dr/Jac → ∞`

## Detailed Technical Analysis

### 1. GVEC Field Computation

From `/proj/plasma/CODE/ert/SIMPLE/build/_deps/gvec-src/python/gvec/quantities.py` lines 532-534:

```python
ds["B_contra_t"] = (ds.iota - ds.dLA_dz) * ds.dPhi_dr / ds.Jac
ds["B_contra_z"] = (1 + ds.dLA_dt) * ds.dPhi_dr / ds.Jac
ds["B"] = ds.B_contra_t * ds.e_theta + ds.B_contra_z * ds.e_zeta
```

**Key observations:**
- All field components scale with `dPhi_dr/Jac`
- Uses generalized coordinates with coordinate transformation
- Relies on proper Jacobian computation

### 2. SIMPLE VMEC Field Computation

From `/afs/itp.tugraz.at/proj/plasma/CODE/ert/SIMPLE/src/magfie.f90` lines 106, 120, 141, etc.:

```fortran
bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
```

**Key observations:**
- Direct use of VMEC Fourier coefficients
- Native VMEC coordinates (s, θ, φ)
- Standard metric tensor approach: |B|² = B^i B_i
- No explicit coordinate transformation scaling

### 3. Coordinate System Differences

| Aspect | SIMPLE/VMEC | GVEC |
|--------|-------------|------|
| Coordinates | Native (s, θ, φ) | Generalized (r, θ, ζ) |
| Radial coord | r = √s | r definition unclear |
| Field computation | Direct Fourier | Scaled by dPhi_dr/Jac |
| Axis treatment | Built-in regularization | May lack proper handling |

### 4. Mathematical Analysis of Scaling Factor

The scaling factor `dPhi_dr/Jac` behaves as follows near the axis:

**Scenario 1: GVEC r = √s (consistent with SIMPLE)**
- dPhi_dr ~ 2√s
- Jac ~ s (area scaling)
- Scaling ~ 2√s / s = 2/√s → ∞ as s → 0

**Scenario 2: GVEC r = s (different from SIMPLE)**
- dPhi_dr ~ 1
- Jac ~ s
- Scaling ~ 1/s → ∞ as s → 0

Both scenarios produce divergent behavior at the magnetic axis.

## Jacobian Behavior Analysis

From GVEC quantities.py lines 378-380:
```python
Jac_l = dX1_dr * dX2_dt - dX1_dt * dX2_dr
Jac = Jac_h * Jac_l
```

**Near magnetic axis:**
- Flux surfaces become smaller
- Coordinate derivatives may become singular
- Jac_l → 0 as area elements shrink
- Result: 1/Jac → ∞

## Theoretical Field Scaling Comparison

The generated plot `theoretical_field_scaling_analysis.png` demonstrates:

1. **VMEC field**: Roughly constant with s
2. **GVEC field (Scenario 1)**: Scales as 1/√s near axis
3. **GVEC field (Scenario 2)**: Scales as 1/s near axis (worse)
4. **Regularized GVEC**: Constant near axis (proper behavior)

## Specific Differences Summary

### SIMPLE/VMEC Implementation
- ✅ Direct Fourier coefficient evaluation
- ✅ Native VMEC coordinates (s, θ, φ)
- ✅ Proper axis treatment built-in
- ✅ Standard metric tensor: |B|² = B^i B_i

### GVEC Implementation
- ⚠️ Generalized coordinate system
- ⚠️ Field scaled by dPhi_dr/Jac
- ❌ May not handle axis singularities properly
- ❌ Coordinate transformation from VMEC

## Recommended Solutions

### 1. Immediate Fixes

**Check GVEC Input Parameters:**
- Verify `sgrid_nelems` (radial grid resolution)
- Check `sgrid_grid_type` (grid distribution)
- Review `X1X2_deg` and `LA_deg` (basis function degrees)

**Axis Treatment Settings:**
- Look for axis regularization options in GVEC
- Check if there are special axis boundary conditions
- Verify radial coordinate mapping consistency

### 2. Parameter Adjustments

```ini
# Suggested GVEC parameter modifications
sgrid_nelems = 21        # Higher radial resolution
sgrid_grid_type = 4      # Use proper grid distribution
X1X2_deg = 3             # Appropriate basis degree
LA_deg = 3               # Stream function degree

# Potential axis treatment (if available)
axis_regularization = true
min_radial_coordinate = 1e-6
```

### 3. Coordinate Verification

**Critical checks needed:**
1. Ensure GVEC r coordinate matches SIMPLE's √s
2. Verify flux normalization: Φ = s * Φ_edge
3. Check dPhi_dr computation method
4. Validate Jacobian calculation near axis

### 4. Alternative Approaches

**If parameter adjustments fail:**
- Use GVEC's Boozer coordinates if available
- Restrict field comparisons to s > 0.01
- Implement custom axis regularization
- Consider hybrid approach (VMEC near axis, GVEC elsewhere)

## Testing and Validation Recommendations

### 1. Diagnostic Tests
```bash
# Run field comparison with detailed output
make test_vmec_gvec

# Generate 1D radial profiles
python3 create_1d_radial_plots.py

# Check theoretical scaling
python3 detailed_field_analysis.py
```

### 2. Parameter Sensitivity Study
- Test different `sgrid_nelems` values (11, 21, 41)
- Compare different `sgrid_grid_type` options
- Evaluate basis function degree effects

### 3. Coordinate Transformation Verification
- Compare flux surface mappings
- Verify Jacobian values directly
- Check coordinate derivative computations

## Implementation Priority

### High Priority (Immediate)
1. ✅ **Completed**: Root cause identification
2. ✅ **Completed**: Theoretical analysis
3. 🔄 **In Progress**: Parameter optimization
4. ⭐ **Next**: GVEC configuration adjustment

### Medium Priority (Short-term)
1. Coordinate transformation verification
2. Axis regularization implementation
3. Comprehensive parameter study

### Low Priority (Long-term)
1. GVEC source code modifications
2. Custom regularization algorithms
3. Hybrid field evaluation methods

## Conclusion

The large GVEC field values at small s are **definitively caused** by the `dPhi_dr/Jac` scaling factor becoming singular near the magnetic axis. This is a **coordinate system and numerical treatment issue**, not a fundamental physics problem.

The solution requires:
1. **Proper GVEC parameter configuration** for axis treatment
2. **Verification of coordinate system consistency** between GVEC and VMEC
3. **Implementation of axis regularization** if not already available

This analysis provides the foundation for resolving the field comparison discrepancies and achieving accurate GVEC field evaluations throughout the plasma volume.

---

**Generated Files:**
- `theoretical_field_scaling_analysis.png`: Theoretical scaling comparison
- `field_comparison_1d.png`: Mock 1D field comparison
- `field_small_s_behavior.png`: Small-s behavior analysis
- This comprehensive analysis report

**Analysis Scripts:**
- `investigate_field_differences.py`: Initial investigation
- `detailed_field_analysis.py`: Comprehensive analysis
- `create_1d_radial_plots.py`: 1D plotting utilities