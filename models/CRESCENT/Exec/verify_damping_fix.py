#!/usr/bin/env python3
"""
Verify that damping is applied correctly on all four boundaries after the fix.
Tests the damping function logic.
"""
import numpy as np

# Domain configuration (from multibox_mpi4_test.log)
L1_lo_x, L1_hi_x = 37.5, 162.5  # L1 box boundaries
L1_lo_y, L1_hi_y = 37.5, 162.5
x_lo_S, x_hi_S = 56.25, 143.75  # Spectral subdomain boundaries
y_lo_S, y_hi_S = 56.25, 143.75
W_damp_x = 18.75  # Damping zone width
W_damp_y = 18.75
alpha = 2.0  # Steepness parameter

def compute_chi_D(x, y):
    """
    Compute chi_D mask using the FIXED algorithm
    """
    # X-direction damping
    if x < x_lo_S:
        # Left damping zone (FIXED: use L1_lo_x, not 0)
        d = x - L1_lo_x
        arg = alpha * (2.0 * d / W_damp_x - 1.0)
        chi_D_x = 0.5 * (1.0 + np.tanh(arg))
    elif x > x_hi_S:
        # Right damping zone
        d = x - x_hi_S
        arg = alpha * (2.0 * d / W_damp_x - 1.0)
        chi_D_x = 0.5 * (1.0 - np.tanh(arg))
    else:
        chi_D_x = 1.0
    
    # Y-direction damping
    if y < y_lo_S:
        # Bottom damping zone (FIXED: use L1_lo_y, not 0)
        d = y - L1_lo_y
        arg = alpha * (2.0 * d / W_damp_y - 1.0)
        chi_D_y = 0.5 * (1.0 + np.tanh(arg))
    elif y > y_hi_S:
        # Top damping zone
        d = y - y_hi_S
        arg = alpha * (2.0 * d / W_damp_y - 1.0)
        chi_D_y = 0.5 * (1.0 - np.tanh(arg))
    else:
        chi_D_y = 1.0
    
    return chi_D_x * chi_D_y

# Test points at boundaries
test_points = [
    ("Left edge", 38.28, 100.0),
    ("Right edge", 161.72, 100.0),
    ("Bottom edge", 100.0, 38.28),
    ("Top edge", 100.0, 161.72),
    ("Center", 100.0, 100.0),
]

print("=" * 60)
print("Damping Mask Verification (AFTER FIX)")
print("=" * 60)
print()

for name, x, y in test_points:
    chi_D = compute_chi_D(x, y)
    print(f"{name:15s} ({x:6.2f}, {y:6.2f}): chi_D = {chi_D:.6f}")
    if "edge" in name.lower() and "center" not in name.lower():
        if chi_D < 0.1:
            print(f"                ✓ Strong damping (expected)")
        elif chi_D < 0.5:
            print(f"                ⚠ Moderate damping")
        else:
            print(f"                ✗ NO damping (BUG!)")

print()
print("Expected behavior:")
print("  - All four edges should show strong damping (chi_D < 0.1)")
print("  - Center should show NO damping (chi_D ≈ 1.0)")
print()

# Create a 2D map of chi_D
nx, ny = 80, 80
dx = (L1_hi_x - L1_lo_x) / nx
dy = (L1_hi_y - L1_lo_y) / ny

chi_D_map = np.zeros((ny, nx))
for j in range(ny):
    for i in range(nx):
        x = L1_lo_x + (i + 0.5) * dx
        y = L1_lo_y + (j + 0.5) * dy
        chi_D_map[j, i] = compute_chi_D(x, y)

print("Summary statistics:")
print(f"  chi_D min: {chi_D_map.min():.6f}")
print(f"  chi_D max: {chi_D_map.max():.6f}")
print(f"  chi_D at corners (should be very small):")
print(f"    Bottom-left:  {chi_D_map[0, 0]:.6f}")
print(f"    Bottom-right: {chi_D_map[0, -1]:.6f}")
print(f"    Top-left:     {chi_D_map[-1, 0]:.6f}")
print(f"    Top-right:    {chi_D_map[-1, -1]:.6f}")

