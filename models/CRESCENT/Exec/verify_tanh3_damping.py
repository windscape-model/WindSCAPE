#!/usr/bin/env python3
"""Test tanh³ damping function"""
import numpy as np

L1_lo_x, L1_hi_x = 37.5, 162.5
x_lo_S, x_hi_S = 56.25, 143.75
W_damp_x = 18.75
alpha = 2.0

def compute_chi_D_tanh3(x):
    """Compute chi_D using tanh³ (improved)"""
    if x < x_lo_S:
        d = x - L1_lo_x
        arg = alpha * (2.0 * d / W_damp_x - 1.0)
        tanh_val = np.tanh(arg)
        chi_D_x = 0.5 * (1.0 + tanh_val**3)
    elif x > x_hi_S:
        d = x - x_hi_S
        arg = alpha * (2.0 * d / W_damp_x - 1.0)
        tanh_val = np.tanh(arg)
        chi_D_x = 0.5 * (1.0 - tanh_val**3)
    else:
        chi_D_x = 1.0
    return chi_D_x

print("Damping function comparison:")
print(f"{'Position':15s} {'tanh':>12s} {'tanh³':>12s} {'Improvement':>15s}")
print("-" * 60)

test_x = [37.5, 38.28, 40.0, 45.0, 50.0, 56.25, 100.0]
for x in test_x:
    # Original tanh
    if x < x_lo_S:
        d = x - L1_lo_x
        arg = alpha * (2.0 * d / W_damp_x - 1.0)
        chi_tanh = 0.5 * (1.0 + np.tanh(arg))
    else:
        chi_tanh = 1.0
    
    chi_tanh3 = compute_chi_D_tanh3(x)
    
    if x < x_lo_S:
        improvement = chi_tanh / chi_tanh3 if chi_tanh3 > 1e-10 else 0
        print(f"{x:6.2f} m:      {chi_tanh:12.6e} {chi_tanh3:12.6e}   {improvement:6.1f}x smaller")
    else:
        print(f"{x:6.2f} m:      {chi_tanh:12.6f} {chi_tanh3:12.6f}   (no damping)")

print()
print(f"At L1 boundary (x=37.5): chi_D = {compute_chi_D_tanh3(37.5):.6e}")
print(f"  → Meets requirement: {compute_chi_D_tanh3(37.5) < 1e-4}")

