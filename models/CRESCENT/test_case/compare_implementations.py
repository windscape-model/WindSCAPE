#!/usr/bin/env python3
"""
This script compares the results of the C++ and Fortran implementations
of the HOS method for ocean wave simulations.

Requirements:
- NumPy
- Matplotlib
- yt (for reading AMReX data)
- h5py (for reading Fortran HDF5 output)

Usage:
1. Run both implementations with the same parameters
2. Adjust file paths below to point to outputs
3. Run this script to generate comparison plots
4. Use the --test-type flag to specify 'jonswap' or 'stokes' test

Example:
    python compare_implementations.py --test-type stokes

For the Stokes wave test, machine precision comparison is expected.
For the JONSWAP test, only statistical similarity should be expected
due to the use of random numbers.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse

try:
    import yt
    from yt.frontends.boxlib.data_structures import AMReXDataset
except ImportError:
    print("Error: yt package with AMReX support is required")
    print("Install with: pip install yt")
    sys.exit(1)

try:
    import h5py
except ImportError:
    print("Error: h5py package is required")
    print("Install with: pip install h5py")
    sys.exit(1)

# Default paths - will be updated based on test type
CPP_OUTPUT_DIR = ""
FORTRAN_OUTPUT_DIR = ""

def read_amrex_data(filename):
    """Read surface elevation data from AMReX output"""
    ds = yt.load(filename)

    # Extract data from level 0 (assuming single-level simulation)
    data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                          dims=ds.domain_dimensions)

    # Extract relevant fields
    eta = data['eta'].d
    u = data['xvel'].d
    v = data['yvel'].d
    w = data['zvel'].d

    # Reshape from 3D to 2D (assuming z dimension is 1)
    eta = eta[:,:,0]
    u = u[:,:,0]
    v = v[:,:,0]
    w = w[:,:,0]

    # Get domain information
    domain_size = ds.domain_right_edge - ds.domain_left_edge
    nx, ny = eta.shape

    # Create coordinate arrays
    x = np.linspace(0, domain_size[0], nx)
    y = np.linspace(0, domain_size[1], ny)

    return {'eta': eta, 'u': u, 'v': v, 'w': w, 'x': x, 'y': y}

def read_fortran_data(filename):
    """Read surface elevation data from Fortran HDF5 output"""
    with h5py.File(filename, 'r') as f:
        # These field names may need to be adjusted based on actual output format
        eta = f['eta'][:]
        u = f['u_hos'][:]
        v = f['v_hos'][:]
        w = f['w_hos'][:]

        # Get domain information
        # Adjust these based on how domain info is stored in the Fortran output
        nx, ny = eta.shape
        lx = f.attrs.get('lx', 1.0)  # domain length in x
        ly = f.attrs.get('ly', 1.0)  # domain length in y

        # Create coordinate arrays
        x = np.linspace(0, lx, nx)
        y = np.linspace(0, ly, ny)

    return {'eta': eta, 'u': u, 'v': v, 'w': w, 'x': x, 'y': y}

def compare_data(cpp_data, fortran_data):
    """Compare data between C++ and Fortran implementations"""
    # Calculate statistics
    cpp_eta_max = np.max(cpp_data['eta'])
    cpp_eta_min = np.min(cpp_data['eta'])
    cpp_eta_mean = np.mean(cpp_data['eta'])
    cpp_eta_std = np.std(cpp_data['eta'])

    fortran_eta_max = np.max(fortran_data['eta'])
    fortran_eta_min = np.min(fortran_data['eta'])
    fortran_eta_mean = np.mean(fortran_data['eta'])
    fortran_eta_std = np.std(fortran_data['eta'])

    # Print comparison
    print("=== Surface Elevation Statistics ===")
    print(f"{'':20} {'C++':>12} {'Fortran':>12} {'Diff (%)':>12}")
    print(f"{'Maximum':20} {cpp_eta_max:12.6f} {fortran_eta_max:12.6f} {100*(cpp_eta_max-fortran_eta_max)/fortran_eta_max:12.6f}")
    print(f"{'Minimum':20} {cpp_eta_min:12.6f} {fortran_eta_min:12.6f} {100*(cpp_eta_min-fortran_eta_min)/fortran_eta_min:12.6f}")
    print(f"{'Mean':20} {cpp_eta_mean:12.6f} {fortran_eta_mean:12.6f} {100*(cpp_eta_mean-fortran_eta_mean)/fortran_eta_mean:12.6f}")
    print(f"{'Standard Deviation':20} {cpp_eta_std:12.6f} {fortran_eta_std:12.6f} {100*(cpp_eta_std-fortran_eta_std)/fortran_eta_std:12.6f}")

    # Calculate wave energy spectrum
    cpp_spectrum = np.abs(np.fft.fft2(cpp_data['eta']))**2
    fortran_spectrum = np.abs(np.fft.fft2(fortran_data['eta']))**2

    # Azimuthal average for 1D spectrum
    def azimuthal_average(data):
        ny, nx = data.shape
        kx = np.fft.fftfreq(nx) * nx
        ky = np.fft.fftfreq(ny) * ny
        kx, ky = np.meshgrid(kx, ky)

        # Radial coordinates
        k_radius = np.sqrt(kx**2 + ky**2)
        k_radius = k_radius.flatten()

        # Bin the data
        max_k = int(np.max(k_radius)) + 1
        k_bins = np.linspace(0, max_k, max_k + 1)
        spectrum_1d = np.zeros(max_k)

        flat_data = data.flatten()
        for i in range(max_k):
            mask = (k_radius >= i) & (k_radius < i+1)
            if np.any(mask):
                spectrum_1d[i] = np.mean(flat_data[mask])

        return spectrum_1d, k_bins[:-1]

    cpp_spec_1d, k = azimuthal_average(cpp_spectrum)
    fortran_spec_1d, _ = azimuthal_average(fortran_spectrum)

    # Create plots directory
    os.makedirs("comparison_plots", exist_ok=True)

    # Plot 1: Surface elevation comparison
    plt.figure(figsize=(12, 10))

    plt.subplot(2,2,1)
    plt.pcolormesh(cpp_data['x'], cpp_data['y'], cpp_data['eta'], shading='auto')
    plt.colorbar(label='Surface Elevation (m)')
    plt.title('C++ Implementation')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')

    plt.subplot(2,2,2)
    plt.pcolormesh(fortran_data['x'], fortran_data['y'], fortran_data['eta'], shading='auto')
    plt.colorbar(label='Surface Elevation (m)')
    plt.title('Fortran Implementation')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')

    plt.subplot(2,2,3)
    diff = cpp_data['eta'] - fortran_data['eta']
    plt.pcolormesh(cpp_data['x'], cpp_data['y'], diff, shading='auto', cmap='RdBu_r')
    plt.colorbar(label='Difference (m)')
    plt.title('Difference (C++ - Fortran)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')

    plt.subplot(2,2,4)
    # Cross section at middle of y
    mid_y = cpp_data['eta'].shape[1] // 2
    plt.plot(cpp_data['x'], cpp_data['eta'][:,mid_y], 'b-', label='C++')
    plt.plot(fortran_data['x'], fortran_data['eta'][:,mid_y], 'r--', label='Fortran')
    plt.title(f'Surface Elevation at y-midpoint')
    plt.xlabel('x (m)')
    plt.ylabel('Surface Elevation (m)')
    plt.legend()

    plt.tight_layout()
    plt.savefig('comparison_plots/surface_elevation_comparison.png', dpi=300)

    # Plot 2: Energy spectrum comparison
    plt.figure(figsize=(10, 6))
    plt.loglog(k[1:], cpp_spec_1d[1:], 'b-', label='C++')
    plt.loglog(k[1:], fortran_spec_1d[1:], 'r--', label='Fortran')
    plt.title('Wave Energy Spectrum Comparison')
    plt.xlabel('Wavenumber')
    plt.ylabel('Spectral Density')
    plt.legend()
    plt.grid(True, which='both', ls='-', alpha=0.7)
    plt.tight_layout()
    plt.savefig('comparison_plots/spectrum_comparison.png', dpi=300)

    # Plot 3: Velocity magnitude comparison
    plt.figure(figsize=(12, 10))

    # Calculate velocity magnitude
    cpp_vel_mag = np.sqrt(cpp_data['u']**2 + cpp_data['v']**2 + cpp_data['w']**2)
    fortran_vel_mag = np.sqrt(fortran_data['u']**2 + fortran_data['v']**2 + fortran_data['w']**2)

    plt.subplot(2,2,1)
    plt.pcolormesh(cpp_data['x'], cpp_data['y'], cpp_vel_mag, shading='auto')
    plt.colorbar(label='Velocity Magnitude (m/s)')
    plt.title('C++ Implementation')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')

    plt.subplot(2,2,2)
    plt.pcolormesh(fortran_data['x'], fortran_data['y'], fortran_vel_mag, shading='auto')
    plt.colorbar(label='Velocity Magnitude (m/s)')
    plt.title('Fortran Implementation')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')

    plt.subplot(2,2,3)
    vel_diff = cpp_vel_mag - fortran_vel_mag
    plt.pcolormesh(cpp_data['x'], cpp_data['y'], vel_diff, shading='auto', cmap='RdBu_r')
    plt.colorbar(label='Difference (m/s)')
    plt.title('Difference (C++ - Fortran)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')

    plt.subplot(2,2,4)
    # Velocity magnitude histogram
    plt.hist(cpp_vel_mag.flatten(), bins=50, alpha=0.5, label='C++')
    plt.hist(fortran_vel_mag.flatten(), bins=50, alpha=0.5, label='Fortran')
    plt.title('Velocity Magnitude Distribution')
    plt.xlabel('Velocity Magnitude (m/s)')
    plt.ylabel('Count')
    plt.legend()

    plt.tight_layout()
    plt.savefig('comparison_plots/velocity_comparison.png', dpi=300)

    # For Stokes waves, add a high-precision comparison plot
    if np.max(np.abs(cpp_data['eta'] - fortran_data['eta'])) < 1e-3:  # Close enough to warrant precision analysis
        plt.figure(figsize=(10, 8))
        plt.subplot(2,1,1)
        # Cross section at middle of y - zoomed in
        mid_y = cpp_data['eta'].shape[1] // 2
        plt.plot(cpp_data['x'], cpp_data['eta'][:,mid_y], 'b-', label='C++')
        plt.plot(fortran_data['x'], fortran_data['eta'][:,mid_y], 'r--', label='Fortran')
        plt.title(f'Surface Elevation at y-midpoint (Zoomed)')
        plt.xlabel('x (m)')
        plt.ylabel('Surface Elevation (m)')
        plt.legend()

        # Zoom in on a section
        x_range = cpp_data['x'].max() - cpp_data['x'].min()
        plt.xlim(cpp_data['x'].min() + 0.4*x_range, cpp_data['x'].min() + 0.6*x_range)

        plt.subplot(2,1,2)
        # Show the absolute difference
        plt.plot(cpp_data['x'], np.abs(cpp_data['eta'][:,mid_y] - fortran_data['eta'][:,mid_y]), 'k-')
        plt.title('Absolute Difference |C++ - Fortran|')
        plt.xlabel('x (m)')
        plt.ylabel('Difference (m)')
        plt.yscale('log')
        plt.xlim(cpp_data['x'].min() + 0.4*x_range, cpp_data['x'].min() + 0.6*x_range)
        plt.grid(True)

        plt.tight_layout()
        plt.savefig('comparison_plots/precision_comparison.png', dpi=300)

    return cpp_spec_1d, fortran_spec_1d, k

def main():
    parser = argparse.ArgumentParser(description="Compare C++ and Fortran HOS implementations")
    parser.add_argument("--test-type", choices=["jonswap", "stokes"], default="jonswap",
                        help="Type of test to compare: jonswap or stokes")
    args = parser.parse_args()

    # Set paths based on test type
    global CPP_OUTPUT_DIR, FORTRAN_OUTPUT_DIR
    if args.test_type == "jonswap":
        CPP_OUTPUT_DIR = "../build/Exec/comparison_test00020"
        FORTRAN_OUTPUT_DIR = "../HOS Fortran Code/output/hos_jonswap_output_00020.h5"
        precision_test = False
    else:  # stokes
        CPP_OUTPUT_DIR = "../build/Exec/stokes_test00020"
        FORTRAN_OUTPUT_DIR = "../HOS Fortran Code/output/hos_stokes_output_00020.h5"
        precision_test = True

    print(f"Running {args.test_type.upper()} comparison test")
    print(f"Precision comparison mode: {'ON' if precision_test else 'OFF'}")

    print("Reading C++ AMReX data...")
    try:
        cpp_data = read_amrex_data(CPP_OUTPUT_DIR)
    except Exception as e:
        print(f"Error reading C++ data: {e}")
        print(f"Looking for file: {CPP_OUTPUT_DIR}")
        sys.exit(1)

    print("Reading Fortran HDF5 data...")
    try:
        fortran_data = read_fortran_data(FORTRAN_OUTPUT_DIR)
    except Exception as e:
        print(f"Error reading Fortran data: {e}")
        print(f"Looking for file: {FORTRAN_OUTPUT_DIR}")
        sys.exit(1)

    print("Comparing data...")
    cpp_spec, fortran_spec, k = compare_data(cpp_data, fortran_data)

    if precision_test:
        # Calculate absolute max difference for precision testing
        abs_diff = np.max(np.abs(cpp_data['eta'] - fortran_data['eta']))
        rel_diff = np.max(np.abs((cpp_data['eta'] - fortran_data['eta']) /
                                 (np.abs(fortran_data['eta']) + 1e-10)))

        print("\n=== PRECISION COMPARISON RESULTS ===")
        print(f"Maximum absolute difference in surface elevation: {abs_diff:.3e}")
        print(f"Maximum relative difference in surface elevation: {rel_diff:.3e}")
        print(f"Machine epsilon for double precision: {np.finfo(float).eps:.3e}")

        if abs_diff < 1e-10:
            print("\nCONCLUSION: Implementations match with high precision!")
        else:
            print("\nCONCLUSION: Implementations show some differences.")

    print("\nComparison complete. Plots saved to 'comparison_plots' directory.")

if __name__ == "__main__":
    main()