#include "CRESCENT_Physics_Helpers.H"

#include <AMReX_GpuLaunch.H>

namespace crescent::physics::helpers {

using amrex::Box;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Real;
using amrex::TilingIfNotGPU;

void computeEtaPowersWithFactorial(
    MultiFab& eta_pow,
    const MultiFab& eta,
    int perturbation_order)
{
    // Compute powers of surface elevation divided by factorial for Taylor expansion
    // eta_pow[k] stores η^(k+1) / (k+1)! needed for vertical coordinate transformation

    // Start with eta_pow[0] = η
    MultiFab::Copy(eta_pow, eta, 0, 0, 1, 0);

    // For each order k=1,...,perturbation_order-2:
    for (int k = 1; k < perturbation_order - 1; ++k)
    {
        // Nonlinear update for this order
        for (MFIter mfi(eta_pow, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto eta_pow_arr = eta_pow.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int z) noexcept
            {
                // Recursively build higher powers: η^(k+1)/(k+1)! = [η^k/k!] * η / (k+1)
                eta_pow_arr(i, j, z, k) = eta_pow_arr(i, j, z, k-1) * eta_pow_arr(i, j, z, 0) / static_cast<Real>(k + 1);
            });
        }
    }
}

} // namespace crescent::physics::helpers
