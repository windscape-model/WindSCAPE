#include "CRESCENT_Dynamics.H"

#include <cmath>

#include "CRESCENT_Constants.H"
#include "CRESCENT_Spectral.H"
#include "shared/CRESCENT_SpectralFilterExecution.H"
#include "shared/CRESCENT_SpectralFilterKernels.H"
#include "shared/CRESCENT_SpectralGridMetrics.H"

namespace crescent::dynamics {

using amrex::Box;
using amrex::GeometryData;
using amrex::IntVect;
using amrex::make_alias;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Real;
using amrex::TilingIfNotGPU;
namespace FFT = amrex::FFT;

void advanceSurfaceRK4(MultiFab& eta, MultiFab& phi,
                       MultiFab& eta_saved, MultiFab& phi_saved,
                       const MultiFab& rhs_eta, const MultiFab& rhs_phi,
                       const Real dt, const int rk_stage)
{
    // Valid region only: the RHS is only ever written on valid cells, so
    // updating ghost cells here would advance them from unwritten RHS ghosts.
    // No consumer reads state ghosts between stages (periodic scratch copies
    // use zero ghosts and all spectral operations act on valid data).
    AMREX_ASSERT(rk_stage < 4);

    if (rk_stage == 0) {
        // first stage: eta_saved = eta, phi_saved = phi
        // Fused kernel for both copies
        for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& eta_arr = eta.const_array(mfi);
            auto const& phi_arr = phi.const_array(mfi);
            auto eta0_arr = eta_saved.array(mfi);
            auto phi0_arr = phi_saved.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                eta0_arr(i,j,k) = eta_arr(i,j,k);
                phi0_arr(i,j,k) = phi_arr(i,j,k);
            });
        }
    } else {
        Real fac = 0.5 * dt;
        if (rk_stage == 3) fac = dt; // fourth stage

        const int comp = rk_stage - 1;

        // Update with weighted RHS: state = saved + factor * derivative
        // Fused kernel for both updates
        for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& eta0_arr = eta_saved.const_array(mfi);
            auto const& phi0_arr = phi_saved.const_array(mfi);
            auto const& feta_arr = rhs_eta.const_array(mfi);
            auto const& fphi_arr = rhs_phi.const_array(mfi);
            auto eta_arr = eta.array(mfi);
            auto phi_arr = phi.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                eta_arr(i,j,k) = eta0_arr(i,j,k) + fac * feta_arr(i,j,k,comp);
                phi_arr(i,j,k) = phi0_arr(i,j,k) + fac * fphi_arr(i,j,k,comp);
            });
        }
    }
}

void computeSurfaceRHS(MultiFab& rhs_eta, MultiFab& rhs_phi,
                       const MultiFab& eta,
                       const MultiFab& eta_x, const MultiFab& eta_y,
                       const MultiFab& phi_x, const MultiFab& phi_y,
                       const MultiFab& phi_z,
                       const MultiFab& pressure,
                       MultiFab& compute_rhs_metric,
                       MultiFab& compute_rhs_phi_z_sq,
                       std::unique_ptr<FFT::R2C<Real>>& fft,
                       const GeometryData& geom,
                       const Real inverse_gravity,
                       const int perturbation_order,
                       const int rk_stage,
                       MultiFab& dealiasing_scratch)
{
    AMREX_ALWAYS_ASSERT(rk_stage >= 0);
    AMREX_ALWAYS_ASSERT(rk_stage < rhs_eta.nComp());
    AMREX_ALWAYS_ASSERT(rk_stage < rhs_phi.nComp());
    AMREX_ALWAYS_ASSERT(phi_z.nComp() >= 1);

    const int d_nrk = rk_stage;

    if (perturbation_order > 1) {
        // Use pre-allocated workspace to avoid memory allocation overhead
        MultiFab& metric = compute_rhs_metric;
        MultiFab& phi_z_sq = compute_rhs_phi_z_sq;

        for (MFIter mfi(eta_x, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto eta_x_arr = eta_x.const_array(mfi);
            auto eta_y_arr = eta_y.const_array(mfi);
            auto phi_z_arr = phi_z.const_array(mfi);
            auto metric_arr = metric.array(mfi);
            auto phi_z_sq_arr = phi_z_sq.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                metric_arr(i,j,k) = 1.0 + eta_x_arr(i,j,k)*eta_x_arr(i,j,k) + eta_y_arr(i,j,k)*eta_y_arr(i,j,k);
                phi_z_sq_arr(i,j,k) = phi_z_arr(i,j,k) * phi_z_arr(i,j,k);
            });
        }

        // Compute metric term: metric = 1 + |∇η|²
        crescent::spectral::dealiasing(metric, fft, geom, dealiasing_scratch);

        // Surface evolution: dη/dt includes advection and vertical velocity
        for (MFIter mfi(rhs_eta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto eta_x_arr = eta_x.const_array(mfi);
            auto eta_y_arr = eta_y.const_array(mfi);
            auto phi_x_arr = phi_x.const_array(mfi);
            auto phi_y_arr = phi_y.const_array(mfi);
            auto phi_z_arr = phi_z.const_array(mfi);
            auto metric_arr = metric.const_array(mfi);

            auto feta_arr = rhs_eta.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                feta_arr(i,j,k,d_nrk) = -(phi_x_arr(i,j,k) * eta_x_arr(i,j,k) + phi_y_arr(i,j,k) * eta_y_arr(i,j,k)) + metric_arr(i,j,k)*phi_z_arr(i,j,k);
            });
        }
        MultiFab feta_nrk(rhs_eta, make_alias, rk_stage, 1);
        crescent::spectral::dealiasing(feta_nrk, fft, geom, dealiasing_scratch);

        // Dealias vertical velocity squared
        crescent::spectral::dealiasing(phi_z_sq, fft, geom, dealiasing_scratch);

        // Bernoulli equation for potential: includes kinetic energy and gravity
        for (MFIter mfi(rhs_phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto phi_x_arr = phi_x.const_array(mfi);
            auto phi_y_arr = phi_y.const_array(mfi);
            auto metric_arr = metric.const_array(mfi);
            auto phi_z_sq_arr = phi_z_sq.const_array(mfi);
            auto eta_arr = eta.const_array(mfi);

            auto fphi_arr = rhs_phi.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real vsq = phi_x_arr(i,j,k)*phi_x_arr(i,j,k) + phi_y_arr(i,j,k)*phi_y_arr(i,j,k);
                fphi_arr(i,j,k,d_nrk) = -0.5 * vsq
                                       + 0.5 * metric_arr(i,j,k) * phi_z_sq_arr(i,j,k)
                                       - eta_arr(i,j,k) / inverse_gravity;
            });
        }
        MultiFab fphi_nrk(rhs_phi, make_alias, rk_stage, 1);
        crescent::spectral::dealiasing(fphi_nrk, fft, geom, dealiasing_scratch);

        // The applied pressure is not band-limited by the projection: it is
        // subtracted after the dealiasing of the nonlinear terms.
        for (MFIter mfi(rhs_phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto fphi_arr = rhs_phi.array(mfi);
            auto pa_arr = pressure.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                fphi_arr(i, j, k, d_nrk) -= pa_arr(i, j, k);
            });
        }
    } else {
        // solve linear equation if perturbation_order=1
        MultiFab::Copy(rhs_eta, phi_z, 0, rk_stage, 1, 0);

        // Linear case: only gravity and pressure terms
        for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto eta_arr = eta.const_array(mfi);
            auto pa_arr = pressure.const_array(mfi);

            auto fphi_arr = rhs_phi.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                fphi_arr(i,j,k,d_nrk) = -eta_arr(i,j,k) / inverse_gravity - pa_arr(i,j,k);
            });
        }
    }

}

Real computeAdjustmentRampWeight(Real stage_time, Real adjustment_time)
{
    if (adjustment_time <= 0.0 || stage_time >= 2.0 * adjustment_time) {
        return 1.0;
    }
    const Real ratio = stage_time / adjustment_time;
    return 1.0 - std::exp(-(ratio * ratio * ratio * ratio));
}

void applyAdjustmentRampToStageRHS(
    MultiFab& rhs_eta, MultiFab& rhs_phi,
    const MultiFab& eta, const MultiFab& phi,
    const MultiFab& pressure,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom,
    const Real stage_time, const Real adjustment_time,
    const int rk_stage,
    MultiFab& w_lin_scratch)
{
    if (adjustment_time <= 0.0 || stage_time >= 2.0 * adjustment_time) {
        return; // W = 1 to round-off: fully nonlinear
    }
    const Real ramp =
        computeAdjustmentRampWeight(stage_time, adjustment_time);

    // w_lin = IFFT(|k| phi_hat): the flat-surface (order-1) DNO of the same
    // phi whose derivatives assembled the RHS.
    MultiFab::Copy(w_lin_scratch, phi, 0, 0, 1, IntVect(0));
    {
        const auto spectral_grid =
            crescent::shared::makeSpectralGridMetrics(w_lin_scratch, geom);
        const int ny = spectral_grid.ny;
        const Real dk_x = spectral_grid.dk_x;
        const Real dk_y = spectral_grid.dk_y;
        const Real norm = Real(spectral_grid.nx) * Real(spectral_grid.ny);
        const auto multiplier =
            [=] AMREX_GPU_DEVICE(int i, int j, int,
                                 amrex::GpuComplex<Real>& c) noexcept {
                const Real k_mag =
                    crescent::spectral::filter_kernels::spectralMagnitude(
                        i, j, ny, dk_x, dk_y);
                c = c * (k_mag / norm);
            };
        crescent::spectral::applyInPlaceSpectralFilter(
            w_lin_scratch, fft, multiplier);
    }

    const int comp = rk_stage;
    const Real g = 1.0 / crescent::constants::INV_GRAVITY;
    for (MFIter mfi(rhs_eta, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto feta = rhs_eta.array(mfi);
        auto fphi = rhs_phi.array(mfi);
        auto wlin = w_lin_scratch.const_array(mfi);
        auto eta_arr = eta.const_array(mfi);
        auto pa = pressure.const_array(mfi);
        const Real w_ramp = ramp;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real feta_lin = wlin(i, j, k);
            const Real fphi_lin = -g * eta_arr(i, j, k) - pa(i, j, k);
            feta(i, j, k, comp) =
                feta_lin + w_ramp * (feta(i, j, k, comp) - feta_lin);
            fphi(i, j, k, comp) =
                fphi_lin + w_ramp * (fphi(i, j, k, comp) - fphi_lin);
        });
    }
}

void updateSurface(MultiFab& eta, MultiFab& phi,
                   const MultiFab& eta_saved, const MultiFab& phi_saved,
                   const MultiFab& rhs_eta, const MultiFab& rhs_phi,
                   const Real dt)
{
    // RK4 final update: combine all stages with standard weights
    // New value = saved + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    const Real rk_weight = dt * (1.0 / 6.0);
    for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto eta_arr = eta.array(mfi);
        auto phi_arr = phi.array(mfi);
        auto eta_saved_arr = eta_saved.const_array(mfi);
        auto phi_saved_arr = phi_saved.const_array(mfi);
        auto rhs_eta_arr = rhs_eta.const_array(mfi);
        auto rhs_phi_arr = rhs_phi.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            eta_arr(i,j,k) = eta_saved_arr(i,j,k) + rk_weight * (
                rhs_eta_arr(i,j,k,0) + 2.0 * rhs_eta_arr(i,j,k,1) +
                2.0 * rhs_eta_arr(i,j,k,2) + rhs_eta_arr(i,j,k,3));
            phi_arr(i,j,k) = phi_saved_arr(i,j,k) + rk_weight * (
                rhs_phi_arr(i,j,k,0) + 2.0 * rhs_phi_arr(i,j,k,1) +
                2.0 * rhs_phi_arr(i,j,k,2) + rhs_phi_arr(i,j,k,3));
        });
    }
}

} // namespace crescent::dynamics
