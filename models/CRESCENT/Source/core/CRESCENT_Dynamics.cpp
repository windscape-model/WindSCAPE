#include "CRESCENT_Dynamics.H"
#include "CRESCENT_Spectral.H"

namespace crescent::dynamics {

using amrex::Box;
using amrex::GeometryData;
using amrex::IntVect;
using amrex::make_alias;
using amrex::makeSlab;
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
    const int ng = eta.nGrow();
    AMREX_ASSERT(rk_stage < 4);

    if (rk_stage == 0) {
        // first stage: eta_saved = eta, phi_saved = phi
        // Fused kernel for both copies
        for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(ng);
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
            const Box& bx = mfi.growntilebox(ng);
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
    IntVect ng(0,0,0);
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
            auto pa_arr = pressure.const_array(mfi);

            auto fphi_arr = rhs_phi.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real vsq = phi_x_arr(i,j,k)*phi_x_arr(i,j,k) + phi_y_arr(i,j,k)*phi_y_arr(i,j,k);
                fphi_arr(i,j,k,d_nrk) = -0.5 * vsq + 0.5 * metric_arr(i,j,k) * phi_z_sq_arr(i,j,k) - eta_arr(i,j,k) / inverse_gravity - pa_arr(i,j,k);
            });
        }
        MultiFab fphi_nrk(rhs_phi, make_alias, rk_stage, 1);
        crescent::spectral::dealiasing(fphi_nrk, fft, geom, dealiasing_scratch);
    } else {
        // solve linear equation if perturbation_order=1
        MultiFab::Copy(rhs_eta, phi_z, rk_stage, 0, 1, ng);

        // Linear case: only gravity and pressure terms
        for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& sbx = makeSlab(bx, 2, rk_stage);

            auto eta_arr = eta.const_array(mfi);
            auto pa_arr = pressure.const_array(mfi);

            auto fphi_arr = rhs_phi.array(mfi);

            ParallelFor(sbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                fphi_arr(i,j,k,d_nrk) = -eta_arr(i,j,k) / inverse_gravity - pa_arr(i,j,k);
            });
        }
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
