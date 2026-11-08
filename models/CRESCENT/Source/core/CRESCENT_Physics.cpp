#include "CRESCENT_Physics.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_GpuComplex.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_Spectral.H"
#include "CRESCENT_SpectralDerivativeExecution.H"
#include "CRESCENT_SpectralFilterKernels.H"
#include "CRESCENT_SpectralGridMetrics.H"

namespace crescent::physics {

using amrex::Box;
using amrex::cMultiFab;
using amrex::Copy;
using amrex::GeometryData;
using amrex::GpuComplex;
using amrex::make_alias;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Real;
using amrex::TilingIfNotGPU;

namespace FFT = amrex::FFT;

namespace {

void computeDirectionalDerivativesNoFilter(
    MultiFab& field_x,
    MultiFab& field_y,
    const MultiFab& field,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const crescent::shared::SpectralGridMetrics& spectral_grid,
    Real dkx,
    Real dky,
    cMultiFab& input_spectral_scratch,
    cMultiFab& x_spectral_scratch,
    cMultiFab& y_spectral_scratch)
{
    const auto unit_spectral_weight =
        [] AMREX_GPU_HOST_DEVICE(int, int) noexcept -> Real {
            return Real(1.0);
        };

    crescent::spectral::executeDirectionalDerivativePair(
        field_x,
        field_y,
        field,
        fft,
        spectral_grid,
        dkx,
        dky,
        unit_spectral_weight,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch);
}

// eta_pow(n) = eta^(n+1)/(n+1)!, n = 0..M-2. All powers are formed first;
// then every power above the first is projected with the 2/3 rule.
void computeEtaPowers(MultiFab& eta_pow,
                      const MultiFab& eta,
                      std::unique_ptr<FFT::R2C<Real>>& fft,
                      const GeometryData& geom,
                      const int perturbation_order,
                      MultiFab& dealiasing_scratch)
{
    MultiFab::Copy(eta_pow, eta, 0, 0, 1, 0);
    for (int k = 1; k < perturbation_order - 1; ++k)
    {
        for (MFIter mfi(eta_pow, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto eta_pow_arr = eta_pow.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int z) noexcept
            {
                eta_pow_arr(i, j, z, k) = eta_pow_arr(i, j, z, k-1) * eta_pow_arr(i, j, z, 0) / static_cast<Real>(k + 1);
            });
        }
    }

    for (int k = 1; k < perturbation_order - 1; ++k)
    {
        MultiFab tmp(eta_pow, make_alias, k, 1);
        crescent::spectral::dealiasing(tmp, fft, geom, dealiasing_scratch);
    }
}

void forwardPhiExpansionComponent(
    cMultiFab& spectral_phi_expansion,
    MultiFab& phi_expansion,
    int component,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom)
{
    MultiFab physical_component(phi_expansion, make_alias, component, 1);
    cMultiFab spectral_component(spectral_phi_expansion, make_alias, component, 1);
    crescent::spectral::forwardFFT(
        fft, geom, physical_component, spectral_component);
}

// deriv_spectral = |k|^(eta_order+1) * spectrum of phi_(phi_component)
void multiplySpectralExpansionTerm(
    cMultiFab& deriv_spectral,
    cMultiFab& phi_expansion_spectrum,
    const MultiFab& wavenumber_pow,
    int eta_order,
    int phi_component)
{
    for (MFIter mfi(phi_expansion_spectrum, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();

        auto phi_spectrum_arr = phi_expansion_spectrum.array(mfi);
        auto wavenumber_arr = wavenumber_pow.const_array(mfi);
        auto dr_c_arr = deriv_spectral.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dr_c_arr(i,j,k) =
                phi_spectrum_arr(i,j,k,phi_component) *
                wavenumber_arr(i,j,k,eta_order);
        });
    }
}

// phi_expansion(kw) -= eta_pow(k1) * expansion_deriv
void subtractWeightedPhysicalContribution(
    MultiFab& phi_expansion,
    const MultiFab& expansion_deriv,
    const MultiFab& eta_pow,
    int kw,
    int k1)
{
    for (MFIter mfi(expansion_deriv, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto phi_arr = phi_expansion.array(mfi);
        auto dr_arr = expansion_deriv.const_array(mfi);
        auto eta_pow_arr = eta_pow.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi_arr(i,j,k,kw) = phi_arr(i,j,k,kw) - eta_pow_arr(i,j,k,k1) * dr_arr(i,j,k);
        });
    }
}

// Two-thirds rule as a per-axis (square) mode-index cutoff, the same rule as
// crescent::spectral::dealiasing.
void applyTwoThirdsSpectralCutoff(cMultiFab& spectral_data, int nx, int ny)
{
    int cutoff_x = static_cast<int>(nx * crescent::constants::ONE_THIRD_CUTOFF);
    int cutoff_y = static_cast<int>(ny * crescent::constants::ONE_THIRD_CUTOFF);

    for (MFIter mfi(spectral_data, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        auto spectral_arr = spectral_data.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const int kx =
                crescent::spectral::filter_kernels::spectralKxIndex(i, nx);
            const int ky =
                crescent::spectral::filter_kernels::spectralKyIndex(j, ny);

            if (crescent::spectral::filter_kernels::spectralModeExceedsCutoff(
                    kx, ky, cutoff_x, cutoff_y)) {
                spectral_arr(i,j,k) = GpuComplex<Real>(0.0, 0.0);
            }
        });
    }
}

// Potential coefficients phi_m (m = 0..M-1) of the HOS expansion at z = 0:
//   phi_0 = phi_s,  phi_m = - sum_{n=0}^{m-1} eta^(n+1)/(n+1)! |k|^(n+1) phi_(m-n-1).
// Each phi_m is kept in spectral space for the later orders; at the end
// every coefficient is truncated with the 2/3 rule and transformed back
// once.
void expandVelocityPotential(
    const MultiFab& surface_phi,
    const MultiFab& wavenumber_pow,
    const VerticalVelocityWorkspace& workspace,
    int perturbation_order,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom)
{
    auto& phi_expansion = workspace.phi_expansion;

    Copy(phi_expansion, surface_phi, 0, 0, 1, 0);
    forwardPhiExpansionComponent(
        workspace.phi_expansion_spectrum, phi_expansion, 0, fft, geom);

    for (int m = 1; m < perturbation_order; ++m)
    {
        phi_expansion.setVal(0.0, m, 1, 0);

        for (int n = 0; n < m; ++n)
        {
            multiplySpectralExpansionTerm(
                workspace.deriv_spectral,
                workspace.phi_expansion_spectrum,
                wavenumber_pow,
                n,
                m - n - 1);
            crescent::spectral::backwardFFT(
                fft, workspace.deriv_spectral, workspace.expansion_deriv);
            subtractWeightedPhysicalContribution(
                phi_expansion, workspace.expansion_deriv, workspace.eta_pow, m, n);
        }

        forwardPhiExpansionComponent(
            workspace.phi_expansion_spectrum, phi_expansion, m, fft, geom);
    }

    const Box& mf_box = phi_expansion.boxArray().minimalBox();
    const int nx = mf_box.length(0);
    const int ny = mf_box.length(1);
    for (int m = 0; m < perturbation_order; ++m)
    {
        cMultiFab spectral_component(
            workspace.phi_expansion_spectrum, make_alias, m, 1);
        applyTwoThirdsSpectralCutoff(spectral_component, nx, ny);

        MultiFab physical_component(phi_expansion, make_alias, m, 1);
        fft->backward(spectral_component, physical_component);
    }
}

// output = |k|^wvn_power applied to input (forward and backward transform);
// kpow holds |k|^(n+1) in component n.
void applyWavenumberPower(
    const MultiFab& input,
    MultiFab& output,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    amrex::Array4<const Real> kpow,
    int wvn_power)
{
    const Box& mf_box = input.boxArray().minimalBox();
    const int nx = mf_box.length(0);
    const int ny = mf_box.length(1);
    const Real normalization = static_cast<Real>(nx) * static_cast<Real>(ny);
    const int component = wvn_power - 1;

    fft->forwardThenBackward(input, output, [=] AMREX_GPU_DEVICE (int i, int j, int k, GpuComplex<Real> &fft_x) noexcept
    {
        const Real wvn_val = kpow(i, j, k, component);

        fft_x = fft_x * wvn_val / normalization;
    });
}

// output += weights(weight_comp) * input
void accumulateWeightedIntoField(
    MultiFab& output,
    const MultiFab& weights,
    const MultiFab& input,
    int weight_comp)
{
    for (MFIter mfi(input, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto output_arr = output.array(mfi);
        auto input_arr = input.const_array(mfi);
        auto weight_arr = weights.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            output_arr(i,j,k) += weight_arr(i,j,k,weight_comp) * input_arr(i,j,k);
        });
    }
}

void accumulateIntoField(MultiFab& output, const MultiFab& input)
{
    for (MFIter mfi(input, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto output_arr = output.array(mfi);
        auto input_arr = input.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            output_arr(i,j,k) += input_arr(i,j,k);
        });
    }
}

// Vertical-velocity series from the potential coefficients: with the
// cumulative sums P_n = sum_{m<=n} phi_m,
//   phi_z = sum_{n=0}^{M-2} eta^(n+1)/(n+1)! |k|^(n+2) P_(M-n-2) + |k| P_(M-1).
void accumulateVerticalVelocity(
    MultiFab& phi_z,
    const VerticalVelocityWorkspace& workspace,
    int perturbation_order,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    Real dk_x,
    Real dk_y)
{
    auto& pk = workspace.expansion_accum;
    const auto& phi_expansion = workspace.phi_expansion;

    const int ny = pk.boxArray().minimalBox().length(1);
    workspace.kpow_table.update(
        *fft,
        perturbation_order,
        {dk_x, dk_y, Real(ny)},
        [=] AMREX_GPU_DEVICE(int i, int j, int n) noexcept {
            const Real k_mag =
                crescent::spectral::filter_kernels::spectralMagnitude(
                    i, j, ny, dk_x, dk_y);
            return std::pow(k_mag, n + 1);
        });
    const auto kpow = workspace.kpow_table.localWeights();

    Copy(pk, phi_expansion, 0, 0, 1, 0);
    for (int n = 0; n < perturbation_order - 1; ++n) {
        for (MFIter mfi(pk, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto pk_arr = pk.array(mfi);
            auto phi_arr = phi_expansion.const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                pk_arr(i,j,k,n+1) = pk_arr(i,j,k,n) + phi_arr(i,j,k,n+1);
            });
        }
    }

    phi_z.setVal(0.0);

    for (int n = 0; n <= perturbation_order - 2; ++n)
    {
        MultiFab pk_component(pk, make_alias, perturbation_order - n - 2, 1);
        applyWavenumberPower(
            pk_component, workspace.spectral_work, fft, kpow, n + 2);
        accumulateWeightedIntoField(
            phi_z, workspace.eta_pow, workspace.spectral_work, n);
    }

    MultiFab pk_last(pk, make_alias, perturbation_order - 1, 1);
    applyWavenumberPower(
        pk_last, workspace.spectral_work, fft, kpow, 1);
    accumulateIntoField(phi_z, workspace.spectral_work);
}

} // namespace

void computeSurfaceDerivatives(MultiFab& eta_x, MultiFab& eta_y,
                                MultiFab& phi_x, MultiFab& phi_y,
                                const MultiFab& eta,
                                const MultiFab& phi,
                                std::unique_ptr<FFT::R2C<Real>>& fft,
                                const GeometryData& geom,
                                const Real dk_x, const Real dk_y,
                                cMultiFab& input_spectral_scratch,
                                cMultiFab& x_spectral_scratch,
                                cMultiFab& y_spectral_scratch)
{
    // Compute each field's x/y derivatives from one forward FFT, instead of
    // separate x and y transforms for the same input field.

    const auto spectral_grid = crescent::shared::makeSpectralGridMetrics(eta, geom);

    computeDirectionalDerivativesNoFilter(
        eta_x,
        eta_y,
        eta,
        fft,
        spectral_grid,
        dk_x,
        dk_y,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch);
    computeDirectionalDerivativesNoFilter(
        phi_x,
        phi_y,
        phi,
        fft,
        spectral_grid,
        dk_x,
        dk_y,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch);
}

void computeVerticalVelocity(MultiFab& phi_z,
                              const MultiFab& eta,
                              const MultiFab& phi,
                              const MultiFab& wavenumber_pow,
                              const VerticalVelocityWorkspace& workspace,
                              std::unique_ptr<FFT::R2C<Real>>& fft,
                              const GeometryData& geom,
                              int perturbation_order,
                              Real dk_x,
                              Real dk_y,
                              bool apply_final_projection)
{
    BL_PROFILE("crescent::physics::computeVerticalVelocity");

    computeEtaPowers(
        workspace.eta_pow,
        eta,
        fft,
        geom,
        perturbation_order,
        workspace.dealiasing_scratch);
    expandVelocityPotential(
        phi, wavenumber_pow, workspace, perturbation_order, fft, geom);
    accumulateVerticalVelocity(
        phi_z, workspace, perturbation_order, fft, dk_x, dk_y);

    if (apply_final_projection) {
        crescent::spectral::dealiasing(
            phi_z, fft, geom, workspace.dealiasing_scratch);
    }
}

void computeHorizontalVelocities(
    MultiFab& u,
    MultiFab& v,
    crescent::shared::DerivativeBundleView derivatives,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom,
    MultiFab& dealiasing_scratch,
    bool apply_spectral_dealiasing)
{
    for (MFIter mfi(derivatives.eta_x, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto us_arr = u.array(mfi);
        auto vs_arr = v.array(mfi);
        auto phi_z_arr = derivatives.phi_z.const_array(mfi);
        auto eta_x_arr = derivatives.eta_x.const_array(mfi);
        auto eta_y_arr = derivatives.eta_y.const_array(mfi);
        auto phi_x_arr = derivatives.phi_x.const_array(mfi);
        auto phi_y_arr = derivatives.phi_y.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            us_arr(i,j,k) = phi_x_arr(i,j,k) - phi_z_arr(i,j,k) * eta_x_arr(i,j,k);
            vs_arr(i,j,k) = phi_y_arr(i,j,k) - phi_z_arr(i,j,k) * eta_y_arr(i,j,k);
        });
    }

    if (apply_spectral_dealiasing) {
        crescent::spectral::dealiasing(u, fft, geom, dealiasing_scratch);
        crescent::spectral::dealiasing(v, fft, geom, dealiasing_scratch);
    }
}

} // namespace crescent::physics
