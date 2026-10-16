#ifndef CRESCENT_SPECTRAL_WEIGHT_TABLE_H
#define CRESCENT_SPECTRAL_WEIGHT_TABLE_H

#include <vector>

#include <AMReX_Array4.H>
#include <AMReX_FFT.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFab.H>

namespace crescent::shared {

// Real per-mode weights on the spectral layout of an R2C transform, built
// once for a parameter set and reused by every transform with the same
// parameters. A filter that reads a weight gets the value that evaluating
// the same expression inside the transform callback would give, without the
// per-mode pow and exp of every call.
class SpectralWeightTable
{
public:
    // Builds the table unless it already holds ncomp components for these
    // parameters. fill(i, j, n) gives component n at spectral mode (i, j).
    template <typename Fill>
    void update(amrex::FFT::R2C<amrex::Real>& fft,
                int ncomp,
                const std::vector<amrex::Real>& parameters,
                Fill fill)
    {
        if (weights_.ok() && weights_.nComp() == ncomp && parameters == parameters_) {
            return;
        }
        auto const& [cba, cdm] = fft.getSpectralDataLayout();
        weights_.clear();
        weights_.define(cba, cdm, ncomp, 0);
        for (amrex::MFIter mfi(weights_); mfi.isValid(); ++mfi) {
            auto weights = weights_.array(mfi);
            amrex::ParallelFor(mfi.validbox(), ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    weights(i, j, k, n) = fill(i, j, n);
                });
        }
        parameters_ = parameters;
    }

    // Weights of this rank's spectral box, indexed like the (i, j, k) of the
    // transform callback; the R2C spectral layout has at most one box per
    // rank.
    [[nodiscard]] amrex::Array4<const amrex::Real> localWeights() const
    {
        amrex::Array4<const amrex::Real> local;
        int nlocal = 0;
        for (amrex::MFIter mfi(weights_); mfi.isValid(); ++mfi) {
            local = weights_.const_array(mfi);
            ++nlocal;
        }
        AMREX_ALWAYS_ASSERT(nlocal <= 1);
        return local;
    }

private:
    amrex::MultiFab weights_;
    std::vector<amrex::Real> parameters_;
};

} // namespace crescent::shared

#endif // CRESCENT_SPECTRAL_WEIGHT_TABLE_H
