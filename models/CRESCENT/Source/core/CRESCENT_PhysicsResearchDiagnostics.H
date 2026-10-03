#ifndef CRESCENT_PHYSICS_RESEARCH_DIAGNOSTICS_H
#define CRESCENT_PHYSICS_RESEARCH_DIAGNOSTICS_H

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <limits>
#include <string>

namespace crescent::physics::research_diagnostics {

inline bool shouldLogDetailedDiagnostics(int verbosity)
{
    return verbosity >= 2;
}

inline void printLeadingComponentValuesOnRow(
    const std::string& label,
    amrex::Array4<const amrex::Real> const& arr,
    int i_begin,
    int j,
    int k,
    int sample_count,
    int perturbation_order)
{
    for (int comp = 0; comp < perturbation_order; ++comp) {
        amrex::Print() << "  " << label << "[" << comp << "](i="
                       << i_begin << ":" << (i_begin + sample_count - 1)
                       << ", j=" << j << ", k=" << k << ") = ";
        for (int n = 0; n < sample_count; ++n) {
            amrex::Print() << arr(i_begin + n, j, k, comp);
            if (n < sample_count - 1) {
                amrex::Print() << ", ";
            }
        }
        amrex::Print() << std::endl;
    }
}

struct SpectralEnergyExtrema
{
    amrex::Real j0_max = 0.0;
    amrex::Real other_max = 0.0;
};

inline SpectralEnergyExtrema measureSpectralEnergyExtrema(
    const amrex::cMultiFab& spectral_data)
{
    SpectralEnergyExtrema extrema;

    for (amrex::MFIter mfi(spectral_data); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto spectral_arr = spectral_data.const_array(mfi);
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        if (lo.y <= 0 && 0 <= hi.y) {
            for (int i = lo.x; i <= hi.x; ++i) {
                const auto coeff = spectral_arr(i, 0, 0);
                const amrex::Real energy =
                    coeff.real() * coeff.real() + coeff.imag() * coeff.imag();
                extrema.j0_max = std::max(extrema.j0_max, energy);
            }
        }

        for (int j = lo.y; j <= hi.y; ++j) {
            if (j != 0) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    const auto coeff = spectral_arr(i, j, 0);
                    const amrex::Real energy =
                        coeff.real() * coeff.real() + coeff.imag() * coeff.imag();
                    extrema.other_max = std::max(extrema.other_max, energy);
                }
            }
        }
    }

    return extrema;
}

struct ColumnExtrema
{
    amrex::Real min = std::numeric_limits<amrex::Real>::max();
    amrex::Real max = std::numeric_limits<amrex::Real>::lowest();
};

struct ColumnProbeExtrema
{
    int i = 0;
    amrex::Real x = 0.0;
    ColumnExtrema extrema;
};

inline int centerColumnIndex(const amrex::MultiFab& field)
{
    const amrex::Box domain = field.boxArray().minimalBox();
    return (domain.smallEnd(0) + domain.bigEnd(0)) / 2;
}

inline amrex::Real cellCenterX(const amrex::GeometryData& geom, int i)
{
    return geom.ProbLo(0) +
           (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * geom.CellSize(0);
}

inline ColumnExtrema measureColumnExtremaAtIndex(
    const amrex::MultiFab& field,
    int component,
    int i_sample)
{
    ColumnExtrema extrema;

    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto field_arr = field.const_array(mfi);
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        if (lo.x <= i_sample && i_sample <= hi.x) {
            for (int j = lo.y; j <= hi.y; ++j) {
                const amrex::Real value = field_arr(i_sample, j, 0, component);
                extrema.min = std::min(extrema.min, value);
                extrema.max = std::max(extrema.max, value);
            }
        }
    }

    return extrema;
}

inline ColumnProbeExtrema measureCenterColumnExtrema(
    const amrex::MultiFab& field,
    const amrex::GeometryData& geom,
    int component)
{
    const int i_sample = centerColumnIndex(field);
    return {
        i_sample,
        cellCenterX(geom, i_sample),
        measureColumnExtremaAtIndex(field, component, i_sample)
    };
}

} // namespace crescent::physics::research_diagnostics

#endif // CRESCENT_PHYSICS_RESEARCH_DIAGNOSTICS_H
