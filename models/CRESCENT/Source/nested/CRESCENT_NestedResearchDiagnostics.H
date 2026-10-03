#ifndef CRESCENT_NESTED_RESEARCH_DIAGNOSTICS_H
#define CRESCENT_NESTED_RESEARCH_DIAGNOSTICS_H

#include <algorithm>
#include <cmath>

#include <AMReX_MultiFab.H>

#include "CRESCENT_NestedTypes.H"

namespace crescent::nested::research_diagnostics {

inline amrex::Real computeMaxAdjacentYDifference(const amrex::MultiFab& mf)
{
    amrex::Real max_dy = 0.0;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& arr = mf.array(mfi);
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        for (int j = lo.y; j < hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                const amrex::Real dy = std::abs(arr(i, j + 1, 0) - arr(i, j, 0));
                max_dy = std::max(max_dy, dy);
            }
        }
    }

    return max_dy;
}

inline crescent::nested::ColumnExtrema sampleColumnExtrema(
    const amrex::MultiFab& mf,
    int column_i,
    int component)
{
    crescent::nested::ColumnExtrema extrema {1e30, -1e30, -1, -1, 0};

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& arr = mf.array(mfi);
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        if (lo.x <= column_i && column_i <= hi.x) {
            for (int j = lo.y; j <= hi.y; ++j) {
                const amrex::Real val = arr(column_i, j, 0, component);
                extrema.sample_count++;
                if (val < extrema.min_value) {
                    extrema.min_value = val;
                    extrema.min_j = j;
                }
                if (val > extrema.max_value) {
                    extrema.max_value = val;
                    extrema.max_j = j;
                }
            }
        }
    }

    return extrema;
}

inline crescent::nested::BoundaryColumnExtrema sampleBoundaryColumnExtrema(
    const amrex::MultiFab& mf,
    int left_column_i,
    int right_column_i,
    int component)
{
    return {sampleColumnExtrema(mf, left_column_i, component),
            sampleColumnExtrema(mf, right_column_i, component)};
}

inline amrex::Real columnVariation(const crescent::nested::ColumnExtrema& extrema)
{
    return (extrema.sample_count > 0) ? (extrema.max_value - extrema.min_value) : 0.0;
}

} // namespace crescent::nested::research_diagnostics

#endif // CRESCENT_NESTED_RESEARCH_DIAGNOSTICS_H
