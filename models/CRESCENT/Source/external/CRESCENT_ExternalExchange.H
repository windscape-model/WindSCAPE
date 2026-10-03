#ifndef CRESCENT_EXTERNAL_EXCHANGE_H
#define CRESCENT_EXTERNAL_EXCHANGE_H

#include <AMReX_REAL.H>

namespace amrex {

class BoxArray;
class DistributionMapping;
class Geometry;
class MultiFab;

} // namespace amrex

namespace crescent::external_exchange {

struct SurfaceExchangeMetadata
{
    int level = -1;
    int finest_level = -1;
    int step = -1;
    amrex::Real time = 0.0;
};

struct FreeSurfaceExchangeView
{
    SurfaceExchangeMetadata metadata;
    // Non-owning AMReX field views on the solver's current level data.
    amrex::MultiFab& eta;
    amrex::MultiFab& phi;
    amrex::MultiFab& pressure;          // non-const: also the external write target (ERF -> HOS forcing)
    amrex::MultiFab& u;                 // surface velocity x-component
    amrex::MultiFab& v;                 // surface velocity y-component
    amrex::MultiFab& vertical_velocity; // phi_z (surface vertical velocity)
    const amrex::Geometry& geometry;
    const amrex::BoxArray& box_array;
    const amrex::DistributionMapping& distribution_mapping;
};

struct ConstFreeSurfaceExchangeView
{
    SurfaceExchangeMetadata metadata;
    const amrex::MultiFab& eta;
    const amrex::MultiFab& phi;
    const amrex::MultiFab& pressure;
    const amrex::MultiFab& u;
    const amrex::MultiFab& v;
    const amrex::MultiFab& vertical_velocity;
    const amrex::Geometry& geometry;
    const amrex::BoxArray& box_array;
    const amrex::DistributionMapping& distribution_mapping;
};

ConstFreeSurfaceExchangeView asConstFreeSurfaceExchange(
    FreeSurfaceExchangeView view);

class SurfaceExchangeProvider
{
public:
    virtual ~SurfaceExchangeProvider() = default;

    [[nodiscard]] virtual int finestSurfaceExchangeLevel() const = 0;
    [[nodiscard]] virtual FreeSurfaceExchangeView surfaceExchangeView(int lev) = 0;
    [[nodiscard]] virtual ConstFreeSurfaceExchangeView surfaceExchangeView(int lev) const = 0;
};

} // namespace crescent::external_exchange

#endif // CRESCENT_EXTERNAL_EXCHANGE_H
