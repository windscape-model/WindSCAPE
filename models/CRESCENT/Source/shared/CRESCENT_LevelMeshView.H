#ifndef CRESCENT_LEVEL_MESH_VIEW_H
#define CRESCENT_LEVEL_MESH_VIEW_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>

namespace crescent::shared {

struct LevelMeshView
{
    const amrex::BoxArray& grids;
    const amrex::DistributionMapping& dmap;
    const amrex::Geometry& geom;

    amrex::GeometryData domain() const
    {
        return geom.data();
    }
};

} // namespace crescent::shared

#endif // CRESCENT_LEVEL_MESH_VIEW_H
