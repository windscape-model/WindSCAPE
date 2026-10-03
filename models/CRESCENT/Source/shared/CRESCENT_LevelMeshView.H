#ifndef CRESCENT_LEVEL_MESH_VIEW_H
#define CRESCENT_LEVEL_MESH_VIEW_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

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

inline LevelMeshView makeLevelMeshView(
    const amrex::Vector<amrex::BoxArray>& grids,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    const amrex::Vector<amrex::Geometry>& geom,
    int lev)
{
    return {grids[lev], dmap[lev], geom[lev]};
}

} // namespace crescent::shared

#endif // CRESCENT_LEVEL_MESH_VIEW_H
