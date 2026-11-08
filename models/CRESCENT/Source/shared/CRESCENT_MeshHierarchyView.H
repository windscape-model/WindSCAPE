#ifndef CRESCENT_MESH_HIERARCHY_VIEW_H
#define CRESCENT_MESH_HIERARCHY_VIEW_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_Vector.H>

#include "CRESCENT_LevelMeshView.H"

namespace crescent::shared {

struct MeshHierarchyView
{
    const amrex::Vector<amrex::BoxArray>& grids;
    const amrex::Vector<amrex::DistributionMapping>& dmap;
    const amrex::Vector<amrex::Geometry>& geom;
    const amrex::Vector<amrex::IntVect>& ref_ratio;

    LevelMeshView levelMesh(int lev) const
    {
        return {grids[lev], dmap[lev], geom[lev]};
    }

    amrex::IntVect refinementRatio(int lev) const
    {
        if (lev == 0 || lev > static_cast<int>(ref_ratio.size())) {
            return amrex::IntVect::TheUnitVector();
        }
        return ref_ratio[lev - 1];
    }
};

} // namespace crescent::shared

#endif // CRESCENT_MESH_HIERARCHY_VIEW_H
