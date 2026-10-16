#ifndef CRESCENT_ROOT_LAYOUTS_H
#define CRESCENT_ROOT_LAYOUTS_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>

namespace crescent::shared {

struct RootLayouts
{
    amrex::Vector<amrex::BoxArray> ba;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<amrex::BoxArray> cba;
    amrex::Vector<amrex::DistributionMapping> cdmap;

    void resize(int nlevs_max)
    {
        ba.resize(nlevs_max);
        dmap.resize(nlevs_max);
        cba.resize(nlevs_max);
        cdmap.resize(nlevs_max);
    }

    void defineLevel(int lev, const amrex::Geometry& geom)
    {
        ba[lev] = amrex::BoxArray(geom.Domain());
        amrex::Vector<int> root_map(1, amrex::ParallelDescriptor::IOProcessorNumber());
        dmap[lev] = amrex::DistributionMapping(ba[lev]);
        dmap[lev].define(root_map);
    }
};

} // namespace crescent::shared

#endif // CRESCENT_ROOT_LAYOUTS_H
