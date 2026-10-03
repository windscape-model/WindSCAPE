#include "CRESCENT_Utils.H"

#include <AMReX_MultiFab.H>

namespace crescent::utils {

using amrex::MultiFab;
using amrex::ParallelCopy;

void gatherToRoot(const MultiFab& mf, MultiFab& root_mf)
{
    // send MF on all ranks to rank 0 - NOTE: root_mf and mf have different BA and DM
    // root_mf and mf are expected to be defined before this call
    ParallelCopy(root_mf, mf, 0, 0, 1);
}

} // namespace crescent::utils
