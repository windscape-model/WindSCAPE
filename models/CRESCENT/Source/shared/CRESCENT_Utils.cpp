#include "CRESCENT_Utils.H"

#include <AMReX_MultiFab.H>

namespace crescent::utils {

using amrex::MultiFab;
using amrex::ParallelCopy;

void gatherToRoot(const MultiFab& mf, MultiFab& root_mf)
{
    // Copy mf from all ranks into root_mf, which lives on rank 0 with its own
    // box array and distribution; both must be defined.
    ParallelCopy(root_mf, mf, 0, 0, 1);
}

} // namespace crescent::utils
