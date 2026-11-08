#pragma once

#include <AMReX_MultiFab.H>
namespace crescent::utils {

/**
 * @brief Gather MultiFab from all processors to root processor
 *
 * @param mf MultiFab distributed across all processors
 * @param root_mf MultiFab on root processor
 */
void gatherToRoot(const amrex::MultiFab& mf, amrex::MultiFab& root_mf);

} // namespace crescent::utils
