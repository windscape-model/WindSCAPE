#include "solver/CRESCENT_RuntimeOutputs.H"

// Production runtime output implementation.

#include <cmath>

#include <AMReX_MultiFab.H>
#include <AMReX_NFiles.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

#include "CRESCENT_IO.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_RuntimeReportingContext.H"

using amrex::Box;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Print;
using amrex::Real;
using amrex::TilingIfNotGPU;
using amrex::Vector;

namespace {

const amrex::Vector<std::string>& plotVariableNames()
{
    static const amrex::Vector<std::string> plot_variable_names = {
        "eta",
        "eta_x",
        "eta_y",
        "phi",
        "phi_x",
        "phi_y",
        "rhs_eta",
        "rhs_phi",
        "pressure",
        "xvel",
        "yvel",
        "zvel"
    };
    return plot_variable_names;
}

void clearTinyPlotfileValues(amrex::MultiFab& mf, amrex::Real zero_tol)
{
    const int ncomp = mf.nComp();
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto arr = mf.array(mfi);
        ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            if (std::abs(arr(i, j, k, n)) < zero_tol) {
                arr(i, j, k, n) = 0.0;
            }
        });
    }
}

void writeCRESCENTMultilevelTecplot(
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    auto& solution = reporting_context.solutionState();
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    const auto& root_layouts = reporting_context.rootLayouts();
    const auto& geom = reporting_context.meshHierarchy().geom;
    const int finest_level = reporting_context.finestLevel();
    const auto& coupling_zones = reporting_context.nestedRuntime().coupling_zones;
    const auto& coupling_derivatives =
        reporting_context.nestedRuntime().coupling_derivatives;

    const crescent::io::NestedTecplotDiagnosticLevels diagnostic_levels = {
        {
            solution.eta_x,
            solution.eta_y,
            solution.phi_x,
            solution.phi_y,
            solution.phi_z
        },
        {
            coupling_derivatives.l1_blended_derivatives.eta_x,
            coupling_derivatives.l1_blended_derivatives.eta_y,
            coupling_derivatives.l1_blended_derivatives.phi_x,
            coupling_derivatives.l1_blended_derivatives.phi_y,
            coupling_derivatives.l1_blended_derivatives.phi_z
        },
        coupling_zones.chi_D_mask,
        coupling_zones.sigma_blend,
        coupling_derivatives.core_derivatives_on_l1.eta_x,
        coupling_derivatives.core_derivatives_on_l1.phi_x
    };

    crescent::io::writeTecplotMultiLevelCRESCENT(
        reporting.plot_file,
        runtime.istep[0],
        finest_level,
        geom,
        solution.eta,
        solution.phi,
        solution.u,
        solution.v,
        root_layouts.ba,
        root_layouts.dmap,
        runtime.t_new[0],
        diagnostic_levels);

    crescent::io::writeTecplotMovie(
        reporting.plot_file,
        runtime.istep[0],
        finest_level,
        geom,
        solution.eta,
        root_layouts.ba,
        root_layouts.dmap,
        runtime.t_new[0],
        coupling_zones.domains);
}

void writeStandardMultilevelTecplot(
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    auto& solution = reporting_context.solutionState();
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    const auto& root_layouts = reporting_context.rootLayouts();
    const auto& geom = reporting_context.meshHierarchy().geom;

    crescent::io::writeTecplotMultiLevel(
        reporting.plot_file,
        runtime.istep[0],
        reporting_context.finestLevel(),
        geom,
        solution.eta,
        solution.u,
        solution.v,
        root_layouts.ba,
        root_layouts.dmap,
        runtime.t_new[0]);
}

void writeSingleLevelTecplot(
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    auto& solution = reporting_context.solutionState();
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    const auto& root_layouts = reporting_context.rootLayouts();
    const auto& geom = reporting_context.meshHierarchy().geom;
    const int lev = 0;

    crescent::io::writeTecplot(
        reporting.plot_file,
        runtime.istep[0],
        geom[lev],
        solution.eta[lev],
        solution.u[lev],
        solution.v[lev],
        root_layouts.ba[lev],
        root_layouts.dmap[lev],
        runtime.t_new[lev]);
}

} // namespace

namespace crescent::solver::runtime_outputs {

void writePlotFile(const crescent::solver::RuntimeReportingContext& reporting_context)
{
    auto& solution = reporting_context.solutionState();
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    const auto& nested = reporting_context.nestedConfig();
    const auto mesh = reporting_context.meshHierarchy();
    const auto& grids = mesh.grids;
    const auto& dmap = mesh.dmap;
    const auto& geom = mesh.geom;
    const auto& ref_ratio = mesh.ref_ratio;
    const auto& varnames = plotVariableNames();
    const int finest_level = reporting_context.finestLevel();
    std::string plotfilename =
        amrex::Concatenate(reporting.plot_file, runtime.istep[0], 5);

    const int mf_out_ncomp = varnames.size();
    if (mf_out_ncomp == 0) {
        return;
    }

    Vector<MultiFab> mf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], mf_out_ncomp, 0);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Plotfiles are direct SolutionState snapshots, not the nested diagnostic view.
        int mf_comp = 0;
        // Plotfile fields follow the fixed order declared by plotVariableNames().
        MultiFab::Copy(mf[lev], solution.eta[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.eta_x[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.eta_y[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.phi[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.phi_x[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.phi_y[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.rhs_eta[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.rhs_phi[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.pressure[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.u[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.v[lev], 0, mf_comp++, 1, 0);
        MultiFab::Copy(mf[lev], solution.phi_z[lev], 0, mf_comp++, 1, 0);

        AMREX_ALWAYS_ASSERT(mf_comp == mf_out_ncomp);

        if (runtime.istep[0] == 0 && nested.use_crescent_algorithm == 1 && lev > 0) {
            // Keep initialization plotfiles bitwise stable without mutating solver state.
            constexpr Real init_plot_zero_tol = 1.0e-27;
            clearTinyPlotfileValues(mf[lev], init_plot_zero_tol);
        }
    }

    Vector<int> istep_vec(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        istep_vec[lev] = runtime.istep[lev];
    }

    crescent::io::writePlotfile(
        reporting.plot_file,
        runtime.istep[0],
        finest_level,
        varnames,
        mf,
        geom,
        runtime.t_new[0],
        istep_vec,
        ref_ratio);
}

void writeTecplotFile(const crescent::solver::RuntimeReportingContext& reporting_context)
{
    const int finest_level = reporting_context.finestLevel();

    if (finest_level == 0) {
        writeSingleLevelTecplot(reporting_context);
        return;
    }

    if (crescent::nested::hasActiveNestedL1(
            reporting_context.nestedConfig(),
            finest_level)) {
        writeCRESCENTMultilevelTecplot(reporting_context);
        return;
    }

    writeStandardMultilevelTecplot(reporting_context);
}

} // namespace crescent::solver::runtime_outputs
