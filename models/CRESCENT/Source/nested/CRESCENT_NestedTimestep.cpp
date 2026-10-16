#include "CRESCENT_NestedTimestep.H"

#include "CRESCENT_Dynamics.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedCorrectionField.H"
#include "CRESCENT_NestedTotalField.H"
#include "CRESCENT_SolverContext.H"

namespace crescent::nested::free_surface_evolution {

void advanceL1Timestep(crescent::shared::LevelContext& level_context,
                       amrex::Real time)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    const bool correction_field =
        level_context.nestedConfig().l1_coupling_mode == 2;

    const amrex::Real c[4] = {0.0, 0.5, 0.5, 1.0};

    for (int irk = 0; irk < 4; ++irk) {
        crescent::dynamics::advanceSurfaceRK4(
            solution.eta[lev],
            solution.phi[lev],
            solution.eta_saved[lev],
            solution.phi_saved[lev],
            solution.rhs_eta[lev],
            solution.rhs_phi[lev],
            runtime.dt[lev],
            irk);

        if (correction_field) {
            // Correction-field coupling: no periodic copies and no C3
            // repair; the whole stage runs on compact fields.
            crescent::nested::correction_field::runL1CorrectionFieldStage(
                level_context,
                irk,
                time + c[irk] * runtime.dt[lev],
                c[irk]);
        } else {
            crescent::nested::total_field::runStage(level_context, irk);
        }
    }

    crescent::dynamics::updateSurface(
        solution.eta[lev],
        solution.phi[lev],
        solution.eta_saved[lev],
        solution.phi_saved[lev],
        solution.rhs_eta[lev],
        solution.rhs_phi[lev],
        runtime.dt[lev]);

    if (correction_field) {
        // The correction taper is the CCF sponge and drain.
        crescent::nested::correction_field::
            applyPostRK4CorrectionFieldTaper(level_context);
    } else {
        crescent::nested::total_field::applyStateRelaxation(level_context);
    }
}

void recomputeL1OutputDerivatives(crescent::shared::LevelContext& level_context)
{
    if (!crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            level_context.level())) {
        return;
    }
    if (level_context.nestedConfig().l1_coupling_mode == 2) {
        crescent::nested::correction_field::
            recomputeCorrectionFieldOutputDerivatives(level_context);
    } else {
        crescent::nested::total_field::prepareOutputDerivatives(level_context);
    }
}

} // namespace crescent::nested::free_surface_evolution
