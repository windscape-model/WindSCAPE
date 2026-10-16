#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <iostream>

#include "CRESCENT_NestedC3Dispatch.H"
#include "CRESCENT_NestedAlgorithmInputs.H"
#include "CRESCENT_NestedMasks.H"
#include "CRESCENT_WaveState.H"

namespace {

bool testRepairOptionForwarding()
{
    crescent::nested::AlgorithmConfig config;
    config.c3_boundary_penalty = 0.25;
    config.c3_coefficient_regularization = 3.0;

    crescent::nested::CouplingZoneState coupling_zones;
    coupling_zones.resize(2);
    coupling_zones.c3_buffer_width_x_computed[1] = 17;
    coupling_zones.c3_buffer_width_y_computed[1] = 11;

    const auto options = crescent::nested::c3::makeInterfaceRepairOptions(
        1, config, coupling_zones);
    return options.buffer_width_x == 17 &&
           options.buffer_width_y == 11 &&
           options.boundary_penalty_weight == 0.25 &&
           options.coefficient_regularization_weight == 3.0;
}

bool testCorrectedConsistencyDefaults()
{
    const crescent::nested::AlgorithmConfig config;
    return config.l1_coupling_mode == 1;
}

bool testBlendingMaskClamping()
{
    using namespace amrex;

    const Box domain_box(IntVect(0, 0, 0), IntVect(15, 15, 0));
    const RealBox real_box({0.0, 0.0, 0.0}, {16.0, 16.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 0};
    const Geometry geom(domain_box, &real_box, 0, periodic);
    const BoxArray grids(domain_box);
    const DistributionMapping mapping(grids);
    MultiFab sigma(grids, mapping, 1, 0);

    crescent::nested::domains::DomainStructure domain{};
    domain.x_lo_A = 4.0;
    domain.x_hi_A = 12.0;
    domain.y_lo_A = 4.0;
    domain.y_hi_A = 12.0;
    domain.x_lo_S = 2.0;
    domain.x_hi_S = 14.0;
    domain.y_lo_S = 2.0;
    domain.y_hi_S = 14.0;

    crescent::nested::masks::generateBlendingFunction(sigma, geom, domain);
    Gpu::streamSynchronize();

    const auto values = sigma.const_arrays()[0];
    bool ok = true;
    for (int j = 0; j < 16; ++j) {
        for (int i = 0; i < 16; ++i) {
            const Real x = i + 0.5;
            const Real y = j + 0.5;
            const Real s = values(i, j, 0);
            const bool in_accurate = x > 4.0 && x < 12.0 && y > 4.0 && y < 12.0;
            const bool outside_spectral = x < 2.0 || x > 14.0 || y < 2.0 || y > 14.0;
            ok = ok && s >= 0.0 && s <= 1.0;
            ok = ok && (!in_accurate || s == 0.0);
            ok = ok && (!outside_spectral || s == 1.0);
        }
    }
    return ok;
}

int runTest()
{
    const bool forwarding_ok = testRepairOptionForwarding();
    const bool defaults_ok = testCorrectedConsistencyDefaults();
    const bool mask_ok = testBlendingMaskClamping();
    if (!forwarding_ok || !defaults_ok || !mask_ok) {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::cerr << "nested control test failed: forwarding="
                      << forwarding_ok << ", corrected defaults="
                      << defaults_ok << ", blending mask="
                      << mask_ok << '\n';
        }
        return 1;
    }
    return 0;
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    const int result = runTest();
    amrex::Finalize();
    return result;
}
