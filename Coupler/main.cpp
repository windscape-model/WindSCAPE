#include <iostream>
#include <cstring>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include "ERF.H"
#include "Coupler.H"

using namespace amrex;

// ERF references a file-global inputs_name (used in plotfile job info).
std::string inputs_name;

// ERF's AMReX amr defaults -- replicated so the coupled run reproduces standalone ERF
// (refine_grid_layout (1,1,0), n_proper=2, large max_grid_size, blocking_factor=1).
static void add_par ()
{
    ParmParse pp("amr");
    pp.add("refine_grid_layout_x",1);
    pp.add("refine_grid_layout_y",1);
    pp.add("refine_grid_layout_z",0);
    pp.add("n_proper",2);
    int max_grid_size = 2048; pp.queryAdd("max_grid_size",max_grid_size);
    int blocking_factor = 1;  pp.queryAdd("blocking_factor",blocking_factor);
    int n_error_buf = 0;      pp.queryAdd("n_error_buf",n_error_buf);
}

/**
 * ERF + CRESCENT coupled driver.
 *
 * The Coupler owns the master time loop and drives ERF and CRESCENT as peers:
 *   for each coarse step, ERF::advanceOneCoarseStep runs ComputeDt, then (via the
 *   pre-advance hook) the Coupler pushes ERF surface pressure to CRESCENT, sub-cycles
 *   CRESCENT, and pulls wave fields back; then ERF takes its dynamics step.
 * Standalone ERF (erf exe) and CRESCENT (crescent exe) remain available separately.
 */
int main (int argc, char* argv[])
{
#if defined(AMREX_USE_MPI)
    MPI_Init(&argc, &argv);
#endif

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <inputs_erf> [args...]\n";
        return 1;
    }

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, add_par);

    if (!std::strchr(argv[1], '=')) { inputs_name = argv[1]; }

    BL_PROFILE_VAR("main()", pmain);
    const Real strt_total = amrex::second();
    {
        // --- ERF (atmosphere) ---
        ERF erf;
        erf.InitData();

        // --- CRESCENT (waves): create + initialize via the coupling layer ---
        void* crescent_solver  = nullptr;
        bool  use_coupling   = true;
        Real  water_density  = 1000.0;
        int   crescent_dt_ratio   = 5;
        int   crescent_grid_ratio = 5;
        coupler::InitializeCRESCENTCoupling(
            crescent_solver, erf.Geom(0), erf.boxArray(0),
            use_coupling, water_density, crescent_dt_ratio, crescent_grid_ratio);

        // Wave fields on ERF's grid (computed from CRESCENT each step). Not fed back into
        // ERF's terrain here -- matches the prior coupling (one-way ERF<-CRESCENT diagnostics
        // plus ERF->CRESCENT pressure forcing). Terrain feedback is a later step.
        MultiFab crescent_eta(erf.boxArray(0), erf.DistributionMap(0), 1, 0);
        MultiFab crescent_u  (erf.boxArray(0), erf.DistributionMap(0), 1, 0);
        MultiFab crescent_v  (erf.boxArray(0), erf.DistributionMap(0), 1, 0);
        MultiFab crescent_w  (erf.boxArray(0), erf.DistributionMap(0), 1, 0);
        crescent_eta.setVal(0.0); crescent_u.setVal(0.0); crescent_v.setVal(0.0); crescent_w.setVal(0.0);

        // --- Coupling exchange, injected into ERF's time loop each coarse step ---
        // ERF::Evolve owns the time loop and calls this hook after ComputeDt, before timeStep
        // (the exact point needed for the wave -> moving-terrain feedback). The Coupler drives
        // the exchange: ERF surface pressure -> CRESCENT, sub-cycle CRESCENT, waves -> ERF terrain.
        erf.setCouplingHook(
            [&] (Real cur_time, Real dt_coarse)
            {
                if (!use_coupling || crescent_solver == nullptr) { return; }
                MultiFab p_surf = erf.extractSurfacePressure(0);
                coupler::TransferPressureToCRESCENT(
                    crescent_solver, erf.Geom(0), erf.boxArray(0), erf.DistributionMap(0),
                    crescent_grid_ratio, p_surf, cur_time, 0);
                coupler::AdvanceCRESCENTWithSubcycling(crescent_solver, dt_coarse, crescent_dt_ratio);
                coupler::TransferWavesFromCRESCENT(
                    crescent_solver, erf.Geom(0), erf.boxArray(0), erf.DistributionMap(0),
                    crescent_grid_ratio, crescent_eta, crescent_u, crescent_v, crescent_w, 0);
                // Two-way feedback: push the wave surface (end-of-coarse-step time) into ERF's
                // moving terrain so the fitted mesh follows the waves.
                erf.setWaveElevation(crescent_eta, cur_time + dt_coarse, 0);
            });

        // Run the coupled simulation (ERF drives the loop; the hook does the exchange each step).
        erf.Evolve();

        coupler::DestroyCRESCENTSolver(crescent_solver);

        Real end_total = amrex::second() - strt_total;
        ParallelDescriptor::ReduceRealMax(end_total, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "\nTotal Time: " << end_total << '\n';
    }
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
#if defined(AMREX_USE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
