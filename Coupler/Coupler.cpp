#include "Coupler.H"


// Include CRESCENT header here, not in Coupler.H, so the Coupler header stays free of solver types
#include "CRESCENT.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

namespace coupler {

void InitializeCRESCENTCoupling(
    void*& crescent_solver,
    const amrex::Geometry& erf_geom,
    const amrex::BoxArray& erf_grids,
    bool& use_crescent_coupling,
    amrex::Real& water_density,
    int& crescent_dt_ratio,
    int& crescent_grid_ratio)
{
    // Read coupling parameters from input file
    amrex::ParmParse pp("erf");
    pp.query("crescent_coupling", use_crescent_coupling);

    if (!use_crescent_coupling) {
        amrex::Print() << "\n[ERF-CRESCENT] Coupling disabled\n" << std::endl;
        crescent_solver = nullptr;
        return;
    }

    // Read coupling configuration
    pp.query("water_density", water_density);
    pp.query("crescent_dt_ratio", crescent_dt_ratio);
    pp.query("crescent_grid_ratio", crescent_grid_ratio);

    std::string crescent_input_file = "inputs_crescent";
    pp.query("crescent_input_file", crescent_input_file);

    // Load CRESCENT parameters from its input file into ParmParse database
    amrex::Print() << "[CRESCENT Coupling] Loading CRESCENT parameters from: " << crescent_input_file << "\n";
    amrex::ParmParse::addfile(crescent_input_file.c_str());

    // Get ERF domain info
    auto prob_lo_gpu = erf_geom.ProbLoArray();
    auto prob_hi_gpu = erf_geom.ProbHiArray();
    amrex::Vector<amrex::Real> prob_lo{prob_lo_gpu[0], prob_lo_gpu[1], prob_lo_gpu[2]};
    amrex::Vector<amrex::Real> prob_hi{prob_hi_gpu[0], prob_hi_gpu[1], prob_hi_gpu[2]};

    // Calculate CRESCENT grid size (finer than ERF)
    // Use full domain, not first box (which may be decomposed across MPI ranks)
    const amrex::Box& domain = erf_geom.Domain();
    int nx_erf = domain.length(0);
    int ny_erf = domain.length(1);
    int nz_erf = domain.length(2);
    int nx_crescent = nx_erf * crescent_grid_ratio;
    int ny_crescent = ny_erf * crescent_grid_ratio;

    amrex::Print() << "\n[CRESCENT Coupling] Setting up CRESCENT parameters:\n";
    amrex::Print() << "  CRESCENT Grid: " << nx_crescent << " × " << ny_crescent << " × 1\n";
    amrex::Print() << "  CRESCENT Domain: [" << prob_lo[0] << ", " << prob_lo[1] << ", 0] to ["
                   << prob_hi[0] << ", " << prob_hi[1] << ", 10]\n\n";

    // Strategy B: Use AMReX with per-instance ParmParse prefix support
    // CRESCENT will read from crescent.amr.* and crescent.geometry.*
    // ERF reads from erf.amr.* and erf.geometry.* (or unprefixed for backward compatibility)

    amrex::Print() << "[CRESCENT Coupling] Creating CRESCENT object with prefix 'crescent'...\n";

    // Create CRESCENT with prefix - this will read from crescent.amr.* and crescent.geometry.*
    CRESCENT* crescent = new CRESCENT("crescent");

    // NOTE: CRESCENT reads its parameters in the constructor (via the "crescent" prefix),
    // so there is no separate ReadParameters() step.
    amrex::Print() << "[CRESCENT Coupling] Initializing CRESCENT solver...\n";
    crescent->Init();

    amrex::Print() << "[CRESCENT Coupling] CRESCENT initialization complete!\n";

    // Assertions: Verify CRESCENT got the correct geometry
    const auto& crescent_geom = crescent->Geom(0);
    const auto& crescent_domain = crescent_geom.Domain();

    amrex::Print() << "\n[CRESCENT Coupling] Verifying CRESCENT geometry:\n";
    amrex::Print() << "  Domain: " << crescent_domain << "\n";
    amrex::Print() << "  ProbLo: [" << crescent_geom.ProbLo(0) << ", " << crescent_geom.ProbLo(1) << ", " << crescent_geom.ProbLo(2) << "]\n";
    amrex::Print() << "  ProbHi: [" << crescent_geom.ProbHi(0) << ", " << crescent_geom.ProbHi(1) << ", " << crescent_geom.ProbHi(2) << "]\n";

    // Critical assertions for CRESCENT grid correctness
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        crescent_domain.length(0) == nx_crescent && crescent_domain.length(1) == ny_crescent && crescent_domain.length(2) == 1,
        "CRESCENT grid dimensions mismatch! Expected (" + std::to_string(nx_crescent) + "," +
        std::to_string(ny_crescent) + ",1), got (" + std::to_string(crescent_domain.length(0)) + "," +
        std::to_string(crescent_domain.length(1)) + "," + std::to_string(crescent_domain.length(2)) + ")"
    );

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::abs(crescent_geom.ProbLo(0) - prob_lo[0]) < 1e-10 &&
        std::abs(crescent_geom.ProbLo(1) - prob_lo[1]) < 1e-10 &&
        std::abs(crescent_geom.ProbHi(0) - prob_hi[0]) < 1e-10 &&
        std::abs(crescent_geom.ProbHi(1) - prob_hi[1]) < 1e-10,
        "CRESCENT domain bounds mismatch!"
    );

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        crescent_geom.isPeriodic(0) && crescent_geom.isPeriodic(1) && !crescent_geom.isPeriodic(2),
        "CRESCENT periodicity mismatch! Expected (1,1,0)"
    );

    amrex::Print() << "[CRESCENT Coupling] ✓ CRESCENT geometry verified successfully\n\n";

    // Store as void pointer
    crescent_solver = static_cast<void*>(crescent);

    // Print coupling information
    amrex::Print() << "\n";
    amrex::Print() << "╔═══════════════════════════════════════════════════════╗\n";
    amrex::Print() << "║  ERF-CRESCENT In-Process Coupling Initialized             ║\n";
    amrex::Print() << "╠═══════════════════════════════════════════════════════╣\n";
    amrex::Print() << "║  ERF grid:   " << nx_erf << " × " << ny_erf
                   << " × " << nz_erf << std::endl;
    amrex::Print() << "║  CRESCENT grid:   " << nx_crescent << " × " << ny_crescent << " × 1" << std::endl;
    amrex::Print() << "║  Grid ratio:  " << crescent_grid_ratio << "x" << std::endl;
    amrex::Print() << "║  dt ratio:    " << crescent_dt_ratio << "x (CRESCENT sub-cycles)" << std::endl;
    amrex::Print() << "║  ρ_water:     " << water_density << " kg/m³" << std::endl;
    amrex::Print() << "╚═══════════════════════════════════════════════════════╝\n";
    amrex::Print() << "\n";
}

void DestroyCRESCENTSolver(void*& crescent_solver)
{
    if (crescent_solver != nullptr) {
        CRESCENT* crescent = static_cast<CRESCENT*>(crescent_solver);
        delete crescent;
        crescent_solver = nullptr;
    }
}

void TransferPressureToCRESCENT(
    void* crescent_solver_ptr,
    const amrex::Geometry& erf_geom,
    const amrex::BoxArray& erf_grids,
    const amrex::DistributionMapping& erf_dmap,
    int crescent_grid_ratio,
    const amrex::MultiFab& erf_p_surface,  // 2D surface pressure (already extracted)
    amrex::Real time,
    int crescent_lev)
{
    CRESCENT* crescent_solver = static_cast<CRESCENT*>(crescent_solver_ptr);

    const amrex::IntVect ratio(crescent_grid_ratio, crescent_grid_ratio, 1);

    // Use the working pattern from TransferWavesFromCRESCENT:
    // Step 1: Create refined CRESCENT BoxArray (fine resolution on CRESCENT's index space)
    const auto& crescent_ba = crescent_solver->boxArray(crescent_lev);
    amrex::BoxArray crescent_ba_coarse = crescent_ba;
    crescent_ba_coarse.coarsen(ratio);

    // Diagnostic: print box domains (no IOProcessor check - amrex::Print already handles it)
    amrex::Print() << "  [Coupling] ERF BoxArray minBox: " << erf_grids.minimalBox() << std::endl;
    amrex::Print() << "  [Coupling] CRESCENT BoxArray (fine) minBox: " << crescent_ba.minimalBox() << std::endl;
    amrex::Print() << "  [Coupling] CRESCENT BoxArray (coarse) minBox: " << crescent_ba_coarse.minimalBox() << std::endl;

    // Step 2: Create temporary coarse MultiFab on CRESCENT's DistributionMapping
    const auto& crescent_dm = crescent_solver->DistributionMap(crescent_lev);
    amrex::MultiFab crescent_p_coarse(crescent_ba_coarse, crescent_dm, 1, 0);

    // Initialize to a sentinel value to detect unfilled cells
    crescent_p_coarse.setVal(0.0);

    // Step 3: ParallelCopy ERF coarse data to CRESCENT coarse (handles DM differences)
    crescent_p_coarse.ParallelCopy(erf_p_surface, 0, 0, 1);

    // Diagnostic: check after ParallelCopy (min/max are collective MPI ops - all ranks must participate)
    amrex::Print() << "  [Coupling] After ParallelCopy: crescent_p_coarse min=" << crescent_p_coarse.min(0)
                   << ", max=" << crescent_p_coarse.max(0) << " Pa" << std::endl;

    // Step 4: Refine locally (same DM, so no communication needed)
    amrex::MultiFab crescent_p_fine(crescent_ba, crescent_dm, 1, 0);

    for (amrex::MFIter mfi(crescent_p_coarse); mfi.isValid(); ++mfi) {
        const amrex::Box& coarse_bx = mfi.validbox();
        amrex::Box fine_bx = coarse_bx;
        fine_bx.refine(ratio);

        auto p_coarse = crescent_p_coarse.const_array(mfi);
        auto p_fine = crescent_p_fine.array(mfi);

        // Piecewise constant interpolation (replicate coarse cell to fine cells)
        // Use ParallelFor with proper index mapping
        amrex::ParallelFor(fine_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // Map fine cell index to coarse cell index using AMReX coarsening
            amrex::IntVect fine_idx(i, j, k);
            amrex::IntVect coarse_idx = amrex::coarsen(fine_idx, ratio);

            // Replicate coarse pressure value to fine grid
            p_fine(i, j, k) = p_coarse(coarse_idx[0], coarse_idx[1], 0);
        });
    }

    // Diagnostic: check after refinement (before division) (collective MPI operation)
    amrex::Print() << "  [Coupling] After refinement: crescent_p_fine min=" << crescent_p_fine.min(0)
                   << ", max=" << crescent_p_fine.max(0) << " Pa" << std::endl;

    // Apply physics: divide by water density (convert to pressure head for CRESCENT)
    // CRESCENT evolution equation uses p_atm/ρ_water
    const amrex::Real water_density = 1000.0;  // kg/m³
    crescent_p_fine.mult(1.0 / water_density);

    // Diagnostic: check values BEFORE passing to CRESCENT (collective MPI operation)
    amrex::Print() << "  [Coupling] Passing pressure to CRESCENT at t=" << time
                   << ", p/ρ_w range=[" << crescent_p_fine.min(0)
                   << ", " << crescent_p_fine.max(0) << "] m" << std::endl;

    // Pass to CRESCENT solver
    crescent_solver->setAtmosphericPressure(crescent_p_fine, time, crescent_lev);
}

void AdvanceCRESCENTWithSubcycling(
    void* crescent_solver_ptr,
    amrex::Real dt_erf,
    int crescent_dt_ratio)
{
    CRESCENT* crescent_solver = static_cast<CRESCENT*>(crescent_solver_ptr);
    amrex::Real dt_crescent = dt_erf / static_cast<amrex::Real>(crescent_dt_ratio);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "  [Coupling] Advancing CRESCENT with " << crescent_dt_ratio
                       << " substeps (dt_crescent=" << dt_crescent << ")" << std::endl;
    }

    // Sub-cycle: CRESCENT takes multiple small steps
    for (int i = 0; i < crescent_dt_ratio; ++i) {
        crescent_solver->advanceSingleTimestep(dt_crescent);
    }

    // Check if CRESCENT needs to write a checkpoint
    if (crescent_solver->shouldWriteCheckpoint()) {
        crescent_solver->WriteCheckpointFile();
    }
}

void TransferWavesFromCRESCENT(
    void* crescent_solver_ptr,
    const amrex::Geometry& erf_geom,
    const amrex::BoxArray& erf_grids,
    const amrex::DistributionMapping& erf_dmap,
    int crescent_grid_ratio,
    amrex::MultiFab& erf_eta,
    amrex::MultiFab& erf_u,
    amrex::MultiFab& erf_v,
    amrex::MultiFab& erf_w,
    int crescent_lev)
{
    amrex::Print() << "\n[TransferWavesFromCRESCENT] START - Entering function\n";
    CRESCENT* crescent_solver = static_cast<CRESCENT*>(crescent_solver_ptr);
    amrex::Print() << "[TransferWavesFromCRESCENT] Cast crescent_solver pointer: " << crescent_solver << "\n";

    // Get CRESCENT wave fields (fine grid) - use simple direct accessors
    amrex::Print() << "[TransferWavesFromCRESCENT] About to call getEta(0)...\n";
    const amrex::MultiFab& crescent_eta = crescent_solver->getEta(crescent_lev);
    amrex::Print() << "[TransferWavesFromCRESCENT] getEta(0) succeeded, eta size = " << crescent_eta.size() << "\n";

    amrex::Print() << "[TransferWavesFromCRESCENT] About to call getU(0)...\n";
    const amrex::MultiFab& crescent_u = crescent_solver->getU(crescent_lev);
    amrex::Print() << "[TransferWavesFromCRESCENT] getU(0) succeeded, u size = " << crescent_u.size() << "\n";

    amrex::Print() << "[TransferWavesFromCRESCENT] About to call getV(0)...\n";
    const amrex::MultiFab& crescent_v = crescent_solver->getV(crescent_lev);
    amrex::Print() << "[TransferWavesFromCRESCENT] getV(0) succeeded, v size = " << crescent_v.size() << "\n";

    amrex::Print() << "[TransferWavesFromCRESCENT] About to call getVerticalVelocity(0)...\n";
    const amrex::MultiFab& crescent_w = crescent_solver->getVerticalVelocity(crescent_lev);
    amrex::Print() << "[TransferWavesFromCRESCENT] getVerticalVelocity(0) succeeded, w size = " << crescent_w.size() << "\n";

    const amrex::IntVect ratio(crescent_grid_ratio, crescent_grid_ratio, 1);
    amrex::Print() << "[TransferWavesFromCRESCENT] Coarsening ratio = " << ratio << "\n";

    // Use a simpler approach: Create fine MultiFabs on ERF's DM (refined grid),
    // ParallelCopy CRESCENT data, then average locally

    // Step 1: Create refined ERF BoxArray (fine resolution on ERF's index space)
    amrex::BoxArray erf_ba_fine = erf_grids;
    erf_ba_fine.refine(ratio);
    amrex::Print() << "[TransferWavesFromCRESCENT] Created refined ERF BoxArray, size = " << erf_ba_fine.size() << "\n";
    amrex::Print() << "[TransferWavesFromCRESCENT] ERF fine BA: " << erf_ba_fine << "\n";
    amrex::Print() << "[TransferWavesFromCRESCENT] CRESCENT BA:      " << crescent_eta.boxArray() << "\n";

    // Check if domains match
    const amrex::Box& erf_domain = erf_ba_fine.minimalBox();
    const amrex::Box& crescent_domain = crescent_eta.boxArray().minimalBox();
    amrex::Print() << "[TransferWavesFromCRESCENT] ERF fine domain: " << erf_domain << "\n";
    amrex::Print() << "[TransferWavesFromCRESCENT] CRESCENT domain:      " << crescent_domain << "\n";

    // Step 2: Create temporary fine MultiFabs on ERF's DistributionMapping
    amrex::MultiFab erf_eta_fine(erf_ba_fine, erf_dmap, 1, 0);
    amrex::MultiFab erf_u_fine(erf_ba_fine, erf_dmap, 1, 0);
    amrex::MultiFab erf_v_fine(erf_ba_fine, erf_dmap, 1, 0);
    amrex::MultiFab erf_w_fine(erf_ba_fine, erf_dmap, 1, 0);
    amrex::Print() << "[TransferWavesFromCRESCENT] Created fine MultiFabs on ERF's DM\n";

    // Set to zero first to check if ParallelCopy works
    erf_eta_fine.setVal(0.0);
    erf_u_fine.setVal(0.0);
    erf_v_fine.setVal(0.0);
    erf_w_fine.setVal(0.0);

    // Step 3: ParallelCopy CRESCENT fine data to ERF fine MultiFabs
    amrex::Print() << "[TransferWavesFromCRESCENT] ParallelCopy CRESCENT data to ERF fine grid...\n";
    erf_eta_fine.ParallelCopy(crescent_eta, 0, 0, 1);
    amrex::Print() << "[TransferWavesFromCRESCENT] eta ParallelCopy done, checking values...\n";
    amrex::Print() << "[TransferWavesFromCRESCENT] erf_eta_fine min/max = " << erf_eta_fine.min(0) << " / " << erf_eta_fine.max(0) << "\n";

    erf_u_fine.ParallelCopy(crescent_u, 0, 0, 1);
    erf_v_fine.ParallelCopy(crescent_v, 0, 0, 1);
    erf_w_fine.ParallelCopy(crescent_w, 0, 0, 1);
    amrex::Print() << "[TransferWavesFromCRESCENT] ParallelCopy complete\n";

    // Step 4: Average down locally (same DM, so no communication needed)
    amrex::Print() << "[TransferWavesFromCRESCENT] Averaging down to coarse grid...\n";

    // Disable vectorization for this loop to avoid SIGILL
    #ifdef __clang__
    #pragma clang loop vectorize(disable)
    #endif
    for (amrex::MFIter mfi(erf_eta); mfi.isValid(); ++mfi) {
        const amrex::Box& coarse_bx = mfi.validbox();
        amrex::Box fine_bx = coarse_bx;
        fine_bx.refine(ratio);

        auto eta_f = erf_eta_fine.const_array(mfi);
        auto u_f = erf_u_fine.const_array(mfi);
        auto v_f = erf_v_fine.const_array(mfi);
        auto w_f = erf_w_fine.const_array(mfi);

        auto eta_c = erf_eta.array(mfi);
        auto u_c = erf_u.array(mfi);
        auto v_c = erf_v.array(mfi);
        auto w_c = erf_w.array(mfi);

        const int rx = ratio[0];
        const int ry = ratio[1];
        const amrex::Real inv_num_cells = 1.0 / static_cast<amrex::Real>(rx * ry);

        // Conservative averaging - manually written loop instead of ParallelFor to avoid vectorization issues
        const auto& lo = amrex::lbound(coarse_bx);
        const auto& hi = amrex::ubound(coarse_bx);

        for (int kc = lo.z; kc <= hi.z; ++kc) {
            for (int jc = lo.y; jc <= hi.y; ++jc) {
                for (int ic = lo.x; ic <= hi.x; ++ic) {
                    amrex::Real eta_sum = 0.0;
                    amrex::Real u_sum = 0.0;
                    amrex::Real v_sum = 0.0;
                    amrex::Real w_sum = 0.0;

                    // Average over fine cells
                    for (int jj = 0; jj < ry; ++jj) {
                        for (int ii = 0; ii < rx; ++ii) {
                            int ifine = ic * rx + ii;
                            int jfine = jc * ry + jj;
                            eta_sum += eta_f(ifine, jfine, kc);
                            u_sum += u_f(ifine, jfine, kc);
                            v_sum += v_f(ifine, jfine, kc);
                            w_sum += w_f(ifine, jfine, kc);
                        }
                    }

                    eta_c(ic, jc, kc) = eta_sum * inv_num_cells;
                    u_c(ic, jc, kc) = u_sum * inv_num_cells;
                    v_c(ic, jc, kc) = v_sum * inv_num_cells;
                    w_c(ic, jc, kc) = w_sum * inv_num_cells;
                }
            }
        }
    }
    amrex::Print() << "[TransferWavesFromCRESCENT] Averaging complete\n";

    amrex::Print() << "[TransferWavesFromCRESCENT] Computing min/max for output...\n";
    // All ranks must participate in min/max (collective MPI ops)
    amrex::Print() << "  [Coupling] Transferred waves from CRESCENT (averaged " << ratio[0] << "×" << ratio[1] << "):" << std::endl;
    amrex::Print() << "    eta range: [" << erf_eta.min(0) << ", " << erf_eta.max(0) << "] m" << std::endl;
    amrex::Print() << "    u range:   [" << erf_u.min(0) << ", " << erf_u.max(0) << "] m/s" << std::endl;
    amrex::Print() << "    v range:   [" << erf_v.min(0) << ", " << erf_v.max(0) << "] m/s" << std::endl;
    amrex::Print() << "    w range:   [" << erf_w.min(0) << ", " << erf_w.max(0) << "] m/s" << std::endl;
    amrex::Print() << "[TransferWavesFromCRESCENT] COMPLETED SUCCESSFULLY\n";
}

void ApplyWaveBoundaryCondition(
    const amrex::MultiFab& erf_eta,
    const amrex::MultiFab& erf_u,
    const amrex::MultiFab& erf_v,
    const amrex::MultiFab& erf_w,
    amrex::MultiFab& z_phys_nd,
    const amrex::Geometry& erf_geom,
    int lev)
{
    // Update ERF's z_phys_nd (terrain height) to reflect moving surface
    // This implements ERF's moving terrain capability for wave-following mesh
    //
    // ERF uses terrain-following coordinates where z_phys_nd defines
    // the physical height of each computational grid point.
    // By updating z_phys_nd with wave elevation, ERF's mesh follows the waves.

    const int k_surface = 0;

    for (amrex::MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto z_phys = z_phys_nd.array(mfi);
        auto eta = erf_eta.const_array(mfi);
        auto u_surf = erf_u.const_array(mfi);
        auto v_surf = erf_v.const_array(mfi);
        auto w_surf = erf_w.const_array(mfi);

        // Update terrain height at all vertical levels
        // For terrain-following coordinates, we need to update the entire column
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // Get wave elevation at this horizontal position
            amrex::Real wave_height = eta(i, j, 0);

            if (k == k_surface) {
                // Bottom surface follows wave exactly
                z_phys(i, j, k) = wave_height;
            } else {
                // Upper levels: shift by wave height (terrain-following)
                // This maintains vertical grid structure while following waves
                // NOTE: In production, may need smooth transition function
                amrex::Real base_height = z_phys(i, j, k);
                z_phys(i, j, k) = base_height + wave_height;
            }
        });
    }

    // Surface velocity boundary conditions
    // In ERF's moving terrain formulation, surface velocities are imposed
    // as Dirichlet boundary conditions at k=0
    //
    // NOTE: This requires access to ERF's state arrays (u_bc, v_bc, w_bc)
    // which are not passed to this function. In actual ERF integration,
    // this would be done in ERF's boundary condition filling routine.
    //
    // The proper implementation location is:
    //   ERF/Source/BoundaryConditions/ERF_PhysBCFunct.cpp
    // where wave velocities would be applied during BC fill:
    //
    //   if (use_crescent_coupling && boundary_type == SURFACE) {
    //       u_bc(i, j, k_surface) = crescent_u(i, j);
    //       v_bc(i, j, k_surface) = crescent_v(i, j);
    //       w_bc(i, j, k_surface) = crescent_w(i, j);
    //   }

    // Recompute terrain metrics for moving mesh
    // ERF's terrain-following coordinates require several derived quantities:
    //   - detJ: Jacobian determinant of coordinate transformation
    //   - z_t: ∂z/∂t (time derivative of terrain height)
    //   - Terrain slopes and curvatures
    //
    // NOTE: This requires calling ERF's terrain metric computation
    // In actual integration, would call:
    //   make_terrain_metrics(lev);
    // from within ERF class where this function is available.

    // All ranks must participate in min/max (collective MPI ops)
    amrex::Print() << "  [Coupling] Applied wave boundary condition:" << std::endl;
    amrex::Print() << "    Updated z_phys_nd with wave elevation" << std::endl;
    amrex::Print() << "    Wave height range: [" << erf_eta.min(0)
                   << ", " << erf_eta.max(0) << "] m" << std::endl;
}

} // namespace coupler


