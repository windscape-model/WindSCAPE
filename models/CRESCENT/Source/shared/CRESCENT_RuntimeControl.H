#ifndef CRESCENT_RUNTIME_CONTROL_H
#define CRESCENT_RUNTIME_CONTROL_H

#include <limits>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace crescent::shared {

struct RuntimeReportingState
{
    std::string plot_file = "";
    int plot_int = -1; // number of steps to write plotfiles
    int tecplot_int = -1; // number of steps to write Tecplot files (-1 = disabled)
    int netcdf_int = -1;  // number of steps to write NetCDF files (-1 = disabled, for future use)
    int hdf5_int = -1;    // number of steps to write HDF5 files (-1 = disabled, for future use)
    int write_omega_a_edge_diagnostics = 0; // 0=off, 1=write y-variation edge diagnostics for L1
    int omega_a_edge_diagnostic_int = 1;    // write every N timesteps when enabled
    amrex::Real energy_initial_total = 0.0; // Baseline total energy for crescent_energy.csv percent change
    bool energy_tracking_initialized = false;
    bool omega_a_stage_diagnostics_initialized = false;
    bool omega_a_edge_diagnostics_initialized = false;
};

struct RuntimeCheckpointState
{
    std::string file_prefix = "crescent_chk";
    int interval = -1; // number of steps between checkpoints; -1 disables checkpointing
    std::string restart_file = "";
    bool restarted = false;
};

struct RuntimeControlState
{
    int max_step = std::numeric_limits<int>::max();
    amrex::Real start_time = 0.0;
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    int verbosity = 0;  // 0=quiet, 1=minimal, 2=diagnostic (default: quiet for performance)

    amrex::Vector<int> istep;
    amrex::Vector<int> nsubsteps;
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;
    amrex::Real dt_override = -1.0;
    // Reporting owns runtime outputs and reporting-only bookkeeping, leaving the
    // top-level state focused on time integration control.
    RuntimeReportingState reporting;
    RuntimeCheckpointState checkpoint;

    void resizeLevelState(int nlevs_max)
    {
        istep.resize(nlevs_max, 0);
        nsubsteps.resize(nlevs_max, 1);
        t_new.resize(nlevs_max, 0.0);
        t_old.resize(nlevs_max, -1.0e100);
        dt.resize(nlevs_max, 1.0e100);
    }
};

} // namespace crescent::shared

#endif // CRESCENT_RUNTIME_CONTROL_H
