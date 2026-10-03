#pragma once

// Primary solver declaration header.

#include <memory>
#include <string>

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include "CRESCENT_ExternalExchange.H"

namespace crescent::ownership {
struct SolverOwnership;
}

class CRESCENT : public amrex::AmrCore,
               public crescent::external_exchange::SurfaceExchangeProvider
{
public:
    CRESCENT();

    // Strategy-B constructor: per-instance ParmParse prefix so an in-process
    // instance (e.g. ERF coupling) reads <prefix>.amr.*, <prefix>.geometry.*
    // and <prefix>.* solver parameters instead of the default "crescent" prefix.
    explicit CRESCENT(std::string const& pp_prefix);

    ~CRESCENT() override;

    CRESCENT(CRESCENT&&) noexcept = delete;
    CRESCENT& operator=(CRESCENT&& other) noexcept = delete;
    CRESCENT(const CRESCENT& other) = delete;
    CRESCENT& operator=(const CRESCENT& other) = delete;

    void Init();
    void Run();

    // ===== ERF-CRESCENT coupling interface =====
    // ERF -> CRESCENT "receive": write atmospheric pressure forcing into the surface pressure
    // field. The no-level form targets level 0 (the level coupled to ERF); the level-aware
    // overload lets a driver target a specific (e.g. refined) level.
    void setAtmosphericPressure(const amrex::MultiFab& p_atm, amrex::Real time);
    void setAtmosphericPressure(const amrex::MultiFab& p_atm, amrex::Real time, int lev);
    // Externally driven single wave timestep (ERF sub-cycles this).
    void advanceSingleTimestep(amrex::Real dt);
    void setManualGeometryFlag(bool flag) { manual_geometry_setup = flag; }
    // HOS -> ERF: surface field accessors (level lev).
    [[nodiscard]] const amrex::MultiFab& getEta(int lev) const;
    [[nodiscard]] const amrex::MultiFab& getU(int lev) const;
    [[nodiscard]] const amrex::MultiFab& getV(int lev) const;
    [[nodiscard]] const amrex::MultiFab& getVerticalVelocity(int lev) const;
    // Checkpoint helpers used by the ERF coupling driver.
    [[nodiscard]] bool shouldWriteCheckpoint() const;
    void WriteCheckpointFile();

    [[nodiscard]] int finestSurfaceExchangeLevel() const override;
    [[nodiscard]] crescent::external_exchange::FreeSurfaceExchangeView surfaceExchangeView(int lev) override;
    [[nodiscard]] crescent::external_exchange::ConstFreeSurfaceExchangeView surfaceExchangeView(int lev) const override;

    // ------------------------------------
    // Members to override from AmrCore:
    // ------------------------------------
    // Tag cells for refinement
    void ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time, int ngrow) override;

    void MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray &ba,
                                const amrex::DistributionMapping &dm) override;
    void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray &ba,
                     const amrex::DistributionMapping &dm) override;
    void ClearLevel(int lev) override;
    void MakeNewLevelFromScratch(int lev, amrex::Real time, const amrex::BoxArray &ba,
                                 const amrex::DistributionMapping &dm) override;
    // ------------------------------------

private:
    [[nodiscard]] crescent::ownership::SolverOwnership& solverOwnership() noexcept;
    [[nodiscard]] const crescent::ownership::SolverOwnership& solverOwnership() const noexcept;

    // Shared construction body used by both the default and prefix constructors.
    void construct();

    std::unique_ptr<crescent::ownership::SolverOwnership> solver_ownership;

    // Solver-parameter ParmParse prefix (default "crescent"; "<prefix>" under Strategy B).
    std::string m_solver_pp_prefix = "crescent";
    // Coupling flags (parity with the legacy coupled HOS; inert for standalone runs).
    bool use_atmospheric_pressure = false;
    bool manual_geometry_setup = false;
};
