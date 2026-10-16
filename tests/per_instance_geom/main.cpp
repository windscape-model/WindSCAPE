// Multi-level per-instance geometry test.
//
// Constructs TWO AmrCore instances (prefixes "a" and "b") with DIFFERENT physical
// domains and periodicity, each with max_level >= 1, and prints each instance's
// physical domain at level 0 AND the refined level(s). With the per-instance
// geometry fix, instance "a" must use a.geometry.* at every level and instance "b"
// must use b.geometry.* -- proving the AMR levels of each model stay within that
// model's own boundaries, independent of the process-wide Geometry singleton.

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_Print.H>

class TestAmr : public amrex::AmrCore
{
public:
    explicit TestAmr (std::string const& prefix) : amrex::AmrCore(prefix) {}
    // Pure-virtual AmrCore hooks -- unused here (we only inspect geom[], set in the ctor).
    void ErrorEst (int, amrex::TagBoxArray&, amrex::Real, int) override {}
    void MakeNewLevelFromScratch (int, amrex::Real, const amrex::BoxArray&, const amrex::DistributionMapping&) override {}
    void MakeNewLevelFromCoarse (int, amrex::Real, const amrex::BoxArray&, const amrex::DistributionMapping&) override {}
    void RemakeLevel (int, amrex::Real, const amrex::BoxArray&, const amrex::DistributionMapping&) override {}
    void ClearLevel (int) override {}
};

static void report (const char* name, amrex::AmrCore& m)
{
    for (int lev = 0; lev <= m.maxLevel(); ++lev) {
        const auto& g = m.Geom(lev);
        amrex::Print() << "  " << name << " Geom(" << lev << "): ProbHi=["
                       << g.ProbHi(0) << ", " << g.ProbHi(1) << ", " << g.ProbHi(2) << "]"
                       << "  isPeriodic=(" << g.isPeriodic(0) << "," << g.isPeriodic(1)
                       << "," << g.isPeriodic(2) << ")"
                       << "  domain=" << g.Domain() << "\n";
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        TestAmr a("a");
        TestAmr b("b");
        amrex::Print() << "\n=== Per-instance multi-level geometry ===\n";
        report("a", a);
        report("b", b);
        amrex::Print() << "=== expect: a -> a.geometry.* at ALL levels; b -> b.geometry.* (different) ===\n";
    }
    amrex::Finalize();
    return 0;
}
