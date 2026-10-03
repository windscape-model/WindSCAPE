#include "Prob.H"

using namespace amrex;

namespace crescent {

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* problo, const amrex::Real* probhi)
{
    amrex::Print() << " Initializing problem... " << std::endl;
}

void Problem::init_data(const amrex::Box &bx,
                        const amrex::GeometryData &geomdata,
                        amrex::Array4<amrex::Real> const &eta,
                        amrex::Array4<amrex::Real> const &phi)
{

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
    {
        // Initialize to zero - actual wave field is set by initializeStokesWave() or generateJonswapWaves()
        eta(i, j, k) = 0.0;
        phi(i, j, k) = 0.0;
    });

    // Note: This is called once per box/tile during grid initialization
    // Actual wave initialization happens later in the solver wave-setup path.
}

} // namespace crescent
