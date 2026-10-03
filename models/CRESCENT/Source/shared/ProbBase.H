#ifndef PROBBASE_H_
#define PROBBASE_H_

#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

// Namespaced under `crescent` so that, when the CRESCENT engine is compiled in-process
// with ERF (which has its own global ProblemBase / amrex_probinit), there is no
// ODR collision or duplicate-symbol clash at link time.
namespace crescent {

class ProblemBase
{
public:
    virtual ~ProblemBase() = default;

    virtual void init_data(const amrex::Box &bx,
                           const amrex::GeometryData &geomdata,
                           amrex::Array4<amrex::Real> const &eta,
                           amrex::Array4<amrex::Real> const &phi) {

                           }
};

extern std::unique_ptr<ProblemBase> amrex_probinit(const amrex_real* problo,
                                                   const amrex_real* probhi) AMREX_ATTRIBUTE_WEAK;

} // namespace crescent

#endif
