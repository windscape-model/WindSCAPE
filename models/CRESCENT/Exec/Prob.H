#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#include "ProbBase.H"

// Namespaced under `crescent` so the concrete problem class does not collide with
// ERF's own global `Problem` when the CRESCENT engine is compiled in-process.
namespace crescent {

class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* problo, const amrex::Real* probhi);

    void init_data(const amrex::Box &bx,
                const amrex::GeometryData &geomdata,
                amrex::Array4<amrex::Real> const &eta,
                amrex::Array4<amrex::Real> const &phi) override;

};

} // namespace crescent

#endif
