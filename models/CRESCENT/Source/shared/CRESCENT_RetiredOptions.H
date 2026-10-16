#ifndef CRESCENT_RETIRED_OPTIONS_H
#define CRESCENT_RETIRED_OPTIONS_H

#include <sstream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>

namespace crescent::shared {

// Options retired by the 2026-10 code review. Each keeps only its former
// default, which a deck may still state; any other value stops the run. The
// legacy and exploratory branches are available at commit f39cacf and
// earlier.
template <typename T>
void requireRetiredValue(amrex::ParmParse& pp, const char* name, T supported)
{
    T value = supported;
    if (pp.query(name, value) && !(value == supported)) {
        std::ostringstream message;
        message << "crescent." << name << " was retired; the only supported "
                << "value is " << supported << " (the other values are "
                << "available at commit f39cacf and earlier)";
        amrex::Abort(message.str());
    }
}

} // namespace crescent::shared

#endif // CRESCENT_RETIRED_OPTIONS_H
