cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# ============================================================================
# ERF + CRESCENT coupled superbuild
#
# Builds, in order, each component with its OWN (unmodified) CMake against a
# single shared, vendored AMReX, then the top-level Coupler that links them:
#
#   shared/amrex   -> install (one shared AMReX)
#   models/ERF     -> erf_api  library (coupling-agnostic; built vs external AMReX)
#   models/CRESCENT  -> crescent_static library (vs external AMReX)
#   Coupler/       -> windscape exe (installs the exchange hook into ERF's time loop)
#
# Standalone ERF and CRESCENT executables remain available from their own builds;
# "coupling on/off" = which executable you run.
#
# Usage:  cmake -S . -B build/superbuild && cmake --build build/superbuild
# ============================================================================
project(ERF_CRESCENT_Superbuild NONE)

include(ExternalProject)

set(SB        "${CMAKE_BINARY_DIR}")
set(AMREX_INSTALL "${SB}/install/amrex")
set(ERF_BUILD     "${SB}/erf-build")
set(CRESCENT_BUILD  "${SB}/crescent-build")

set(_build_type "Release" CACHE STRING "Build type for all components")
# Debug aids for the shared AMReX (propagate to all components via the installed
# AMReX_Config.H): assertions and Array4 bound checking. OFF for production.
set(_amrex_assertions  "OFF" CACHE STRING "AMReX_ASSERTIONS for the shared AMReX")
set(_amrex_bound_check "OFF" CACHE STRING "AMReX_BOUND_CHECK for the shared AMReX")

# 1) Shared AMReX (our vendored fork) -> installed once.
ExternalProject_Add(amrex
    SOURCE_DIR  "${CMAKE_SOURCE_DIR}/shared/amrex"
    INSTALL_DIR "${AMREX_INSTALL}"
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${_build_type}
        -DCMAKE_INSTALL_PREFIX=${AMREX_INSTALL}
        -DAMReX_SPACEDIM=3 -DAMReX_PRECISION=DOUBLE
        -DAMReX_MPI=ON -DAMReX_OMP=OFF
        -DAMReX_EB=ON -DAMReX_LINEAR_SOLVERS=ON -DAMReX_LINEAR_SOLVERS_INCFLO=ON
        -DAMReX_PARTICLES=ON -DAMReX_PARTICLES_PRECISION=DOUBLE -DAMReX_PIC=ON
        -DAMReX_FFT=ON -DAMReX_PLOTFILE_TOOLS=ON
        -DAMReX_FORTRAN=OFF -DAMReX_ENABLE_TESTS=OFF
        -DAMReX_ASSERTIONS=${_amrex_assertions} -DAMReX_BOUND_CHECK=${_amrex_bound_check}
    BUILD_ALWAYS OFF
)

# 2) ERF library (coupling-agnostic), built against the shared AMReX. No install
#    (ERF's install pulls in all Exec problem exes); the Coupler links the build tree.
ExternalProject_Add(erf
    DEPENDS amrex
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/models/ERF"
    BINARY_DIR "${ERF_BUILD}"
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${_build_type}
        -DERF_USE_INTERNAL_AMREX=OFF
        -DAMReX_DIR=${AMREX_INSTALL}/lib/cmake/AMReX
        -DERF_ENABLE_MPI=ON -DERF_ENABLE_FFT=ON
        -DERF_DIM=3 -DERF_PRECISION=DOUBLE
        -DERF_ENABLE_TERRAIN_VELOCITY=ON
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target erf_api
    INSTALL_COMMAND ""
    BUILD_ALWAYS OFF
)

# 3) CRESCENT static library, built against the shared AMReX. No install.
ExternalProject_Add(crescent
    DEPENDS amrex
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/models/CRESCENT"
    BINARY_DIR "${CRESCENT_BUILD}"
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${_build_type}
        -DCRESCENT_USE_INTERNAL_AMREX=OFF
        -DAMReX_DIR=${AMREX_INSTALL}/lib/cmake/AMReX
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target crescent_static
    INSTALL_COMMAND ""
    BUILD_ALWAYS OFF
)

# 4) The Coupler executable -- links the prebuilt ERF + CRESCENT libs + shared AMReX.
ExternalProject_Add(coupler
    DEPENDS erf crescent
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/Coupler"
    BINARY_DIR "${SB}/coupler-build"
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${_build_type}
        -DAMREX_INSTALL=${AMREX_INSTALL}
        -DERF_BUILD=${ERF_BUILD}
        -DCRESCENT_BUILD=${CRESCENT_BUILD}
    INSTALL_COMMAND ""
    BUILD_ALWAYS ON
)
