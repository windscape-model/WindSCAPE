cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(HOS CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(CMakePackageConfigHelpers)

option(HOS_USE_INTERNAL_AMREX "Add AMReX as submodule" ON)

# CRESCENT engine requires C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (${HOS_USE_INTERNAL_AMREX})
    # AMREX_LOCATION is overridable from the command line (e.g. to build the
    # bitwise A/B gate against the same vanilla AMReX commit CRESCENT used).
    if (NOT DEFINED AMREX_LOCATION)
        set(AMREX_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX")
    endif()

    include(${CMAKE_SOURCE_DIR}/cmake/SetAmrexOptions.cmake)
    list(APPEND CMAKE_MODULE_PATH "${AMREX_LOCATION}/Tools/CMake")

    add_subdirectory(${AMREX_LOCATION} ${CMAKE_BINARY_DIR}/amrex)
else()
    # Superbuild / external AMReX: link a pre-installed AMReX
    # (configure with -DHOS_USE_INTERNAL_AMREX=OFF -DAMReX_DIR=<install>/lib/cmake/AMReX).
    if (NOT TARGET AMReX::amrex)
        find_package(AMReX CONFIG REQUIRED)
    endif()
endif()

find_package(MPI REQUIRED)

add_subdirectory(Exec)

set(HOS_PUBLIC_HEADERS
    Source/shared/ProbBase.H
    Source/CRESCENT.H
    Source/external/CRESCENT_ExternalExchange.H
    Source/core/CRESCENT_Spectral.H
    Source/shared/CRESCENT_IO.H
    Source/shared/CRESCENT_Utils.H
    Source/init/CRESCENT_Init.H
    Source/core/CRESCENT_Physics.H
    Source/core/CRESCENT_Dynamics.H
    Source/shared/CRESCENT_Interpolation.H)
set_target_properties(crescent_srclib PROPERTIES PUBLIC_HEADER "${HOS_PUBLIC_HEADERS}")
