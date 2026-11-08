cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(ERFCRESCENTCoupler CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Locations of the prebuilt peer components (overridable by the superbuild) ---
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
if(NOT DEFINED AMREX_INSTALL)
    set(AMREX_INSTALL "${_root}/build/install/amrex")
endif()
if(NOT DEFINED ERF_BUILD)
    set(ERF_BUILD "${_root}/build/erf")
endif()
if(NOT DEFINED CRESCENT_BUILD)
    set(CRESCENT_BUILD "${_root}/build/crescent")
endif()
if(NOT DEFINED ERF_SRC)
    set(ERF_SRC "${_root}/models/ERF/Source")
endif()
if(NOT DEFINED CRESCENT_SRC)
    set(CRESCENT_SRC "${_root}/models/CRESCENT/Source")
endif()
# The ERF "problem" (provides the global amrex_probinit + custom terrain). Selectable
# like an ERF Exec case; defaults to the ERF-CRESCENT coupled wave test (CRESCENTWaves).
if(NOT DEFINED ERF_PROBLEM_DIR)
    set(ERF_PROBLEM_DIR "${_root}/models/ERF/Exec/DevTests/ERF_CRESCENT_Test")
endif()

find_package(AMReX REQUIRED PATHS "${AMREX_INSTALL}/lib/cmake/AMReX" NO_DEFAULT_PATH)
find_package(MPI REQUIRED)

# --- ERF: prebuilt static library (built coupling-free vs the shared AMReX) ---
add_library(erf_imported STATIC IMPORTED)
set_target_properties(erf_imported PROPERTIES IMPORTED_LOCATION "${ERF_BUILD}/liberf_api.a")
# ERF include dirs = every directory under Source that holds a header.
file(GLOB_RECURSE _erf_headers "${ERF_SRC}/*.H")
set(ERF_INC_DIRS "")
foreach(_h ${_erf_headers})
    get_filename_component(_d "${_h}" DIRECTORY)
    list(APPEND ERF_INC_DIRS "${_d}")
endforeach()
list(REMOVE_DUPLICATES ERF_INC_DIRS)

# --- CRESCENT: prebuilt static library (engine + coupling interface) ---
add_library(crescent_imported STATIC IMPORTED)
set_target_properties(crescent_imported PROPERTIES IMPORTED_LOCATION "${CRESCENT_BUILD}/Exec/libcrescent_static.a")

# --- The coupler executable ---
add_executable(coupler_app main.cpp Coupler.cpp
    "${ERF_PROBLEM_DIR}/ERF_Prob.cpp")
set_target_properties(coupler_app PROPERTIES OUTPUT_NAME windscape)
# Match ERF's compile defines so ERF.H has the same layout (AMREX_SPACEDIM/NDEBUG come from AMReX::amrex).
target_compile_definitions(coupler_app PRIVATE ERF_USE_FFT ERF_USE_MOISTURE)
target_include_directories(coupler_app PRIVATE
    ${ERF_INC_DIRS}
    "${ERF_PROBLEM_DIR}"
    "${CRESCENT_SRC}" "${CRESCENT_SRC}/core" "${CRESCENT_SRC}/external" "${CRESCENT_SRC}/init"
    "${CRESCENT_SRC}/nested" "${CRESCENT_SRC}/shared" "${CRESCENT_SRC}/solver" "${CRESCENT_SRC}/coupling")
# Link order: ERF + CRESCENT (depend on AMReX) then AMReX (leaf).
target_link_libraries(coupler_app PRIVATE erf_imported crescent_imported AMReX::amrex MPI::MPI_CXX)
# ERF/AMReX reference amrex::buildInfoGet*(); generate the build-info TU for this exe.
# Must come AFTER target_link_libraries: GNU ld scans static archives left to right, so the
# buildinfo archive has to follow liberf_api.a on the link line to satisfy the references its
# objects introduce (macOS ld64 is order-insensitive, which masked this locally).
include(AMReXBuildInfo)
generate_buildinfo(coupler_app "${_root}")
