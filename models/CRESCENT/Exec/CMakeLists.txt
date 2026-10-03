set(crescent_lib_name crescent_srclib)
add_library(${crescent_lib_name} OBJECT)

include(${CMAKE_SOURCE_DIR}/cmake/BuildCRESCENT.cmake)
build_crescent_lib(${crescent_lib_name})

add_executable(crescent "")
target_sources(crescent PRIVATE Prob.cpp)
target_include_directories(crescent PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
build_crescent_exe(crescent ${crescent_lib_name})

# Static library for external linking (e.g. the top-level Coupler): the CRESCENT
# engine objects plus the CRESCENT problem (crescent::amrex_probinit from Prob.cpp).
add_library(crescent_static STATIC $<TARGET_OBJECTS:${crescent_lib_name}>)
target_sources(crescent_static PRIVATE Prob.cpp)
target_include_directories(crescent_static PUBLIC
    ${CMAKE_SOURCE_DIR}/Source
    ${CMAKE_SOURCE_DIR}/Source/core
    ${CMAKE_SOURCE_DIR}/Source/external
    ${CMAKE_SOURCE_DIR}/Source/init
    ${CMAKE_SOURCE_DIR}/Source/nested
    ${CMAKE_SOURCE_DIR}/Source/shared
    ${CMAKE_SOURCE_DIR}/Source/solver
    ${CMAKE_SOURCE_DIR}/Source/coupling)
target_link_libraries(crescent_static PUBLIC AMReX::amrex MPI::MPI_CXX)
