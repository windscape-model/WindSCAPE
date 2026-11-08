cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(PerInstanceGeomTest CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
if(NOT DEFINED AMREX_INSTALL)
    set(AMREX_INSTALL "${_root}/build/install/amrex")
endif()

find_package(AMReX REQUIRED PATHS "${AMREX_INSTALL}/lib/cmake/AMReX" NO_DEFAULT_PATH)
find_package(MPI REQUIRED)

add_executable(pig_app main.cpp)
set_target_properties(pig_app PROPERTIES OUTPUT_NAME per_instance_geom)
include(AMReXBuildInfo)
generate_buildinfo(pig_app "${_root}")
target_link_libraries(pig_app PRIVATE AMReX::amrex MPI::MPI_CXX)
