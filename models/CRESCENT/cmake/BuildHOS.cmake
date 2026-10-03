function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(build_hos_lib hos_lib_name)

    set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
    target_sources(${hos_lib_name} PRIVATE
        ${SRC_DIR}/HOS.cpp
        ${SRC_DIR}/HOS_Spectral.cpp
        ${SRC_DIR}/HOS_IO.cpp
        ${SRC_DIR}/HOS_Utils.cpp
        ${SRC_DIR}/HOS_Init.cpp
        ${SRC_DIR}/HOS_Physics.cpp
        ${SRC_DIR}/HOS_Dynamics.cpp
    )

    include(AMReXBuildInfo)
    generate_buildinfo(${hos_lib_name} ${CMAKE_SOURCE_DIR})
if (${HOS_USE_INTERNAL_AMREX})
    target_include_directories(${hos_lib_name} PUBLIC ${AMREX_LOCATION}/Tools/C_scripts)
endif()

    target_include_directories(${hos_lib_name} PUBLIC ${SRC_DIR})
    target_include_directories(${hos_lib_name} PUBLIC ${CMAKE_BINARY_DIR})

    target_link_libraries_system(${hos_lib_name} PUBLIC AMReX::amrex)

    target_link_libraries(${hos_lib_name} PUBLIC MPI::MPI_CXX)

endfunction()

function(build_hos_exe hos_exe_name)

    target_sources(${hos_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Source/main.cpp)

    target_link_libraries(${hos_exe_name} PUBLIC ${hos_lib_name})

    install(TARGETS ${hos_exe_name}
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)

endfunction()
