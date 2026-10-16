function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(build_crescent_lib crescent_lib_name)

    set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
    target_sources(${crescent_lib_name} PRIVATE
        ${SRC_DIR}/solver/CRESCENT.cpp
        ${SRC_DIR}/external/CRESCENT_ExternalExchange.cpp
        ${SRC_DIR}/core/CRESCENT_CoreFreeSurfaceEvolution.cpp
        ${SRC_DIR}/shared/CRESCENT_DNOExecution.cpp
        ${SRC_DIR}/core/CRESCENT_Spectral.cpp
        ${SRC_DIR}/shared/CRESCENT_IO.cpp
        ${SRC_DIR}/solver/CRESCENT_RunLoop.cpp
        ${SRC_DIR}/solver/CRESCENT_RuntimeDiagnostics.cpp
        ${SRC_DIR}/solver/CRESCENT_RuntimeOutputs.cpp
        ${SRC_DIR}/solver/CRESCENT_Checkpoint.cpp
        ${SRC_DIR}/solver/CRESCENT_CheckpointConfigIO.cpp
        ${SRC_DIR}/solver/CRESCENT_StepDiagnostics.cpp
        ${SRC_DIR}/solver/CRESCENT_PostStep.cpp
        ${SRC_DIR}/shared/CRESCENT_Utils.cpp
        ${SRC_DIR}/init/CRESCENT_Init.cpp
        ${SRC_DIR}/solver/CRESCENT_SetupParameters.cpp
        ${SRC_DIR}/solver/CRESCENT_Setup.cpp
        ${SRC_DIR}/solver/CRESCENT_Constructor.cpp
        ${SRC_DIR}/init/CRESCENT_WaveInit.cpp
        ${SRC_DIR}/init/CRESCENT_JonswapInit.cpp
        ${SRC_DIR}/solver/CRESCENT_AMR.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedSetup.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedAlgorithmInputs.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedTimestep.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedTotalField.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedCoreSnapshots.cpp
        ${SRC_DIR}/solver/CRESCENT_TimestepSupport.cpp
        ${SRC_DIR}/core/CRESCENT_Physics.cpp
        ${SRC_DIR}/core/CRESCENT_Dynamics.cpp
        ${SRC_DIR}/shared/CRESCENT_Interpolation.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedDomains.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedMasks.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedRK4.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedCorrectionField.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedEnergy.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedC3InterfaceRepair.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedSurfaceDerivatives.cpp
        ${SRC_DIR}/nested/CRESCENT_NestedFilter.cpp
        ${SRC_DIR}/coupling/CRESCENT_Coupling.cpp
    )

    include(AMReXBuildInfo)
    generate_buildinfo(${crescent_lib_name} ${CMAKE_SOURCE_DIR})
if (${CRESCENT_USE_INTERNAL_AMREX})
    target_include_directories(${crescent_lib_name} PUBLIC ${AMREX_LOCATION}/Tools/C_scripts)
endif()

    target_include_directories(${crescent_lib_name} PUBLIC ${SRC_DIR})
    target_include_directories(${crescent_lib_name} PUBLIC ${SRC_DIR}/core)
    target_include_directories(${crescent_lib_name} PUBLIC ${SRC_DIR}/external)
    target_include_directories(${crescent_lib_name} PUBLIC ${SRC_DIR}/init)
    target_include_directories(${crescent_lib_name} PUBLIC ${SRC_DIR}/nested)
    target_include_directories(${crescent_lib_name} PUBLIC ${SRC_DIR}/shared)
    target_include_directories(${crescent_lib_name} PUBLIC ${SRC_DIR}/solver)
    target_include_directories(${crescent_lib_name} PUBLIC ${CMAKE_BINARY_DIR})

    target_link_libraries_system(${crescent_lib_name} PUBLIC AMReX::amrex)

    target_link_libraries(${crescent_lib_name} PUBLIC MPI::MPI_CXX)

    if(AMReX_GPU_BACKEND STREQUAL "CUDA")
        setup_target_for_cuda_compilation(${crescent_lib_name})
    endif()

endfunction()

function(build_crescent_exe crescent_exe_name crescent_lib_target)

    target_sources(${crescent_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Source/main.cpp)

    target_link_libraries(${crescent_exe_name} PUBLIC ${crescent_lib_target})

    if(AMReX_GPU_BACKEND STREQUAL "CUDA")
        setup_target_for_cuda_compilation(${crescent_exe_name})
    endif()

    install(TARGETS ${crescent_exe_name}
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)

endfunction()
