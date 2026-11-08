add_executable(
    crescent_test_linear_dynamics
    test_linear_dynamics.cpp)
target_link_libraries(crescent_test_linear_dynamics PRIVATE crescent_srclib)

add_test(
    NAME crescent.linear_dynamics_rhs
    COMMAND crescent_test_linear_dynamics)

add_executable(
    crescent_test_energy
    test_energy.cpp)
target_link_libraries(crescent_test_energy PRIVATE crescent_srclib)

add_test(
    NAME crescent.zakharov_energy
    COMMAND crescent_test_energy)

add_test(
    NAME crescent.rejects_ptf_without_c3_repair
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_ptf_no_c3.inputs)
set_tests_properties(
    crescent.rejects_ptf_without_c3_repair
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.use_c3_repair was retired")

add_test(
    NAME crescent.rejects_unimplemented_cell_center_c3
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_ptf_cell_center_c3.inputs)
set_tests_properties(
    crescent.rejects_unimplemented_cell_center_c3
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.c3_use_interface was retired")

add_test(
    NAME crescent.rejects_invalid_stokes_timestep
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_stokes_zero_dt.inputs)
set_tests_properties(
    crescent.rejects_invalid_stokes_timestep
    PROPERTIES
        PASS_REGULAR_EXPRESSION "non-JONSWAP runs require crescent.dt to be a positive finite number")

add_test(
    NAME crescent.rejects_invalid_spectral_filter
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_spectral_filter.inputs)
set_tests_properties(
    crescent.rejects_invalid_spectral_filter
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.beta1_filter must be a positive finite number")

add_test(
    NAME crescent.rejects_derived_wave_input
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_derived_wave_input.inputs)
set_tests_properties(
    crescent.rejects_derived_wave_input
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.alpha is derived")

add_test(
    NAME crescent.rejects_invalid_nested_roi
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_nested_roi.inputs)
set_tests_properties(
    crescent.rejects_invalid_nested_roi
    PROPERTIES
        PASS_REGULAR_EXPRESSION "0 < roi_spectral_fraction < roi_accurate_fraction < 0.5")

add_test(
    NAME crescent.rejects_invalid_jonswap_dt_override
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_jonswap_zero_dt_override.inputs)
set_tests_properties(
    crescent.rejects_invalid_jonswap_dt_override
    PROPERTIES
        PASS_REGULAR_EXPRESSION "explicit crescent.dt_override.*must be a positive finite number")

add_test(
    NAME crescent.rejects_quiet_nonperiodic_geometry
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_quiet_nonperiodic_geometry.inputs)
set_tests_properties(
    crescent.rejects_quiet_nonperiodic_geometry
    PROPERTIES
        PASS_REGULAR_EXPRESSION "requires periodic geometry in x and y")

add_test(
    NAME crescent.rejects_quiet_multilayer_geometry
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_quiet_multilayer_geometry.inputs)
set_tests_properties(
    crescent.rejects_quiet_multilayer_geometry
    PROPERTIES
        PASS_REGULAR_EXPRESSION "requires exactly one z cell")

add_test(
    NAME crescent.rejects_unsafe_c3_buffer
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_small_ptf_c3_buffer.inputs)
set_tests_properties(
    crescent.rejects_unsafe_c3_buffer
    PROPERTIES
        PASS_REGULAR_EXPRESSION "C3 repair buffers.*overlap Omega_A or exceed the child line")

add_test(
    NAME crescent.rejects_c3_harmonics_above_nyquist
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_c3_nyquist_harmonics.inputs)
set_tests_properties(
    crescent.rejects_c3_harmonics_above_nyquist
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.c3_x_transverse_harmonics was retired")

add_test(
    NAME crescent.rejects_unqualified_ccf_adjustment
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_ccf_adjustment_without_consistency.inputs)
set_tests_properties(
    crescent.rejects_unqualified_ccf_adjustment
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.ccf_adjustment_consistency was retired")

add_test(
    NAME crescent.rejects_invalid_ccf_transform_admissibility
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_ccf_transform_admissibility.inputs)
set_tests_properties(
    crescent.rejects_invalid_ccf_transform_admissibility
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.ccf_transform_admissibility was retired")

add_test(
    NAME crescent.rejects_retired_ccf_band_restrict
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_ccf_band_restrict.inputs)
set_tests_properties(
    crescent.rejects_retired_ccf_band_restrict
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.ccf_band_restrict was retired")

add_test(
    NAME crescent.rejects_invalid_ccf_w1_projection_consistency
    COMMAND $<TARGET_FILE:crescent>
            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/invalid_ccf_w1_projection_consistency.inputs)
set_tests_properties(
    crescent.rejects_invalid_ccf_w1_projection_consistency
    PROPERTIES
        PASS_REGULAR_EXPRESSION "crescent.ccf_w1_projection_consistency was retired")

add_executable(
    crescent_test_checkpoint_config
    test_checkpoint_config.cpp)
target_link_libraries(crescent_test_checkpoint_config PRIVATE crescent_srclib)

add_test(
    NAME crescent.checkpoint_config_roundtrip
    COMMAND crescent_test_checkpoint_config)

add_executable(
    crescent_test_bilinear_interpolation
    test_bilinear_interpolation.cpp)
target_link_libraries(crescent_test_bilinear_interpolation PRIVATE crescent_srclib)

add_test(
    NAME crescent.bilinear_cell_centering
    COMMAND crescent_test_bilinear_interpolation)

add_executable(
    crescent_test_c3_least_squares
    test_c3_least_squares.cpp)
target_link_libraries(crescent_test_c3_least_squares PRIVATE crescent_srclib)

add_test(
    NAME crescent.c3_least_squares_residual
    COMMAND crescent_test_c3_least_squares)

add_executable(
    crescent_test_nested_control_forwarding
    test_nested_control_forwarding.cpp)
target_link_libraries(crescent_test_nested_control_forwarding PRIVATE crescent_srclib)

add_test(
    NAME crescent.nested_control_forwarding
    COMMAND crescent_test_nested_control_forwarding)

add_executable(
    crescent_test_jonswap_spectrum
    test_jonswap_spectrum.cpp)
target_link_libraries(crescent_test_jonswap_spectrum PRIVATE crescent_srclib)

add_test(
    NAME crescent.jonswap_spectrum_construction
    COMMAND crescent_test_jonswap_spectrum)

add_executable(
    crescent_test_wave_initialization_controls
    test_wave_initialization_controls.cpp)
target_link_libraries(crescent_test_wave_initialization_controls PRIVATE crescent_srclib)

add_test(
    NAME crescent.wave_initialization_controls
    COMMAND crescent_test_wave_initialization_controls)
