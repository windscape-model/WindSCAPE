#!/usr/bin/env bash

# WindSCAPE CI test script (runs inside the per-machine node allocation; see
# .gitlab-ci.yml and .gitlab/LC/runners/<machine>.yml). Implements the tiers of
# Documentation/CI.md:
#   Tier 0: superbuild + coupled / nested / orbital / sub-stepping / standalone-CRESCENT
#           smokes, CRESCENT's CTest suite, and the per-instance geometry unit test
#   Tier 1: 20-step coupled regression vs golden plotfile (only if
#           WINDSCAPE_GOLD_REPO is set; otherwise skipped with a notice)
# Style follows ERF's .gitlab/LC/gitlab_test.sh (env-driven, srun-aware).

set -o errexit
set -o nounset
set -o pipefail

echo "Start: $(date)"

echo "==========="
echo "WindSCAPE CI"
echo "==========="

modules=${MODULE_LIST:-""}
mpiexec_executable=${MPIEXEC_EXECUTABLE:-"srun"}
if [[ "${mpiexec_executable}" == "flux" ]]
then
    mpiexec_executable="$(which ${mpiexec_executable}) run"
else
    mpiexec_executable="$(which ${mpiexec_executable})"
fi
mpiexec_preflags=${MPIEXEC_PREFLAGS:-"--overlap"}
mpirun_1="${mpiexec_executable} ${mpiexec_preflags} -n 1"

host=$(hostname)
build_type=${BUILD_TYPE:-"Release"}

echo "HOST: ${host}"
src_dir="${PWD}"
echo "Source directory: ${src_dir}"
build_dir="$(realpath -- "${src_dir}/../windscape_build_${host}_${CI_PIPELINE_ID:-local}_$(date +%F_%H_%M_%S)")"
echo "Build directory: ${build_dir}"

echo "============="
echo "Setup modules"
echo "============="
if [[ -n ${modules} ]]
then
    module load ${modules}
fi
module list 2>&1 || true

echo "==============="
echo "Build WindSCAPE"
echo "==============="
time cmake -S "${src_dir}" -B "${build_dir}" -D_build_type="${build_type}"
time cmake --build "${build_dir}" -j "${CMAKE_BUILD_PARALLEL_LEVEL:-16}"

windscape_exe="${build_dir}/coupler-build/windscape"
crescent_exe="${build_dir}/crescent-build/Exec/crescent"
fcompare_exe="${build_dir}/install/amrex/bin/amrex_fcompare"
test -x "${windscape_exe}"
test -x "${fcompare_exe}"
# The superbuild builds crescent_static; build the standalone exe and CRESCENT's test
# executables on top of it.
time cmake --build "${build_dir}/crescent-build" -j "${CMAKE_BUILD_PARALLEL_LEVEL:-16}"
test -x "${crescent_exe}"

deck_dir="${src_dir}/models/ERF/Exec/DevTests/ERF_CRESCENT_Test"
run_root="${build_dir}/ci-runs"
mkdir -p "${run_root}"

echo "======================"
echo "Tier 0: coupled smoke"
echo "======================"
mkdir -p "${run_root}/coupled-smoke" && cd "${run_root}/coupled-smoke"
cp "${deck_dir}/inputs_erf" "${deck_dir}/inputs_crescent" .
${mpirun_1} "${windscape_exe}" inputs_erf max_step=3 erf.plot_int_1=-1 | tee smoke.log
grep -q "Coarse STEP 3 ends" smoke.log
[[ $(grep -c "Transferred waves from CRESCENT" smoke.log) -eq 3 ]]
! grep -qiE "nan|abort" smoke.log
echo "coupled smoke OK"

echo "============================="
echo "Tier 0: coupled nested smoke"
echo "============================="
mkdir -p "${run_root}/nested-smoke" && cd "${run_root}/nested-smoke"
cp "${deck_dir}/inputs_erf" "${deck_dir}/inputs_crescent_nested" .
${mpirun_1} "${windscape_exe}" inputs_erf erf.crescent_input_file=inputs_crescent_nested \
    max_step=5 erf.plot_int_1=-1 | tee nested.log
grep -q "Coarse STEP 5 ends" nested.log
echo "nested smoke OK"

echo "=============================="
echo "Tier 0: orbital moving-wall BC"
echo "=============================="
mkdir -p "${run_root}/orbital-smoke" && cd "${run_root}/orbital-smoke"
cp "${deck_dir}/inputs_erf_orbital" "${deck_dir}/inputs_crescent" .
${mpirun_1} "${windscape_exe}" inputs_erf_orbital max_step=20 erf.plot_int_1=-1 | tee orbital.log
grep -q "Coarse STEP 20 ends" orbital.log
echo "orbital smoke OK"

echo "================================"
echo "Tier 0: coupled sub-stepping smoke"
echo "================================"
mkdir -p "${run_root}/substep-smoke" && cd "${run_root}/substep-smoke"
cp "${deck_dir}/inputs_erf" "${deck_dir}/inputs_crescent" .
${mpirun_1} "${windscape_exe}" inputs_erf max_step=5 erf.plot_int_1=-1 \
    erf.substepping_type=Implicit erf.fixed_mri_dt_ratio=2 erf.use_lagged_delta_rt=false | tee substep.log
grep -q "Coarse STEP 5 ends" substep.log
echo "sub-stepping smoke OK"

echo "=================================="
echo "Tier 0: standalone CRESCENT smoke"
echo "=================================="
mkdir -p "${run_root}/crescent-smoke" && cd "${run_root}/crescent-smoke"
${mpirun_1} "${crescent_exe}" "${src_dir}/models/CRESCENT/test_case/inputs_stokes_test" \
    max_step=3 crescent.plot_int=-1 crescent.tecplot_int=-1 | tee crescent.log
grep -q "Done timestep 3" crescent.log
echo "standalone CRESCENT smoke OK"

echo "==========================="
echo "Tier 0: CRESCENT CTest suite"
echo "==========================="
# The unit tests and the input-contract rejection decks of the standalone repo
# (models/CRESCENT/tests), built against the shared AMReX. The rejection decks run the
# crescent executable serially through ctest's own launcher.
(cd "${build_dir}/crescent-build" && ctest --output-on-failure)
echo "CRESCENT CTest suite OK"

echo "===================================="
echo "Tier 0: per-instance geometry test"
echo "===================================="
pig_build="${build_dir}/per-instance-geom"
cmake -S "${src_dir}/tests/per_instance_geom" -B "${pig_build}" \
    -DCMAKE_BUILD_TYPE="${build_type}" \
    -DAMReX_DIR="${build_dir}/install/amrex/lib/cmake/AMReX"
cmake --build "${pig_build}" -j "${CMAKE_BUILD_PARALLEL_LEVEL:-16}"
mkdir -p "${run_root}/pig" && cd "${run_root}/pig"
${mpirun_1} "${pig_build}/per_instance_geom" "${src_dir}/tests/per_instance_geom/inputs" | tee pig.log
echo "per-instance geometry OK"

echo "==============================="
echo "Tier 1: coupled regression gate"
echo "==============================="
if [[ -n "${WINDSCAPE_GOLD_REPO:-}" ]]
then
    cd "${run_root}"
    rm -rf windscape-gold-files
    # Job-token clone (grant this repo job-token permission to the gold repo,
    # mirroring ERF's erf-llnl-gold-files setup)
    git clone --branch "${WINDSCAPE_GOLD_GIT_REF:-main}" --depth 1 \
        "https://gitlab-ci-token:${CI_JOB_TOKEN}@${WINDSCAPE_GOLD_REPO}" windscape-gold-files
    mkdir -p "${run_root}/regression" && cd "${run_root}/regression"
    cp "${deck_dir}/inputs_erf" "${deck_dir}/inputs_crescent" .
    ${mpirun_1} "${windscape_exe}" inputs_erf max_step=20 erf.plot_int_1=20 | tee reg.log
    "${fcompare_exe}" "${run_root}/windscape-gold-files/${CI_MACHINE:-dane}/plt00020" plt00020
    echo "coupled regression OK"
else
    echo "WINDSCAPE_GOLD_REPO not set: skipping the golden-plotfile gate"
fi

echo "All WindSCAPE CI tiers passed."
echo "End: $(date)"
