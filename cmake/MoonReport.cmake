# Runs only via the explicitly requested spp_moon_report target, never during
# ordinary configure, library builds, CLI builds or CTest discovery.
cmake_minimum_required(VERSION 3.28)

foreach(required SPP_LUNAR_CHECK SPP_SKY_CHECK SPP_CACHE_CHECK SPP_BENCH SPP_EARTH_BENCH SPP_EVENT_BENCH SPP_PHASE_BENCH SPP_ROTATION_BENCH SPP_INTEGRATOR_BENCH SPP_SPK_BENCH SPP_REFERENCE_DIR)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
        message(FATAL_ERROR "Moon report: missing ${required}")
    endif()
endforeach()

set(report_failed FALSE)
function(report_stage label executable)
    message(STATUS "\n${label}")
    execute_process(COMMAND "${executable}" ${ARGN} RESULT_VARIABLE result)
    if(NOT "${result}" STREQUAL "0")
        message(STATUS "${label}: FAILED (${result})")
        set(report_failed TRUE PARENT_SCOPE)
    endif()
endfunction()

message(STATUS "Moon report: accuracy and throughput are separate measurements.")
set(stages 9)
if(SPP_MOON_REFERENCE_KERNEL)
    set(stages 10)
endif()
report_stage("1/${stages} ELP implementation accuracy, both parameter fits" "${SPP_LUNAR_CHECK}"
    "${SPP_REFERENCE_DIR}/elp_mpp02_states.csv")
report_stage("2/${stages} Interpolation residuals and cache contracts" "${SPP_CACHE_CHECK}")
report_stage("3/${stages} Independent JPL sky and USNO event accuracy" "${SPP_SKY_CHECK}"
    "${SPP_REFERENCE_DIR}/moon_horizons.csv"
    "${SPP_REFERENCE_DIR}/moon_usno_phases.csv")
report_stage("4/${stages} Earth cache residency, eviction and speed" "${SPP_EARTH_BENCH}")
report_stage("5/${stages} Adaptive horizon events: numerical gates and speed" "${SPP_EVENT_BENCH}")
report_stage("6/${stages} Matched-TT independent phase accuracy and search cost" "${SPP_PHASE_BENCH}"
    "${SPP_REFERENCE_DIR}/moon_horizons_phase_tt.csv")
report_stage("7/${stages} IAU orientation, rotating states and independent DE441 comparison" "${SPP_ROTATION_BENCH}"
    "${SPP_REFERENCE_DIR}/moon_iau_rotation.csv" "${SPP_REFERENCE_DIR}/moon_rotation_horizons.csv")
report_stage("8/${stages} Direct, cached and cache-miss Moon throughput" "${SPP_BENCH}"
    --moon-only --samples=128 --rounds=5 --fit=both)
if(SPP_MOON_REFERENCE_KERNEL)
    report_stage("9/${stages} SPK reader vs independent CSPICE and throughput" "${SPP_SPK_BENCH}"
        "${SPP_MOON_REFERENCE_KERNEL}" "${SPP_REFERENCE_DIR}/spk_de440s.csv")
    report_stage("10/${stages} Adaptive integrator: analytic controls and lunar propagation vs DE440s" "${SPP_INTEGRATOR_BENCH}"
        "${SPP_MOON_REFERENCE_KERNEL}")
else()
    report_stage("9/${stages} Adaptive integrator: analytic controls and dense-output throughput" "${SPP_INTEGRATOR_BENCH}")
    message(STATUS "SPK/lunar propagation comparisons omitted: SPP_MOON_REFERENCE_KERNEL is unset.")
endif()

if(report_failed)
    message(FATAL_ERROR "Moon report FAILED; all available stages were attempted. Read their output above.")
endif()
message(STATUS "Moon accuracy checks PASS; speed measurements have no arbitrary pass/fail threshold.")
