# Explicit opt-in CPU comparison; never included by an ordinary build.
cmake_minimum_required(VERSION 3.28)

get_filename_component(source_dir "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(kernel "${source_dir}/context/Spp Moon/reference-data/de440s.bsp")
if(NOT EXISTS "${kernel}")
    message(FATAL_ERROR "The ten-stage comparison requires the existing kernel: ${kernel}. No download was attempted.")
endif()

if(CMAKE_HOST_WIN32)
    if("$ENV{LOCALAPPDATA}" STREQUAL "")
        message(FATAL_ERROR "LOCALAPPDATA is unset; cannot select the local Windows build directory.")
    endif()
    file(TO_CMAKE_PATH "$ENV{LOCALAPPDATA}/SPP/moon-build" build_dir)
    set(machine windows)
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
    set(build_dir "${source_dir}/build/moon-Linux")
    set(machine nas)
else()
    message(FATAL_ERROR "This comparison currently supports Windows and Linux.")
endif()

set(log_dir "${source_dir}/context/Spp Moon/performance/comparison")
file(MAKE_DIRECTORY "${log_dir}")
string(TIMESTAMP stamp "%Y%m%d-%H%M%S" UTC)
set(log_file "${log_dir}/${machine}-${stamp}.txt")
cmake_host_system_information(RESULT host_name QUERY HOSTNAME)
cmake_host_system_information(RESULT processor QUERY PROCESSOR_NAME)
file(WRITE "${log_file}"
    "Machine: ${host_name}\nOS: ${CMAKE_HOST_SYSTEM_NAME}\nProcessor: ${processor}\nBuild: ${build_dir}\nRelease; SPP_NATIVE=OFF; report uses one benchmark thread.\n")
message(STATUS "Build: ${build_dir}")
message(STATUS "Log: ${log_file}")

function(run_stage label)
    message(STATUS "${label}")
    file(APPEND "${log_file}" "\n${label}\n")
    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY "${source_dir}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output ERROR_VARIABLE output
        ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE
    )
    file(APPEND "${log_file}" "${output}\nExit status: ${result}\n")
    if(NOT "${result}" STREQUAL "0")
        message(FATAL_ERROR "${label} failed (${result}). Log: ${log_file}")
    endif()
endfunction()

run_stage("Configure identical portable Release settings"
    "${CMAKE_COMMAND}" --preset moon -S "${source_dir}" -B "${build_dir}"
    -DCMAKE_BUILD_TYPE=Release -DSPP_BUILD_CAPI=OFF -DSPP_NATIVE=OFF
    "-DSPP_MOON_REFERENCE_KERNEL=${kernel}")
run_stage("Build benchmark executables"
    "${CMAKE_COMMAND}" --build "${build_dir}" --target spp_bench
    spp_lunar_check spp_lunar_cache_check spp_moon_sky_check
    spp_earth_cache_bench spp_moon_event_bench spp_moon_phase_bench
    spp_moon_rotation_bench spp_spk_bench spp_integrator_bench --parallel 4)
run_stage("Run the same ten-stage accuracy and CPU performance report"
    "${CMAKE_COMMAND}" --build "${build_dir}" --target spp_moon_report)
message(STATUS "Comparison report saved: ${log_file}")
