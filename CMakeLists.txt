cmake_minimum_required(VERSION 3.28)

project(Spp
    VERSION 1.2.0
    DESCRIPTION "High-precision solar positional astronomy library"
    LANGUAGES CXX
)

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------
# Apps, tools and install rules default to ON only when Spp is the top-level
# project, so add_subdirectory / FetchContent users get just the library.
option(SPP_NATIVE           "Tune for the build machine (-march=native)"                 OFF)
option(SPP_BUILD_APPS       "Build the CLI and available OpenGL viewers"                  ${PROJECT_IS_TOP_LEVEL})
option(SPP_BUILD_VALIDATION "Build the spp_cache_check accuracy tool"                    ${PROJECT_IS_TOP_LEVEL})
option(SPP_BUILD_BENCH      "Build the spp_bench benchmark"                              OFF)
option(SPP_BUILD_FIGURES    "Build spp_figures, which renders the README figures"        OFF)
option(SPP_INSTALL          "Generate install rules and the CMake package"               ${PROJECT_IS_TOP_LEVEL})
option(SPP_BUILD_CAPI       "Build spp_capi, the C interface as a shared library (spp.dll / spp_capi.so)" OFF)
option(SPP_BUILD_MOON       "Build the ELP/MPP02 lunar ephemeris, sky coordinates and phases" ON)
set(SPP_MOON_REFERENCE_KERNEL "" CACHE FILEPATH "Optional existing DE440s SPK for numerical Moon performance/accuracy gates; never downloaded by CMake")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)
include(cmake/SppFlags.cmake)

# ------------------------------------------------------------------------------
# VSOP87D tables: data/VSOP87D.ear is compiled into the library at build time
# ------------------------------------------------------------------------------
# vsop_gen runs on the build machine. Natively, and when cross-compiling with
# CMAKE_CROSSCOMPILING_EMULATOR (which CMake prepends to the command), it is
# built with the project's compiler. When cross-compiling without an emulator,
# it is either given (SPP_VSOP_GEN_EXECUTABLE) or built with a host compiler
# that accepts GCC/Clang options (SPP_HOST_CXX). Its output does not depend on
# the compiler: decimals are converted with correctly rounded std::from_chars.
set(SPP_VSOP_GEN_EXECUTABLE "" CACHE FILEPATH "Prebuilt host vsop_gen, for cross-compiling without an emulator")
set(SPP_HOST_CXX "c++" CACHE STRING "Host C++20 compiler that builds vsop_gen when cross-compiling without an emulator")

if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    if(SPP_VSOP_GEN_EXECUTABLE)
        set(SPP_VSOP_GEN ${SPP_VSOP_GEN_EXECUTABLE})
        set(SPP_VSOP_GEN_DEPENDS ${SPP_VSOP_GEN_EXECUTABLE})
    else()
        set(SPP_VSOP_GEN ${PROJECT_BINARY_DIR}/host/vsop_gen${CMAKE_HOST_EXECUTABLE_SUFFIX})
        set(SPP_VSOP_GEN_DEPENDS ${SPP_VSOP_GEN})
        add_custom_command(
            OUTPUT ${SPP_VSOP_GEN}
            COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/host
            COMMAND ${SPP_HOST_CXX} -std=c++20 -O2 ${PROJECT_SOURCE_DIR}/tools/vsop_gen.cpp -o ${SPP_VSOP_GEN}
            DEPENDS ${PROJECT_SOURCE_DIR}/tools/vsop_gen.cpp
            COMMENT "Building vsop_gen for the host with ${SPP_HOST_CXX}"
            VERBATIM
        )
    endif()
else()
    add_executable(vsop_gen tools/vsop_gen.cpp)
    spp_apply_flags(vsop_gen)
    set(SPP_VSOP_GEN vsop_gen)
    set(SPP_VSOP_GEN_DEPENDS vsop_gen)
endif()

set(SPP_GENERATED_DIR ${PROJECT_BINARY_DIR}/generated)
set(SPP_VSOP_INC ${SPP_GENERATED_DIR}/spp/vsop87d_earth_terms.inc)

add_custom_command(
    OUTPUT ${SPP_VSOP_INC}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${SPP_GENERATED_DIR}/spp
    COMMAND ${SPP_VSOP_GEN} ${PROJECT_SOURCE_DIR}/data/VSOP87D.ear ${SPP_VSOP_INC}
    DEPENDS ${SPP_VSOP_GEN_DEPENDS} ${PROJECT_SOURCE_DIR}/data/VSOP87D.ear
    COMMENT "Generating VSOP87D Earth tables"
    VERBATIM
)

# ------------------------------------------------------------------------------
# Library
# ------------------------------------------------------------------------------
add_library(spp STATIC
    src/core/calendar.cpp
    src/core/frames.cpp
    src/core/time.cpp
    src/earth.cpp
    src/earth_state.cpp
    src/events.cpp
    src/observe.cpp
    src/sky.cpp
    src/sun.cpp
    src/model/espenak_meeus.cpp
    src/model/iau1980.cpp
    src/model/laskar.cpp
    src/model/meeus.cpp
    src/model/saemundsson.cpp
    src/model/spk.cpp
    src/model/vsop87d.cpp
    ${SPP_VSOP_INC}
)
add_library(spp::spp ALIAS spp)

target_include_directories(spp
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${SPP_GENERATED_DIR}
)
target_compile_features(spp PUBLIC cxx_std_20)
# Debug and release archives can share one install prefix: libsppd.a / libspp.a.
set_target_properties(spp PROPERTIES DEBUG_POSTFIX d)
spp_apply_flags(spp)

# ELP/MPP02 uses compiled, untruncated coefficient tables. A solar-only build
# can omit the model, its generator and its data entirely.
if(SPP_BUILD_MOON)
    set(SPP_ELP_GEN_EXECUTABLE "" CACHE FILEPATH "Prebuilt host elp_gen for cross-compiling without an emulator")
    if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
        if(SPP_ELP_GEN_EXECUTABLE)
            set(SPP_ELP_GEN ${SPP_ELP_GEN_EXECUTABLE})
        else()
            set(SPP_ELP_GEN ${PROJECT_BINARY_DIR}/host/elp_gen${CMAKE_HOST_EXECUTABLE_SUFFIX})
            add_custom_command(
                OUTPUT ${SPP_ELP_GEN}
                COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/host
                COMMAND ${SPP_HOST_CXX} -std=c++20 -O2 ${PROJECT_SOURCE_DIR}/tools/elp_gen.cpp -o ${SPP_ELP_GEN}
                DEPENDS ${PROJECT_SOURCE_DIR}/tools/elp_gen.cpp
                COMMENT "Building elp_gen for the host with ${SPP_HOST_CXX}"
                VERBATIM
            )
        endif()
        set(SPP_ELP_GEN_DEPENDS ${SPP_ELP_GEN})
    else()
        add_executable(elp_gen tools/elp_gen.cpp)
        spp_apply_flags(elp_gen)
        set(SPP_ELP_GEN elp_gen)
        set(SPP_ELP_GEN_DEPENDS elp_gen)
    endif()
    file(GLOB SPP_ELP_DATA CONFIGURE_DEPENDS
        ${PROJECT_SOURCE_DIR}/data/elp_mpp02/elp_main.*
        ${PROJECT_SOURCE_DIR}/data/elp_mpp02/elp_pert.*)
    set(SPP_ELP_INC ${SPP_GENERATED_DIR}/spp/elp_mpp02_terms.inc)
    add_custom_command(
        OUTPUT ${SPP_ELP_INC}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${SPP_GENERATED_DIR}/spp
        COMMAND ${SPP_ELP_GEN} ${PROJECT_SOURCE_DIR}/data/elp_mpp02 ${SPP_ELP_INC}
        DEPENDS ${SPP_ELP_GEN_DEPENDS} ${SPP_ELP_DATA}
        COMMENT "Generating full ELP/MPP02 lunar tables"
        VERBATIM
    )
    target_sources(spp PRIVATE src/moon.cpp src/moon_observe.cpp src/moon_events.cpp
        src/model/elp_mpp02.cpp src/model/iau_moon.cpp src/model/lunar_numerical.cpp ${SPP_ELP_INC})
endif()

# ------------------------------------------------------------------------------
# Executables
# ------------------------------------------------------------------------------
# The library is C++20; executables use C++23 for std::print. Some standard
# libraries keep part of std::print in a separate library (libstdc++exp for
# libstdc++, on some platforms and versions), whatever the compiler. Test-link
# a std::print program, and add libstdc++exp only if it is what makes it link.
include(CheckCXXSourceCompiles)
set(SPP_PRINT_TEST "#include <print>\nint main() { std::print(\"{}\", 1); }")
set(SPP_LIBRARY_CXX_STANDARD ${CMAKE_CXX_STANDARD})
set(CMAKE_CXX_STANDARD 23)   # the checks compile with the executables' standard
check_cxx_source_compiles("${SPP_PRINT_TEST}" SPP_PRINT_LINKS)
set(SPP_PRINT_LIBRARIES "")
if(NOT SPP_PRINT_LINKS)
    set(CMAKE_REQUIRED_LIBRARIES stdc++exp)
    check_cxx_source_compiles("${SPP_PRINT_TEST}" SPP_PRINT_LINKS_WITH_STDCXXEXP)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(SPP_PRINT_LINKS_WITH_STDCXXEXP)
        set(SPP_PRINT_LIBRARIES stdc++exp)
    else()
        message(WARNING "std::print does not link with this standard library; the executables will not build")
    endif()
endif()
set(CMAKE_CXX_STANDARD ${SPP_LIBRARY_CXX_STANDARD})

function(spp_add_executable target)
    add_executable(${target} ${ARGN})
    set_target_properties(${target} PROPERTIES CXX_STANDARD 23)
    target_link_libraries(${target} PRIVATE spp::spp ${SPP_PRINT_LIBRARIES})
    spp_apply_flags(${target})
endfunction()

if(SPP_BUILD_APPS)
    spp_add_executable(spp_cli apps/spp/main.cpp)
    set_target_properties(spp_cli PROPERTIES OUTPUT_NAME spp)
    target_compile_definitions(spp_cli PRIVATE SPP_CLI_MOON=$<BOOL:${SPP_BUILD_MOON}>)

    # Sunrise/sunset azimuth spirals: built when GLFW, GLEW and GLM are available.
    find_package(glfw3 CONFIG QUIET)
    find_package(GLEW QUIET)
    if(NOT TARGET GLEW::GLEW)
        find_package(glew CONFIG QUIET)
        if(TARGET GLEW::glew)
            add_library(GLEW::GLEW ALIAS GLEW::glew)
        endif()
    endif()
    find_package(glm CONFIG QUIET)
    find_package(OpenGL QUIET)
    if(TARGET glfw AND TARGET GLEW::GLEW AND TARGET glm::glm AND TARGET OpenGL::GL)
        spp_add_executable(spp_spiral apps/spiral/main.cpp)
        target_link_libraries(spp_spiral PRIVATE glfw GLEW::GLEW glm::glm OpenGL::GL)
        target_include_directories(spp_spiral PRIVATE ${PROJECT_SOURCE_DIR}/apps/common)
    else()
        message(STATUS "spp_spiral not built: GLFW, GLEW, GLM or OpenGL not found")
    endif()
    if(SPP_BUILD_MOON AND TARGET glfw AND TARGET GLEW::GLEW AND TARGET OpenGL::GL)
        find_package(Threads REQUIRED)
        spp_add_executable(spp_moon apps/moon/main.cpp)
        target_include_directories(spp_moon PRIVATE ${PROJECT_SOURCE_DIR}/apps/common)
        target_link_libraries(spp_moon PRIVATE glfw GLEW::GLEW OpenGL::GL Threads::Threads)
        if(SPP_INSTALL)
            install(TARGETS spp_moon RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
            install(FILES data/licenses/DejaVu.txt DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/Spp)
        endif()
        if(UNIX AND NOT APPLE)
            find_package(OpenGL QUIET COMPONENTS EGL)
        endif()
        if(TARGET OpenGL::EGL)
            target_link_libraries(spp_moon PRIVATE OpenGL::EGL)
            target_compile_definitions(spp_moon PRIVATE SPP_MOON_EGL=1)
        else()
            target_compile_definitions(spp_moon PRIVATE SPP_MOON_EGL=0)
        endif()
    elseif(SPP_BUILD_MOON)
        message(STATUS "spp_moon not built: GLFW, GLEW or OpenGL not found")
    endif()
endif()

if(SPP_BUILD_VALIDATION)
    # Uses internal headers (src/) and the generated tables.
    spp_add_executable(spp_cache_check tools/cache_check.cpp)
    target_include_directories(spp_cache_check PRIVATE ${PROJECT_SOURCE_DIR}/src ${SPP_GENERATED_DIR})

    # The whole pipeline against an extended-precision evaluation of the same models.
    spp_add_executable(spp_pipeline_check tools/pipeline_check.cpp)
    target_include_directories(spp_pipeline_check PRIVATE ${PROJECT_SOURCE_DIR}/src ${SPP_GENERATED_DIR})

    # events::day against a dense sampling of its own condition, 60°-90°.
    spp_add_executable(spp_day_scan tools/day_scan.cpp)

    # The figures the README quotes about events, sidereal time and time scales.
    spp_add_executable(spp_readme_numbers tools/readme_numbers.cpp)

    # Library side of tools/reference_check.py.
    spp_add_executable(spp_dump tools/spp_dump.cpp)

    if(SPP_BUILD_MOON)
        spp_add_executable(spp_lunar_check tools/lunar_check.cpp)
        enable_testing()
        add_test(NAME spp_lunar_reference COMMAND spp_lunar_check
            ${PROJECT_SOURCE_DIR}/tools/reference/elp_mpp02_states.csv)
        spp_add_executable(spp_moon_sky_check tools/moon_sky_check.cpp)
        target_include_directories(spp_moon_sky_check PRIVATE ${PROJECT_SOURCE_DIR}/src)
        add_test(NAME spp_moon_sky_reference COMMAND spp_moon_sky_check
            ${PROJECT_SOURCE_DIR}/tools/reference/moon_horizons.csv
            ${PROJECT_SOURCE_DIR}/tools/reference/moon_usno_phases.csv)
        spp_add_executable(spp_lunar_cache_check tools/lunar_cache_check.cpp)
        add_test(NAME spp_lunar_cache COMMAND spp_lunar_cache_check)
    endif()
endif()

if(SPP_BUILD_BENCH)
    spp_add_executable(spp_spk_bench tools/spk_bench.cpp)
    target_include_directories(spp_spk_bench PRIVATE ${PROJECT_SOURCE_DIR}/src)
    spp_add_executable(spp_earth_cache_bench tools/earth_cache_bench.cpp)
    target_include_directories(spp_earth_cache_bench PRIVATE ${PROJECT_SOURCE_DIR}/src)
    spp_add_executable(spp_bench tools/bench.cpp)
    target_compile_definitions(spp_bench PRIVATE SPP_BENCH_MOON=$<BOOL:${SPP_BUILD_MOON}>
        SPP_BENCH_CONFIG="$<CONFIG>")
    if(SPP_BUILD_MOON AND SPP_BUILD_VALIDATION)
        find_package(Threads REQUIRED)
        spp_add_executable(spp_moon_worker_bench tools/moon_worker_bench.cpp)
        target_include_directories(spp_moon_worker_bench PRIVATE ${PROJECT_SOURCE_DIR}/apps/common)
        target_link_libraries(spp_moon_worker_bench PRIVATE Threads::Threads)
        add_test(NAME spp_moon_worker_performance COMMAND spp_moon_worker_bench)
        set_tests_properties(spp_moon_worker_performance PROPERTIES TIMEOUT 30)
        spp_add_executable(spp_integrator_bench tools/integrator_bench.cpp)
        target_include_directories(spp_integrator_bench PRIVATE ${PROJECT_SOURCE_DIR}/src)
        target_compile_definitions(spp_integrator_bench PRIVATE SPP_GRAVITY_REFERENCE="${PROJECT_SOURCE_DIR}/tools/reference/gravity_de440.csv")
        add_test(NAME spp_integrator_performance COMMAND spp_integrator_bench)
        if(SPP_MOON_REFERENCE_KERNEL)
            add_test(NAME spp_spk_performance COMMAND spp_spk_bench "${SPP_MOON_REFERENCE_KERNEL}"
                ${PROJECT_SOURCE_DIR}/tools/reference/spk_de440s.csv)
            add_test(NAME spp_integrator_lunar_performance COMMAND spp_integrator_bench "${SPP_MOON_REFERENCE_KERNEL}")
        endif()
        spp_add_executable(spp_moon_event_bench tools/moon_event_bench.cpp)
        target_include_directories(spp_moon_event_bench PRIVATE ${PROJECT_SOURCE_DIR}/src)
        add_test(NAME spp_moon_event_performance COMMAND spp_moon_event_bench)
        spp_add_executable(spp_moon_phase_bench tools/moon_phase_bench.cpp)
        target_include_directories(spp_moon_phase_bench PRIVATE ${PROJECT_SOURCE_DIR}/src)
        add_test(NAME spp_moon_phase_performance COMMAND spp_moon_phase_bench
            ${PROJECT_SOURCE_DIR}/tools/reference/moon_horizons_phase_tt.csv)
        spp_add_executable(spp_moon_rotation_bench tools/moon_rotation_bench.cpp)
        target_include_directories(spp_moon_rotation_bench PRIVATE ${PROJECT_SOURCE_DIR}/src)
        add_test(NAME spp_moon_rotation_performance COMMAND spp_moon_rotation_bench
            ${PROJECT_SOURCE_DIR}/tools/reference/moon_iau_rotation.csv
            ${PROJECT_SOURCE_DIR}/tools/reference/moon_rotation_horizons.csv)
        # Explicit opt-in report target: build prerequisites, then run every
        # accuracy stage and the Moon-only benchmark, including after failures.
        add_custom_target(spp_moon_report
            COMMAND ${CMAKE_COMMAND}
                "-DSPP_LUNAR_CHECK=$<TARGET_FILE:spp_lunar_check>"
                "-DSPP_SKY_CHECK=$<TARGET_FILE:spp_moon_sky_check>"
                "-DSPP_CACHE_CHECK=$<TARGET_FILE:spp_lunar_cache_check>"
                "-DSPP_BENCH=$<TARGET_FILE:spp_bench>"
                "-DSPP_EARTH_BENCH=$<TARGET_FILE:spp_earth_cache_bench>"
                "-DSPP_EVENT_BENCH=$<TARGET_FILE:spp_moon_event_bench>"
                "-DSPP_PHASE_BENCH=$<TARGET_FILE:spp_moon_phase_bench>"
                "-DSPP_ROTATION_BENCH=$<TARGET_FILE:spp_moon_rotation_bench>"
                "-DSPP_INTEGRATOR_BENCH=$<TARGET_FILE:spp_integrator_bench>"
                "-DSPP_SPK_BENCH=$<TARGET_FILE:spp_spk_bench>"
                "-DSPP_MOON_REFERENCE_KERNEL=${SPP_MOON_REFERENCE_KERNEL}"
                "-DSPP_REFERENCE_DIR=${PROJECT_SOURCE_DIR}/tools/reference"
                -P ${PROJECT_SOURCE_DIR}/cmake/MoonReport.cmake
            DEPENDS spp_lunar_check spp_moon_sky_check spp_lunar_cache_check spp_bench
                spp_earth_cache_bench spp_moon_event_bench spp_moon_phase_bench spp_moon_rotation_bench
                spp_integrator_bench spp_spk_bench
                ${PROJECT_SOURCE_DIR}/cmake/MoonReport.cmake
            COMMENT "Moon accuracy and speed report"
            USES_TERMINAL
            VERBATIM)
    endif()
endif()

if(SPP_BUILD_FIGURES)
    spp_add_executable(spp_figures tools/figures.cpp)
    target_include_directories(spp_figures PRIVATE ${PROJECT_SOURCE_DIR}/apps/common)
endif()

# ------------------------------------------------------------------------------
# C interface: include/spp/spp_c.h, built as a shared library that contains the
# whole engine (spp is linked in statically) and exports only the spp_c_* and
# spp_* functions of the header.
# ------------------------------------------------------------------------------
if(SPP_BUILD_CAPI)
    # The static library is linked into a shared one.
    set_target_properties(spp PROPERTIES POSITION_INDEPENDENT_CODE ON)

    # Version string: project version + source revision ("-dirty" with uncommitted changes).
    set(SPP_SOURCE_REVISION "unknown")
    find_package(Git QUIET)
    if(GIT_FOUND)
        execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short=7 HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE git_result
            OUTPUT_VARIABLE git_revision OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        if(git_result EQUAL 0)
            set(SPP_SOURCE_REVISION ${git_revision})
            execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=no
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE git_changes ERROR_QUIET)
            if(NOT git_changes STREQUAL "")
                string(APPEND SPP_SOURCE_REVISION "-dirty")
            endif()
        endif()
    endif()

    add_library(spp_capi SHARED capi/spp_c.cpp)
    target_link_libraries(spp_capi PRIVATE spp::spp)
    target_include_directories(spp_capi PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
    target_compile_definitions(spp_capi PRIVATE SPP_C_BUILD
        SPP_C_VERSION_STRING="${PROJECT_VERSION}+${SPP_SOURCE_REVISION}")
    # spp.dll on Windows, spp_capi.so on Linux, spp_capi.dylib on macOS (no "lib" prefix).
    set_target_properties(spp_capi PROPERTIES PREFIX "" RUNTIME_OUTPUT_NAME spp CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
    spp_apply_flags(spp_capi)

    if(MINGW)
        # Link libgcc, libstdc++ and winpthread into the DLL, so it depends only
        # on Windows system DLLs (KERNEL32, the UCRT). Release builds strip the
        # symbol table and debug sections.
        target_link_options(spp_capi PRIVATE -static $<$<CONFIG:Release>:-s>)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT APPLE)
        # Keep the symbols of the static library out of the export table.
        target_link_options(spp_capi PRIVATE -Wl,--exclude-libs,ALL)
    endif()
    if(MSVC)
        # The C runtime linked in statically, for the same reason.
        set_target_properties(spp spp_capi PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()

    # A C program that checks the interface: layout, error codes, reference
    # results and the event definitions. Prints PASS or FAIL, exits with 0 or 1.
    enable_language(C)
    add_executable(spp_capi_test capi/spp_capi_test.c)
    set_target_properties(spp_capi_test PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED ON)
    target_link_libraries(spp_capi_test PRIVATE spp_capi)
    if(MINGW)
        target_link_options(spp_capi_test PRIVATE -static-libgcc)
    endif()

    if(SPP_INSTALL)
        install(TARGETS spp_capi
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
endif()

# ------------------------------------------------------------------------------
# Install: the static library, its public headers and a CMake package.
# Consumers:
#     find_package(Spp 1.0 CONFIG REQUIRED)
#     target_link_libraries(app PRIVATE spp::spp)
# ------------------------------------------------------------------------------
if(SPP_INSTALL)
    include(CMakePackageConfigHelpers)
    set(SPP_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/Spp)

    install(TARGETS spp EXPORT SppTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(DIRECTORY include/spp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(EXPORT SppTargets NAMESPACE spp:: DESTINATION ${SPP_CMAKE_DIR})

    configure_package_config_file(cmake/SppConfig.cmake.in
        ${PROJECT_BINARY_DIR}/SppConfig.cmake
        INSTALL_DESTINATION ${SPP_CMAKE_DIR})
    write_basic_package_version_file(${PROJECT_BINARY_DIR}/SppConfigVersion.cmake
        COMPATIBILITY SameMajorVersion)
    install(FILES
        ${PROJECT_BINARY_DIR}/SppConfig.cmake
        ${PROJECT_BINARY_DIR}/SppConfigVersion.cmake
        DESTINATION ${SPP_CMAKE_DIR})
    install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/Spp)
    install(FILES data/licenses/ERFA.txt DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/Spp)
    if(SPP_BUILD_MOON)
        install(FILES data/elp_mpp02/README.md data/elp_mpp02/LICENSE data/elp_mpp02/manifest.json
            DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/Spp/elp_mpp02)
    endif()
endif()
