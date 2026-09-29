cmake_minimum_required(VERSION 3.28)

project(Spp
    VERSION 0.1.0
    DESCRIPTION "High-precision solar positional astronomy library"
    LANGUAGES CXX
)

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------
# Apps, tools and install rules default to ON only when Spp is the top-level
# project, so add_subdirectory / FetchContent users get just the library.
option(SPP_NATIVE           "Tune for the build machine (-march=native)"                 OFF)
option(SPP_BUILD_APPS       "Build the spp command-line tool and the spp_spiral viewer"  ${PROJECT_IS_TOP_LEVEL})
option(SPP_BUILD_VALIDATION "Build the spp_cache_check accuracy tool"                    ${PROJECT_IS_TOP_LEVEL})
option(SPP_BUILD_BENCH      "Build the spp_bench benchmark"                              OFF)
option(SPP_BUILD_FIGURES    "Build spp_figures, which renders the README figures"        OFF)
option(SPP_INSTALL          "Generate install rules and the CMake package"               ${PROJECT_IS_TOP_LEVEL})

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
add_executable(vsop_gen tools/vsop_gen.cpp)
spp_apply_flags(vsop_gen)

set(SPP_GENERATED_DIR ${PROJECT_BINARY_DIR}/generated)
set(SPP_VSOP_INC ${SPP_GENERATED_DIR}/spp/vsop87d_earth_terms.inc)

add_custom_command(
    OUTPUT ${SPP_VSOP_INC}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${SPP_GENERATED_DIR}/spp
    COMMAND vsop_gen ${PROJECT_SOURCE_DIR}/data/VSOP87D.ear ${SPP_VSOP_INC}
    DEPENDS vsop_gen ${PROJECT_SOURCE_DIR}/data/VSOP87D.ear
    COMMENT "Generating VSOP87D Earth tables"
    VERBATIM
)

# ------------------------------------------------------------------------------
# Library
# ------------------------------------------------------------------------------
add_library(spp STATIC
    src/core/calendar.cpp
    src/core/time.cpp
    src/earth.cpp
    src/earth_state.cpp
    src/events.cpp
    src/sky.cpp
    src/sun.cpp
    src/model/espenak_meeus.cpp
    src/model/iau1980.cpp
    src/model/laskar.cpp
    src/model/meeus.cpp
    src/model/saemundsson.cpp
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

# ------------------------------------------------------------------------------
# Executables
# ------------------------------------------------------------------------------
# The library is C++20; executables use C++23 for std::print, which needs
# libstdc++exp with GCC.
function(spp_add_executable target)
    add_executable(${target} ${ARGN})
    set_target_properties(${target} PROPERTIES CXX_STANDARD 23)
    target_link_libraries(${target} PRIVATE spp::spp)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_link_libraries(${target} PRIVATE stdc++exp)
    endif()
    spp_apply_flags(${target})
endfunction()

if(SPP_BUILD_APPS)
    spp_add_executable(spp_cli apps/spp/main.cpp)
    set_target_properties(spp_cli PROPERTIES OUTPUT_NAME spp)

    # Sunrise/sunset azimuth spirals: built when GLFW, GLEW and GLM are available.
    find_package(glfw3 CONFIG QUIET)
    find_package(glew CONFIG QUIET)
    find_package(glm CONFIG QUIET)
    find_package(OpenGL QUIET)
    if(TARGET glfw AND TARGET GLEW::glew AND TARGET glm::glm AND TARGET OpenGL::GL)
        spp_add_executable(spp_spiral apps/spiral/main.cpp)
        target_link_libraries(spp_spiral PRIVATE glfw GLEW::glew glm::glm OpenGL::GL)
        target_include_directories(spp_spiral PRIVATE ${PROJECT_SOURCE_DIR}/apps/common)
    else()
        message(STATUS "spp_spiral not built: GLFW, GLEW, GLM or OpenGL not found")
    endif()
endif()

if(SPP_BUILD_VALIDATION)
    # Uses internal headers (src/) and the generated tables.
    spp_add_executable(spp_cache_check tools/cache_check.cpp)
    target_include_directories(spp_cache_check PRIVATE ${PROJECT_SOURCE_DIR}/src ${SPP_GENERATED_DIR})
endif()

if(SPP_BUILD_BENCH)
    spp_add_executable(spp_bench tools/bench.cpp)
endif()

if(SPP_BUILD_FIGURES)
    spp_add_executable(spp_figures tools/figures.cpp)
    target_include_directories(spp_figures PRIVATE ${PROJECT_SOURCE_DIR}/apps/common)
endif()

# ------------------------------------------------------------------------------
# Install: the static library, its public headers and a CMake package.
# Consumers:
#     find_package(Spp 0.1 CONFIG REQUIRED)
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
        COMPATIBILITY SameMinorVersion)
    install(FILES
        ${PROJECT_BINARY_DIR}/SppConfig.cmake
        ${PROJECT_BINARY_DIR}/SppConfigVersion.cmake
        DESTINATION ${SPP_CMAKE_DIR})
    install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/Spp)
endif()
